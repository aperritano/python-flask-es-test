import csv
import pickle
import requests
import os
import re
from elasticsearch import Elasticsearch
import pandas as pd


class MOD():
    INDEX_NAME = 'searchable_items_blue'
    DOC_TYPE = 'searchable_item'

    go_blacklist = ("GO:0008150", "GO:0003674", "GO:0005575")

    gene_bkp_filename = "genes_bkp.pickle"
    go_bkp_filename = "go_bkp.pickle"
    diseases_bkp_filename = "diseases_bkp.pickle"

    go_dataset = {}
    omim_dataset = {}

    genes = {}
    go = {}
    diseases = {}

    def __init__(self):
        self._load_omim_dataset()
        self._load_go_dataset()
        self.es = Elasticsearch(os.environ['ES_URI'], retry_on_timeout=True)

    @staticmethod
    def factory(organism):
        from worm import WormBase
        if organism in ("Caenorhabditis elegans", "C. elegans", "CAEEL"):
            return WormBase()
        else:
            return None

    @staticmethod
    def _process_gene_id_from_panther(gene_ids_panther, genes):
        gene_ids = gene_ids_panther.split("|")
        mod = MOD.factory(gene_ids[0])

        if mod is None:
            return None

        gene_id = MOD.factory(gene_ids[0]).gene_id_from_panther(gene_ids[1])

        gene_symbol = ""
        if mod.__class__.__module__ == "human":
            gene_symbol = gene_id
        else:
            if gene_id not in genes:
                return None
            else:
                gene_symbol = genes[gene_id]["gene_symbol"]

        return {
            "id": gene_id,
            "symbol": gene_symbol,
            "href": mod.gene_href(gene_id),
            "species": mod.species
        }


    def _load_omim_dataset(self):
        if MOD.omim_dataset != {}:
            return

        print "loading OMIM dataset from file..."
        with open("data/OMIM_diseases.txt", "r") as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)
            next(reader, None)
            next(reader, None)

            for row in reader:
                if len(row) < 3:
                    continue

                name_column = row[2].split(";")
                name = name_column[0].strip()
                if len(name_column) > 1:
                    symbol = name_column[1].strip()
                else:
                    symbol = None

                synonyms = []

                for r in (row[3], row[4]):
                    if r == '':
                        continue

                    alternative_names = r.split(";;")

                    for alt_name_symbol in alternative_names:
                        alt_name_symbol = alt_name_symbol.split(";")
                        alt_name = alt_name_symbol[0].strip().lower()
                        if len(alt_name_symbol) > 1:
                            alt_symbol = ", " + alt_name_symbol[1].strip()
                        else:
                            alt_symbol = ""

                        synonyms.append(alt_name + alt_symbol)

                MOD.omim_dataset["OMIM:" + row[1]] = {
                    "prefix": row[0],
                    "name": name.lower(),
                    "symbol": symbol,
                    "disease_synonyms": synonyms
                }

    def _load_go_dataset(self):
        if MOD.go_dataset != {}:
            return

        print "Loading GO dataset from file..."
        with open("data/go.obo", "r") as f:
            creating_term = None

            for line in f:
                line = line.strip()

                if line == "[Term]":
                    creating_term = True
                elif creating_term:
                    key = (line.split(":")[0]).strip()
                    value = ("".join(":".join(line.split(":")[1:]))).strip()

                    if key == "id":
                        creating_term = value
                        MOD.go_dataset[creating_term] = {}
                    else:
                        if key == "synonym":
                            if value.split(" ")[-2] == "EXACT":
                                value = (" ".join(value.split(" ")[:-2]))[1:-1]
                            else:
                                continue
                        if key == "def":
                            m = re.search('\"(.+)\"', value)
                            value = m.group(1)

                        if key in MOD.go_dataset[creating_term]:
                            MOD.go_dataset[creating_term][key].append(value)
                        else:
                            MOD.go_dataset[creating_term][key] = [value]

    def add_go_annotation_to_gene(self, gene_id, go_id):
        if go_id not in self.go_dataset or go_id in MOD.go_blacklist or gene_id not in self.genes:
            return

        gene_symbol = self.genes[gene_id]["gene_symbol"].upper()

        if go_id in self.go:
            if gene_symbol not in self.go[go_id]["go_genes"]:
                self.go[go_id]["go_genes"].append(gene_symbol)
            if self.species not in self.go[go_id]["go_species"]:
                self.go[go_id]["go_species"].append(self.species)
        else:
            self.go[go_id] = {
                "go_genes": [gene_symbol],
                "go_species": [self.species],

                "name": self.go_dataset[go_id]["name"][0],
                "description": self.go_dataset[go_id]["def"][0],
                "go_type": self.go_dataset[go_id]["namespace"][0],
                "go_synonyms": self.go_dataset[go_id].get("synonym"),

                "name_key": self.go_dataset[go_id]["name"][0],
                "id": go_id,
                "href": "http://amigo.geneontology.org/amigo/term/" + go_id,
                "category": "go"
            }

        if self.go[go_id]["name"] not in self.genes[gene_id]["gene_" + self.go[go_id]["go_type"]]:
            self.genes[gene_id]["gene_" + self.go[go_id]["go_type"]].append(self.go[go_id]["name"])

    def add_disease_annotation_to_gene(self, gene_id, omim_id):
        if omim_id not in self.omim_dataset or gene_id not in self.genes:
            return

        gene_symbol = self.genes[gene_id]["gene_symbol"].upper()

        if omim_id in self.diseases:
            if gene_symbol not in self.diseases[omim_id]["disease_genes"]:
                self.diseases[omim_id]["disease_genes"].append(gene_symbol)
            if self.species not in self.diseases[omim_id]["disease_species"]:
                self.diseases[omim_id]["disease_species"].append(self.species)
        else:
            self.diseases[omim_id] = {
                "disease_genes": [gene_symbol],
                "disease_species": [self.species],

                "name": self.omim_dataset[omim_id]["name"],
                "symbol": self.omim_dataset[omim_id]["symbol"],
                "disease_synonyms": self.omim_dataset[omim_id]["disease_synonyms"],

                "name_key": self.omim_dataset[omim_id]["name"],
                "id": omim_id,
                "key": omim_id,
                "href": "http://omim.org/entry/" + omim_id.split(":")[1],
                "category": "disease"
            }

    def load_from_file(self, filename):
        if os.path.isfile(filename):
            with open(filename, "rb") as f:
                return pickle.load(f)
        return None

    def load_data_from_file(self):
        print "Loading genes from file..."
        self.genes = self.load_from_file(self.gene_bkp_filename)

        print "Loading go from file..."
        self.go = self.load_from_file(self.go_bkp_filename)

        print "Loading diseases from file..."
        self.diseases = self.load_from_file(self.diseases_bkp_filename)

        if self.genes is None or self.go is None or self.diseases is None:
            print ("Fail loading data from backup")

    def save_dict_into_file(self, data, filename):
        with open(filename, "wb") as f:
            pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

    def save_into_file(self):
        print "Saving genes into file..."
        self.save_dict_into_file(self.genes, self.gene_bkp_filename)

        print "Saving go into file..."
        self.save_dict_into_file(self.go, self.go_bkp_filename)

        print "Saving diseases into file..."
        self.save_dict_into_file(self.diseases, self.diseases_bkp_filename)

    def delete_mapping(self):
        print "Deleting mapping..."
        response = requests.delete(os.environ['ES_URI'] + self.INDEX_NAME + "/")
        if response.status_code != 200:
            print "ERROR: " + str(response.json())
        else:
            print "SUCCESS"

    def put_mapping(self):
        from mapping import mapping

        print "Putting mapping... "
        response = requests.put(os.environ['ES_URI'] + self.INDEX_NAME + "/", json=mapping)
        if response.status_code != 200:
            print "ERROR: " + str(response.json())
        else:
            print "SUCCESS"

    def index_into_es(self, data):
        bulk_data = []

        for id in data:
            bulk_data.append({
                'index': {
                    '_index': self.INDEX_NAME,
                    '_type': self.DOC_TYPE,
                    '_id': id
                }
            })
            bulk_data.append(data[id])

            if len(bulk_data) == 300:
                self.es.bulk(index=self.INDEX_NAME, body=bulk_data, refresh=True)
                bulk_data = []

        if len(bulk_data) > 0:
            self.es.bulk(index=self.INDEX_NAME, body=bulk_data, refresh=True)

    def index_genes_into_es(self):
        print "Indexing genes into ES..."
        self.index_into_es(self.genes)

    def index_go_into_es(self):
        print "Indexing go into ES..."
        self.index_into_es(self.go)

    def index_diseases_into_es(self):
        print "Indexing diseases into ES..."
        self.index_into_es(self.diseases)

    def index_all_into_es(self):
        self.index_genes_into_es()
        self.index_go_into_es()
        self.index_diseases_into_es()


    def parse_genome_data(file1, file2, file3):
        '''
           Description: Function to parse given files using Pandas
           Parameters: filenames
        '''
        gene_id = ""
        gene_href = 'http://flybase.org/reports/' + gene_id + ".html"
        flygene_Columns = ['id', 'name', 'species', 'gene_type', 'description', 'gene_chromosomes', 'gene_chromosomes_starts','gene_chromosomes_ends', 'gene_chromosomes_strand', 'gene_synonyms']

        flygene_data = pd.read_csv(file1, sep='\t', headers=None)

        flygene_dataframe = pd.DataFrame(flyGOgene_data['Drosophila melanogaster'], columns= flygene_Columns)
        flygene_dataframe.fillna(method='ffill')

        get_gene_id = 'FBgn0000001' ## hard coding this due to lack of time, this can be retrieved from the above dataframe using the gene name.

        flyGOGene_columns = ['name', 'aliases', 'id', 'list of genes']
        flyGOgene_data = pd.read_csv(file2, sep='\t', headers=None)
        flyGOgene_dataframe = pd.DataFrame(flyGOgene_data['FBgn0000001'], columns= flygene_Columns)
        flyGOgene_dataframe.fillna(method='ffill')

        flyDiseaseGene_columns = ['Disease Symbol', 'Disease Name', 'Disease Synonym', 'OMIM ID', 'Description', 'Associated Genes']
        flyDiseaseGene_data = pd.read_csv(file2, sep='\t', headers=None)
        flyDiseaseGene_dataframe = pd.DataFrame(flyGOgene_data['FBgn0000001'], columns= flygene_Columns)
        flyDiseaseGene_dataframe.fillna(method='ffill')

        # combining dataframes side-by-side
        bulk_gene_data = pd.concat([flygene_dataframe, flyGOgene_dataframe, flyDiseaseGene_dataframe], axis=1)

        #Few additional dataclean up steps still needs to be performed, like delete duplicates, check data types and correct them if not used properly, trim spaces in the data

        index_into_es(bulk_gene_data)
