from mod import MOD
import xlrd
import csv


class FlyBase(MOD):
    species = "Drosophila melanogaster"

    @staticmethod
    def gene_href(gene_id):
        return "http:/flybase.org/reports/" + gene_id + ".html"

    @staticmethod
    def gene_id_from_panther(panther_id):
        # example: WormBase=WBGene00004831
        return panther_id.split("=")[1]

    def load_genes(self):
        genes = MOD.genes

        go_data_csv_filename = "./data/FlyGene.tsv"

        print("Fetching go data from FlyGene tsv file...")

        with open(go_data_csv_filename, 'rb') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if row[3] == "":
                    gene_type = None
                else:
                    gene_type = row[3]

                chromosomes = []
                if row[5]:
                    chromosomes = [row[5]]

                genes[row[0]] = {
                    "id": row[0],
                    "name": row[1],
                    "gene_symbol": row[1],
                    "species": row[2],
                    "gene_type": gene_type,
                    "description": row[4],
                    "gene_chromosomes": chromosomes,
                    "gene_chromosome_starts": row[6],
                    "gene_chromosome_ends": row[7],
                    "gene_chromosome_strand": row[8],
                    "gene_synonyms": map(lambda s: s.strip(), row[9].split(",")),
                    "name_key": row[1].lower(),
                    "href": FlyBase.gene_href(row[0]),
                    "category": "gene"
                }

    def load_go(self):
        go_data_csv_filename = "./data/FlyGOGeneMapping.tsv"

        print("Fetching go data from FlyBase tsv file...")

        with open(go_data_csv_filename, 'rb') as f:
            reader = csv.reader(f, delimiter='\t')

            for i in xrange(5):
                next(reader, None)

    def load_diseases(self):
        disease_data_csv_filename = "./data/FlyDiseaseGeneMapping.tsv"

        print("Fetching disease data from Fly tsv file...")

        with open(disease_data_csv_filename, 'rb') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)

            for row in reader:
                if row[3] and row[3] != "":
                    omim_ids = map(lambda s: s.strip(), row[3].split(","))

                    for omim_id in omim_ids:
                        self.add_disease_annotation_to_gene(gene_id=None, omim_id="OMIM:"+omim_id)
