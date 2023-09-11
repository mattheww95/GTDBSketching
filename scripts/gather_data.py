"""Script for generating a mash sketch for GTDB

2023-09-11: Matthew Wells
"""
import os
import sys
import subprocess
import argparse
from tqdm import tqdm


class BacterialTaxonomy:
    """Intake the GTDB taxonomy list
    """
    __prefixes = frozenset(["RS_", "GB_"])
    
    def __init__(self, file_path) -> None:
        self.file_path = file_path
        self.data = dict()
        self.read_taxonomy(self.file_path)
    
    def read_taxonomy(self, file_path):
        sys.stderr.write(f"Reading in taxonomy file.\n")
        strip_prefix = lambda x: x[x.index("_")+1:]
        with open(file_path, 'r', encoding="utf8") as taxonomy:
            for i in taxonomy:
                values = i.strip().split("\t")
                self.data[strip_prefix(values[0])] = values[1]

class IndexFiles:
    """Create a list of paths to each assembly the GTDB data dump.

    The assemblies can all be found in the genomic_files_reps/gtdb_genomes.../database folder
    """

    def __init__(self, db_path, suffix_expr=".fna.gz"):
        self.suffix_expr = suffix_expr
        self.file_index = dict()
        self.db_path = db_path
        self.recurse_directory(self.db_path)

    def recurse_directory(self, path):
        sys.stderr.write(f"Indexing files ending with {self.suffix_expr}.\n")
        entries = [*os.scandir(path)]
        clean_name = lambda x: x[:x.rindex("_")]
        try:
            while entry := entries.pop():
                if entry.is_dir(follow_symlinks=False):
                    entries.extend(os.scandir(entry))
                else:
                    if entry.path.endswith(self.suffix_expr):
                        self.file_index[clean_name(entry.name)] = entry.path
        except IndexError:
            pass


class SketchData:

    def __init__(self, database, taxonomy, kmer_size=21, sketch_size=1000, prob_low_kmer=0.01) -> None:
        if kmer_size > 31 or kmer_size < 1:
            sys.stderr.write(f"Inaccurate kmer sizes must be between 1-32.\n")
            sys.exit()
        if sketch_size < 10:
            sys.stderr.write(f"sketch size is too small, must be greater than 10.\n")
            sys.exit()
        if prob_low_kmer > 1.0 or prob_low_kmer == 0.0:
            sys.stderr.write(f"warning low kmer size, must be a float between 0 and 1\n")
            sys.exit()
        
        #self.amino_acid_alph = "-a" if amino_acid_alph else ""
        self.prob_low_kmer = str(prob_low_kmer)
        self.sketch_size = str(sketch_size)
        self.kmer_size = str(kmer_size)
        self.taxonomy = BacterialTaxonomy(taxonomy)
        self.database = IndexFiles(database)
        self.create_commands()

    def create_commands(self):
        """Create mash commands from merged data
        """
        key_errors = 0
        sys.stderr.write(f"Creating mash sketches.\n")
        file_idx_len = len(self.database.file_index)
        with open("error_log.txt", 'w') as errors:
            for key, value in tqdm(self.database.file_index.items()):
                try:
                    comment = self.taxonomy.data[key]
                except KeyError:
                    key_errors +=1
                    errors.write(f"The assembly {key} does not have any taxonomic information attached to it.\n")
                else:
                    #print(["mash", "sketch", "-I", key, "-C",  f'"{comment}"', "-k", self.kmer_size, "-s", self.sketch_size, "-w", self.prob_low_kmer, value]) 
                    proc = subprocess.Popen(["mash", "sketch", "-I", key, "-C",  f'"{comment}"', "-k", self.kmer_size, "-s", self.sketch_size, "-w", self.prob_low_kmer, value], stdout=errors, stderr=errors) 
                    proc.wait()
        print(f"Entry Errors: {key_errors} | This corresponds to {round(key_errors/float(file_idx_len)*100, 2)}% of assemblies missing taxomonomic data and have been excluded from the sketch.")
        print(f"Please check the error log for assemblies missing taxonomic data.")


class BigPaste:
    """Find all sketches to prepare for a pasting together with mash
    """

    __fofn_name = "SketchedAssemblies.txt"
    def __init__(self, fp, sketch_name) -> None:
        self.file_path = fp
        self.sketches = IndexFiles(self.file_path, ".fna.gz.msh")
        self.sketch_name = sketch_name
        self.create_file_list(self.sketches)
        self.mega_paste()
    
    def mega_paste(self):
        """This is where the magic happens.
        """
        sys.stderr.write("Pasting sketches!\n")
        big_proc = subprocess.Popen(["mash", "paste", self.sketch_name, "-l", self.__fofn_name])
        big_proc.wait()


    def create_file_list(self, files_obj):
        with open(self.__fofn_name, 'w') as output_f:
            for path in files_obj.file_index.values():
                output_f.write(f"{path}\n")


def main():
    parser = argparse.ArgumentParser(prog="GTDB Sketching", 
                                    description="Generate a mash sketch with a Kraken2 like taxonomy comment from a GTDB download.")
    parser.add_argument("-t", "--taxonomy", help="A file in the GTDB download containing the taxonomy \
                        strings for the various assemblies. It will be called something along the lines of bac120_taxonomy.tsv. The assembly ID will be prefaced by either RS_ or GB_")
    parser.add_argument("-d", "--database", help="Directory containing the assembled zipped genomes (ending in .fna.gz). \
                        Specify the database to the directory you wish to use to generate your sketch.\
                        A path to the database could like 'data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps_r214/database/GCF/' \
                        to create a sketch of the refseq genomes")    
    parser.add_argument("-n", "--sketch-name", help="The name you would like for the final sketch (ending in '.msh')")
    parser.add_argument("-k", "--kmer-size", help="Kmer size, an integer between 3 and 31, default: 21", type=int, default=21)
    parser.add_argument("-s", "--sketch-size", help="sketch size, an integer greater than 10: default 1000", type=int, default=1000)
    args = parser.parse_args()
    SketchData(database=args.database, taxonomy=args.taxonomy, kmer_size=args.kmer_size, sketch_size=args.sketch_size)
    BigPaste(args.database, args.sketch_name)

if __name__=="__main__":
    taxonomy = "../GTDB_data/data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv"
    database = "../GTDB_data/data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps_r214/database/GCF/"
    #SketchData(database=database, taxonomy=taxonomy)
    main()
    #BigPaste(database, "GTDB_out_20230911.msh")