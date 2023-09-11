# GTDB Sketching

## Requirments
Python >= 3.5  
tqdm  
[mash](https://github.com/marbl/Mash)

## Usage
You must first download and extract a GTDB data download, and verify that you have unzipped the taxonomy.tsv file, and the folder containing the representative genomes.

```
Generate a mash sketch with a Kraken2 like taxonomy comment from a GTDB download.

optional arguments:
  -h, --help            show this help message and exit
  -t TAXONOMY, --taxonomy TAXONOMY
                        A file in the GTDB download containing the taxonomy strings for the various assemblies. It will be called something along the lines of bac120_taxonomy.tsv. The assembly ID will be prefaced by either RS_ or GB_
  -d DATABASE, --database DATABASE
                        Directory containing the assembled zipped genomes (ending in .fna.gz). Specify the database to the directory you wish to use to generate your sketch. A path to the database could like
                        'data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps_r214/database/GCF/' to create a sketch of the refseq genomes
  -n SKETCH_NAME, --sketch-name SKETCH_NAME
                        The name you would like for the final sketch (ending in '.msh')
  -k KMER_SIZE, --kmer-size KMER_SIZE
                        Kmer size, an integer between 3 and 31, default: 21
  -s SKETCH_SIZE, --sketch-size SKETCH_SIZE
                        sketch size, an integer greater than 10: default 1000

```