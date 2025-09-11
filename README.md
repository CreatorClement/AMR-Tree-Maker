# AMR-Tree-Maker
AMR Tree Maker is meant to be ran after AMRFinderPlus [AMRFinderPlus](https://github.com/ncbi/amr). 


Getting Started

While running AMRFinderPlus, please generate the following: 

The --nucleotide_output flag in AMRFinderPlus causes the program to produce FASTA files of the identified antimicrobial resistance (AMR) genes at the nucleotide level. By default, AMRFinderPlus only produces a report, so using this flag allows you to extract the actual nucleotide sequences that were matched. 


Input/output:
  -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        path to query signature
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        location of an output directory
  -p PREFIX, --prefix PREFIX
                        prefix to describe the input sample read files
