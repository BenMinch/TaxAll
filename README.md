# TaxAll
A fast program for assigning taxonomy to any sequence (basically running diamond blast in bulk with a custom database). 

## Setup
1. Installing the custom diamond database: The custom diamond database is unfortunately too large to host publicly. The only difference between the custom database and the RefSeq database available from NCBI is the inclusion of many more NCLDV (giant virus) genomes. If you aren't interested in categorizing GV genomes, you can go ahead and just download refseq along with the accompanying names.dmp and nodes.dmp files. Taxall should work fine with this database. If you want the custom database, reach out via email bsm122@miami.edu.
2. Install dependencies: Pandas and diamon blast are the only two.

## Running Taxall

`python taxall.py -i proteins.faa -o taxall_output.csv -d refseq.dmnd -names names.dmp -nodes nodes.dmp -n 5`

### Flags
1. -i : Input protein file (amino acid sequences).
2. -o : What you want the output to be called.
3. -d: Path to your database (refseq or custom)
4. -names: Names file (downloaded from NCBI)
5. -nodes: Nodes file (downloaded from NCBI)
6. -n: Number of hits per protein to show. (Default is 3).

## Interpreting the output
Taxall will assign taxonomy down to the species level for every protein in your fasta file. It will also show percent identity and e-value which can be helpful in knowing if it is a good hit or not. Typically, we only use kingdom classifications as these tend to be pretty accurate at distinguishing virus proteins from bacteria or eukarya.
