# TaxAll
A fast program for assigning taxonomy to any sequence (basically running diamond blast in bulk with a custom database). 

## Setup
1. Installing the custom diamond database: The custom diamond database is unfortunately too large to host publicly. The only difference between the custom database and the RefSeq database available from NCBI is the inclusion of many more NCLDV (giant virus) genomes. If you aren't interested in categorizing GV genomes, you can go ahead and just download refseq along with the accompanying names.dmp and nodes.dmp files. Taxall should work fine with this database. If you want the custom database, reach out via email bsm122@miami.edu.
2. Install dependencies: Pandas and diamon blast are the only two.

## Running Taxall

`python taxall.py -i proteins.faa -o taxall_output.csv -d refseq.dmnd -names names.dmp -nodes nodes.dmp -n 5`
