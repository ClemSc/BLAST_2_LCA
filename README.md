# LCA4BLAST

LCA4BLAST performs Last Common Ancestor (LCA) estimation from Blast search results.

## Installation
LCA4BLAST is a stand-alone module that requires Python3 and the Pandas library, which can be installed with pip:
pip install pandas
or conda
conda install pandas
## Use
### Prerequisite.
LCA4BLAST use as inputs a local blast search result file and NCBI taxonomy information files: nodes.dmp and names.dmp which are packaged in taxdump:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
Input files can be provided uncompressed or gzipped.

Typically, the blast search command should look like this:

```console
blastn -query reads.fasta -task megablast -db nt -out blast_otus.tsv -outfmt "6 qseqid saccver pident qcovs length evalue bitscore staxid" -num_threads 8 -evalue 1e-05
```
Mind that your blast database must be built with taxid references. Typically with:
```console
makeblastdb -in nt.fasta -dbtype nucl -out nt -taxid_map acc2taxid_map_file -parse_seqids
```

Adapt to your liking, but outfmt needs to be '6' and the following fields are mandatory: "qseqid pident length evalue staxid". Any extra fields will be shown in the final result, with values corresponding to the best hit for each surviving queries.

### Example command:
```console
LCA4BLAST/LCA4BLAST.py -t nodes.dmp.gz -n names.dmp.gz -i blast_otus.tsv.gz -o RESULTS_LCA.tsv -L 80 -H 95 -p 90 -l 350  -f "qseqid saccver pident qcovs length evalue bitscore staxid"
```

### Mandatory parameters:
- Input and output filename are mandatory and are set with `-t/--nodes`, `-n/--names`, `-i/--input`, `-o/--output` (check the example command).
`-f/--fields` = Blast output fields. Use the same values and same order as in your Blast search command `-outfmt` without the leading "6", e.g. "qseqid saccver pident qcovs length evalue bitscore staxid"

### Optional parameters (a default value will be used when not specified):
- `-H/-high_pident` = high similary threshold (percentage), default is 80
- `-L/--low_pident` = low similary threshold (percentage), default is 95
  LCA4BLAST uses a high similarity and a low similarity threshold to decide whether a species level identification is relevant or not. Only queries for which at least one hit was found with percentage identity (pident) > H will get a chance to be identified to species. Otherwise, genus will be used as the lower rank possible. Hits with pident < L are not used. In general, H should be set where the user estimate the species gap to be.

- `-p/--p_hits` = percentage of hits thresholds, default is 90.
Blast database is not perfect, and contains misidentified sequences. This parameter specify over which majority threshold a taxa should be accepted. E.g. if p = 90 and 95% of the hits are assigned to Megalothorax minimus, 5% to Megalothorax willemi, then the query will be assigned to Megalothorax minimus. If only 85% of the hits are assigned to M. minimus and 15 to M. willemi, then the query will be assigned to "Megalothorax".

- `-l/--length` = minimum alignment length, default is 350.
The minimum alignment length query/subject to be accepted.




