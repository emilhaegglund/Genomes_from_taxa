# Python script to download bacterial genomes from a taxonomic id
Still under development!

## Usage
Example to download all refseq genomes of planctomycetes:
```
python download_genomes --taxid 203682 --output planctomycetes --filter complete,scaffold,contig --workers 12
```
```
Download genomes from NCBI

optional arguments:
  -h, --help            show this help message and exit
  --taxid TAXID, -t TAXID
                        Taxonomical ID
  --update UPDATE, -u UPDATE
                        Update
  --output OUTPUT, -o OUTPUT
                        Output directory to store downloaded genomes
  --filter FILTER, -f FILTER
                        Filter on assembly status. Allowed values are
                        complete, scaffold, contig. Example use -f
                        complete,contig will download genomes with Complete
                        genome or contig as assembly status
  --workers WORKERS, -w WORKERS
                        Number of parallel downloads
```

## Todo
[] Add switch for genbank or refseq
[] Add an all argument for the filter switch
[] Fix the update function

