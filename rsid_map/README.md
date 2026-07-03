# RSID MAPPING

Helper scripts to move from rsid to chrom_pos_ref_alt notation.

## rsid_map.py

```
usage: rsid_map.py [-h] -o OUT --bim BIM [--rsids RSIDS] [--prefix PREFIX]

Process a file

optional arguments:
  -h, --help         show this help message and exit
  -o OUT, --out OUT  folder in which to save the results
  --bim BIM          bim file
  --rsids RSIDS      optional list of rsids
  --prefix PREFIX    prefix of output files for rsid filtering
```

This script downloads (if not provided in the output folder) a huge vcf `ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz" with SNP info and builds a mapping from rsid to chrom_pos.
One can subfilter the mapping to only position one is interested in with the `--rsids` flag. 

## convert_rsids.py

```
usage: convert_rsids.py [-h] -o OUT --file FILE --map MAP [--gz] [--no-header]
                        --metadata [METADATA [METADATA ...]]
                        [--columns [COLUMNS [COLUMNS ...]]]
                        (--to-rsid | --to-chrompos)

Process a file

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --out OUT     folder in which to save the results
  --file FILE, -f FILE  File to map
  --map MAP             Mapping file from rsid to chrompos
  --gz                  Compress output file to gz
  --no-header           Flag to use when no header is present
  --metadata [METADATA [METADATA ...]], -m [METADATA [METADATA ...]]
                        columns required for parsing. Should be SNPID,A1,A2 columns
                        for to-chrompos and SNPID for to-rsid
  --columns [COLUMNS [COLUMNS ...]]
                        column that need to be kept, either numerical integers
                        or column names
  --to-rsid
  --to-chrompos
```

This script takes `FILE` and outputs a new file, keeping the required `COLUMNS`. 
