# LDSC

Wrapper pipeline based on https://github.com/bulik/ldsc.

Most of the work is done by the wdl itself, but some preprocessing steps are needed, mainly due to the fact that the nature of the input sumstats can be different.

## Preprocressing
### Step1 - Munge data

In this step we need to make sure that input sumstats are coherent with the requirements by ldsc for its own munging step.

The required input format is as follows:

```
SNP	A1	A2	BETA	P
rs74337086	A	G	0.0923	0.5059
rs76388980	A	G	0.1227	0.2945
rs562172865	T	C	-0.0262	0.8142
rs780596509	A	G	-0.2202	0.1545
rs778009914	A	G	-0.3938	0.3044
rs564223368	T	C	0.2195	0.03913
rs71628921	C	A	0.1763	0.3682
rs577189614	A	G	0.0845	0.5341
rs77357188	T	C	-0.0414	0.3383
```

For FG sumstats the command
```
gunzip -c $SUMSTATS | \
awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++) a[$i]=i; print "SNP","A1","A2","BETA","P"}   NR>1 {if($a["rsids"]!="") print $a["rsids"],$a["alt"],$a["ref"],$a["beta"],$a["pval"]}' | \
awk 'BEGIN{FS=OFS="\t"} {n=split($1,a,","); for(i=1;i<=n;i++) print a[i],$0}' | \
cut -f1,3- | gzip > $PHENO.premunge.gz
```

will do the trick.

To simplify, one can simply run the `munge_fg.wdl` provided updating the ```munge_fg.meta_fg``` input file, which is a tsv that contains `PHENO\tPATH`:
```
C3_BREAST_EXALLC	gs://finngen-production-library-green/finngen_R5/finngen_R5_analysis_data/summary_stats/release/finngen_R5_C3_BREAST_EXALLC.gz
C3_BRONCHUS_LUNG_EXALLC	gs://finngen-production-library-green/finngen_R5/finngen_R5_analysis_data/summary_stats/release/finngen_R5_C3_BRONCHUS_LUNG_EXALLC.gz
C3_PROSTATE_EXALLC	gs://finngen-production-library-green/finngen_R5/finngen_R5_analysis_data/summary_stats/release/finngen_R5_C3_PROSTATE_EXALLC.gz
G6_PARKINSON	gs://finngen-production-library-green/finngen_R5/finngen_R5_analysis_data/summary_stats/release/finngen_R5_G6_PARKINSON.gz
```

### Step2 - Metadata table

The pipeline requres a metadata table for each input list of sumstats to be processed. The input file needs to be structured as `PHENO\tPATH\tN` where `N` is the total number of valid cases+controls of each pheno.
```
C3_BREAST_EXALLC	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-0/C3_BREAST_EXALLC.premunge.gz	110611
C3_BRONCHUS_LUNG_EXALLC	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-1/C3_BRONCHUS_LUNG_EXALLC.premunge.gz	180418
C3_PROSTATE_EXALLC	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-2/C3_PROSTATE_EXALLC.premunge.gz	83146
G6_PARKINSON	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-3/G6_PARKINSON.premunge.gz	224566
H7_AMD	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-4/H7_AMD.premunge.gz	214660
```

If needed, the script `build_tables.py` helps in this matter.
By running `python3 build_tables.py --endpoints GZIPPED_ENDPOINT_FILE --sumstats LIST_OF_PHENOS -o OUTPATH` a file called `count.txt` is generated where the structure is `PHENO\tN\CASES\CONTROLS\NAs`
```
C3_BREAST_EXALLC	110611	8597	102014	113955
C3_BRONCHUS_LUNG_EXALLC	180418	1735	178683	44148
C3_PROSTATE_EXALLC	83146	6464	76682	141420
G6_PARKINSON	224566	2207	222359	0
H7_AMD	214660	3867	210793	9906
H7_CATARACTSENILE	222085	27421	194664	2481
H7_GLAUCPRIMOPEN	220275	4558	215717	4291
I9_AF	142885	22559	120326	81681
I9_MI_STRICT	204766	11909	192857	19800
J10_ASTHMA	160321	21128	139193	64245
```

One can then easily use the `join` command to build the table.E.g.
`while read f; do PHENO=$(basename $f .premunge.gz) && echo -e "$PHENO\t$f" | join - -t $'\t' <(cut -f 1,2 counts.txt); done < flagship_paths_munged.txt > input_meta.txt`


## WDL time!

Now all the data is ready to run the big run. A brief summary of the logic of the wdl.

As input one cane have two lists of sumstats or just one. If two are provided then all the cross scores are computed between the two lists. If only one is passed instead the first list is duplicated as a second, thus running an inner product on itself. The total number of comparison that needs to run is `N*L/2` jobs to be run where `N` and `L` are the lengths of the two lists. In case of inner product, the number is `N(N-1)/2` instead. This means that the growht is quadratic and thus I recommed first testing the pipeline with a smaller set of sumstats.

In the `ldsc_sandbox.json` file the comparison list is by default a test (10 phenotypes) version with the latest (r9) data. Please run first in test mode and once sure that everything works, remove the suffix `_test` and you can run the full version.

`filter_meta` splits the input lists into chunks for munging. The number of chunks in this step is given by `  "ldsc_rg.filter_meta.filter_chunks": Int`. This step is quite fast anyways and, in principle, should only run once since its output is cached.

Each input chunk list is passed to `munge_ldsc`. In this step, the input sumstats are munged by ldsc directly and, while we're at it, the heritability is calculated. The outputs of the heritability calculation are passed to `gather_h2` which builds an output table and json with all the summaries, as well as merging all the logs.

`return_couples` is the step that harmonizes the two input lists. It checks for duplicates in the two lists and makes sure that the smallest subset of pairwise comparisons are actually ran. The total list of couple pheno comparisons is then split into chunks (defined by  the input ` "ldsc_rg.return_couples.chunks": Int `) producing also a list of sumstats, so that only the requried sumstats are localized.

The list of pairwise couples and list of paths is pased to `multi_rg` where a wrapper script runs the `ldsc.py --rg` command in a parallelized way. Therefore, increasing the number of CPUs for the job will lead to faster running times.

Finally, the output and logs of the step are passed to `gather_summaries` which is the analagous to the previous gather_h2 step:
```
p1	p2	rg	se	z	p	h2_obs	h2_obs_se	h2_int	h2_int_se	gcov_int	gcov_int_se
AB1_AMOEBIASIS	AB1_AMOEBIASIS	1.0006	0.0009	1139.1313	0.0	0.0009	0.0013	0.9682	0.0061	0.9682	0.0061
AB1_AMOEBIASIS	AB1_ANOGENITAL_HERPES_SIMPLEX	0.5613	0.7023	0.7992	0.4242	0.0039	0.0015	0.9872	0.0069	0.0007	0.0044
AB1_AMOEBIASIS	AB1_ARTHROPOD	-1.0562	1.0868	-0.9718	0.3311	0.0024	0.0015	0.9977	0.0071	0.0089	0.005
AB1_AMOEBIASIS	AB1_ASPERGILLOSIS	0.7488	0.9701	0.7719	0.4402	0.0019	0.0016	0.981	0.0065	-0.0074	0.0048
```

### Other inputs

The choice of LD scores to be used is between the default European 1kg and the custom built Finngen one. This can be chosen in the json inputs by selection ``` "ldsc_rg.population" ``` as being either ```eur|fin``` as the ```ld_path``` root path will be  automatically updated to download the desired files.

FInally, ```dsc_rg.multi_rg.args``` is an optional input that allows to add flags to the default ldsc command.