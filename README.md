# LDSC

Wrapper pipeline based on https://github.com/bulik/ldsc.

Most of the work is done by the wdl itself, but some preprocessing steps are needed, mainly due to the fact that the nature of the input sumstats can be different.

## WDL
The wdl takes a list of sumstats and generates heritabilites and (optional) genetic correlation between all N(N-1)/2 pairs or iff two separate lists are passed then only between cross N*M pairs.

### Inputs


| Parameter | Description |
|-----------|-------------|
| `ldsc_rg.docker` | Docker image used for pipeline tasks. |
| `ldsc_rg.only_het` | If `true`, computes only heritabilities (not genetic correlations). |
| `ldsc_rg.meta_fg` | Metadata table for primary summary statistics (TSV). |
| `ldsc_rg.comparison_fg` | (Optional) Metadata table for secondary sumstats (for cross-trait analysis). |
| `ldsc_rg.name` | Output prefix for result files. |
| `ldsc_rg.population` | Population label for LD score reference (e.g., "EUR", "FIN"). |
| `ldsc_rg.ld_path` | A map linking population codes to LD score files. |
| `ldsc_rg.filter_meta.filter_chunks` | Number of chunks to split input tables for parallel processing (increase for large datasets). |
| `ldsc_rg.premunge_ss.p_col` | Column name for p-value in your sumstats. |
| `ldsc_rg.premunge_ss.a1_effect_col` | Column name for effect allele. |
| `ldsc_rg.premunge_ss.a2_ne_col` | Column name for non-effect/reference allele. |
| `ldsc_rg.premunge_ss.beta_col` | Column name for effect size (beta). |
| `ldsc_rg.premunge_ss.rsid_col` | Column name for rsIDs. |
| `ldsc_rg.premunge_ss.chrom_col` | (Optional) Column name for chromosome (if not using rsIDs). |
| `ldsc_rg.premunge_ss.pos_col` | (Optional) Column name for position (if not using rsIDs). |
| `ldsc_rg.munge_ldsc.snplist` | Path to snplist file for LD score regression. |
| `ldsc_rg.return_couples.chunks` | Number of parallel batches for genetic correlation computations. |
| `ldsc_rg.multi_rg.cpus` | Number of CPUs to use for multi_rg step. |
| `ldsc_rg.multi_rg.args` | (Optional) Extra arguments passed to ldsc.py. |


The metadata tables should be  structured as `PHENO\tPATH\tN` where `N` is the total number of valid cases+controls of each pheno.
```
C3_BREAST_EXALLC	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-0/C3_BREAST_EXALLC.premunge.gz	110611
C3_BRONCHUS_LUNG_EXALLC	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-1/C3_BRONCHUS_LUNG_EXALLC.premunge.gz	180418
C3_PROSTATE_EXALLC	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-2/C3_PROSTATE_EXALLC.premunge.gz	83146
G6_PARKINSON	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-3/G6_PARKINSON.premunge.gz	224566
H7_AMD	gs://fg-cromwell_fresh/munge_fg/d17c3b71-2510-4d89-8bfb-3f788b50bd59/call-munge/shard-4/H7_AMD.premunge.gz	214660
```

### Munging

The wdl now contains internally a `premunge_ss` step where input sumstats are processed to match the LDSC notation, which is 

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


Therefore now the the inputs also require to pass the relevant column names for the munging. In case the data is not in rsid format, the script will automatically map chrom/pos --> rsid if needed. `chrom_col` and `pos_col` are required *only* if the `rsid_col` is missing

### Long description

Now all the data is ready to run the big run. A brief summary of the logic of the wdl.

The boolean flag `het_only` if set to `True` only produces heritabilities and does not compute genetic correlations. It's set to false by default, thus automatically calculating correlations. Please make sure it's your intention to do so.

As input one can have two lists of sumstats or just one. If two are provided then all the cross scores are computed between the two lists. If only one is passed instead the first list is duplicated as a second, thus running an inner product on itself. The total number of comparison that needs to run is `N*L/2` jobs to be run where `N` and `L` are the lengths of the two lists. In case of inner product, the number is `N(N-1)/2` instead. This means that the growht is quadratic and thus I recommed first testing the pipeline with a smaller set of sumstats.

`filter_meta` splits the input lists into chunks for munging. The number of chunks in this step is given by `  "ldsc_rg.filter_meta.filter_chunks": Int`. This step is quite fast anyways and, in principle, should only run once since its output is cached.

Each input chunk list is passed to `premunge_ss` that prepares the input sumstats for the ldsc pipeline as described above.

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
