workflow ldsc_rg {

    File meta_fg
    File meta_prs
    File snplist

    String docker
    Boolean test


    Array[Array[String]] sumstats_fg = read_tsv(meta_fg)
    Array[Array[String]] sumstats_prs = read_tsv(meta_prs)

    Array[Array[String]] fg_meta  = if test then [sumstats_fg[0],sumstats_fg[1]]   else sumstats_fg
    Array[Array[String]] prs_meta = if test then [sumstats_prs[0],sumstats_prs[1]] else sumstats_prs

    scatter (a in prs_meta) {
        call munge_prs {
            input: docker=docker, study=a[0], sumstats=a[1], n=a[2],snplist = snplist
        }
    }
    scatter (a in fg_meta) {
        call munge_fg {
            input: docker=docker, pheno=a[0], sumstats=a[1], n=a[2],snplist = snplist
            }
        call rg {
            input: docker=docker, file1=munge_fg.out, files2=munge_prs.out
            }
      }
}



task rg {

    File file1
    Array[File] files2
    String docker


    command <<<

        f1=`basename ${file1}`
        mv ${file1} $f1
        for file2 in ${sep=" " files2}; do
            f2=`basename $file2`
            mv $file2 $f2
            ldsc.py \
            --rg $f1,$f2 \
            --ref-ld-chr /eur_w_ld_chr/ \
            --w-ld-chr /eur_w_ld_chr/ \
            --out ldsc_$f1-$f2
        done

    >>>

    output {
        Array[File] logs = glob("ldsc_*.log")
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "6 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

task munge_prs {

    File sumstats
    String study
    String docker
    Int n
    File snplist

    command <<<

        munge_sumstats.py \
        --sumstats ${sumstats} \
        --N ${n} \
        --out ${study}.ldsc \
        --merge-alleles ${snplist}

    >>>

    output {
        File out = study + ".ldsc.sumstats.gz"
        File log = study + ".ldsc.log"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}
task munge_fg {

  String pheno
  File sumstats
  Int n
  String docker
  File snplist

  command <<<

    gunzip -c ${sumstats} | \
    awk 'BEGIN{FS=OFS="\t"} \
    NR==1{for(i=1;i<=NF;i++) a[$i]=i; print "SNP","A1","A2","BETA","P"} \
    NR>1 {if($a["rsids"]!="") print $a["rsids"],$a["alt"],$a["ref"],$a["beta"],$a["pval"]}' | \
    awk 'BEGIN{FS=OFS="\t"} {n=split($1,a,","); for(i=1;i<=n;i++) print a[i],$0}' | \
    cut -f1,3- | gzip > ${pheno}.premunge.gz

    munge_sumstats.py \
    --sumstats ${pheno}.premunge.gz \
    --N ${n} \
    --out ${pheno}.ldsc \
    --merge-alleles ${snplist}

    >>>

    output {
      File out = pheno + ".ldsc.sumstats.gz"
      File log = pheno + ".ldsc.log"
    }

    runtime {
      docker: "${docker}"
      cpu: 1
      memory: "4 GB"
      disks: "local-disk 200 HDD"
      zones: "europe-west1-b"
      preemptible: 2
      noAddress: true
      }
}
