task munge_fg {

    String pheno
    File sumstats
    Int n
    String docker

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
        --merge-alleles /w_hm3.snplist

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

task munge_premunged {

    File sumstats
    String sumstats_base = basename(sumstats)
    String docker

    command <<<

        munge_sumstats.py \
        --sumstats ${sumstats} \
        --out ${sumstats_base}.ldsc \
        --merge-alleles /w_hm3.snplist

    >>>

    output {
        File out = sumstats_base + ".ldsc.sumstats.gz"
        File log = sumstats_base + ".ldsc.log"
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

task gather {

    Array[Array[File]] logs2D
    Array[File] logs = flatten(logs2D)
    String outfile
    String docker

    command <<<

        grep -A 1 -E "^Summary" ${logs[0]} | tail -1 | awk 'BEGIN{OFS="\t"} {$1=$1; print $0}' > ${outfile}
        for file in ${sep=" " logs}; do
            grep -A 2 -E "^Summary" $file | tail -1 | awk 'BEGIN{OFS="\t"} {gsub(".ldsc.sumstats.gz", ""); $1=$1; print $0}' >> ${outfile}
        done

    >>>

    output {
        File out = outfile
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

workflow ldsc_rg {

    File meta_fg
    File meta_premunged
    String docker

    Array[Array[String]] sumstats_fg = read_tsv(meta_fg)
    Array[String] sumstats_premunged = read_lines(meta_premunged)

    scatter (sumstats in sumstats_premunged) {
        call munge_premunged {
            input: docker=docker, sumstats=sumstats
        }
    }

    scatter (a in sumstats_fg) {
        call munge_fg {
            input: docker=docker, pheno=a[0], sumstats=a[1], n=a[2]
        }
        call rg {
            input: docker=docker, file1=munge_fg.out, files2=munge_premunged.out
        }
    }

    call gather {
        input: docker=docker, logs2D=rg.logs
    }
}
