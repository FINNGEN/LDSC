version 1.0

workflow munge_fg {

    input {
        File sumstats_loc
        String docker
    }

    Array[String] sumstats = read_lines(sumstats_loc)

    scatter (a in sumstats) {
        call munge {
            input:
                sumstats = a,
                docker = docker
        }
    }
}


task munge {

    input {
        File sumstats
        String rsid_col = "rsid"
        String a1_col = "ALT"
        String a2_col = "REF"
        Array[Array[String]] study_specific_cols
        String file_suffix = "_meta_out.tsv.gz"
        String docker
    }

    Int file_size = 2 * ceil(size(sumstats, "GB")) + 2
    String pheno = basename(sumstats, file_suffix)

    command <<<
    
        gunzip -c ~{sumstats} | awk '
        BEGIN {
            FS=OFS="\t"
        }
        FNR==NR {
            for(i=1; i<=NF; i++) a[FNR, i]=$i
            n_studies = FNR
            next
        }
        FNR==1 {
            for(i=1;i<=NF;i++) h[$i]=i
        }
        FNR>1 && $h["~{rsid_col}"] ~ /^rs/ {
            rsid = $h["~{rsid_col}"]
            a1 = $h["~{a1_col}"]
            a2 = $h["~{a2_col}"]
            for(i=1; i<=n_studies; i++) {
                if ($h[a[i, 3]] > 0 && $h[a[i, 3]] < 1 && $h[a[i, 2]] < 1e6 && $h[a[i, 2]] > -1e6) {
                    if (!seen[a[i, 1]]++) print "SNP", "A1", "A2", "BETA", "P" > "~{pheno}_"a[i, 1]".premunge.tsv"
                    print rsid, a1, a2, $h[a[i, 2]], $h[a[i, 3]] > "~{pheno}_"a[i, 1]".premunge.tsv"
                }
            }
        }' ~{write_tsv(study_specific_cols)} -

        bgzip ~{pheno}*.premunge.tsv
    
    >>>
    
    output {
        Array[File] out = glob("~{pheno}*.premunge.tsv.gz")
    }
        
    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ${file_size} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}
