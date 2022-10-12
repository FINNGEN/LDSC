version 1.0

workflow munge_raw {
    input {
        File sumstats_loc
        String docker
    }
    Array[String] sumstats_raw = read_lines(sumstats_loc)

    scatter (ss in sumstats_raw) {
        call munge {
            input: docker=docker, sumstats=ss
        }
    }
}


task munge {

  input {
    File sumstats
    File annotation
    String docker
  }

  command <<<

    join -t $'\t' -1 1 -2 1 \
    <(zcat ~{sumstats} | awk '
    BEGIN{OFS="\t"}
    NR==1{print "#variant",$0}
    NR>1 {print $1":"$2":"$3":"$4,$0}
    ' | sort -k1,1) \
    <(zcat ~{annotation} | sort -k1,1) | \
    awk '
    BEGIN{OFS="\t"}
    NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "SNP","A1","A2","BETA","P"}
    NR>1&&$h["rsid"]!="NA"&&$h["pval"]>0 {
      print $h["rsid"],$h["alt"],$h["ref"],$h["beta"],$h["pval"]
    }' | gzip > ~{basename(sumstats, ".gz")}.premunged.gz

    >>>

    output {
      File out = basename(sumstats, ".gz") + ".premunged.gz"
    }

    runtime {
      docker: "${docker}"
      cpu: 1
      memory: "3 GB"
      disks: "local-disk 200 HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      preemptible: 2
      noAddress: true
    }
}
