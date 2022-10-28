version 1.0


workflow munge_fg {
  input {
    File meta_fg
    Boolean test
    String docker
  }
    Array[Array[String]] sumstats_fg = read_tsv(meta_fg)
    Array[Array[String]] fg_meta  = if test then [sumstats_fg[0],sumstats_fg[1]]   else sumstats_fg


    scatter (a in fg_meta) {
        call munge {
            input: docker=docker, pheno=a[0], sumstats=a[1]
            }
      }
}


task munge {
  input {
    String pheno
    File sumstats
    String docker
  }
  Int file_size = 2*ceil(size(sumstats,"GB")) + 2

  command <<<
  
  gunzip -c ~{sumstats} | \
      awk 'BEGIN{FS=OFS="\t"} \
    NR==1{for(i=1;i<=NF;i++) a[$i]=i; print "SNP","A1","A2","BETA","P"} \
    NR>1 {if($a["rsids"]!="") print $a["rsids"],$a["alt"],$a["ref"],$a["beta"],$a["pval"]}' | \
      awk 'BEGIN{FS=OFS="\t"} {n=split($1,a,","); for(i=1;i<=n;i++) print a[i],$0}' | \
      cut -f1,3- | gzip > ~{pheno}.premunge.gz
  
  >>>
  
  output {
    File out = pheno + ".premunge.gz"
  }
    
  runtime {
    docker: "${docker}"
    cpu: 1
    memory: "2 GB"
    disks: "local-disk ${file_size} HDD"
    zones: "europe-west1-b"
    preemptible: 1
    noAddress: true
  }
}
