workflow ldsc_rg {

    File meta_fg
    File meta_comparison
    File snplist

    String docker
    Boolean test


    Array[Array[String]] sumstats_fg = read_tsv(meta_fg)
    Array[Array[String]] sumstats_comparison = read_tsv(meta_comparison)

    Array[Array[String]] fg_meta  = if test then [sumstats_fg[0],sumstats_fg[1]]   else sumstats_fg
    Array[Array[String]] comparison_meta = if test then [sumstats_comparison[0],sumstats_comparison[1]] else sumstats_comparison

    scatter (a in comparison_meta) {
        call munge_ldsc as munge_comparison {
            input: docker=docker, pheno=a[0], sumstats=a[1], n=a[2],snplist = snplist
        }
    }
    scatter (fg in fg_meta) {
        call munge_ldsc {
            input: docker=docker, pheno=fg[0], sumstats=fg[1], n=fg[2],snplist = snplist
            }
      }

    call multi_rg {
      input: fg_files = munge_ldsc.out,comparison_files=munge_comparison.out,docker = docker

    }
}


task multi_rg {

  Array[File] fg_files
  Array[File] comparison_files

  String args
  Int cpus

  Int disk_size = ceil(size(fg_files[0],"MB"))*length(fg_files) + ceil(size(comparison_files[0],'MB'))*length(comparison_files)
  Int final_disk_size = ceil(disk_size/1000) + 20

  String docker
  String? multi_docker
  String? final_docker = if defined(multi_docker) then multi_docker else docker


  command <<<
  echo ${disk_size} ${final_disk_size}

  python3 /scripts/ldsc_mult.py \
  --ldsc-path "ldsc.py" \
  --list-1  ${write_lines(fg_files)} \
  --list-2  ${write_lines(comparison_files)} \
  -o /cromwell_root/results/ \
  --args ${args}

  >>>

  output {
      Array[File] out = glob("/cromwell_root/results/*")
  }

  runtime {
      docker: "${final_docker}"
      cpu: "${cpus}"
      memory: "${cpus} GB"
      disks: "local-disk ${final_disk_size} HDD"
      zones: "europe-west1-b"
      preemptible: 2
      noAddress: true
  }


}


task munge_ldsc {

    File sumstats
    String pheno
    String docker
    Int n
    File snplist

    Int file_size = ceil(size(sumstats,"GB")) + 1

    command <<<

        munge_sumstats.py \
        --sumstats ${sumstats} \
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
        disks: "local-disk ${file_size} HDD"
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
