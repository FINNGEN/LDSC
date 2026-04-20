version 1.0

workflow ldsc_rg {

  input {
    File meta_fg
    File? comparison_fg
    File meta_other = if defined(comparison_fg) then comparison_fg else meta_fg
    Boolean only_het
    
    String docker
    String name
    
    String population
    Map[String,File] ld_path

  }

  File ld_list = ld_path[population]
  String final_name = name + "_" + population
   
  # merge sumstats keeping unique values across them and creating munge chunks
  call filter_meta {input: meta_fg = meta_fg,meta_other = meta_other, docker = docker}

  #scatter over chunks and run heritability
  scatter (chunk in filter_meta.chunk_list)
  {
    call premunge_ss {input:chunk = chunk,docker = docker }
    call munge_ldsc{input:docker = docker, chunk = chunk,  ld_list = ld_list,premunged=premunge_ss.premunged}
  }
  
  # gather h2 data from munge step
  call gather_h2{input: name = final_name,docker = docker,het_jsons = munge_ldsc.het_json,het_log= munge_ldsc.het_log}

  if (!only_het)
  {
    # returns all valid couples of phenocodes that need to be run in terms of absolute paths from the output of the previous task
    call return_couples { input : file1 = meta_fg,file2 = meta_other, munged_chunks = munge_ldsc.munged, docker = docker}
    scatter (i in range(length(return_couples.couples))) {
      # where the actual work is done
      call multi_rg {input: couples=return_couples.couples[i], paths_list = return_couples.paths_list[i],jobs = return_couples.jobs, name = i,ld_list =ld_list,docker=docker}
    }
    # gather all the outputs of multi_rg in order to create a final table
    call gather_summaries {  input:   docker = docker,   name = final_name,   summaries = multi_rg.summary, logs = multi_rg.log   }
  }
}


task premunge_ss {
  input {
    File chunk              # e.g. chunk file, a tsv with header "pheno, path, N, ..."
    String docker
    String beta_col 
    String p_col
    String a1_effect_col
    String a2_ne_col
    String? chrom_col
    String? pos_col
    String? rsid_col
  }
  File rsid_map = "gs://finngen-production-library-green/rsids/convert/finngen.rsid.map.tsv.pickle"
  Array[File] sumstats = transpose(read_tsv(chunk))[1]
  Int mem = if defined(rsid_col) then 4 else 8
  Int disk_size = 10 + 2*ceil(size(sumstats[0],'GB')) * length(sumstats) + ceil(size(rsid_map,'GB'))
  command <<<
  set -euo pipefail
  cat ~{write_lines(sumstats)} > sumstats.txt
  head sumstats.txt
  while read -r f; do
      bash /scripts/ldsc_rsid_munge.sh --convert-script /rsid_map/scripts/convert_rsids.py \
           --input "$f" \
           --outdir "$(pwd)" \
           --rsid-map ~{rsid_map} \
           --beta-col '~{beta_col}' \
           --p-col '~{p_col}' \
           --a1-col '~{a1_effect_col}' \
           --a2-col '~{a2_ne_col}' \
           ~{if defined(chrom_col) then "--chrom-col '" + chrom_col + "'" else ""}  ~{if defined(pos_col) then "--pos-col '" + pos_col + "'" else ""}  ~{if defined(rsid_col) then "--rsid '" + rsid_col +"'" else ""} 
  done < sumstats.txt
  ls -l *.premunge.gz
  >>>
  output {
    Array[File] premunged = glob("*.premunge.gz")
  }
  runtime {
    docker: docker
    cpu: 1
    memory: "~{mem} GB"
    disks: "local-disk ~{disk_size} HDD"
  }
}

task multi_rg {

  input {
    File couples
    File paths_list
    File ld_list
    String name
    Int cpus
    Int jobs
    String docker
    String? args
  }
  Array[File] sumstats = read_lines(paths_list)
  Array[File] ld_files = read_lines(ld_list)
  
  Int final_cpus = if jobs > cpus then cpus else jobs
  Int mem = 2*cpus
  Int disk_size = 30 + ceil(size(sumstats[0],"MB")*length(sumstats)/1000)

  
  command <<<
  echo ~{disk_size} ~{final_cpus} ~{jobs}
  # get ld_path from first file in ld file list
  ld_path="$(dirname ~{ld_files[0]})/"
  cat ~{write_lines(sumstats)} > sumstats.txt &&  wc -l sumstats.txt
  cat ~{couples} > couples.txt && wc -l couples.txt

  python3 /scripts/ldsc_mult.py --ldsc-path "python3 /ldsc-2-to-3/ldsc.py "   --list sumstats.txt --couples couples.txt -o ./results/ --ld-path $ld_path  ~{if defined(args) then "--args " + args else ""}
    
  echo -e "\nDONE"
  # write to file list of ldsc log files
  for f in ./results/*log; do echo $f &&  echo $f >> summaries.txt ; done
  # copy content of each log into single log file
  while read f; do cat $f >> ~{name}.log ; done < summaries.txt
  # extract metadata from each log file
  python3 /scripts/extract_metadata.py  --summaries summaries.txt  --name ~{name}

  >>>

  output {
      File log = "${name}.log"
      File summary = "${name}.ldsc.summary.log"
  }

  runtime {
      docker: "${docker}"
      cpu: "${final_cpus}"
      memory: "${mem} GB"
      disks: "local-disk ${disk_size} HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      preemptible: 1
      noAddress: true
  }
}

task return_couples {

  input {
    File file1
    File file2
    Int chunks
    String docker
    Array[Array[String]] munged_chunks
  }
  
  Array[Array[String]] list1 = read_tsv(file1)
  Array[Array[String]] list2 = read_tsv(file2)
  Array[String] munged_sumstats = flatten(munged_chunks)

  command <<<

  # Initalize variables/inputs
  cat ~{write_lines(munged_sumstats)} > path_list.txt
  cat ~{write_tsv(list1)} > phenos1.txt
  cat ~{write_tsv(list2)} > phenos2.txt
  num_chunks=~{chunks}
  
  # 1. Generate unique, sorted pairs
  awk 'NR==FNR{a[$1];next} {for(i in a) print ($1 < i ? $1"\t"i : i"\t"$1)}' phenos1.txt phenos2.txt | sort -u > all_pairs.tmp

  # 2. Calculate chunk size
  total_pairs=$(wc -l < all_pairs.tmp)
  n=$(( (total_pairs + num_chunks - 1) / num_chunks ))

  # 3. Split into chunk files
  split -l "$n" -d -a 2 all_pairs.tmp chunk_

  # 4. Builds list of required sumstats for each chunk
  python3 -c "import os, sys; d={os.path.basename(f).replace('.premunge.gz',''): (f, n) for f, n in zip(open('fnames.txt').read().splitlines(), open('ns.txt').read().splitlines())}; [print(f'{p}\t{d[p][0]}\t{d[p][1]}') if p in d else sys.exit(f'Missing: {p}') for p in open('phenos.txt').read().splitlines()]" | nl --number-format=rn --number-width=2 > meta.txt

  # 5. Cleanup and Metadata
  echo "$n" > jobs.txt
  rm all_pairs.tmp
  >>>

  output {
      Array[File] couples = glob("./chunk*")
      Array[File] paths_list = glob("./paths*")
      Int jobs = read_int("jobs.txt")
  }

  runtime {
      docker: "${docker}"
      cpu: 2
      memory: "4 GB"
      disks: "local-disk 10 HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      preemptible: 2
      noAddress: true
  }
}

task munge_ldsc{

  input {
    File chunk
    String docker
    File snplist
    String? args    
    File ld_list
    Array[File] premunged
  }
  
  Array[Array[String]] by_type = transpose(read_tsv(chunk))
  Array[String] phenos = by_type[0]
  Array[File] fnames = premunged
  Array[String] ns = by_type[2]
  Array[File] ld_files = read_lines(ld_list)
  Int disk_size = 30 + 2*ceil(size(fnames[0],'GB')) * length(fnames) + ceil(size(ld_files,'GB'))
  
  command <<<

  df -h .
  cat ~{write_lines(phenos)} > phenos.txt &&   cat ~{write_lines(fnames)} > fnames.txt &&  cat ~{write_lines(ns)} > ns.txt
  # update the table with the premunged paths
  python3 -c "import sys,os; phenos=[l.strip() for l in open('phenos.txt')]; fnames=[l.strip() for l in open('fnames.txt')]; ns=[l.strip() for l in open('ns.txt')]; [(os.path.basename(f).replace('.premunge.gz','')==p or sys.exit(f'mismatch: {p} {f}')) for p,f in zip(phenos,fnames)]; [print(f'{p}\t{f}\t{n}') for p,f,n in zip(phenos,fnames,ns)]" | nl --number-format=rn --number-width=2 > meta.txt
  head meta.txt

  # get ld_path from first file in ld file list
  ld_path="$(dirname ~{ld_files[0]})/"
  echo $ld_path
  # read through file and munge it + calculate heritability
  while read line ; do
      arr=($line)
      echo -ne "\r${arr[0]}/~{length(phenos)} ${arr[1]}                  "
      zcat ${arr[2]} > ${arr[1]}.tmp.txt
      python3 /ldsc-2-to-3/munge_sumstats.py  --sumstats ${arr[1]}.tmp.txt    --N ${arr[3]} --out ${arr[1]}.ldsc  --merge-alleles ~{snplist} 1> /dev/null
      python3 /scripts/het.py  --ldsc-path "python3 /ldsc-2-to-3/ldsc.py"   --sumstats ${arr[1]}.ldsc.sumstats.gz --ld-path $ld_path ~{if defined(args) then "--args " + args else ""} -o . 1> /dev/null ;
  done < meta.txt
  # merge log files
  cat *ldsc.log >> munge.log &&  cat *ldsc.h2.log >> het.log
  
  # merge jsons into one
  jq -s 'reduce .[] as $item ({}; . * $item)' ./*json > het.json 

  >>>
  output {
    Array[File] munged = glob("./*sumstats.gz")
    File munge_log = "./munge.log"
    File het_json = "./het.json"
    File het_log = " ./het.log"
  }
  runtime {
      docker: "${docker}"
      cpu: 1
      memory: "4 GB"
      disks: "local-disk ${disk_size} HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      preemptible: 2
  }
}

task filter_meta {

  input {
    File meta_fg
    File meta_other
    String docker
    Int filter_chunks
  }
  command <<<
  cat ~{meta_fg} > tmp.txt
  cat ~{meta_other} >> tmp.txt
  sort tmp.txt | uniq >> meta.txt
  split -en r/~{filter_chunks} -d --additional-suffix=.txt meta.txt chunk
  >>>

  output {Array[File] chunk_list =  glob("./chunk*")}
  runtime {
      docker: "${docker}"
      cpu: 1
      memory: "4 GB"
      disks: "local-disk 10 HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      preemptible: 2
      noAddress: true
  }
}

task gather_summaries {

  input {
    Array[File] summaries
    Array[File] logs
    String name
    String docker
  }

  Int disk_size = ceil(size(summaries[0],'MB')*length(summaries)) + ceil(size(logs[0],'MB')*length(logs))
  Int final_disk_size = ceil(disk_size/1000) + 50

  command <<<

  df -h .
  while read f; do cat $f >> ~{name}.ldsc.logs.txt ; done < ~{write_lines(logs)}

  head -n1 ~{summaries[0]} > ~{name}.ldsc.summary.tsv
  while read f; do cat $f | sed -E 1d >> ~{name}.ldsc.summary.tsv ; done < ~{write_lines(summaries)}
  >>>

  output {
    File log = "${name}.ldsc.logs.txt"
    File summary = "${name}.ldsc.summary.tsv"
  }

  runtime {
      docker: "${docker}"
      cpu: 2
      memory: "4 GB"
      disks: "local-disk ${final_disk_size} HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      preemptible: 2
      noAddress: true
  }
}

task gather_h2{

  input {
    Array[File] het_jsons
    Array[File] het_log
    String name
    String docker
  }
  

  command <<<
  cat ~{write_lines(het_jsons)} >> h2.txt
  python3 /scripts/extract_metadata.py \
  --het h2.txt \
  --name ~{name}
  while read f; do cat $f >> ~{name}.ldsc.heritability.log; done <  ~{write_lines(het_log)}
  python3 /scripts/plot_summary.py  --het ~{name}.ldsc.heritability.tsv  --columns INT RATIO
  >>>

  output {
    File herit = "${name}.ldsc.heritability.json"
    File herit_tsv = "${name}.ldsc.heritability.tsv"
    File fig = "${name}.ldsc.heritability.pdf"
    File log = "${name}.ldsc.heritability.log"
  }

  runtime {
    docker: "${docker}"
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 20 HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    noAddress: true
  }
}
