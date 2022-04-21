version 1.0

workflow ldsc_rg {


  input {
    File meta_fg
    File? comparison_fg
    File meta_other = if defined(comparison_fg) then comparison_fg else meta_fg

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
  scatter (chunk in filter_meta.chunk_list){
    call munge_ldsc{input:docker = docker, chunk = chunk,ld_list = ld_list}
  }
  # returns all valid couples of phenocodes that need to be run in terms of absolute paths from the output of the previous task
  call return_couples { input : file1 = meta_fg,file2 = meta_other, munged_chunks = munge_ldsc.munged, docker = docker}

  scatter (i in range(length(return_couples.couples))) {
    # where the actual work is done
    call multi_rg {input: couples=return_couples.couples[i], paths_list = return_couples.paths_list[i],jobs = return_couples.jobs, name = i, docker = docker,ld_list =ld_list}
  }


  # gather all the outputs of multi_rg in order to create a final table
  call gather_summaries {  input:   docker = docker,   name = final_name,   summaries = multi_rg.summary, logs = multi_rg.log   }
  # gather h2 data from munge step
  call gather_h2{input: name = final_name,docker = docker,het_jsons = munge_ldsc.het_json,het_log= munge_ldsc.het_log}

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
  Int disk_size = 20 + ceil(size(sumstats[0],"MB")*length(sumstats)/1000)

  

  command <<<
    echo ~{disk_size} ~{final_cpus} ~{jobs}
    # get ld_path from first file in ld file list
    ld_path="$(dirname ~{ld_files[0]})/"
    cat ~{write_lines(sumstats)} > sumstats.txt &&  wc -l sumstats.txt
    cat ~{couples} > couples.txt

    python3 /scripts/ldsc_mult.py \
    --ldsc-path "ldsc.py" \
    --list sumstats.txt   \
    --couples couples.txt \
    -o /cromwell_root/results/ \
    --ld-path $ld_path \
    ~{if defined(args) then "--args " + args else ""}
    

   for f in /cromwell_root/results/*log; do echo $f >> summaries.txt ; done
   while read f; do cat $f >> ~{name}.log ; done < summaries.txt
   cat summaries.txt

   python3 /scripts/extract_metadata.py \
   --summaries summaries.txt \
   --name ~{name}

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
      zones: "europe-west1-b"
      preemptible: 2
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

  # write to file all absolute paths of sumstats
  cat ~{write_lines(munged_sumstats)} > path_list.txt

  python3 <<CODE

  import itertools,os

  # open input pheno lists
  with open('~{write_tsv(list1)}') as i: first_list = [elem.strip().split()[0] for elem in i.readlines()]
  with open('~{write_tsv(list2)}') as i: second_list = [elem.strip().split()[0] for elem in i.readlines()]

  #create all possible combinations
  couples = itertools.product(list(set(first_list)),list(set(second_list)))

  #sort them based on first element and return unique elements
  sorted_couples = sorted([list(elem) for elem in set([tuple(sorted(list(couple))) for couple in couples])])
  print(len(first_list))
  print(len(second_list))
  print(len(sorted_couples))
  # return subchunks
  n = 1 + len(sorted_couples)//~{chunks}
  sublists = [sorted_couples[i:i + n] for i in range(0, len(sorted_couples), n)]

  # create pheno to file mapping
  pheno_file_dict ={}
  with open("./path_list.txt",'rt') as i:
    for line in i:
      f = line.strip()
      pheno_file_dict[os.path.basename(f).split('.ldsc.sumstats.gz')[0]] = f

  # convert pheno couples to file couples and create chunks with absolute paths to munged sumstats
  for i,sub in enumerate(sublists):
      with open(f"./chunk_{i}",'wt') as o,open(f"./paths_{i}",'wt') as p:
        paths = []
        for c in sub:
          o.write('\t'.join(c) + '\n')
          for pheno in c:paths.append(pheno_file_dict[pheno])
        p.write('\n'.join(set(paths)))


  with open('jobs.txt','wt') as o:
      o.write(str(n))

  CODE
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
      zones: "europe-west1-b"
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

    }
    
  Array[Array[String]] by_type = transpose(read_tsv(chunk))
  Array[String] phenos = by_type[0]
  Array[File] fnames = by_type[1]
  Array[String] ns = by_type[2]

  Int disk_size = 2 + ceil(size(fnames[0],'GB')) * length(fnames)
  
  Array[File] ld_files = read_lines(ld_list)
  
  command <<<

    df -h .
    # rebuild the original tsv with metadata but now localized filepaths and add line number
    paste -d '\t' ~{write_lines(phenos)} ~{write_lines(fnames)} ~{write_lines(ns)} | nl --number-format=rn --number-width=2  > meta.txt
    wc -l meta.txt


    # get ld_path from first file in ld file list
    ld_path="$(dirname ~{ld_files[0]})/"
    echo $ld_path
    
    # read through file and munge it + calculate heritability
    while read line ; do arr=($line) && echo -ne "\r${arr[0]}/~{length(phenos)} ${arr[1]}                  "
    munge_sumstats.py  --sumstats ${arr[2]}  \
    --N ${arr[3]} --out ${arr[1]}.ldsc  --merge-alleles ~{snplist} 1> /dev/null && \
    python3 /scripts/het.py  --ldsc-path "ldsc.py" \
    --sumstats ${arr[1]}.ldsc.sumstats.gz --ld-path $ld_path ~{if defined(args) then "--args " + args else ""} -o . 1> /dev/null ;
    done < meta.txt

    ls 
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
      zones: "europe-west1-b"
      preemptible: 2
      noAddress: true
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
  split -edn r/~{filter_chunks} --additional-suffix=.txt meta.txt chunk

  >>>

  output {Array[File] chunk_list =  glob("./chunk*")}
  runtime {
      docker: "${docker}"
      cpu: 1
      memory: "4 GB"
      disks: "local-disk 10 HDD"
      zones: "europe-west1-b"
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
      zones: "europe-west1-b"
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
      zones: "europe-west1-b"
      preemptible: 2
      noAddress: true
  }
}
