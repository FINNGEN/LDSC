version 1.0

workflow ldsc_rg {

  input {
    File meta_fg
    File meta_other
    Boolean only_het
    String name
    String population
    String ld_root = "gs://finngen-production-library-green/ldsc/POP_ld.txt"
    String docker =  "eu.gcr.io/finngen-sandbox-v3-containers/ldsc:rsid_munge.2"
    File snplist = "gs://finngen-production-library-green/ldsc/w_hm3.snplist"
    Int filter_chunk_size = 30
    Int couples_chunk_size = 1000
  }

  File ld_list =sub(ld_root,"POP",population)
  String final_name = name + "_" + population

  # Split each input file into chunks; other-only entries are deduplicated against fg
  call filter_meta { input: meta_fg = meta_fg, meta_other = meta_other, docker = docker, chunk_size = filter_chunk_size }

  # Scatter over fg chunks
  scatter (chunk in filter_meta.chunk_fg) {
    call premunge_ss as premunge_fg {
      input: chunk = chunk, docker = docker
    }
    call munge_ldsc as munge_fg {
      input: docker = docker, chunk = chunk, ld_list = ld_list, premunged = premunge_fg.premunged, snplist = snplist
    }
  }

  # Only scatter over other chunks when the two input files differ
  if (meta_fg != meta_other) {
    scatter (chunk in filter_meta.chunk_other) {
      call premunge_ss as premunge_other {
        input: chunk = chunk, docker = docker
      }
      call munge_ldsc as munge_other {
        input: docker = docker, chunk = chunk, ld_list = ld_list, premunged = premunge_other.premunged, snplist = snplist
      }
    }
  }

  # Combine outputs from both scatter arms; other arm is empty when files are identical
  Array[File] all_het_jsons       = flatten([munge_fg.het_json, select_first([munge_other.het_json, []])])
  Array[File] all_het_logs        = flatten([munge_fg.het_log,  select_first([munge_other.het_log,  []])])
  # Wraps the two scattered arrays into one, then flattens twice
  Array[String] all_munged = flatten(flatten([munge_fg.munged, select_first([munge_other.munged, []])]))


  call gather_h2 {
    input: name = final_name, docker = docker, het_jsons = all_het_jsons, het_log = all_het_logs
  }

  if (!only_het) {
    call return_couples {
      input: file1 = meta_fg, file2 = meta_other, munged_sumstats = all_munged, docker = docker, chunk_size = couples_chunk_size

    }
    scatter (i in range(length(return_couples.couples))) {
      call multi_rg {
        input: couples = return_couples.couples[i], paths_list = return_couples.paths_list[i],
               jobs = return_couples.jobs, name = i, ld_list = ld_list, docker = docker
      }
    }
    call gather_summaries {
      input: docker = docker, name = final_name, summaries = multi_rg.summary, logs = multi_rg.log
    }
  }

  output {
    File herit_tsv    = gather_h2.herit_tsv
    File herit_log    = gather_h2.log
    Array[String] munged_ss = all_munged
    File corr_summary = select_first([gather_summaries.summary,gather_h2.herit_tsv])
    File corr_log     = select_first([gather_summaries.log,gather_h2.log])
  } 
}


task premunge_ss {
  input {
    File chunk
    String docker
    String beta_col
    String p_col
    String a1_effect_col
    String a2_ne_col
    String chrom_col = ""
    String pos_col = ""
    String rsid_col = ""
  }
  File rsid_map = "gs://finngen-production-library-green/rsids/convert/finngen.rsid.map.tsv.pickle"
  Array[String] phenos = transpose(read_tsv(chunk))[0]
  Array[File] sumstats = transpose(read_tsv(chunk))[1]
  Int mem = if rsid_col == "" then 16 else 4
  Int disk_size = 10 + 2*ceil(size(sumstats[0],'GB')) * length(sumstats) + ceil(size(rsid_map,'GB'))
  command <<<
  set -euo pipefail
  cat ~{write_lines(phenos)} > phenos.txt
  cat ~{write_lines(sumstats)} > sumstats.txt
  head sumstats.txt
  mapfile -t pheno_arr < phenos.txt
  mapfile -t ss_arr < sumstats.txt
  for i in "${!ss_arr[@]}"; do
      bash /scripts/ldsc_rsid_munge.sh --convert-script /rsid_map/scripts/convert_rsids.py \
           --input "${ss_arr[$i]}" \
           --outdir "$(pwd)" \
           --pheno "${pheno_arr[$i]}" \
           --rsid-map ~{rsid_map} \
           --beta-col '~{beta_col}' \
           --p-col '~{p_col}' \
           --a1-col '~{a1_effect_col}' \
           --a2-col '~{a2_ne_col}' \
           ~{if chrom_col != "" then "--chrom-col '" + chrom_col + "'" else ""} \
           ~{if pos_col != "" then "--pos-col '" + pos_col + "'" else ""} \
           ~{if rsid_col != "" then "--rsid '" + rsid_col + "'" else ""}
      rm -f *pre_convert*
  done

  ls -l *.premunge.gz
  >>>
  output {
    Array[File] premunged = glob("./*.premunge.gz")
  }
  runtime {
    docker: docker
    cpu: 1
    memory: "~{mem} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 2
  }
}


task munge_ldsc {

  input {
    File chunk
    String docker
    String args = ""
    File ld_list
    Array[File] premunged
    File snplist
  }

  Array[String] phenos = transpose(read_tsv(chunk))[0]
  Array[String] sumstats = transpose(read_tsv(chunk))[1]
  Array[String] ns = transpose(read_tsv(chunk))[2]
  Array[File] ld_files = read_lines(ld_list)
  Int disk_size = 30 + 2*ceil(size(premunged[0],'GB')) * length(premunged) + ceil(size(ld_files,'GB'))

  command <<<

  # build PHENO\tPATH\tN table to loop over. premunge_ss names files BASE.PHENO.premunge.gz,
  # where BASE is the input basename without extension; match on the full name since phenos can contain dots
  cat ~{write_lines(premunged)} > premunged.txt
  paste <(cat ~{write_lines(phenos)}) <(cat ~{write_lines(sumstats)}) <(cat ~{write_lines(ns)}) > chunk_meta.txt
  while IFS=$'\t' read pheno ss n; do
      base=$(basename "$ss")
      for ext in .txt.gz .tsv.gz .gz .txt .tsv; do base="${base%$ext}"; done
      f=$(awk -v t="${base}.${pheno}.premunge.gz" '{n = split($0, p, "/")} p[n] == t' premunged.txt)
      if [[ -z "$f" ]]; then echo "ERROR: no premunged file ${base}.${pheno}.premunge.gz for $pheno" >&2; exit 1; fi
      printf '%s\t%s\t%s\n' "$pheno" "$f" "$n"
  done < chunk_meta.txt > pheno_meta.txt
  nl --number-format=rn --number-width=2 pheno_meta.txt > meta.txt
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
      python3 /scripts/het.py  --ldsc-path "python3 /ldsc-2-to-3/ldsc.py"   --sumstats ${arr[1]}.ldsc.sumstats.gz --ld-path $ld_path ~{if args != "" then "--args " + args else ""} -o . 1> /dev/null ;
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
    File het_log = "./het.log"
  }
  runtime {
      docker: "${docker}"
      cpu: 1
      memory: "8 GB"
      disks: "local-disk ${disk_size} HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      preemptible: 2
  }
}


task return_couples {

  input {
    File file1
    File file2
    Int chunk_size
    String docker
    Array[String] munged_sumstats
  }

  Array[Array[String]] list1 = read_tsv(file1)
  Array[Array[String]] list2 = read_tsv(file2)

  command <<<

  # Initalize variables/inputs
  cat ~{write_lines(munged_sumstats)} > path_list.txt
  cat ~{write_tsv(list1)} > phenos1.txt
  cat ~{write_tsv(list2)} > phenos2.txt

  # 1. Generate unique, sorted pairs
  awk 'NR==FNR{a[$1];next} {for(i in a) print ($1 < i ? $1"\t"i : i"\t"$1)}' phenos1.txt phenos2.txt | sort -u > all_pairs.tmp

  # 2. Split into chunk files of chunk_size pairs each
  split -l ~{chunk_size} -d -a 5 all_pairs.tmp chunk_

  # 3. Builds list of required sumstats for each chunk
  # Files are PHENO.ldsc.sumstats.gz; key is the full basename before .ldsc.sumstats.gz (phenos can contain dots)
  python3 -c "import os, glob; d={os.path.basename(f.strip()).removesuffix('.ldsc.sumstats.gz'): f.strip() for f in open('path_list.txt')};[open(f.replace('chunk_', 'paths_'), 'w').write('\n'.join(set(d[p] for line in open(f) for p in line.strip().split('\t') if p in d))) for f in glob.glob('chunk_*')]"

  # 4. Cleanup and Metadata
  echo "~{chunk_size}" > jobs.txt
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


task filter_meta {

  input {
    File meta_fg
    File meta_other
    String docker
    Int chunk_size
  }
  command <<<
  set -euo pipefail
  # Drop exact duplicate rows, then fail if a pheno name still has more than one distinct row
  dedup_meta() {
      awk '!seen[$0]++' "$1" > "$2"
      dups=$(cut -f1 "$2" | sort | uniq -d)
      if [[ -n "$dups" ]]; then
          echo "ERROR: pheno names with more than one distinct row in $1:" >&2
          echo "$dups" >&2
          exit 1
      fi
  }
  dedup_meta ~{meta_fg} meta_fg.tsv
  dedup_meta ~{meta_other} meta_other.tsv
  # Chunk fg file
  split -l ~{chunk_size} -d --additional-suffix=.txt meta_fg.tsv chunk_fg
  # Chunk other file
  split -l ~{chunk_size} -d --additional-suffix=.txt meta_other.tsv chunk_other
  >>>
  output {
    Array[File] chunk_fg    = glob("./chunk_fg*")
    Array[File] chunk_other = glob("./chunk_other*")
  }
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

task gather_h2 {

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


task multi_rg {

  input {
    File couples
    File paths_list
    File ld_list
    String name
    Int cpus = 16
    Int jobs
    String docker
    String args = ""
  }
  Array[File] sumstats = read_lines(paths_list)
  Array[File] ld_files = read_lines(ld_list)

  Int final_cpus = if jobs > cpus then cpus else jobs
  Int mem = 8*cpus
  Int disk_size = 30 + ceil(size(sumstats[0],"MB")*length(sumstats)/1000)

  command <<<
  echo ~{disk_size} ~{final_cpus} ~{jobs}
  # get ld_path from first file in ld file list
  ld_path="$(dirname ~{ld_files[0]})/"
  cat ~{write_lines(sumstats)} > sumstats.txt &&  wc -l sumstats.txt
  cat ~{couples} > couples.txt && wc -l couples.txt

  python3 /scripts/ldsc_mult.py --ldsc-path "python3 /ldsc-2-to-3/ldsc.py "   --list sumstats.txt --couples couples.txt -o ./results/ --ld-path $ld_path  ~{if args != "" then "--args " + args else ""}

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
