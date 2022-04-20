from utils import tmp_bash,progressBar
import argparse,gzip,os,multiprocessing,time
from collections import Counter,defaultdict
import pandas as pd

def return_shared_phenos(endpoints,sumstats):
    """
    Gets shared phenos between sumstats list and endpoints
    """
    #with open(sumstats) as i:paths = set([os.path.splitext(os.path.basename(path.strip()))[0] for path in i])
    with open(sumstats) as i:paths = [elem.strip() for elem in i]
    print(f"{len(paths)} phenos present")

    with gzip.open(endpoints) as i:
        header = set([elem.decode("utf-8") for elem in next(i).strip().split()])
    print(f"{len(header)} columns present in endpoints")

    shared = list(header.intersection(paths))
    print(f"{len(shared)} phenos shared between the two files.")
    return sorted(shared)


def read_pheno(pheno,endpoints,test):
    nrows = 100 if test else None
    # read in data, casting NAs as 2s so we can reognize them
    df = pd.read_csv(endpoints,sep='\t',usecols = [pheno],nrows=nrows,dtype=pd.Int64Dtype()).fillna(2)
    
    count_dict = defaultdict(int)
    counter = Counter(df[pheno])
    count_dict.update(dict(counter))                      
    return '\t'.join(map(str,[pheno,count_dict[0],count_dict[1],count_dict[2]])) + '\n'

def wrapper_read(iter):
    return read_pheno(*iter)

def write_results(out_path,results):
    with open(os.path.join( out_path,"counts.txt"),'wt') as o:
        for entry in results:
            o.write(entry)


def main(args):
    shared_phenos = return_shared_phenos(args.endpoints,args.sumstats)
    params = [(pheno,args.endpoints,args.test) for pheno in shared_phenos]
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    results = pool.map_async(wrapper_read,params,chunksize=1)
    while not results.ready():
        time.sleep(2)
        progressBar(len(params) - results._number_left,len(params))

    progressBar(len(params) - results._number_left,len(params))
    print('\ndone.')
    results = results.get()
    pool.close()
    write_results(args.out_path,results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manipulation of sumstats for building tables required by LDSC.")
    parser.add_argument("--endpoints",type = str, help = "Endpoint file.",required = True)
    parser.add_argument("--sumstats",type = str, help = "File with sumstats",required = True)
    parser.add_argument('-o',"--out_path",type = str, help = "Folder in which to save the results", required = True)
    parser.add_argument("--test", action = 'store_true')


    args = parser.parse_args()
    main(args)
