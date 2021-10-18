from utils import tmp_bash,progressBar
import argparse,gzip,os
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


def filter_dataframe(endpoints,shared_phenos,out_path):

    out_path = out_path + 'shared.csv'
    if not os.path.isfile(out_path):
        df = pd.read_csv(endpoints,sep='\t',usecols = shared_phenos)
        df.to_csv(out_path,sep='\t',index = False)
        print("data imported.")
    else:
        print(f"{out_path} already exists")
    return out_path

def count(endpoints,shared_phenos,test,out_path):
    """
    Function that creates tables of cases/controls/nas for phenotypes
    """
    nrows = 1000 if test else None
    cols = shared_phenos[:10] if test else shared_phenos
    df = pd.read_csv(endpoints,nrows = nrows,sep='\t',usecols = cols)
    print("data imported.")
    samples = int(len(df))
    with open(out_path + "counts.txt",'wt') as o:
        for i,pheno in enumerate(cols):
            progressBar(i,len(cols))
            p_data = df[pheno]
            nas = int(samples - p_data.count())
            cases = int(p_data.sum())
            controls = int(samples - nas - cases)
            o.write('\t'.join(map(str,[pheno,cases+controls,cases,controls,nas])) + '\n')
    progressBar(1,1)
    print('\ndone.')

def main(args):

    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Manipulation of sumstats for building tables required by LDSC.")
    parser.add_argument("--endpoints",type = str, help = "Endpoint file.",required = True)
    parser.add_argument("--sumstats",type = str, help = "File with sumstats",required = True)
    parser.add_argument('-o',"--out_path",type = str, help = "Folder in which to save the results", required = True)
    parser.add_argument("--test", action = 'store_true')


    args = parser.parse_args()
    shared_phenos = return_shared_phenos(args.endpoints,args.sumstats)
    new_endpoints = filter_dataframe(args.endpoints,shared_phenos,args.out_path)
    count(new_endpoints,shared_phenos,args.test,args.out_path)
