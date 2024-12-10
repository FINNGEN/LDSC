import gzip,argparse,os,multiprocessing,itertools,time
from utils import make_sure_path_exists,progressBar,tmp_bash
from functools import partial

def get_unique_pairs(file_list,out_path):


    # read lists and keep only one instance of each file
    with open(file_list) as i: files = [elem.strip() for elem in i.readlines()]


    # get all couples possible as the product of both lists
    couples = itertools.product(set(files),set(files))
    # sort all couples into tuples --> create set in order to remove duplicate entries
    sorted_couples = sorted(list(set([tuple(sorted(list(couple))) for couple in couples])))

    with open(os.path.join(out_path,'couples.txt'),'wt') as o:
        for couple in sorted_couples:
            phenos = [os.path.basename(elem).split(".ldsc.sumstats")[0] for elem in  couple]
            o.write('\t'.join(phenos) + '\n')
    return sorted_couples


def return_file_couples(couples,file_list):

    """
    From the pheno couples, create file lists
    """
    with open(file_list) as i: files = [elem.strip() for elem in i.readlines()]


    # create mapping of basenames to sumstats
    pheno_file_dict ={}
    for f in files:
        base = os.path.basename(f).split('.ldsc.sumstats.gz')[0]
        pheno_file_dict[base] = f

    # recreate couples based on file names
    file_couples = []
    with open(couples) as i:
        for line in i:
            pheno1,pheno2 = line.strip().split()
            if all([os.path.isfile(elem) for elem in [pheno1,pheno2]]):
                file_couples.append([pheno1,pheno2])
            else:
                file_couples.append([pheno_file_dict[pheno1],pheno_file_dict[pheno2]])

    return file_couples



def multiproc(couples,ld_path,args,cpus,out_path,ldsc_path,debug):

    cmd = f"{ldsc_path} --ref-ld-chr {ld_path} --w-ld-chr {ld_path} {args} --rg FILE1,FILE2 --out {os.path.join(out_path,'ldsc_')}"
    print(cmd)

    params = [[elem[0],elem[1],cmd,debug] for elem in couples]

    print(f"{len(params)} sumstats couples to loop over")

    
    pool = multiprocessing.Pool(cpus)
    results = pool.map_async(multi_wrapper,params,chunksize=1)
    while not results.ready():
        time.sleep(2)
        progressBar(len(params) - results._number_left,len(params))

    progressBar(len(params) - results._number_left,len(params))
    results = results.get()
    pool.close()


def multi_wrapper(arguments):
    multi_func(*arguments)
    return

def multi_func(file1,file2,cmd,check=False):
    cmd = cmd.replace("FILE1",file1).replace("FILE2",file2) + f"{os.path.basename(file1).split('.ldsc')[0]}_{os.path.basename(file2).split('.ldsc')[0]}"
    tmp_bash(cmd,check=check)


def main(args):

    if not args.couples:
        couples = get_unique_pairs(args.list,args.o)

    else:
        couples = return_file_couples(args.couples,args.list)

    #for elem in couples:print('\t'.join([os.path.basename(f) for f in elem]))
    multiproc(couples,args.ld_path,args.args,args.cpus,args.o,args.ldsc_path,args.debug)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Parallelize ldsc generation.")
    parser.add_argument('--list',help ='Input files 1',required = True)
    parser.add_argument('--couples',help ='List of couples to be analyzed using basenames of files')

    parser.add_argument('--ldsc-path',help ='Path to ldsc path', required = False, default = 'ldsc.py')
    parser.add_argument('--ld-path',help ='Path to ld ref scores', required = True)
    parser.add_argument("--debug",action="store_true")
    parser.add_argument("-o",help ="Out path")
    parser.add_argument("--args",help = "ldsc args",type = str,default ="")
    parser.add_argument("--cpus",type = int, help = "number of cpus to use", default =  multiprocessing.cpu_count())

    args = parser.parse_args()
    print(args)
    make_sure_path_exists(args.o)

    main(args)
