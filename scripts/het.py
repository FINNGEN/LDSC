import os,argparse,re,json
from utils import make_sure_path_exists,tmp_bash


def h2(ldsc_path,args,sumstats,out_path,debug):

    pheno = os.path.basename(sumstats).split(".ldsc")[0]
    out_file = f"{os.path.join(out_path,pheno + '.ldsc.h2')}"
    if debug:print(sumstats,pheno)
    cmd = f"{ldsc_path} {args} --h2  {sumstats} --out {out_file} "
    if debug:print(cmd)
    tmp_bash(cmd,True)

    log_file = out_file+ ".log"
    h2 = get_het(log_file,pheno)
    with open(os.path.join(out_path, pheno +'.ldsc.h2.json'), 'w') as fp:json.dump(h2, fp)


def get_het(f,pheno):
    h2_dict = {pheno:["NA",'NA']}
    with open(f) as i:
        for line in i:
            if "Total Observed scale h2" in line:#fins line with h2 data
                h2 = re.findall("\d+\.\d+", line)#returns floats in line
            if "Intercept:" in line:
                int = re.findall("\d+\.\d+", line)
            if "Ratio" in line:
                ratio = re.findall("\d+\.\d+", line)
                ratio =  ["NA",'NA'] if len(ratio) != 2 else ratio



    h2_dict[pheno] = h2 + int + ratio

    return h2_dict

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Run het or return NA.")
    parser.add_argument('--ldsc-path',help ='Path to ldsc path', required = False, default = 'ldsc.py')
    parser.add_argument("--sumstats",help = "Sumstats file",required = True,type = str)
    parser.add_argument("--args",help = "ldsc args",required = True,type = str)
    parser.add_argument("-o",help ="Out path")
    parser.add_argument("--debug",action="store_true")
    args = parser.parse_args()
    make_sure_path_exists(args.o)
    if args.debug:print(args)
    h2(args.ldsc_path,args.args,args.sumstats,args.o,args.debug)
