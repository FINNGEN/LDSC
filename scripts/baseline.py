import argparse,os,logging,re,gzip,subprocess,shlex
from utils import tmp_bash,make_sure_path_exists,pretty_print,file_exists,log_levels,progressBar
import numpy as np


def filter_maf(bed_file,maf,bed_out,chroms,force=False):
    """
    Creates plink data set filtered for maf.
    """
    #check that all files exists
    paths = [bed_out +  f'.{c}.bed' for c in chroms]
    if not all([os.path.isfile(elem) for elem in paths])  or force:
        print('One or more files missing or force flag passed. Generating new plink data set.')
        force = True
        plink_cmd = f"plink2   --bfile {bed_file.replace('.bed','')} --maf {maf} --max-maf {1-float(maf)} --chr CHROM --make-bed --out  {bed_out}.CHROM "
        logging.debug(plink_cmd)
        for i,chrom in enumerate(chroms):
            progressBar(i+1,len(chroms))
            tmp_bash(plink_cmd.replace("CHROM",chrom))
            tmp_bash(f"plink2 --bfile {bed_out}.{chrom} --freq --out {bed_out}.{chrom}")

        print("\ndone.")
    else:
        print("All chrom files already generated.")

    return force
    
def build_annotation(bed_chrom_root,out_root,chroms,af_bins,force = False):

    """
    Builds annotation files for each chrom with binary annotation (low/high af).
    """
    annot_root = f"{out_root}.CHROM.annot.gz"
    paths = [os.path.isfile(annot_root.replace("CHROM",chrom)) for chrom in chroms]
    if not all(paths) or  force:       
        # builds bins 
        bin_names,bins = [],[]
        with open(af_bins) as i:
            for line in i:
                line = line.strip().split()
                bin_names.append(line[0])
                bins.append(line[1])
        logging.debug(bins)
    
        print(f"Generating at  {annot_root}")
        for i,chrom in enumerate(chroms):
            progressBar(i+1,len(chroms))
            # define inputs and outputs based on chrom
            freq_file = f"{bed_chrom_root}.{chrom}.afreq"
            bim_file = f"{bed_chrom_root}.{chrom}.bim"
            out_file = annot_root.replace("CHROM",chrom)
            # extract data (bp and freq) from .afreq file
            with open(freq_file) as i,open(bim_file) as b:
                next(i)
                variant_data,freqs = [],[]
                for line, bim_line in zip(i,b):
                    c,variant,_,pos,*_ = bim_line.strip().split() # get chrom,pos and variant id from bim file
                    _,freq_variant,_,_,freq,_ = line.strip().split() # get matching freq from
                    assert variant==freq_variant
                    variant_data.append('\t'.join([c,pos,variant,'0','1']))
                    freqs.append(min(float(freq),1-float(freq)))
            # bin data according to freq
            zeros = np.zeros((len(freqs),len(bin_names)),dtype =int)
            logging.info(zeros)
            freqs = np.array(freqs)
            logging.info(freqs.shape)
            b = np.digitize(freqs,np.array(bins,dtype=float)) -1
            logging.info(b.shape)
            zeros[np.arange(len(b)), b] = 1
            logging.info(zeros)

            # check that it's all ok and write to file
            assert len(variant_data) == len(zeros)
            with gzip.open(out_file,'wt') as o:
                # for each chromsome, calculate freq
                o.write( '\t'.join(["CHR","BP","SNP","CM","base"] + bin_names)+ '\n')
                for elem in zip(variant_data,zeros):o.write(elem[0] + '\t' +  '\t'.join(map(str,elem[1])) + '\n')
                
            print("\ndone.")
    else:
        print('All files already generated.')       

    return annot_root,force
        
def ld_scores(bed_chrom_root,annot_root,out_root,chroms,ldsc_path,force=False):

    score_root = f"{out_root}.CHROM.l2.ldscore.gz"
    paths = [os.path.isfile(score_root.replace("CHROM",chrom)) for chrom in chroms]
    if not  all(paths) or force:
        force = True
        
        # CREATE BASIC COMMAND REPLACING ALL CHROM NUMBERS WITH PLACEHOLDER
        ldsc_cmd = f"{ldsc_path} --l2 --bfile {bed_chrom_root}.CHROM  --ld-wind-kb 1000 --annot {annot_root} --out {out_root}.CHROM"
        logging.debug(ldsc_cmd.replace("CHROM","1"))
    
        for i,chrom in enumerate(chroms):
            progressBar(i+1,len(chroms))
            tmp_bash(ldsc_cmd.replace("CHROM",chrom))
        print("\ndone.")

    else:
        print('All files already generated.')
        
    return score_root,force


def ld_weights(bed_chrom_root,out_root,chroms,ldsc_path,force=False):

    weights_root = f"{out_root}.CHROM.l2.ldscore.gz"
    paths = [os.path.isfile(weights_root.replace("CHROM",chrom)) for chrom in chroms]
    if not all(paths) or force:        
        force = True
        # CREATE BASIC COMMAND REPLACING ALL CHROM NUMBERS WITH PLACEHOLDER
        weights_cmd = f"{ldsc_path} --l2 --bfile {bed_chrom_root}.CHROM  --ld-wind-kb 1000  --out {out_root}.CHROM"
        logging.debug(weights_cmd.replace("CHROM","1"))
        for i,chrom in enumerate(chroms):
            progressBar(i+1,len(chroms))
            tmp_bash(weights_cmd.replace("CHROM",chrom))
        print("\ndone.")

    else:
        print('All files already generated.')
    
    return weights_root,force


def ld_tau(bed_chrom_root,out_root,sumstats,trait,scores_root,weights_root,chroms,ldsc_path,force=False):

    tau_file = f"{out_root}.log"
    if not os.path.isfile(tau_file) or force:
        force = True
        # CREATE BASIC COMMAND REPLACING ALL CHROM NUMBERS WITH PLACEHOLDER
        tau_cmd = f"{ldsc_path} --h2 {sumstats} --ref-ld-chr {scores_root.split('CHROM')[0]} --w-ld-chr {weights_root.split('CHROM')[0]} --print-coefficients --print-delete-vals --not-M-5-50 --out {out_root} "
        logging.debug(tau_cmd.replace("CHROM","22"))
        tmp_bash(tau_cmd)
        print("\ndone.")

    else:
        print('Tau log file already generated')
    return tau_file,force


def create_file_rscript(annot_root,out_root,force=False):

    m_file = f"{out_root}.M"
    if not os.path.isfile(m_file) or force:
        force = True
        print(f"{m_file} missing or force.")
        cmd = f"Rscript {os.path.join(os.path.split(os.path.abspath(__file__))[0],'create_M_file.r')} {annot_root.split('.CHROM')[0]} {m_file}"
        logging.debug(cmd)
        tmp_bash(cmd)
    else:
        print(f"{m_file} already generated")
    return m_file,force

def compute_enrichments(tau_root,m_file,out_root):


    cmd = f"Rscript {os.path.join(os.path.split(os.path.abspath(__file__))[0],'compute_enrichments.r')} {tau_root} {m_file}"
    logging.debug(cmd)
    out_log = f"{out_root}.log"
    with  open(out_log, "w") as f:
        subprocess.call(shlex.split(cmd),stdout =f,stderr=f)


def main(args):

    #paths manipulation
    prefix = args.prefix + f".maf_{args.maf}"
    out_root = os.path.join(args.out,prefix)
    logging.debug(out_root)
    
    pretty_print("MAF FILTERING")
    plink_path = os.path.join(args.out,"plink")
    make_sure_path_exists(plink_path)
    bed_chrom_root = os.path.join(plink_path,prefix )
    force =filter_maf(args.bed,args.maf,bed_chrom_root,args.chroms,args.force)
    
    pretty_print("ANNOTATION")
    annot_path = os.path.join(args.out,"annot")
    make_sure_path_exists(annot_path)
    annot_chrom_root = os.path.join(annot_path,prefix )
    annot_root,force = build_annotation(bed_chrom_root,annot_chrom_root,args.chroms,args.maf_bins,force)

    pretty_print("LD SCORES")
    scores_path = os.path.join(args.out,"scores")
    make_sure_path_exists(scores_path)
    scores_chrom_root = os.path.join(scores_path,prefix)
    scores_root,force = ld_scores(bed_chrom_root,annot_root,scores_chrom_root,args.chroms,args.ldsc_path,force)

    pretty_print("WEIGHTS")
    weights_path = os.path.join(args.out,"weights")
    make_sure_path_exists(weights_path)
    weights_chrom_root = os.path.join(weights_path,prefix)
    weights_root,force = ld_weights(bed_chrom_root,weights_chrom_root,args.chroms,args.ldsc_path,force)

    pretty_print("TAU")
    tau_path = os.path.join(args.out,"tau")
    make_sure_path_exists(tau_path)
    tau_root = os.path.join(tau_path,f"{prefix}.{args.trait}")
    tau_file,force = ld_tau(bed_chrom_root,tau_root,args.sumstats,args.trait,scores_root,weights_root,args.chroms,args.ldsc_path,force)

    pretty_print("RSCRIPT")
    m_file,force = create_file_rscript(annot_root,out_root,force)
    compute_enrichments(tau_root,m_file,out_root)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Pipeline for LDSC baseline version.")
    parser.add_argument('--ldsc-path',help ='Path to ldsc path', required = False, default = '/home/pete/Tools/ldsc/ldsc.py',type = file_exists)
    parser.add_argument("-out",help ="Out path")
    parser.add_argument("--prefix",help ="Prefix output",default = 'finngen_test')

    parser.add_argument('--bed',help ='Path to plink dataset', required = True,type = file_exists)
    parser.add_argument('--maf',help ='Maf filter', required = False, default = '0.005')
    parser.add_argument('--rare-af',help ='Filter frequency to be considered rare or not', required = False, default = '0.05')


    parser.add_argument('--maf-bins',help ='Path to  tsv file with bins (name,start,end)', required = True,type = file_exists)
    parser.add_argument('--sumstats',help ='Path to sumstats', required = True,type = file_exists)

    parser.add_argument('--trait',help ='Trait name')
    parser.add_argument('--chroms',help ='Chroms to use.', nargs="+",required = False, default = list(map(str,range(1,23))))

    parser.add_argument("--args",help = "ldsc args",type = str,default ="")
    parser.add_argument('--force',action = 'store_true',help = 'Flag for forcing re-run.')
    parser.add_argument( "-log",  "--log",  default="warning", choices = log_levels, help=(  "Provide logging level. " "Example --log debug"))

    args = parser.parse_args()

    if  args.trait is  None: args.trait = os.path.basename(args.sumstats).split(".ldsc")[0]
    print(logging.getLevelName(10))
    level = log_levels[args.log]
    logging.basicConfig(level=level,format="%(levelname)s: %(message)s")

    logging.info(args)
    make_sure_path_exists(args.out)   
    main(args)
    
