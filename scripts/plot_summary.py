import numpy as np
import pandas as pd
import seaborn as sns
import pylab,argparse,os
from matplotlib import pyplot as plt
import matplotlib as mpl
from pathlib import Path
from utils import make_sure_path_exists


def plot_data(het,columns,fig_path):
    pylab.ioff()

    # SET UP FIG
    df = pd.read_csv(het,sep='\t',index_col =[0],usecols =["PHENO"] + columns).dropna()
    print(df.head())

    g = sns.JointGrid(data=df,x = columns[0],y=columns[1])
    g.plot_joint(sns.scatterplot,s=20, alpha=.5)
    g.plot_marginals(sns.histplot, kde=True)

    g.savefig(fig_path)
    print('done')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ="Plot het.")
    parser.add_argument('--het',help ='het tsv file', required = True)
    parser.add_argument('--columns',nargs = 2,help ='what to plot', default =["INT","RATIO"])
    parser.add_argument("-o",help ="Out path",default = ".")
    args = parser.parse_args()
    make_sure_path_exists(args.o)

    fig_path = os.path.join(args.o,Path(args.het).stem + ".pdf")
    print(fig_path)
    plot_data(args.het,args.columns,fig_path)
