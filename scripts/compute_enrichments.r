param <- commandArgs(trailingOnly=T)

file  = eval(paste(text=param[1]))
Mfile = eval(paste(text=param[2]))

if (length(param)!=2) {
  cat("\nError: two arguments must be supplied.\n")
  cat("ARG1 : the file header from the .log and\n")
  cat("       the .part_delete files\n")
  cat("ARG2 : the .M file.\n\n")
  stop("two arguments must be supplied.\n", call.=FALSE)
}


se_jacknife <- function(theta,theta_j){
	theta_hatJ=sum(theta-theta_j)+sum(theta_j/200)
	tau=200*theta-199*theta_j
	sqrt((1/200)*sum((tau-theta_hatJ)**2/(199)))
}

# read file.log, file.part_delete files
baselineLF.log = read.table(paste(file,".log",sep=""),h=F,fill=T);
rg1            = which(baselineLF.log[,1]=="Coefficients:");
rg2            = which(baselineLF.log[,1]=="Coefficient")-1;
tau            = NULL
for(i in c(rg1:rg2)){
	for (j in 1:ncol(baselineLF.log)){
        if (baselineLF.log[i,j]!="" & baselineLF.log[i,j]!="Coefficients:") {
			tau = c(tau,as.numeric(as.character(baselineLF.log[i,j])))
        }
	}
}
baselineLF.JN  = read.table(paste(file,".part_delete",sep=""),h=F); 

# read Mfile file, and find binnary annotations
sumannots   = read.table(Mfile,h=T)
list_lowfrq = gsub(".lowfrq","",sumannots[grep("\\.lowfrq",sumannots[,1]),1])
list_common = gsub(".common","",sumannots[grep("\\.common",sumannots[,1]),1])
categories  = intersect(list_lowfrq,list_common) 
row.names(sumannots)=paste(sumannots[,1],"",sep="")
sumannots=sumannots[,-1]
#check number of annotations
if (length(tau)!=ncol(sumannots)){
    stop(paste("\nNumber of annotations in the .log file = ",length(tau),"\nNumber of annotations in the .M file   = ",ncol(sumannots),"\n",sep=""))
}
cat("Number of binary annotations in the .M file :",length(categories),"\n\n")

#enr lowfrq vs common
cat.lf          = paste(categories,".lowfrq",sep="")
cat.cm          = paste(categories,".common",sep="")
Vcm             = apply(sumannots[which(rownames(sumannots)=="Common_MAF_bin_1"):which(rownames(sumannots)=="Common_MAF_bin_10"),],2,sum)
Vlf             = apply(sumannots[which(rownames(sumannots)=="Low.frequency_MAF_bin_1"):which(rownames(sumannots)=="Low.frequency_MAF_bin_5"),],2,sum)
Ncm             = sum(Vcm[which(colnames(sumannots)=="Common_MAF_bin_1"):which(colnames(sumannots)=="Common_MAF_bin_10")])
Nlf             = sum(Vlf[which(colnames(sumannots)=="Low.frequency_MAF_bin_1"):which(colnames(sumannots)=="Low.frequency_MAF_bin_5" )])
Ncm.annot       = apply(sumannots[cat.cm,which(colnames(sumannots)=="Common_MAF_bin_1"):which(colnames(sumannots)=="Common_MAF_bin_10")],1,sum)
Nlf.annot       = apply(sumannots[cat.lf,which(colnames(sumannots)=="Low.frequency_MAF_bin_1"):which(colnames(sumannots)=="Low.frequency_MAF_bin_5") ],1,sum)
N               = Ncm+Nlf

#compute h2c
baselineLF.h2cm     = sum(Vcm*tau)
#compute h2lf
baselineLF.h2lf     = sum(Vlf*tau)

#compute LFVE
baselineLF.enr.lf   = as.numeric(((apply(tau*t(sumannots[as.character(cat.lf),]),2,sum))/(baselineLF.h2lf))/(Nlf.annot/Nlf))
#compute CVE
baselineLF.enr.cm   = as.numeric(((apply(tau*t(sumannots[as.character(cat.cm),]),2,sum))/(baselineLF.h2cm))/(Ncm.annot/Ncm))
#compute LFVE/CVE ratio
baselineLF.ratio    = baselineLF.enr.lf / baselineLF.enr.cm

#compute LFVE statistic
h2cat=apply(tau*t(sumannots[as.character(cat.lf),]),2,sum); C=Nlf.annot; h2=baselineLF.h2lf; M=Nlf
baselineLF.stat.lf  = (h2cat/C)-((h2-h2cat)/(M-C))
#compute CVE statistic
h2cat=apply(tau*t(sumannots[as.character(cat.cm),]),2,sum); C=Ncm.annot; h2=baselineLF.h2cm; M=Ncm
baselineLF.stat.co  = (h2cat/C)-((h2-h2cat)/(M-C))

#compute Jacknife estimates of h2lf, h2c, LFVE, CVE, LFVE/CVE, LFVE stat and CVE stat
baselineLF.h2lfJN    = NULL 
baselineLF.h2cmJN    = NULL
baselineLF.enrJN.lf  = matrix(0,200,length(cat.lf))
baselineLF.enrJN.co  = matrix(0,200,length(cat.cm))
baselineLF.ratioJN   = matrix(0,200,length(cat.cm))
baselineLF.statJN.lf = matrix(0,200,length(cat.lf))
baselineLF.statJN.co = matrix(0,200,length(cat.cm))
for (i in 1:200){
	coeffJN                  = as.numeric(baselineLF.JN[i,])
	baselineLF.h2lfJN        = c(baselineLF.h2lfJN, sum(Vlf*coeffJN))
	baselineLF.h2cmJN        = c(baselineLF.h2cmJN, sum(Vcm*coeffJN))
	baselineLF.enrJN.lf[i,]  = as.numeric(((apply(coeffJN*t(sumannots[as.character(cat.lf),]),2,sum))/(baselineLF.h2lfJN[i]))/(Nlf.annot/Nlf))
	baselineLF.enrJN.co[i,]  = as.numeric(((apply(coeffJN*t(sumannots[as.character(cat.cm),]),2,sum))/(baselineLF.h2cmJN[i]))/(Ncm.annot/Ncm))
	baselineLF.ratioJN[i,]   = baselineLF.enrJN.lf[i,] / baselineLF.enrJN.co[i,]
    h2cat=apply(coeffJN*t(sumannots[as.character(cat.lf),]),2,sum); C=Nlf.annot; h2=baselineLF.h2lfJN[i]; M=Nlf
	baselineLF.statJN.lf[i,] = (h2cat/C)-((h2-h2cat)/(M-C))
	h2cat=apply(coeffJN*t(sumannots[as.character(cat.cm),]),2,sum); C=Ncm.annot; h2=baselineLF.h2cmJN[i]; M=Ncm
	baselineLF.statJN.co[i,] = (h2cat/C)-((h2-h2cat)/(M-C))

}
#write.table(baselineLF.enrJN.co,file=paste(file,".JNCVE" ,sep=""),quote=F,sep="\t",col.names=F,row.names=F)
#write.table(baselineLF.enrJN.lf,file=paste(file,".JNLFVE",sep=""),quote=F,sep="\t",col.names=F,row.names=F)

#compute sd estimates of h2lf and h2c
baselineLF.h2lf.sd  = se_jacknife(baselineLF.h2lf  , baselineLF.h2lfJN )
baselineLF.h2cm.sd  = se_jacknife(baselineLF.h2cm  , baselineLF.h2cmJN )

#compute sd estimates of LFVE, CVE, LFVE/CVE, LFVE stat and CVE stat
baselineLF.enrsd.lf  = NULL
baselineLF.enrsd.cm  = NULL
baselineLF.ratiosd   = NULL
baselineLF.statsd.lf = NULL
baselineLF.statsd.cm = NULL
for (j in 1:length(cat.lf)){
	baselineLF.enrsd.lf  = c(baselineLF.enrsd.lf  , se_jacknife(baselineLF.enr.lf[j]  , baselineLF.enrJN.lf[,j] ) )
	baselineLF.enrsd.cm  = c(baselineLF.enrsd.cm  , se_jacknife(baselineLF.enr.cm[j]  , baselineLF.enrJN.co[,j] ) )
	baselineLF.ratiosd   = c(baselineLF.ratiosd   , se_jacknife(baselineLF.ratio[j]   , baselineLF.ratioJN[,j]  ) )
	baselineLF.statsd.lf = c(baselineLF.statsd.lf , se_jacknife(baselineLF.stat.lf[j] , baselineLF.statJN.lf[,j]) )
	baselineLF.statsd.cm = c(baselineLF.statsd.cm , se_jacknife(baselineLF.stat.co[j] , baselineLF.statJN.co[,j]) )
}

# print h2c and h2lf
cat("        Common variant heritability (h2c)  =",round(baselineLF.h2cm,4),", se =",round(baselineLF.h2cm.sd,4),"\n") 
cat(" Low-frequency variant heritability (h2lf) =",round(baselineLF.h2lf,4),", se =",round(baselineLF.h2lf.sd,4),"\n") 
cat("\nWARNING: If you are using BOLT-LMM summary statistics,\nplease divide h2c and h2lf by Neff/N to interpret them\n(see our manuscript for justification)\n\n")


# print CVE and LFVE in $file.enrichments
out=cbind(
as.character(categories),
Ncm.annot/Ncm,baselineLF.enr.cm,baselineLF.enrsd.cm,
baselineLF.stat.co,baselineLF.statsd.cm,2*pnorm(-abs(baselineLF.stat.co/baselineLF.statsd.cm)),
Nlf.annot/Nlf,baselineLF.enr.lf,baselineLF.enrsd.lf,
baselineLF.stat.lf,baselineLF.statsd.lf,2*pnorm(-abs(baselineLF.stat.lf/baselineLF.statsd.lf)),
baselineLF.ratio,baselineLF.ratiosd,
2*pnorm(-abs((baselineLF.enr.lf-baselineLF.enr.cm)/sqrt((baselineLF.enrsd.cm)**2+(baselineLF.enrsd.lf)**2))),
2*pnorm(-abs(((baselineLF.enr.lf*Nlf.annot/Nlf)-(baselineLF.enr.cm*Ncm.annot/Ncm))/sqrt((baselineLF.enrsd.cm*Ncm.annot/Ncm)**2+(baselineLF.enrsd.lf*Nlf.annot/Nlf)**2))
))

colnames(out) = c(
"annotation" ,
"prop_common_snps", "CVE" ,"CVE_se" ,
"CVE_stat" ,"CVE_stat_se" ,"CVE_p",
"prop_lowfrq_snps", "LFVE","LFVE_se",
"LFVE_stat","LFVE_stat_se","LFVE_p",
"LFVE/CVE_ratio","LFVE/CVE_ratio_se",
"LFVE_vs_CVE_p","%h2lf_vs_%h2c_p")

write.table(out,file=paste(file,".enrichments",sep=""),quote=F,sep="\t",col.names=T,row.names=F)

# Outpouts columns are:
# * annotation: name of the binary annotation
# * prop_common_snps: proportion of common snps in the annotation
# * CVE: common variant enrichment (CVE)
# * CVE_se: common variant enrichment standard error
# * CVE_stat: common variant enrichment statistic
# * CVE_stat_se: common variant enrichment statistic standard error
# * CVE_p: common variant enrichment P value
# * prop_lowfrq_snps: proportion of low-frequency snps in the annotation
# * LFVE: low-frequency variant enrichment (LFVE)
# * LFVE_se: low-frequency variant enrichment standard error
# * LFVE_stat: low-frequency variant enrichment statistic
# * LFVE_stat_se: low-frequency variant enrichment statistic standard error
# * LFVE_p: low-frequency variant enrichment P value
# * LFVE/CVE_ratio: LFVE/CVE ratio
# * LFVE/CVE_ratio_se: LFVE/CVE ratio standard error
# * LFVE_vs_CVE_p: statistical difference between CVE and LFVE
# * %h2lf_vs_%h2c_p: statistical difference between % of common variant heritability 
#                    explained by the annotation and % of low-frequency variant  
#                    heritability explained by the annotation 

