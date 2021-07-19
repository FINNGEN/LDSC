from ubuntu:18.04

RUN apt-get update && apt-get install git wget python-minimal python-pip -y && git clone https://github.com/bulik/ldsc.git &&
ENV PATH /ldsc:$PATH

RUN wget http://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 && tar -xjf eur_w_ld_chr.tar.bz2 && \
    rm eur_w_ld_chr.tar.bz2 && wget http://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 && \
     bzip2 -d w_hm3.snplist.bz2

RUN pip install pandas scipy bitarray
