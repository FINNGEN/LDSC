from ubuntu:18.04

ENV PATH /ldsc:$PATH
RUN apt-get update && apt-get install git wget -y

# https://github.com/kkellysci/munge_small_pvals.git (version that can handle low pvals)
RUN git clone https://github.com/bulik/ldsc.git

RUN wget http://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 && tar -xjf eur_w_ld_chr.tar.bz2 && \
    rm eur_w_ld_chr.tar.bz2 && wget http://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 && \
     bzip2 -d w_hm3.snplist.bz2

# ptyhon2 for ldsc
RUN apt-get install libjpeg-dev zlib1g-dev python-minimal python-pip -y
RUN pip install pandas scipy bitarray

# python3 for custom scripts
RUN apt-get install python3  python3-pip -y
RUN pip3 install Pillow
RUN pip3 install seaborn
RUN apt-get install jq -y

ADD scripts ../scripts/
