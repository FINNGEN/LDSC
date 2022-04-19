from ubuntu:18.04

<<<<<<< HEAD
RUN apt-get update && apt-get install git wget python-minimal python-pip -y && git clone https://github.com/bulik/ldsc.git &&
ENV PATH /ldsc:$PATH
=======
ENV PATH /ldsc:$PATH
RUN apt-get update && apt-get install git wget -y

RUN git clone https://github.com/bulik/ldsc.git
>>>>>>> 5930f7d3dfe462bdac2bb99cdfc30f86553b6033

RUN wget http://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 && tar -xjf eur_w_ld_chr.tar.bz2 && \
    rm eur_w_ld_chr.tar.bz2 && wget http://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 && \
     bzip2 -d w_hm3.snplist.bz2

<<<<<<< HEAD
RUN pip install pandas scipy bitarray
=======
# ptyhon2 for ldsc
RUN apt-get install libjpeg-dev zlib1g-dev python-minimal python-pip -y
RUN pip install pandas scipy bitarray

# python3 for custom scripts
RUN apt-get install python3  python3-pip -y
RUN pip3 install Pillow
RUN pip3 install seaborn
RUN apt-get install jq -y

ADD scripts ../scripts/
>>>>>>> 5930f7d3dfe462bdac2bb99cdfc30f86553b6033
