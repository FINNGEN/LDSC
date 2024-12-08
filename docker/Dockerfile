FROM eu.gcr.io/finngen-refinery-dev/bioinformatics:0.8


RUN apt-get update && apt-get upgrade --yes

# GET LDSC
RUN apt-get install wget -y


# PYTHON3 VERSION!
RUN wget https://github.com/belowlab/ldsc/archive/refs/heads/2-to-3.zip  -O master.zip
RUN unzip master.zip
RUN pip install -r /ldsc-2-to-3/requirements.txt 
RUN apt-get install python3  python3-pip -y
RUN pip3 install Pillow seaborn bitarray
RUN apt-get install jq -y

ADD scripts ../scripts/
