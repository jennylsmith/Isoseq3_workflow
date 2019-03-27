FROM ubuntu:16.04
## aws, curl, git, unzip, python and python3   are already installed in fredhutch/ubuntu
# https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda
# http://www.science.smith.edu/dftwiki/index.php/Tutorial:_Docker_Anaconda_Python_--_4


#Updating Ubuntu packages
RUN apt-get update && yes|apt-get upgrade
RUN apt-get install -y bzip2 curl unzip

#Create working directory
WORKDIR /home/ubuntu/

#Installing Anaconda
RUN curl -LO https://repo.continuum.io/archive/Anaconda2-2018.12-Linux-x86_64.sh
RUN bash ./Anaconda2-2018.12-Linux-x86_64.sh -b -p /home/ubuntu/anaconda2
RUN rm Anaconda2-2018.12-Linux-x86_64.sh

#export path to find executable
ENV PATH /home/ubuntu/anaconda2/bin:$PATH

#Update the conda installation
RUN conda update conda

#Install Isoseq3 and its dependencies
RUN conda install -p /home/ubuntu/anaconda2 -c conda-forge awscli
RUN conda install -p /home/ubuntu/anaconda2 biopython
RUN conda install -p /home/ubuntu/anaconda2 -c bioconda bx-python
RUN conda install -p /home/ubuntu/anaconda2 -c bioconda isoseq3
RUN conda install -p /home/ubuntu/anaconda2 -c bioconda pbcoretools # for manipulating PacBio datasets
RUN conda install -p /home/ubuntu/anaconda2 -c bioconda bamtools    # for converting BAM to fasta
RUN conda install -p /home/ubuntu/anaconda2 -c bioconda pysam       # for making CSV

#Update Isoseq3 version
RUN conda update -p /home/ubuntu/anaconda2 -c bioconda isoseq3

#Install Fetch and Run
RUN curl -LO https://github.com/jennylsmith/aws-batch-helpers/archive/master.zip
RUN unzip master.zip
RUN cp aws-batch-helpers-master/fetch-and-run/fetch_and_run.sh /usr/local/bin/
RUN rm -rf master.zip aws-batch-helpers-master

#Set entry point
ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]
