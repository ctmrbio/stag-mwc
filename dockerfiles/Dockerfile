FROM ubuntu:16.04
LABEL author="Fredrik Boulund"\
      maintainer="fredrik.boulund@ki.se"

RUN export DEBIAN_FRONTEND=noninteractive && \
    apt update && \
    apt-get install -yq build-essential wget unzip && \
    apt-get install -yq mailutils

WORKDIR /mocat2

RUN wget http://vm-lux.embl.de/~kultima/share/MOCAT/v2.0/MOCAT2.zip && \
    unzip MOCAT2.zip && \
    rm MOCAT2.zip && \
    cd MOCAT && \
    yes | perl setup.MOCAT.pl && \
    sed -i "s/EMAIL=''/EMAIL='\$EMAIL'/" /mocat2/MOCAT/src/runMOCAT.sh
# Set env variable EMAIL before running runMOCAT.sh

ENV PATH="${PATH}:/mocat2/MOCAT/src"
