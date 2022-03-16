FROM google/cloud-sdk:slim

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y wget make zlib1g-dev unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV HTSLIB_VERSION 1.14
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 -P ~ && \
    tar xjvf ~/htslib-${HTSLIB_VERSION}.tar.bz2 -C ~ && \
    cd ~/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr --disable-bz2 --disable-lzma && \
    make && make install && \
    rm -r ~/htslib-${HTSLIB_VERSION}*

# Stable (beta 6.24, 6 Jun)
ENV PLINK_VERSION 20210606
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK_VERSION}.zip -P ~ && \
    mkdir -p ~/plink && \
    unzip ~/plink_linux_x86_64_${PLINK_VERSION}.zip -d ~/plink && \
    cp ~/plink/plink /usr/bin/plink && \
    rm -r ~/plink*

# Development (11 Oct)
ENV PLINK2_VERSION 20211011
RUN wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_${PLINK2_VERSION}.zip -P ~ && \
    mkdir -p ~/plink2 && \
    unzip ~/plink2_linux_x86_64_${PLINK2_VERSION}.zip -d ~/plink2 && \
    cp ~/plink2/plink2 /usr/bin/plink2 && \
    rm -r ~/plink2*

VOLUME [ "/root/.config/gcloud" ]
CMD [ "/bin/bash" ]
