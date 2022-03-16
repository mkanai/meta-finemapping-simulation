FROM nobodyxu/intel-mkl:2020-debian-buster as builder
ENV commit=6

RUN apt-get update --fix-missing && \
    apt-get install -y build-essential git curl && \
    git clone https://github.com/mkanai/plink-ng && \
    cd plink-ng && \
    git checkout add_se_meta && \
    curl -O http://zlib.net/zlib-1.2.11.tar.gz && \
    tar -xzvf zlib-1.2.11.tar.gz && \
    cd zlib-1.2.11 && \
    ./configure && \
    make && \
    cd ../1.9/build_static64 && \
    sed -i 's/libz-64.a/libz.a/' Makefile && \
    make && \
    cp plink /usr/bin/

FROM gcr.io/encode-uk-biobank-restrict/plink-htslib:20211119 as htslib
FROM google/cloud-sdk:slim as final

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install R
RUN apt-get update --fix-missing && \
    apt-get install -y software-properties-common && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' && \
    add-apt-repository 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' && \
    apt-get update --fix-missing && \
    apt-get install -y r-base && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("argparse", "dplyr", "stringr", "data.table", "here", "R.utils"), noCache = TRUE, repos="https://cloud.r-project.org")'

COPY --from=builder /usr/bin/plink /usr/bin/plink
COPY --from=htslib /usr/bin/bgzip /usr/bin/bgzip

VOLUME [ "/root/.config/gcloud" ]
CMD [ "/bin/bash" ]
