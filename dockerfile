FROM centos:7 as build

WORKDIR /opt
RUN yum -y update \
    && yum -y install centos-release-scl scl-utils wget \
    && yum -y install gzip devtoolset-8-gcc-c++ git \
RUN git clone https://github.com/WGLab/LongGF \
    && wget https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz \
    && gzip -d sambamba-1.0.1-linux-amd64-static.gz
WORKDIR LongGF/bin
RUN scl enable devtoolset-8 -- bash
RUN scl enable devtoolset-8 -- bash
RUN g++ -g -O3 -std=c++11 -I ./include -L ./lib -Wl,--enable-new-dtags,-rpath,"\$ORIGIN"/lib \
    -lhts -o LongGF _com_fun_.c _gtf_struct_.c get_gfFrombam.c -Wl,--no-as-needed