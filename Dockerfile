FROM ubuntu:17.10

MAINTAINER Rob Syme <rob.syme@gmail.com>

RUN apt-get update \
&& apt-get install -qqy \
git \
wget \
unzip \
build-essential

# Install Augustus
WORKDIR /usr/local
RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.3.tar.gz \
&& tar -xvf augustus*.tar.gz

# Install ProgressiveCactus
RUN apt-get install -qqy \
 python \
 zlib1g-dev \
 python-dev \
 python-numpy

RUN ln -s /usr/lib/python2.7/plat-*/_sysconfigdata_nd.py /usr/lib/python2.7/

RUN git clone git://github.com/glennhickey/progressiveCactus.git \
&& cd progressiveCactus \
&& git checkout tags/0.1 -b 0.1 \
&& git submodule update --init

COPY patches/ktremote.patch progressiveCactus/submodules/kyototycoon/
COPY patches/ktulog.patch progressiveCactus/submodules/kyototycoon/
RUN cd progressiveCactus/submodules/kyototycoon && patch < ktulog.patch && patch < ktremote.patch
COPY patches/kyotocabinet-1.2.76-gcc6.patch progressiveCactus/submodules/kyotocabinet/
RUN cd progressiveCactus/submodules/kyotocabinet && patch < kyotocabinet-1.2.76-gcc6.patch
RUN cd progressiveCactus && make

# Install Hisat2
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip \
&& unzip hisat2-2.1.0-Linux_x86_64.zip \
&& mv hisat2-2.1.0 hisat2 \
&& rm *.zip

# Install MASH v2.0
RUN wget https://github.com/marbl/Mash/releases/download/v2.0/mash-Linux64-v2.0.tar \
&& tar -xvf mash*.tar \
&& rm mash*.tar \
&& mv mash-* mash

# Install R and packages
RUN apt-get install -qqy r-base \
 r-cran-phangorn \
 r-cran-reshape2 \
 r-cran-magrittr \
 r-cran-dplyr

ENV PYTHONPATH /usr/local/progressiveCactus/submodules
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/augustus:/usr/local/progressiveCactus/bin:/usr/local/progressiveCactus/submodules/kentToolBinaries:/usr/local/hisat2:/usr/local/mash
