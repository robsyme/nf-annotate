FROM ubuntu:17.10

MAINTAINER Rob Syme <rob.syme@gmail.com>

RUN apt-get update \
&& apt-get install -qqy \
git \
wget \
unzip \
build-essential

WORKDIR /usr/local

# Install Augustus
RUN apt-get install -qqy \
libboost-iostreams-dev \
libboost-graph-dev \
libsuitesparse-dev \
liblpsolve55-dev \
libbamtools-dev \
libgsl-dev \ 
zlib1g-dev

RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.3.tar.gz \
&& tar -xvf augustus*.tar.gz \
&& rm augustus*.tar.gz
COPY patches/augustusmk.patch /usr/local/augustus/
RUN cd augustus \
&& patch < augustusmk.patch \
&& make \
&& make install


# Install ProgressiveCactus
RUN apt-get install -qqy \
 python \
 zlib1g-dev \
 python-dev \
 python-numpy \
 iputils-ping \
 time

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

# Install newick utils
RUN wget http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz \
&& tar -xvf newick-utils*.tar.gz \
&& rm newick-utils*.tar.gz \
&& mv newick-utils-* newick-utils \
&& cd newick-utils \
&& ./configure \
&& make

# Install aragorn (tRNA prediction)
RUN apt-get install -qqy aragorn

# Install lua
RUN apt-get install -qqy lua5.1 

# Install exonerate
RUN apt-get install -qqy exonerate

# Install genometools
RUN apt-get install -qqy genometools

# Install Cufflinks 
RUN apt-get install -qqy cufflinks

# Install CodingQuarry
RUN wget https://downloads.sourceforge.net/project/codingquarry/CodingQuarry_v2.0.tar.gz \
&& tar -xvf CodingQuarry*.tar.gz \
&& rm CodingQuarry*.tar.gz \
&& mv CodingQuarry_v2.0 codingquarry \
&& cd codingquarry \
&& make

ENV AUGUSTUS_CONFIG_PATH /usr/local/augustus/config
ENV PYTHONPATH /usr/local/progressiveCactus/submodules
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/augustus/bin:/usr/local/augustus/scripts:/usr/local/progressiveCactus/bin:/usr/local/progressiveCactus/submodules/kentToolBinaries:/usr/local/hisat2:/usr/local/mash:/usr/local/progressiveCactus/submodules/hal/bin:/usr/local/newick-utils/src:/usr/local/codinquarry
