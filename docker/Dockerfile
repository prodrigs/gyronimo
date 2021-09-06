FROM        ubuntu:latest
MAINTAINER  jsferreira 

# update and install dependencies
RUN    apt-get update && \
       apt-get install -y \
               software-properties-common \
               wget && \
       add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
       apt-get update && \
       apt-get install -y \
               make \
               git \
               curl \
               vim \
               cmake \
               gsl-bin \
               libgsl23 \
               libgsl-dev \
               libboost1.71-all-dev \
               doxygen \
               gcc-10 g++-10 gcc-10-base && \
       update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 && \
       update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100 && \
       rm -rf /var/lib/apt/lists/*

CMD    bash
