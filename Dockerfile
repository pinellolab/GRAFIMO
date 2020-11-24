# Use the VG Docker image and build over it
FROM quay.io/vgteam/vg:v1.27.1

# Update Ubuntu Software repository
RUN apt-get update

# Upgrade Ubuntu Software repository
RUN apt-get upgrade -y

# make sure having gcc, make and git
RUN apt-get install build-essential -y
RUN apt-get install make -y
RUN apt-get install sudo -y
RUN apt-get install git -y

# Install curl
RUN apt-get install curl -y

# Build VG from source (static binaries no more available)
RUN sudo apt-get -y install build-essential git cmake pkg-config libncurses-dev libbz2-dev  \
                     protobuf-compiler libprotoc-dev libprotobuf-dev libjansson-dev \
                     automake libtool jq bc rs curl unzip redland-utils \
                     librdf-dev bison flex gawk lzma-dev liblzma-dev liblz4-dev \
                     libffi-dev libcairo-dev

# Download and install python3.7
WORKDIR /grafimo_wd
RUN apt-get install python3.7 -y
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 50
RUN apt-get install python3-distutils -y
RUN apt-get install manpages-dev -y
RUN apt-get install python-dev -y
RUN apt-get install python3-dev -y
RUN apt-get install libpython3.7-dev -y

# Install GRAFIMO external dependencies
RUN apt-get install tabix -y
RUN apt-get install graphviz -y

# Making sure that pip is available
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python get-pip.py
RUN apt-get install python3-pip -y
RUN pip3 install pandas
RUN pip3 install numpy
RUN pip3 install statsmodels
RUN pip3 install numba
RUN pip3 install Cython
RUN pip3 install setuptools
RUN pip3 install wheel
RUN pip3 install sphinx
RUN pip3 install -U colorama 

# remove the pip installation file
RUN rm get-pip.py

# Download and build GRAFIMO
RUN pip3 install grafimo

# check we have GRAFIMO, vg, tabix and graphviz available and in PATH
RUN which grafimo && which vg && which tabix && which dot
