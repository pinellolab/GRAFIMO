# Use the VG Docker image and build over it
#FROM quay.io/vgteam/vg:v1.27.1
FROM ubuntu:18.04 AS run

# Update Ubuntu Software repository
RUN apt-get update

# Upgrade Ubuntu Software repository
RUN apt-get upgrade -y

# Install wget and curl
RUN apt-get install wget -y
RUN apt-get install curl -y

# Download vg binaries and set them available in PATH
RUN mkdir /vg_dir
RUN cd /vg_dir
RUN wget https://github.com/vgteam/vg/releases/download/v1.27.1/vg
RUN chmod +x ./vg
RUN mv ./vg /usr/bin/
RUN cd $HOME

# Download and install python3.7
RUN apt-get install python3.7 -y
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 50
RUN apt-get install python3-distutils -y
RUN apt-get install build-essential -y
RUN apt-get install manpages-dev -y
RUN apt-get install python-dev -y
RUN apt-get install python3-dev -y
RUN apt-get install libpython3.7-dev -y

# Make sure to have git available
RUN apt-get install git -y

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

# Set expose to port 80 and 443
EXPOSE 80 443
