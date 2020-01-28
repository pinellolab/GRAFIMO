# Use the VG Docker image and build over it
FROM quay.io/vgteam/vg:v1.21.0

# Update Ubuntu Software repository
RUN apt-get update

# Upgrade Ubuntu Software repository
RUN apt-get upgrade -y

# Download and install python3.7
RUN apt-get install python3.7 -y
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 50
RUN apt-get install python3-distutils -y
RUN apt-get install build-essential -y
RUN apt-get install manpages-dev -y
RUN apt-get install python-dev -y
RUN apt-get install python3-dev -y
RUN apt-get install libpython3.7-dev -y

# Making sure that pip is available
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python get-pip.py
RUN apt-get install python3-pip -y
RUN pip3 install pandas
RUN pip3 install numpy
RUN pip3 install statsmodels
RUN pip3 install numba
RUN pip3 install Cython

# remove the pip installation file
RUN rm get-pip.py

# Make sure to have git available
RUN apt-get install git -y

# Download and build GRAFIMO
RUN git clone https://github.com/pinellolab/GRAFIMO.git
RUN pip3 install GRAFIMO/

# Set expose to port 80 and 443
EXPOSE 80 443
