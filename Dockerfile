# Use the VG Docker image and build over it
FROM quay.io/vgteam/vg:v1.19.0

# Update Ubuntu Software repository
RUN apt-get update

# Upgrade Ubuntu Software repository
RUN apt-get upgrade -y

# Download and install python3.7
RUN apt-get install python3.7 -y
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 50
RUN apt-get install python3-distutils -y

# Making sure that pip is available
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
RUN python get-pip.py
RUN pip install pandas
RUN pip install numpy

# remove the pip installation file
RUN rm get-pip.py

# Make sure to have git available
RUN apt-get install git -y

# Download and build GRAFIMO
RUN git clone https://github.com/InfOmics/GRAFIMO.git
RUN pip install GRAFIMO/

# Set expose to port 80 and 443
EXPOSE 80 443
