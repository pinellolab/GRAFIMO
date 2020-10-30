#!/bin/bash

# Install external dependendecies required to run GRAFIMO:
# - vg
# - tabix
# - graphviz

mkdir vg
cd vg
wget https://github.com/vgteam/vg/releases/download/v1.27.1/vg
chmod +x ./vg
vgpath=`pwd`
export PATH=$vgpath:$PATH 
cd ..

#sudo apt-get install tabix
#sudo apt-get install graphviz
