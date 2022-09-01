#!/usr/bin/bash
### Setting a working directory structure for tools ======================================
workdir=$(pwd)
mkdir raw_data processed_data scripts results software
mkdir processed_data/fpkm
mkdir software/bin

softwarepath=$workdir/software/bin
echo $softwarepath
cd software

### setting up the PATH for installed softwares
export PATH=$PATH:$softwarepath

