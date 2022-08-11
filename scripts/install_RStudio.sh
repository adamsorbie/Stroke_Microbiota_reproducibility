#!/bin/bash 

# install prerequisites 	
apt-get install gdebi-core
wget https://download2.rstudio.org/server/bionic/amd64/rstudio-2022.07.0-548-amd64.deb
gdebi rstudio-2022.07.0-548-amd64.deb
