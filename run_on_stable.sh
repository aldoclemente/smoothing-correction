#!/bin/bash

IMAGE=rhub/debian-gcc-release

docker pull $IMAGE 

docker run --name=tmp-debian-gcc -v $(pwd)/../:/root/smoothing-correction --rm -ti $IMAGE bin/bash -c '

export RGL_USE_NULL=TRUE
export DISPLAY=99.0

apt-get install -y libgl1-mesa-dev libglu1-mesa-dev
apt-get install -y libxml2-dev
apt-get install -y libssl-dev

cd /root/smoothing-correction/ 
Rscript dependencies.R

cd DensityEstimation/
Rscript run_on_stable.R

cd ../SpatialRegression/
Rscript run_on_stable.R

exit
'
