#!/bin/bash

#singularity build ~/storage/MeerKATpol.simg docker://tpasini/pol_meerkat:latest
#singularity shell --pid --writable-tmpfs --cleanenv -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg/,/iranet/groups/lofar/containers/ ~/lofar/MeerKATpol.simg
singularity shell --pid --writable-tmpfs --cleanenv -B/home/baq1889,/localwork/fdg ~/lofar/MeerKATpol.simg
