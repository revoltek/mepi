# mepi

A container file with all required softwares and packages can be found at docker://tpasini/pol_meerkat:latest

To run the singularity one can use (for instance):

`singularity build ~/MeerKATpol.simg docker://tpasini/pol_meerkat:latest`

`singularity shell --pid --writable-tmpfs --cleanenv -B/homes/fdg,/local/work/fdg,/iranet/groups/ulu/fdg/,/iranet/groups/lofar/containers/ ~/MeerKATpol.simg`

For L and S band the data must have X and Y flipped using the script "correct_parang.py" by Ben Hugo
https://github.com/bennahugo/LunaticPolarimetry/blob/master/correct_parang.py
`python correct_parang.py -f xxx --noparang  --applyantidiag`
the script writes the output in the CORRECTED_DATA column, so the corrected data needs to be split after the correction

plots are done with shadems https://github.com/ratt-ru/shadeMS
