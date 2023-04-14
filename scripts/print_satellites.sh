#!/bin/bash

listfile=fits_list
if [[ -n "$1" && "$1" != "-" ]]; then
   listfile=$1
fi

tle=sat.tle
if [[ -n "$2" && "$2" != "-" ]]; then
   tle=$2
fi

outdir=satellites/

if [[ ! -s ${tle} ]]; then
   echo "ERROR : TLE file $tle not found -> please provide proper value in 2nd parameter"
   exit
fi


mkdir -p ${outdir}

for fitsname in `cat $listfile`
do
   # chan_204_20220913T181307_I_diff.fits
   utc=`echo $fitsname | cut -b 10-24 | awk '{gsub("T","_");print $0;}'`
   ux=`date2date -ut2ux=${utc} | awk '{print $3;}'`
   echo "$fitsname -> $utc UTC -> $ux"
   
   outregfile=${fitsname%%fits}reg
   outfile=${fitsname%%fits}txt
   stdoutfile=${fitsname%%fits}out
   
   # sattest %.4f -tle=%s -all -mwa -qth=mro.qth -ra=%.8f  -dec=%.8f -outregfile=%.4f.reg  -outfile=%.4f.txt -print_header -no_spaces -min_elevation=-1000.0000 -radius=360 > %.4f.out 2>&1
   echo "sattest $ux -tle=${tle} -all -mwa -qth=mro.qth -outregfile=${outdir}/${outregfile} -outfile=${outdir}/${outfile} -print_header -no_spaces -min_elevation=0.00 > ${outdir}/${stdoutfile}"
   sattest $ux -tle=${tle} -all -mwa -qth=mro.qth -outregfile=${outdir}/${outregfile} -outfile=${outdir}/${outfile} -print_header -no_spaces -min_elevation=0.00 > ${outdir}/${stdoutfile}
done
