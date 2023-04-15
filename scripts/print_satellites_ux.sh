#!/bin/bash

start_ux=`date +%s`
if [[ -n "$1" && "$1" != "-" ]]; then
   start_ux=$1
fi

duration=600
if [[ -n "$2" && "$2" != "-" ]]; then
   duration=$2
fi

tle=satelitesdb.tle
if [[ -n "$3" && "$3" != "-" ]]; then
   tle=$3
fi

outdir=satellites/

if [[ ! -s ${tle} ]]; then
   echo "ERROR : TLE file $tle not found -> please provide proper value in 2nd parameter"
   exit
fi

if [[ ! -s mro.qth ]]; then
   echo "ERROR : geo. location config file mro.qth not found (get it from https://github.com/marcinsokolowski/transientlib/tree/main/config )"
   exit
fi

mkdir -p ${outdir}

ux=${start_ux}
end_ux=$(($ux+$duration))

echo "#############################################"
echo "PARAMETERS:"
echo "#############################################"
echo "Start UX = $start_ux"
echo "Duration = $duration"
echo "End UX   = $end_ux"
echo "#############################################"

while [[ $ux -le $end_ux ]];
do
   # chan_204_20220913T181307_I_diff.fits
   utc=`ux2ut! $ux`
   echo "UX TIME =  $ux -> $utc UTC"
   
   outregfile=${ux}.reg
   outfile=${ux}.txt
   stdoutfile=${ux}.out
   
   # sattest %.4f -tle=%s -all -mwa -qth=mro.qth -ra=%.8f  -dec=%.8f -outregfile=%.4f.reg  -outfile=%.4f.txt -print_header -no_spaces -min_elevation=-1000.0000 -radius=360 > %.4f.out 2>&1
   echo "sattest $ux -tle=${tle} -all -mwa -qth=mro.qth -outregfile=${outdir}/${outregfile} -outfile=${outdir}/${outfile} -print_header -no_spaces -min_elevation=0.00 > ${outdir}/${stdoutfile}"
   sattest $ux -tle=${tle} -all -mwa -qth=mro.qth -outregfile=${outdir}/${outregfile} -outfile=${outdir}/${outfile} -print_header -no_spaces -min_elevation=0.00 > ${outdir}/${stdoutfile}

   ux=$(($ux+1))
done

