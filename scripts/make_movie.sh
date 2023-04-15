#!/bin/bash

movie_png_rate=1
if [[ -n "$1" && "$1" != "-" ]]; then
   movie_png_rate=$1
fi

dir=images/movie/
if [[ -n "$2" && "$2" != "-" ]]; then
   dir=$2
fi

publish=1
if [[ -n "$3" && "$3" != "-" ]]; then
   publish=$3
fi

image_type="jpg"
if [[ -n "$4" && "$4" != "-" ]]; then
   image_type="$4"
fi

template="1*${image_type}"
if [[ -n "$5" && "$5" != "-" ]]; then
   template="$5"
fi


echo "#############################################"
echo "PARAMETERS "
echo "#############################################"
echo "movie_png_rate = $movie_png_rate"
echo "dir = $dir"
echo "publish = $publish"
echo "#############################################"


mkdir -p ${dir}
cd ${dir}

echo "make_movie.sh started at :"
date


i=0
for png in `ls ../${template}`; 
do      
   if [[ $(($i % $movie_png_rate)) == 0 ]]; then
      out_i=$(($i / $movie_png_rate))
      i_str=`echo $out_i | awk '{printf("%06d\n",$1);}'`;    
      ln -s ${png} image_${i_str}.${image_type}; 
   fi
   i=$(($i+1));
done      

#      mencoder mf://*.png -mf w=800:h=600:fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o 20200226_eda2.mpg
#      mencoder mf://*.png -mf w=800:h=600:fps=10:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o 20200226_eda2.avi      
echo "rm -f sky_last_24h.mp4"
rm -f sky_last_24h.mp4
echo "ffmpeg -framerate 25 -i image_%6d.${image_type} -c:v libx264 -vf \"drawtext=fontfile=/usr/share/fonts/truetype/freefont/FreeSans.ttf:text=%{n}:x=(w-tw)/2:y=h-(2*lh): fontcolor=white:box=1:boxcolor=0x00000099\" -pix_fmt yuv420p sky_last_24h.mp4"
ffmpeg -framerate 25 -i image_%6d.${image_type} -c:v libx264 -vf "drawtext=fontfile=/usr/share/fonts/truetype/freefont/FreeSans.ttf:text=%{n}:x=(w-tw)/2:y=h-(2*lh): fontcolor=white:box=1:boxcolor=0x00000099" -pix_fmt yuv420p sky_last_24h.mp4

# if [[ $publish -gt 0 ]]; then
#   echo "scp sky_last_24h.mp4 aavs1-server:/exports/eda/${station_name}/tv/"
#   scp sky_last_24h.mp4 aavs1-server:/exports/eda/${station_name}/tv/
#
#   echo "scp channel.txt aavs1-server:/exports/eda/${station_name}/tv/"
#   scp channel.txt aavs1-server:/exports/eda/${station_name}/tv/
#else
#   echo "WARNING : publishing of results on the WWW server is not required"
#fi
cd ../../

echo "make_movie.sh finished at :"
date
