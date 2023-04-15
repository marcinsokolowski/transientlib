#!/bin/bash

ls 1*.txt > list
split --lines=100 list # split to avoid memory leaks in python !!!

ls x?? | awk '{print "python ~/github/transientlib/scripts/plotevents.py "$1" --informat=sattest > "$1".out 2>&1";}' > plot_all!

chmod +x plot_all!

echo "./plot_all!"
./plot_all!

echo "~/github/transientlib/scripts/make_movie.sh"
~/github/transientlib/scripts/make_movie.sh
