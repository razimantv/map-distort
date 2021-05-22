#!/bin/bash

mkdir -p images data/processed

cd data
echo "Generating map data"
echo "" original.dat
echo ../states_boundary.dat population.dat 1 | ../../weighted_distort \
  > processed/original.dat 2>processed/log_original.dat

for i in *.dat
do
  echo "" $i
  echo ../states_boundary.dat $i 10000 | ../../weighted_distort \
    > processed/$i 2>processed/log_$i
done
echo done
echo ""

cd ..
echo "Generating svg plots"
gnuplot plotter
echo done
echo ""

echo "Generating png images"
cd images
for i in *.svg
do
  echo "" "${i%.svg}.png"
  convert $i "${i%.svg}.png"
done
echo done
