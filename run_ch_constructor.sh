data=ChainPairD  #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ChainPairD
echo "data: " "$data"

#chconstructor
in=$data"FMAXSPEED.txt"
out=$data"_ch_out.graph"
chconstructor/build/ch_constructor -i OsmGraphCreator/build/creator/$in -o chconstructor/build/$out  -p DP -f FMI_DIST -g FMI_CH
