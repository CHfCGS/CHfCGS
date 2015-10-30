data=bremen   #basicLine  #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ChainPair
echo "data: " "$data"

#chconstructor
in=$data"FMAXSPEED.txt"
out=$data"_ch_out.graph"
valgrind --leak-check=yes --log-file="logfile.txt" --read-var-info=yes chconstructor/build/ch_constructor -i OsmGraphCreator/build/creator/$in -o chconstructor/build/$out  -p DP -f FMI_DIST -g FMI_CH
