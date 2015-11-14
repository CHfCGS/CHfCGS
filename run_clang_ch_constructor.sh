data=ZipOrderingTest   #basicLine  #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ChainPair #ZipOrderingTest
echo "data: " "$data"

export GLIBCXX_FORCE_NEW

#chconstructor
in=$data"FMAXSPEED.txt"
out=$data"_ch_out.graph"
chconstructor/clang-build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -f FMI_DIST -g FMI_CH

