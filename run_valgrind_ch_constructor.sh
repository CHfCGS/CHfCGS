data=ZipOrderingTest   #basicLine  #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ChainPair #ZipOrderingTest
echo "data: " "$data"

export GLIBCXX_FORCE_NEW

#chconstructor
in=$data"FMAXSPEED.txt"
out=$data"_ch_out.graph"
valgrind -v --num-callers=20 --leak-resolution=high --leak-check=yes --log-file="logfile.txt" --read-var-info=yes chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -f FMI_DIST -g FMI_CH
# --leak-check=full --show-leak-kinds=all
