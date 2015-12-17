data=bremen   #basicLine  #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ChainPair #ZipOrderingTest
echo "data: " "$data"

export GLIBCXX_FORCE_NEW

#chconstructor
in=$data"FMAXSPEED.txt"
out=$data"_ch_out.graph"
valgrind -v --num-callers=20 --main-stacksize=1000000000 --leak-resolution=high --leak-check=yes --log-file="logfile.txt" --read-var-info=yes chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s DP -e VE -w P -d EDE -c -v -f FMI_DIST -g FMI_CH
# --leak-check=full --show-leak-kinds=all
#--track-origins=yes 
