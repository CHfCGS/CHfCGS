data=bremen   #basicLine  #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ChainPair #ZipOrderingTest
echo "data: " "$data"

export GLIBCXX_FORCE_NEW

#chmeasurer
in=$data"_ch_out.graph"

valgrind -v --max-stackframe=640020000  --num-callers=20 --leak-resolution=high --leak-check=yes --track-origins=yes --log-file="logfile.txt" --read-var-info=yes chmeasurer/build/ch_measurer -i data/CH/$in -p
# --leak-check=full --show-leak-kinds=all
