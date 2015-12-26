data=bremen  #ChainPair2 # #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ZipOrderingTest #ProP_ILPChainPair
echo "data: " "$data"

in=$data"-latest.osm.pbf"
out=$data"FMAXSPEED.txt"

#graph_creator
#Osm/creator -g fmimaxspeedtext -t time -c ../../data/configs/car.cfg -o bremenFMAXSPEED_Ways.txt ~/Downloads/bremen-latest.osm.pbf
#OsmGraphCreator/build/creator/creator -g fmimaxspeedtext -t time -c OsmGraphCreator/data/configs/car.cfg -o data/Graph/$out data/OSM/$in

#chconstructor
in=$out
out=$data"_ch_out.graph"
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -d EDE -p DP -s DP -e VE -w P -f FMI_DIST -g FMI_CH
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -d EDE -p DP -s DP -e VE -f FMI_DIST -g FMI_CH
chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -p EDGE_DIFF -s DP -e VE -f FMI_DIST -g FMI_CH
 

#measurer
in=$out
cd chmeasurer/build
NOW=$(date +"%Y_%m_%d_%H_%M_%S")
LOGFILE="$NOW.log"


./ch_measurer -i ../../data/CH/$in -l  2>&1 | tee ../logs/$LOGFILE    #-d #-p -l 

#echo ./ch_measurer -i ../../data/CH/$in > ../logs/$LOGFILE -p | tee /dev/tty | foo
