data=bremen #ChainPair2 # #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ZipOrderingTest #ProP_ILPChainPair #PairsToIdentify
echo "data: " "$data"

in=$data"-latest.osm.pbf"
out=$data"FMAXSPEED.txt"

#graph_creator
#Osm/creator -g fmimaxspeedtext -t time -c ../../data/configs/car.cfg -o bremenFMAXSPEED_Ways.txt ~/Downloads/bremen-latest.osm.pbf
#OsmGraphCreator/build/creator/creator -g fmimaxspeedtext -t time -c OsmGraphCreator/data/configs/car.cfg -o data/Graph/$out data/OSM/$in

#chconstructor
in=$out
out=$data"_ch_out.graph"
NOW=$(date +"%Y_%m_%d_%H_%M_%S")
LOGFILE="$NOW.log"
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s NTH -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s DP -e VE -d EDE -w ZO -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s BU -e KE -w CD -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE #-d EDE

#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s NTH -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE
chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s DP -e VE -c -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE #-d EDE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s DP -e VE -d EDE -w P -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE

#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -d EDE -p DP -s DP -c -e VE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -p EDGE_DIFF -v -s DP -e VE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -s DP -e VE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE
 

#measurer
in=$out
cd chmeasurer/build
NOW=$(date +"%Y_%m_%d_%H_%M_%S")
LOGFILE="$NOW.log"


./ch_measurer -i ../../data/CH/$in -e  2>&1 | tee ../logs/$LOGFILE    #-e -c -d -l -p -v 

#echo ./ch_measurer -i ../../data/CH/$in > ../logs/$LOGFILE -p | tee /dev/tty | foo

#visualizer
in=$out
cd ../..
cd simplestGraphRendering-master/build
#./simple -gf ../../data/CH/$in
