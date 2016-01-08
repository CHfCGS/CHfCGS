data=bremen #unbreakedDPVEPunfiltered #ChainPair2 # #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ZipOrderingTest #ProP_ILPChainPair #PairsToIdentify #IntersectionDiff #heapTest
echo "data: " "$data"

in=$data"-latest.osm.pbf"
out=$data"FMAXSPEED.txt"

#graph_creator
#Osm/creator -g fmimaxspeedtext -t time -c ../../data/configs/car.cfg -o bremenFMAXSPEED_Ways.txt ~/Downloads/bremen-latest.osm.pbf
#OsmGraphCreator/build/creator/creator -g fmimaxspeedtext -t time -c OsmGraphCreator/data/configs/car.cfg -o data/Graph/$out data/OSM/$in

#chconstructor
options1=" -p DP -d EDE "
options2=" -s DP -e VE -w P "
#out_options2="${$options2// /_}"
#out_options2=$options2 | sed -e 's/ /_/g'
#new_db_name=`echo "$new_db_name" | sed "s/$replace_string/$replace_with/"`
out_options2=`echo "$options2" | sed "s/ /_/g"`
#echo $options2
#echo $out_options2
in=$out
out=$data$out_options2"ch_out.graph"
NOW=$(date +"%d_%H")
LOGFILE="$out$NOW.log"
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s NTH -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s DP -e VE -d EDE -w ZO -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s BU -e KE -w CD -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE #-d EDE

#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s NTH -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s DP -e VE -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE #-d EDE

#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s DP -e VE -w P -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE
chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  $options1 $options2 -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE

#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -d EDE -p DP -s DP -c -e VE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -p EDGE_DIFF -v -s DP -e VE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE
#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -s DP -e VE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE
 

#measurer
in=$out
cd chmeasurer/build
NOW=$(date +"%d_%H")
LOGFILE="$in$NOW.log"


./ch_measurer -i ../../data/CH/$in -p  2>&1 | tee ../logs/$LOGFILE    #-e -c -d -l -p -v 


#visualizer
in=$out
cd ../..
cd simplestGraphRendering-master/build
#./simple -gf ../../data/CH/$in
