data=bremen #unbreakedDPVEPunfiltered #ChainPair2 # #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez #ZipOrderingTest #ProP_ILPChainPair #PairsToIdentify #IntersectionDiff #heapTest
echo "data: " "$data"

in=$data"-latest.osm.pbf"
out=$data"FMAXSPEED.txt"

#graph_creator
#OsmGraphCreator/build/creator/creator -g fmimaxspeedtext -t time -c OsmGraphCreator/data/configs/car.cfg -o data/Graph/$out data/OSM/$in

#chconstructor

options1=" -p DP -d EDE "
#options1="  "

_options2=" -s BU -e KE "
options2=${1-$_options2}

#options2=" -p NONE " 
echo "options2: " "$options2"

out_options2=`echo "$options2" | sed "s/ /_/g"`
in=$out
out=$data$out_options2"ch_out.graph"
NOW=$(date +"%d_%H")
LOGFILE="$out$NOW.log"

#chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  -p DP -s DP -e VE -w P -d EDE -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE

chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out  $options1 $options2 -f FMI_DIST -g FMI_CH 2>&1 | tee chconstructor/log/$LOGFILE # -d EDE


#measurer
in=$out
cd chmeasurer/build
NOW=$(date +"%d_%H")
LOGFILE="$in$NOW.log"


./ch_measurer -i ../../data/CH/$in -l 2>&1 | tee ../logs/$LOGFILE    #-e -c -d -l -p -v 


#visualizer
in=$out
cd ../..
cd simplestGraphRendering-master/build
#./simple -gf ../../data/CH/$in
