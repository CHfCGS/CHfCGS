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
chconstructor/build/ch_constructor -i data/Graph/$in -o data/CH/$out -p DP -f FMI_DIST -g FMI_CH
# 

#visualizer
in=$out
cd chmeasurer/build
./ch_measurer -gf ../../data/CH/$in