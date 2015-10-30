data=bremen  #ChainPair2 # #bremen #ProDPExample #Intersection #baden-wuerttemberg #stuttgart-regbez
echo "data: " "$data"

in=$data"-latest.osm.pbf"
out=$data"FMAXSPEED.txt"

#graph_creator
#Osm/creator -g fmimaxspeedtext -t time -c ../../data/configs/car.cfg -o bremenFMAXSPEED_Ways.txt ~/Downloads/bremen-latest.osm.pbf
OsmGraphCreator/build/creator/creator -g fmimaxspeedtext -t time -c OsmGraphCreator/data/configs/car.cfg -o OsmGraphCreator/build/creator/$out osm_data/$in

#chconstructor
in=$out
out=$data"_ch_out.graph"
chconstructor/build/ch_constructor -i OsmGraphCreator/build/creator/$in -o chconstructor/build/$out  -p DP -f FMI_DIST -g FMI_CH

#visualizer
in=$out
cd simplestGraphRendering-master/build
./simple -gf ../../chconstructor/build/$in
