cd OsmGraphCreator

mkdir build
cd build
cmake ../
make

cd ..

cd chconstructor
sh create-build.sh
cd ..

cd simplestGraphRendering-master
sh build.sh
