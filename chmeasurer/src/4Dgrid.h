#pragma once

/* 
 * File:   4Dgrid.h
 * Author: tobias
 *
 * Created on 13. October 2015, 16:14
 */

#include <vector>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "nodes_and_edges.h"
#include "chains.h"
#include "grid.h"


using namespace std;

//typedef vector<vector<int> > twoDvector;
//typedef vector<vector<twoDvector> > fourDvector;



template <class GraphT>
class FourDGrid {
    //struct Gridpoint;
    
    struct Gridpoint {
        ChainsOfType::iterator chain;
        Gridpoint *partnerGridpoint;
        bool isMatched;
        
        Gridpoint(ChainsOfType::iterator chain, Gridpoint *partnerGridpoint, bool isMatched):
            chain(chain), partnerGridpoint(partnerGridpoint), isMatched(isMatched) {};
    };
    
    typedef list<Gridpoint> GridCell;

    //typedef vector<vector<GridCell> > twoDvector;
    typedef vector<vector<GridCell> > twoDvector;
    typedef vector<vector<twoDvector> > fourDvector;
private:
    //const std::vector<NodeT> &nodes;
    //constexpr auto M_PI = 3.14159265358979323846;
    double R = 6371.009; //Erdradius
    const double epsilon = 0.0001;
    //Longitude ~ x Latitude ~ y
    double MINLONGITUDE;
    double MAXLONGITUDE;
    double MINLATITUDE;
    double MAXLATITUDE;

    double gridwidth;
    double gridheight;
    double cellsizex;
    double cellsizey;

    int gridsidesize;
    const GraphT &base_graph;
    //ChainsOfType* chainsptr;

    fourDvector Grid;
    std::vector<int> testvector;
    //vector<int> NextChainInSameCell;
    //vector<int> gridtest;

    double pythagoras(double a, double b) {
        return sqrt(pow(a, 2) + pow(b, 2));
    }

    void calculateBorders() {
        for (u_int32_t k = 0; k < base_graph.getNrOfNodes(); k++) {
            if (base_graph.getLon(k) > MAXLONGITUDE) {
                MAXLONGITUDE = base_graph.getLon(k) + epsilon;
            }
            if (base_graph.getLon(k) < MINLONGITUDE) {
                MINLONGITUDE = base_graph.getLon(k) - epsilon;
            }
            if (base_graph.getLat(k) > MAXLATITUDE) {
                MAXLATITUDE = base_graph.getLat(k) + epsilon;
            }
            if (base_graph.getLat(k) < MINLATITUDE) {
                MINLATITUDE = base_graph.getLat(k) - epsilon;
            }
        }
    }

    bool indexInRange(int i, int d) {
        return (i + d >= 0 && i + d < gridsidesize);
    }

    double geoDistComparison(double lon1, double lat1, double lon2, double lat2) {
        double latmed, rdlon, rdlat;
        latmed = (lat1 - lat2) / 2;
        rdlon = (lon1 - lon2) * M_PI / 180;
        rdlat = (lat1 - lat2) * M_PI / 180;
        return R * sqrt((pow(rdlat, 2) + pow(cos(latmed) * rdlon, 2)));
        //return R * (pow(rdlat, 2) + pow(cos(latmed) * rdlon, 2)); //Vergleich ist ok ohne sqrt
    }
    
    double geoDistComparison(NodeID nodeID1, NodeID nodeID2) {
        double lon1 = base_graph.getLon(nodeID1);
        double lat1 = base_graph.getLat(nodeID1);               
        double lon2 = base_graph.getLon(nodeID2);
        double lat2 = base_graph.getLat(nodeID2);
        return geoDistComparison(lon1, lat1, lon2, lat2);
    }
    
    void clearGrid() {
        
        for (int x1 = 0; x1 < gridsidesize; x1++) {
            for (int x2 = 0; x2 < gridsidesize; x2++) {
                for (int y1 = 0; y1 < gridsidesize; y1++) {
                    for (int y2 = 0; y2 < gridsidesize; y2++) {
                        Grid[x1][x2][y1][y2].clear();
                    }

                }
            }
        }
        
        /*
        Grid.resize(gridsidesize);
        Grid.assign (Grid.begin(), Grid.end(), vector<vector<vector<ItList>>>(gridsidesize, vector<vector<ItList>>(gridsidesize, vector<ItList>(gridsidesize, ItList()))));        
         * */        
    }

public:

    FourDGrid(const int size, const GraphT &base_graph) :    
        gridsidesize(size), base_graph(base_graph),
        Grid(size, vector<vector<vector<GridCell>>>(size, vector<vector<GridCell>>(size, vector<GridCell>(size, GridCell())))),
                testvector(size, 0)
        //,NextChainInSameCell(base_graph.getNrOfNodes(), -1)
    {          
        testvector.resize(size *50);
        testvector.clear();
        MINLONGITUDE = numeric_limits<double>::max();
        MAXLONGITUDE = 0;
        MINLATITUDE = numeric_limits<double>::max();
        MAXLATITUDE = 0;
        calculateBorders();
        gridwidth = MAXLONGITUDE - MINLONGITUDE;
        gridheight = MAXLATITUDE - MINLATITUDE;
        assert(gridsidesize > 0);
        cellsizex = gridwidth / gridsidesize;
        cellsizey = gridheight / gridsidesize;  
        
        /*
        for (int x1 = 0; x1 < gridsidesize; x1++) {
            for (int x2 = 0; x2 < gridsidesize; x2++) {
                for (int y1 = 0; y1 < gridsidesize; y1++) {
                    for (int y2 = 0; y2 < gridsidesize; y2++) {
                        Print("cellsize: " << Grid.at(x1).at(x2).at(y1).at(y2).size());
                    }

                }
            }
        }*/
        assert(cellsizex != 0);
        assert(cellsizey != 0);
    }
        
    
    void identifyPairs(Chains_and_Remainder &CaR) {        
        for (int i = 1; i < 11; ++i) { //range of streettypes where two laning can occur
        //int i = 0;
        //for (ChainsOfType &chainsOfType : CaR.chainsAccordingToType) {
            
            ChainsOfType &chainsOfType = CaR.oneWayChainsAccordingToType.at(i);
            Print("Number of chainsOfType: " << chainsOfType.size());
            //int j=i;
            //i++;
            //if (1 <= j && j < 11) {
            if (chainsOfType.size()>0) {
            
                
                loadChains(chainsOfType);            
                //find opposites for all points in grid
                for (int x1 = 0; x1 < gridsidesize; x1++) {
                    for (int x2 = 0; x2 < gridsidesize; x2++) {
                        for (int y1 = 0; y1 < gridsidesize; y1++) {
                            for (int y2 = 0; y2 < gridsidesize; y2++) {
                                
                                GridCell &gridCell = Grid[x1][x2][y1][y2];
                                for (Gridpoint &gridpoint : gridCell) {
                                    gridpoint.partnerGridpoint = findOpposite(*(gridpoint.chain), chainsOfType);
                                }
                                
                            }
                        }
                    }
                }

                //collect mutual partners
                for (int x1 = 0; x1 < gridsidesize; x1++) {
                    for (int x2 = 0; x2 < gridsidesize; x2++) {
                        for (int y1 = 0; y1 < gridsidesize; y1++) {
                            for (int y2 = 0; y2 < gridsidesize; y2++) {
                                
                                for (Gridpoint &gridpoint : Grid[x1][x2][y1][y2]) {                                    
                                    if (gridpoint.isMatched == false && gridpoint.partnerGridpoint != nullptr) {
                                        if (gridpoint.partnerGridpoint->isMatched == false && gridpoint.partnerGridpoint->partnerGridpoint != nullptr) {
                                            if (gridpoint.partnerGridpoint->chain != gridpoint.chain) { //may not be its own partner
                                                if (gridpoint.partnerGridpoint->partnerGridpoint->chain == gridpoint.chain) {
                                                    ChainPair chainpair;
                                                    chainpair.chainTo.splice(chainpair.chainTo.end(), *gridpoint.chain);
                                                    chainpair.chainFrom.splice(chainpair.chainFrom.end(), *gridpoint.partnerGridpoint->chain);
                                                    CaR.chainPairs.push_back(chainpair);
                                                    chainsOfType.erase(gridpoint.chain);
                                                    chainsOfType.erase(gridpoint.partnerGridpoint->chain);
                                                    gridpoint.isMatched = true;
                                                    gridpoint.partnerGridpoint->isMatched = true;
                                                }
                                            }
                                        }
                                    }                                
                                }
                                
                            }
                        }
                    }
                }                        
            }
        }          
    }
         
    
    void loadChains (ChainsOfType &chains) {        
        clearGrid();                
        //chains in Grid einordnen
        //chainsptr = &chains;
        int x1, x2, y1, y2;        
        for (auto it = chains.begin(); it != chains.end(); ++it) {
            const NodeID nodeIdFront = it->front();
            const NodeID nodeIdBack = it->back();
            x1 = (int) ((base_graph.getLon(nodeIdFront) - MINLONGITUDE) / cellsizex);
            y1 = (int) ((base_graph.getLat(nodeIdFront) - MINLATITUDE) / cellsizey);
            x2 = (int) ((base_graph.getLon(nodeIdBack) - MINLONGITUDE) / cellsizex);
            y2 = (int) ((base_graph.getLat(nodeIdBack) - MINLATITUDE) / cellsizey);

            Grid[x1][y1][x2][y2].push_back(Gridpoint(it, nullptr, false));
        }
    }

    Gridpoint* findOpposite(Chain &chainToMatch, ChainsOfType &chains) {
        int x1, x2, y1, y2;
        double dist1, dist2, dist;

        NodeID nodeIdFront = chainToMatch.front();
        NodeID nodeIdBack = chainToMatch.back();
        x1 = (int) ((base_graph.getLon(nodeIdFront) - MINLONGITUDE) / cellsizex);
        y1 = (int) ((base_graph.getLat(nodeIdFront) - MINLATITUDE) / cellsizey);
        x2 = (int) ((base_graph.getLon(nodeIdBack) - MINLONGITUDE) / cellsizex);
        y2 = (int) ((base_graph.getLat(nodeIdBack) - MINLATITUDE) / cellsizey);



        //Suche in 3*3*3*3 Feld        
        //ChainsOfType::iterator closestChain = chainsptr->end();
        //Gridpoint closestChain(chains.end(), nullptr, false);
        Gridpoint* closestChain = nullptr;
        //double waylength = geoDistComparison(nodeIdFront, nodeIdBack);                
        //double closestDistance = numeric_limits<double>::max();
        double closestDistance = 20000;//200; //min(waylength / 10, 1.0);

        for (int dx1 = -1; dx1 <= 1; dx1++) {
            if (indexInRange(x1, dx1)) {
                for (int dx2 = -1; dx2 <= 1; dx2++) {
                    if (indexInRange(x2, dx2)) {
                        for (int dy1 = -1; dy1 <= 1; dy1++) {
                            if (indexInRange(y1, dy1)) {
                                for (int dy2 = -1; dy2 <= 1; dy2++) {
                                    if (indexInRange(y2, dy2)) {

                                        for (Gridpoint &current : Grid[x1 + dx1][y1 + dy1][x2 + dx2][y2 + dy2]) {
                                            NodeID currentNodeIdFront = current.chain->front();
                                            NodeID currentNodeIdBack = current.chain->back();
                                            dist1 = geoDistComparison(nodeIdFront, currentNodeIdBack);
                                            dist2 = geoDistComparison(nodeIdBack, currentNodeIdFront);
                                            dist = pythagoras(dist1, dist2);
                                            if (dist < closestDistance) {
                                                closestDistance = dist;
                                                closestChain = &current;
                                            }
                                            
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return closestChain;

    }

};
