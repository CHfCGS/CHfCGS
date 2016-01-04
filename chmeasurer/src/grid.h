/* 
 * File:   grid.h
 * Author: tobias
 *
 * Created on 28. September 2015, 10:10
 */

#pragma once

#include <vector>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <array>

#include "nodes_and_edges.h"
#include "geoFunctions.h"
#include "window.h"


using namespace std;



template <class GraphT>
class Grid {
    typedef vector<vector<int> > twoDvector;
    //typedef vector<vector<twoDvector> > fourDvector;
    
private:
    //const std::vector<NodeT> &nodes;
    //constexpr auto M_PI = 3.14159265358979323846;
    double R = 6371.009; //Erdradius
    
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
    const GraphT &graph;

    twoDvector FirstNodeGrid;
    vector<int> NextNodeInSameCell;
    //vector<int> gridtest;

    double pythagoras(double a, double b) const {
        return sqrt(pow(a, 2) + pow(b, 2));
    }

    void calculateBorders() {
        for (u_int32_t k = 0; k < graph.getNrOfNodes(); k++) {
            if (graph.getLon(k) > MAXLONGITUDE) {
                MAXLONGITUDE = graph.getLon(k);
            }
            if (graph.getLon(k) < MINLONGITUDE) {
                MINLONGITUDE = graph.getLon(k);
            }
            if (graph.getLat(k) > MAXLATITUDE) {
                MAXLATITUDE = graph.getLat(k);
            }
            if (graph.getLat(k) < MINLATITUDE) {
                MINLATITUDE = graph.getLat(k);
            }
        }
    }
    
    int xCoordinate (double lon) const {
        int x = (int) ((lon - MINLONGITUDE) / cellsizex);
        if (!indexInRange(x)) {
            if (x == -1) {
                x = 0;
            } else if (x == gridsidesize) {
                x = gridsidesize - 1;
            } else {
                assert(false);
            }
        }
        return x;
    }
    int yCoordinate (double lat) const {
        int y = (int) ((lat - MINLATITUDE) / cellsizey);
        if (!indexInRange(y)) {
            if (y == -1) {
                y = 0;
            } else if (y == gridsidesize) {
                y = gridsidesize - 1;
            } else {
                assert(false);
            }
        }
        return y;
    }
    
    bool indexInRange(const int index) const {
        return (0 <= index && index < gridsidesize);
    }
    
    bool indexInRange(const int i, const int d) const {
        return (i + d >= 0 && i + d < gridsidesize);
    }
    /*
    double geoDistComparison(double lon1, double lat1, double lon2, double lat2) {
        double latmed, rdlon, rdlat;
        latmed = (lat1 - lat2) / 2;
        rdlon = (lon1 - lon2) * M_PI / 180;
        rdlat = (lat1 - lat2) * M_PI / 180;
        return R * sqrt((pow(rdlat, 2) + pow(cos(latmed) * rdlon, 2)));
        //return R * (pow(rdlat, 2) + pow(cos(latmed) * rdlon, 2)); //Vergleich ist ok ohne sqrt
    }
     * */

public:

    Grid(const int size, const GraphT &base_graph) : gridsidesize(size), graph(base_graph) {
    //FirstNodeGrid(size, vector<int>(size, -1)),
    //NextNodeInSameCell(base_graph.getNrOfNodes(), -1) {
        //vector<int> assignVector(size, -1);
        
        //FirstNodeGrid.reserve(size);
        FirstNodeGrid.assign(size, vector<int>(size, -1));
        NextNodeInSameCell.assign(base_graph.getNrOfNodes(), -1);
        MINLONGITUDE = numeric_limits<double>::max();
        MAXLONGITUDE = 0;
        MINLATITUDE = numeric_limits<double>::max();
        MAXLATITUDE = 0;
        calculateBorders();
        gridwidth = MAXLONGITUDE - MINLONGITUDE;
        gridheight = MAXLATITUDE - MINLATITUDE;
        cellsizex = gridwidth / gridsidesize;
        cellsizey = gridheight / gridsidesize;

        //wird nur für Konstruktion gebraucht
        twoDvector LastNodeGrid(size, vector<int>(size, -1));
        
        //Knoten in Grid einordnen
        for (u_int32_t k = 0; k < base_graph.getNrOfNodes(); k++) {
            //int x = (int) ((base_graph.getLon(k) - MINLONGITUDE) / cellsizex);
            //int y = (int) ((base_graph.getLat(k) - MINLATITUDE) / cellsizey);
            uint x = xCoordinate(base_graph.getLon(k));
            uint y = yCoordinate(base_graph.getLat(k));

            //Zelle bisher leer? Aufbau einer Liste
            if (FirstNodeGrid[x][y] == -1) {
                FirstNodeGrid[x][y] = k;

            } else {
                NextNodeInSameCell[LastNodeGrid[x][y]] = k;
            }
            LastNodeGrid[x][y] = k;

        }
        
    }
       
    std::list<NodeID> nodeGeoNeighbours(chm::NodeID node_id) const {
        //std::vector<chc::NodeID> neighbours;
        //int x = (int) ((base_graph.getLon(node_id) - MINLONGITUDE) / cellsizex);
        //int y = (int) ((base_graph.getLat(node_id) - MINLATITUDE) / cellsizey);
        
        /*
        for (int dx = -1; dx <= 1; dx++) {
            if (indexInRange(x, dx)) {
                for (int dy = -1; dy <= 1; dy++) {
                    if (indexInRange(y, dy)) {
                        int currentNode_id = FirstNodeGrid[x + dx][y + dy];
                        while (currentNode_id != -1) {
                            neighbours.push_back(currentNode_id);
                            currentNode_id = NextNodeInSameCell[currentNode_id];
                            //numberOfComparedWays++;
                        }
                    }
                    
                }
            }
        }
         * */
        return nodesInNeigbourhood(graph.getLat(node_id), graph.getLon(node_id));
    }
    
    std::list<NodeID> nodesInNeigbourhood(double lat, double lon) const {
        std::list<chm::NodeID> neighbours;
        
        //int x = (int) ((lon - MINLONGITUDE) / cellsizex);
        //int y = (int) ((lat - MINLATITUDE) / cellsizey);
        int x = xCoordinate(lon);
        int y = yCoordinate(lat);
        
        for (int dx = -1; dx <= 1; dx++) {
            if (indexInRange(x, dx)) {
                for (int dy = -1; dy <= 1; dy++) {
                    if (indexInRange(y, dy)) {
                        int currentNode_id = FirstNodeGrid[x + dx][y + dy];
                        while (currentNode_id != -1) {
                            neighbours.push_back(currentNode_id);
                            currentNode_id = NextNodeInSameCell[currentNode_id];                            
                        }
                    }                    
                }
            }
        }
        
        return neighbours;
    }
    
    std::list<NodeID> nodesInWindow(Window window) const {
        std::list<chm::NodeID> windowNodes;
        //get cell coordinates of the window corners
        
        int min_x = xCoordinate(window.MINLONGITUDE); //(int) ((window.MINLONGITUDE - MINLONGITUDE) / cellsizex);
        int min_y = yCoordinate(window.MINLATITUDE); //(int) ((window.MINLATITUDE - MINLATITUDE) / cellsizey);
        int max_x = xCoordinate(window.MAXLONGITUDE); //(int) ((window.MAXLONGITUDE - MINLONGITUDE) / cellsizex);
        int max_y = yCoordinate(window.MAXLATITUDE); //(int) ((window.MAXLATITUDE - MINLATITUDE) / cellsizey);
        assert(min_x <= max_x);
        assert(min_y <= max_y);
        for (int x = min_x; x <= max_x; x++) {
            assert(indexInRange(x, 0));
                for (int y = min_y; y <= max_y; y++) {
                    assert(indexInRange(y, 0));
                        int currentNode_id = FirstNodeGrid[x][y];
                        while (currentNode_id != -1) {                            
                            if (window.isIn(graph, currentNode_id)) {
                                if(graph.isVisibleNode(currentNode_id)) {
                                    windowNodes.push_back(currentNode_id);   
                                }                                
                                //option: let only visible nodes count
                                /*
                                if (graph.isVisibleNode(currentNode_id))  {
                                    
                                }*/
                            }
                            
                            currentNode_id = NextNodeInSameCell[currentNode_id];                            
                        }
                                        
                }
            
        }
        return windowNodes;
    }

};

