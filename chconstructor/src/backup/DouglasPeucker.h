/* 
 * File:   DouglasPeucker.h
 * Author: tobias
 *
 * Created on 4. September 2015, 14:14
 */

#include "nodes_and_edges.h"
#include "grid.h"

#include <limits>
#include <math.h>

#ifndef DOUGLASPEUCKER_H
#define	DOUGLASPEUCKER_H

namespace DP {


template <class GraphT>
class DouglasPeucker{
    private:       
        PrioList nodePrios;
        std::list<Intervall> intervalls;
        //vector<Node> &nodes;
        const GraphT &base_graph;
        Grid<GraphT> grid;

        //naiv TODO: heap implementation;

        PrioNode extractMax(PrioList &nodePrios) {

            
            PrioList::iterator maxIt;
            double max = std::numeric_limits<double>::lowest();

            for (PrioList::iterator it = nodePrios.begin(); it != nodePrios.end(); it++) {
                if (it->prio > max) {
                    maxIt = it;
                    max = it->prio;
                }
            }
            PrioNode r(*maxIt);
            nodePrios.erase(maxIt);
            return r;
        }
    
    
                
    PrioList split() {        
        if (nodePrios.empty()) {
            return PrioList(); //empty list
        }
        else {
            
            PrioNode max = extractMax(nodePrios);
            
            
            //cout << "extracting Node:";
            //cout << nodes[max.node_id].coord.lon << " " <<  nodes[max.node_id].coord.lat << "\n";
            
            std::list<Intervall>::iterator splitIt = max.intervallIt;            
            std::list<PrioNodeRef>::iterator splitPos = max.posInIntervall;                                                                            
            
            Intervall leftIntervall(splitIt->start, splitPos->get().node_id);
            Intervall rightIntervall(splitPos->get().node_id, splitIt->end);                        
            
            std::list<Intervall>::iterator leftIt = intervalls.insert(splitIt, leftIntervall);           
            std::list<Intervall>::iterator rightIt = intervalls.insert(splitIt, rightIntervall);                                                                                                
            
            //left Intervall
            leftIt->prioNodes.splice(leftIt->prioNodes.begin(), splitIt->prioNodes, splitIt->prioNodes.begin(), splitPos);
            //leftIt->prioNodes = std::list<PrioNodeRef>(splitIt->prioNodes.begin(), splitPos);
            //cout << "left: \n";            
            for (std::list<PrioNodeRef>::iterator it = leftIt->prioNodes.begin(); it!= leftIt->prioNodes.end(); it++) {
                //cout << nodes[it->get().node_id].coord.lon << " " << nodes[it->get().node_id].coord.lat << "\n";
                it->get().intervallIt = leftIt;
                it->get().posInIntervall = it;                
            }

            //right Intervall  
            rightIt->prioNodes.splice(rightIt->prioNodes.begin(), splitIt->prioNodes, ++splitPos, splitIt->prioNodes.end());
            //rightIt->prioNodes = std::list<PrioNodeRef>(++splitPos, splitIt->prioNodes.end());
            //cout << "right: \n";            
            for (std::list<PrioNodeRef>::iterator it = rightIt->prioNodes.begin(); it!= rightIt->prioNodes.end(); it++) {                
                //cout << nodes[it->get().node_id].coord.lon << " " << nodes[it->get().node_id].coord.lat << "\n";
                it->get().intervallIt = rightIt;
                it->get().posInIntervall = it;                 
            }                                    
            
            intervalls.erase(splitIt);
                        
            update (*leftIt);
            update (*rightIt);
            
            //cout << "\n";
            
            PrioList r = split();//(max.intervallIt); //, max.posInIntervall);
            r.push_back(max);
            return r;
        }           
    }
    
    void update(Intervall &intervall){
        
        struct vec{
            double x;
            double y;            
        };
        
        //as if the earth would be flat, i.e lat and lon orthogonal, fails at the poles
        double x1 = base_graph.getLon(intervall.start); //prioNodes.front().get().node.coord.lon;
        double y1 = base_graph.getLat(intervall.start); //prioNodes.front().get().node.coord.lat;
        double x2 = base_graph.getLon(intervall.end); //prioNodes.back().get().node.coord.lon;
        double y2 = base_graph.getLat(intervall.end); //prioNodes.back().get().node.coord.lat;
        
        vec parallelVec = {x2-x1, y2-y1};
        //vec orthogonalVec = {parallelVec.y, -parallelVec.x};
        
        
        /*
        //line of the intervall endpoints
        double gradient_Intervall = (y2-y1)/(x2-x1);
        double c_Intervall = y1-gradient_Intervall*x1;
        */
        for (PrioNode &prioNode: intervall.prioNodes) {
          
            //outlier
            double x_o = base_graph.getLon(prioNode.node_id);
            double y_o = base_graph.getLat(prioNode.node_id);
            
            double dist_x;
            double dist_y;
            if (parallelVec.x == 0 && parallelVec.y ==0) { //exception case: divide by zero
                dist_x = x_o - x1;
                dist_y = y_o - y1;
            }else{ //regular case
                //solve linear system for intersection point
                //i.e. (x1, y1) + c * parallelVec = (x_o, y_o) + d * orthogonalVec        
                double c= -((y1-y_o)*parallelVec.y + (x1-x_o)*parallelVec.x)/(pow(parallelVec.x, 2)+pow(parallelVec.y, 2));
                vec intersection = {x1+c*parallelVec.x, y1+c*parallelVec.y};
                dist_x = x_o - intersection.x;
                dist_y = y_o - intersection.y;
            }
                                    
            //calculate distance between outlier and intersection point                    
            double dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
            
            prioNode.prio = dist;            
        }
        return;
    }
    
    public:
        DouglasPeucker(const GraphT &base_graph, const Grid<GraphT> &grid):
        nodePrios(), intervalls(), base_graph(base_graph), grid(grid) {}
        
        PrioList process(const Chain &chain) {                

            assert(chain.node_ids.size() >= 3);

            nodePrios.clear();
            intervalls.clear();
                    
            Intervall initialIntervall(chain.node_ids.front(), chain.node_ids.back());
            intervalls.emplace_back(initialIntervall);
            const auto intervallIt = intervalls.begin();

            for (auto it = ++chain.node_ids.begin(); it != --chain.node_ids.end(); it++) {
                nodePrios.emplace_back(PrioNode(*it, intervallIt));

                //reference each other
                intervallIt->prioNodes.emplace_back(nodePrios.back());
                std::list<PrioNodeRef>::iterator pos = --((intervallIt->prioNodes).end());
                nodePrios.back().posInIntervall = pos;


            }
            update(intervalls.front());     
            
            return split();
        }                     
};

}

#endif	/* DOUGLASPEUCKER_H */

