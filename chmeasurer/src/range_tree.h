#pragma once
#include <CGAL/Cartesian.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Range_tree_k.h>
#include <vector>
#include "nodes_and_edges.h"
#include "chgraph.h"
#include "window.h"

typedef CGAL::Cartesian<double> K;
typedef CGAL::Range_tree_map_traits_2<K, chm::NodeID> Traits;
typedef CGAL::Range_tree_2<Traits> Range_tree_2_type;
typedef Traits::Key Key;
typedef Traits::Interval Interval;


class RangeTree {
    Range_tree_2_type range_tree_2;
    //std::vector<Key> InputList;
public:
    RangeTree(chm::CHGraph<chm::CHNode,chm::CHEdge> &graph) {
        
        std::vector<Key> InputList;
        for (uint node_id = 0; node_id < graph.getNrOfNodes(); node_id++) {
            chm::CHNode node = graph.getNode(node_id);
            Key key = Key(K::Point_2(node.lat, node.lon), node_id);
            InputList.push_back(key);            
        }
        ;
        std::vector<Key> OutputList;
        range_tree_2.make_tree (InputList.begin(), InputList.end());
        //Range_tree_2_type testRange_tree_2(InputList.begin(), InputList.end());
        Interval win(Interval(K::Point_2(0, 3), K::Point_2(4, 5)));
        range_tree_2.window_query(win, std::back_inserter(OutputList));
        /*
        for (auto it = OutputList.begin(); it != OutputList.end(); it++) {
            std::cout << it->first.x() << "," << it->first.y()
                << ":" << it->second << std::endl;
        }*/
    }
    
    std::list<NodeID> rectangleQuery(Window window) {
        return rectangleQuery(window.MINLATITUDE, window.MINLONGITUDE, window.MAXLATITUDE, window.MAXLONGITUDE);
    }
    
    std::list<NodeID> rectangleQuery(double lat1, double lon1, double lat2, double lon2) {
        Interval win(Interval(K::Point_2(lat1, lon1), K::Point_2(lat2, lon2)));
        //Interval win(Interval(K::Point_2(0, 3), K::Point_2(4, 5)));
        std::vector<Key> OutputList;
        range_tree_2.window_query(win, std::back_inserter(OutputList));
        std::list<NodeID> toRet;
        for (auto it = OutputList.begin(); it != OutputList.end(); it++) {
            toRet.push_back(it->second);
        }
        return toRet;
        
    }
};
