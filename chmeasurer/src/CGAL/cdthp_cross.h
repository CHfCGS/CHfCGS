#pragma once

#include "../window.h"
#include "../nodes_and_edges.h"
#include "../chgraph.h"
#include "range_tree.h"
#include "../grid.h"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>



namespace chm {

class CDTHPCross { 

    
    struct Node {
        NodeID node_id;
        double lat;
        double lon;
        Node (NodeID node_id, double lat, double lon): node_id(node_id), lat(lat), lon(lon) {}    
    };



    struct VertexInfo2 {
        //Chain::iterator node_it;
        uint posInChain;
        bool side;

        //std::reference_wrapper<Node> node;
        VertexInfo2() {};
    };
    
    struct FaceInfo2
    {
      FaceInfo2(){}
      int nesting_level;
      bool in_inner() {
          return nesting_level >=1;
      }
      /*
      bool in_domain(){ 
        return nesting_level%2 == 1;
      }*/
    };
    typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel                       K;
    typedef CGAL::Cartesian<NT> K;

    typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo2, K>                      Vbb;
    typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;

    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>    Fbb;
    typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;

    typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
    typedef CGAL::Exact_intersections_tag                                Itag;
    //typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
    typedef CGAL::Constrained_triangulation_2<K, TDS, Itag>  CDT;
    typedef CGAL::Triangulation_hierarchy_2<CDT> CDTH;
    typedef CGAL::Constrained_triangulation_plus_2<CDTH> CDTHP;
    typedef CDTHP::Point                                                Point;
    //typedef CGAL::Polygon_2<Vb>                                        Polygon_2;
    typedef CGAL::Polygon_2<K>                                        Polygon_2;
    typedef K::Point_2                                 Point_2;

    
    CHGraph<CHNode, CHEdge> &graph;
    Grid<CHGraph<CHNode, CHEdge> > &grid;
    
    RangeTree &rangeTree;
    
    
    void 
    mark_domains(CDTHP& ct, 
                 CDTHP::Face_handle start, 
                 int index, 
                 std::list<CDTHP::Edge>& border )
    {
      if(start->info().nesting_level != -1){
        return;
      }
      std::list<CDTHP::Face_handle> queue;
      queue.push_back(start);
      while(! queue.empty()){
        CDTHP::Face_handle fh = queue.front();
        queue.pop_front();
        if(fh->info().nesting_level == -1){
          fh->info().nesting_level = index;
          for(int i = 0; i < 3; i++){
            CDTHP::Edge e(fh,i);
            CDTHP::Face_handle n = fh->neighbor(i);
            if(n->info().nesting_level == -1){
              if(ct.is_constrained(e)) border.push_back(e);
              else queue.push_back(n);
            }
          }
        }
      }
    }
    //explore set of facets connected with non constrained edges,
    //and attribute to each such set a nesting level.
    //We start from facets incident to the infinite vertex, with a nesting
    //level of 0. Then we recursively consider the non-explored facets incident 
    //to constrained edges bounding the former set and increase the nesting level by 1.
    //Facets in the domain are those with an odd nesting level.
    void
    mark_domains(CDTHP& cdt)
    {
      for(CDTHP::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
        it->info().nesting_level = -1;
      }
      std::list<CDTHP::Edge> border;
      mark_domains(cdt, cdt.infinite_face(), 0, border);
      while(! border.empty()){
        CDTHP::Edge e = border.front();
        border.pop_front();
        CDTHP::Face_handle n = e.first->neighbor(e.second);
        if(n->info().nesting_level == -1){
          mark_domains(cdt, n, e.first->info().nesting_level+1, border);
        }
      }
    }
    
    //adds chain and its shortcuts as constraints
    void insertConstraints(CDTHP &cdthp, const Chain &chain) {
        std::list<CDTHP::Vertex_handle> vertexHandles;
    
        for (auto it = chain.begin(); it != chain.end(); it++) {
            CHNode node = graph.getNode(*it);
            CDTHP::Vertex_handle vh = cdthp.push_back(Point(node.lon, node.lat));
            vertexHandles.push_back(vh);
            //CDT::Vertex_handle vh = cdt.insert(Point(it->lon, it->lon));
            // --cdt.vertices_end();
            
            //vh->info().node_it = it;
            vh->info().posInChain = 5;
            vh->info().side = true;
        }    
        for (auto it = vertexHandles.begin(); it != --vertexHandles.end(); it++) {
            auto next = it;
            next++;
            cdthp.insert_constraint(*it, *next);
        }
        //shortcut
        cdthp.insert_constraint(vertexHandles.back(), vertexHandles.front());
    }
    
    uint countInFacets(const CDTHP &cdthp, const std::list<NodeID> &outliers) {
        uint counter = 0;
        
        for (NodeID node_id: outliers) {
            CHNode node = graph.getNode(node_id);
            CDTHP::Face_handle fh = cdthp.locate(Point_2(node.lon, node.lat));
            if (fh->info().in_inner()) {
                counter++;
            }
        }
        return counter;
    }


    void subtractLists(std::list<NodeID> &minuend, std::list<NodeID> subtrahend) {
        minuend.sort();
        minuend.unique();
        subtrahend.sort();
        subtrahend.unique();
        
        auto mIt = minuend.begin();
        auto sIt = subtrahend.begin();        
        while (mIt != minuend.end() && sIt != subtrahend.end()) {
            if (*mIt == *sIt) {
                mIt = minuend.erase(mIt);                
                sIt++;
            } else if (*mIt > *sIt) {
                sIt++;
            } else { //(*mIt < *sIt)
                mIt++;
            }
        }        
    }
    
    std::list<NodeID> getOutliers(const Chain &chain) {
        
        Window window(graph, chain);
        //std::list<NodeID> outliers = rangeTree.rectangleQuery(window);
        std::list<NodeID> outliers = grid.nodesInWindow(window);
        /*
        Print("rangeTree");
        outliers.sort();
        for (NodeID node_id: outliers) {
            Print("node_id" << node_id);
        }
        
        
        Print("grid");
        //outliers2.sort();
        for (NodeID node_id: outliers2) {
            Print("node_id" << node_id);
        }
        Print("chain");
        //outliers2.sort();
        for (NodeID node_id: chain) {
            Print("node_id" << node_id);
        }
         * */
        //assert(outliers.size() == outliers2.size());
        subtractLists(outliers, chain);
        
        return outliers;
    }
    
public:
    CDTHPCross(CHGraph<CHNode, CHEdge> &graph, Grid<CHGraph<CHNode, CHEdge> > &grid, RangeTree &rangeTree):
        graph(graph), grid(grid), rangeTree(rangeTree){                
    }
    
    uint getNofCrossings (const Chain &chain) {
        std::list<NodeID> outliers = getOutliers(chain);
        CDTHP cdthp;
        insertConstraints(cdthp, chain);
        return countInFacets(cdthp, outliers);        
    }
    
};

}