#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <limits>

#include "../geoFunctions.h"
#include "../chains.h"
#include "../chgraph.h"
#include "../nodes_and_edges.h"
#include "prio_nodes.h"
#include "match_functions.h"

namespace chc {
    namespace unit_tests {
        void testCDTMatching();         
    }
}

namespace ls {

namespace mc {
    
struct VertexInfo2 {
    //std::list<PrioNodeHandle>::iterator it;
    PrioNodeHandle pnh;
    //uint posInChain;
    bool side;
    //Chain::iterator node_it;    
    //std::reference_wrapper<Node> node;
    //VertexInfo2() {};
};

struct FaceInfo2
{
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){ 
    return nesting_level%2 == 1;
  }
};





template <class GraphT>
class CDTMatching2 {
private:
    
    typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
    typedef CGAL::Cartesian<NT> K;
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel                       K;
    //typedef CGAL::Simple_cartesian<double> K;
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo2, K>                      Vb;
    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>    Fbb;
    typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
    typedef CGAL::Exact_intersections_tag                                Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
    typedef CDT::Point                                                Point;
    //typedef CGAL::Polygon_2<K>                                        Polygon_2;
    //typedef CGAL::Polygon_2<K>                                        Polygon_2;
    //typedef K::Point_2                                 Point_2;
    
    const GraphT &graph;
    
    void mark_domains(CDT& ct,
            CDT::Face_handle start,
            int index,
            std::list<CDT::Edge>& border) {
        if (start->info().nesting_level != -1) {
            return;
        }
        std::list<CDT::Face_handle> queue;
        queue.push_back(start);
        while (!queue.empty()) {
            CDT::Face_handle fh = queue.front();
            queue.pop_front();
            if (fh->info().nesting_level == -1) {
                fh->info().nesting_level = index;
                for (int i = 0; i < 3; i++) {
                    CDT::Edge e(fh, i);
                    CDT::Face_handle n = fh->neighbor(i);
                    if (n->info().nesting_level == -1) {
                        if (ct.is_constrained(e)) border.push_back(e);
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

    
    void mark_domains(CDT& cdt) {
        for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) {
            it->info().nesting_level = -1;
        }
        std::list<CDT::Edge> border;
        mark_domains(cdt, cdt.infinite_face(), 0, border);
        while (!border.empty()) {
            CDT::Edge e = border.front();
            border.pop_front();
            CDT::Face_handle n = e.first->neighbor(e.second);
            if (n->info().nesting_level == -1) {
                mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
            }
        }
    }

    bool belongsToInner(CDT &cdt, CDT::Edge_circulator ec) {
        CDT::Face_handle fh1 = ec->first;
        CDT::Face_handle fh2 = fh1->neighbor(ec->second);
        //at least one incident facet must lie in the inner area,
        //that is because constraints between start and end are also valid
        return (fh1->info().in_domain() || fh2->info().in_domain());
    }

    CDT::Vertex_handle getOtherNode(CDT &cdt, CDT::Vertex_handle vc, CDT::Edge_circulator ec) {

        CDT::Face_handle fh1 = ec->first;
        CDT::Face_handle fh2 = fh1->neighbor(ec->second);
        std::vector<CDT::Vertex_handle> fh1vertices;
        std::vector<CDT::Vertex_handle> fh2vertices;
        for (int i = 0; i < 3; i++) {
            fh1vertices.push_back(fh1->vertex(i));
            fh2vertices.push_back(fh2->vertex(i));
        }
        //find doubled element unequal to vc
        for (int i = 0; i < 3; i++) {
            CDT::Vertex_handle candidate = fh1vertices[i];
            //we dont want our src vc
            if (candidate == vc) continue;
            //the tgt must also be in the other facet
            for (int j = 0; j < 3; j++) {
                CDT::Vertex_handle candidate2 = fh2vertices[j];
                if (candidate == candidate2) return candidate;
            }
        }
        //should never come until here
        assert(false);

        CDT::Vertex_handle dummyNode;
        return dummyNode;
    }
    
        
    bool arePossiblePartners(VertexInfo2 vi1, VertexInfo2 vi2) {
        //possible partners must lie on different sides
        if (vi1.side != vi2.side) {
            /*
            VertexInfo2 viTo, viFrom;
            if(vi1.side == true) {
                viTo = vi1;
                viFrom = vi2;
            } else {
                viTo = vi2;
                viFrom = vi1;
            }
            //it could happen that the triangulation connects start and end
            if (!((viTo.posInChain == 0 && viFrom.posInChain == size2)
                    || (viTo.posInChain == size1 && viFrom.posInChain == 0))) {
                return true;
            }
             * */
            return true;
        }
        return false;
    }
    /*
    void _setFollower (PrioNodeHandle node_h,
                      PrioNodeHandle follower_h) {
        PrioNode2 &node = *node_h;
        PrioNode2 &follower = *follower_h;
        node.followerValid = true;
        node.follower_h = follower_h;            
        
        const PrioNodeHandle pnh = node_h;
        std::list<PrioNodeHandle> &guidelist = follower.guides;
        guidelist.push_back(pnh);
    }    
    
    double geoDist(const PrioNodeHandle pnh1, const PrioNodeHandle pnh2, const GraphT &graph) {
        const PrioNode2 &pn1 = *pnh1;
        const PrioNode2 &pn2 = *pnh2;
        return geo::geoDist(graph.getNode(pn1.node_id), graph.getNode(pn2.node_id));
    }
    
    void getAndSetNearestFollower(PrioNodeHandle pnh,
                std::list<PrioNodeHandle> possibleFollowers,
                const GraphT &graph) {
        auto nearestFollowerIt = possibleFollowers.end();
        double minDist = std::numeric_limits<double>::max();
        for (auto it = possibleFollowers.begin(); it != possibleFollowers.end(); it++) {
            double dist = geoDist(pnh, *it, graph);
            if (dist < minDist) {
                minDist = dist;
                nearestFollowerIt = it;
            }            
        }
        if(nearestFollowerIt != possibleFollowers.end()) {
            _setFollower(pnh, *nearestFollowerIt);
        }

    }
    */
    
    void insertChainConstraints(const std::list<PrioNodeHandle> &list,
                                    std::list<CDT::Vertex_handle> &vertexHandles,
                                    bool side,
                                    CDT &cdt) {
        //uint posInChain = 0;
        //insert points in cdt
        for (auto it = list.begin(); it != list.end(); it++) {
            const PrioNode2 &pn = *(*it);
            NodeID node_id = pn.node_id;
            auto node = graph.getNode(node_id);
            CDT::Vertex_handle vh = cdt.push_back(Point(node.lon, node.lat));            
            vh->info().pnh = *it;
            //vh->info().posInChain = posInChain;
            vh->info().side = side;
            vertexHandles.push_back(vh);
            //posInChain++;
        }
        //insert constraints
        for (auto it = vertexHandles.begin(); it != --vertexHandles.end(); it++) {
            auto next = it;
            next++;            
            cdt.insert_constraint(*it, *next);
        }
    }
    
    void insertConstraints(const std::list<PrioNodeHandle> &toList,
                            const std::list<PrioNodeHandle> &fromList,
                            CDT &cdt) {
        
        std::list<CDT::Vertex_handle> vertexHandlesTo;
        insertChainConstraints(toList, vertexHandlesTo, true, cdt);        
        
        std::list<CDT::Vertex_handle> vertexHandlesFrom;
        insertChainConstraints(fromList, vertexHandlesFrom, false, cdt);        
        
        //set constraints at begin and end
        //cdt.insert_constraint(vertexHandlesTo.front(), vertexHandlesFrom.front());
        
        cdt.insert_constraint(vertexHandlesTo.front(), vertexHandlesFrom.front());
        cdt.insert_constraint(vertexHandlesTo.back(), vertexHandlesFrom.back());
    }
    
public:
    CDTMatching2(const GraphT &graph) : graph(graph) {        
    }

    void match(const std::list<PrioNodeHandle> &toList,
                const std::list<PrioNodeHandle> &fromList) {
        assert(toList.size() >= 1 && fromList.size() >= 1);
        CDT cdt;                

        /*
        std::list<CDT::Vertex_handle> vertexHandlesTo;
        uint posInChain = 0;
        for (auto it = toList.begin(); it != toList.end(); it++) {
            NodeID node_id = (*it)->node_id;
            auto node = graph.getNode(node_id);
            CDT::Vertex_handle vh = cdt.push_back(Point(node.lon, node.lat));            
            vh->info().it = it;
            vh->info().posInChain = posInChain;
            vh->info().side = true;
            vertexHandlesTo.push_back(vh);
            posInChain++;
        }
        for (auto it = vertexHandlesTo.begin(); it != --vertexHandlesTo.end(); it++) {
            auto next = it;
            next++;
            cdt.insert_constraint(*it, *next);
        }
        
        std::list<CDT::Vertex_handle> vertexHandlesFrom;
        posInChain = 0;
        for (auto it = toList.begin(); it != toList.end(); it++) {
            NodeID node_id = (*it)->node_id;
            auto node = graph.getNode(node_id);
            CDT::Vertex_handle vh = cdt.push_back(Point(node.lon, node.lat));             
            vh->info().it = it;
            vh->info().posInChain = posInChain;
            vh->info().side = false;
            vertexHandlesFrom.push_back(vh);
            posInChain++;
        }
        for (auto it = vertexHandlesFrom.begin(); it != --vertexHandlesFrom.end(); it++) {
            auto next = it;
            next++;
            cdt.insert_constraint(*it, *next);
        }
        
        //set constraints at begin and end
        //cdt.insert_constraint(vertexHandlesTo.front(), vertexHandlesFrom.front());
        cdt.insert_constraint(vertexHandlesTo.front(), vertexHandlesFrom.back());
        cdt.insert_constraint(vertexHandlesTo.back(), vertexHandlesFrom.front());
        */
        
        insertConstraints(toList, fromList, cdt);
        
        mark_domains(cdt);

        uint finiteVerticesCount = 0;
        for (CDT::Finite_vertices_iterator it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); it++) {
            PrioNodeHandle pnh = it->info().pnh;
            PrioNode2 pn = *pnh;
            finiteVerticesCount++;
        }
        
        assert(finiteVerticesCount == toList.size()+ fromList.size());
        

        for (CDT::Finite_vertices_iterator srcIt = cdt.finite_vertices_begin(); srcIt != cdt.finite_vertices_end(); srcIt++) {
            /*
            Point_2 point_2 = it->point();
            CDT::Vertex_handle vh = it;
            auto info = vh->info();
            std::cout << it->info().node_it->node_id << " " << point_2.x() << " " << point_2.y() << std::endl;
            std::cout << std::endl;
            //it->incident_vertices()
            */
            /*
            std::list<std::list<std::list<PPrioNode>::iterator>::iterator> possibleFollowersIts;
            CDT::Vertex_circulator c = it->incident_vertices();
            auto d = c;
            do {
                //Point_2 point_2 = c->point();
                CDT::Vertex_handle vhn = c->handle();

                if (!cdt.is_infinite(vhn)) {
                    if (arePossiblePartners(it->info(), c->info(), toList.size(), fromList.size())) {
                        possibleFollowersIts.push_back(c->info().it);
                    }
                    
                    //std::cout << c->info().node_it->node_id << " " << point_2.x() << " " << point_2.y() << std::endl;
                }
                
            } while (++c != d);
            getAndSetNearestFollower(it->info().it, possibleFollowersIts);
            */
            std::list<PrioNodeHandle> possibleFollowers;
            
            CDT::Vertex_handle src_vh = srcIt->handle();
            CDT::Edge_circulator c = src_vh->incident_edges();
            auto d = c;
            do {
                if (!cdt.is_infinite(c)) {
                    //Point_2 point_2 = c->point();  
                    //CDT::Vertex_handle start_handle = c->handle();
                    if (belongsToInner(cdt, c)) {
                        CDT::Vertex_handle tgt_vh = getOtherNode(cdt, src_vh, c);
                        //std::cout << vh->info().node_it->node_id << " " << vh->point().x() << " " << vh->point().y() << std::endl;
                        
                        if (arePossiblePartners(src_vh->info(), tgt_vh->info())) {
                            possibleFollowers.push_back(tgt_vh->info().pnh);
                        }
                    }
                }
            } while (++c != d);
                        
            getAndSetNearestFollower(src_vh->info().pnh, possibleFollowers, graph);            
            //for (CDT::Finite_vertices_iterator it2 = it->; it != cdt.finite_vertices_end(); it++) {
            //std::cout << std::endl;
        }
        
    }    
};

}

}

