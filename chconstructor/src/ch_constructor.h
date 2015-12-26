#pragma once

#include "defs.h"
#include "nodes_and_edges.h"
#include "graph.h"
#include "chgraph.h"
#include "prioritizer.h"
#include "geoFunctions.h"
#include "simplification/lineSimplifierType.h"
#include "edgeWeight.h"

#include <chrono>
#include <queue>
#include <mutex>
#include <vector>
#include <omp.h>
#include <algorithm>

namespace chc {

    namespace unit_tests {
        void testCHConstructor();        
    }

    namespace {
        uint MAX_UINT(std::numeric_limits<uint>::max());
    }
       
    template <typename NodeT, typename EdgeT>
    class CHConstructor {
    private:
        // typedef CHNode<NodeT> LvlNode;
        typedef CHEdge<EdgeT> Shortcut;
        typedef CHGraph<NodeT, EdgeT> CHGraphT;

        struct CompInOutProduct;
        struct PQElement;
        typedef std::priority_queue<
        PQElement, std::vector<PQElement>, std::greater<PQElement> > PQ;

        CHGraphT& _base_graph;

        struct ThreadData {
            PQ pq;
            std::vector<uint> dists;
            std::vector<uint> reset_dists;
        };
        std::vector<ThreadData> _thread_data;

        uint _num_threads;
        ThreadData& _myThreadData();
        
        //template <typename NodeT, typename EdgeT>
        struct SHTagInfo {
            const Shortcut * shortcut_p; //(const) references can't be sorted
            //Shortcut &shortcut;
            double weight;
            
            SHTagInfo(const Shortcut * shortcut_p, const CHGraphT& _base_graph): shortcut_p(shortcut_p),
                    weight(geo::getTriangleProportion(_base_graph.getNode(shortcut_p->src),
                                                      _base_graph.getNode(shortcut_p->tgt),
                                                      _base_graph.getNode(shortcut_p->center_node))) {}
            
            //double factor = sh_tag_info.weight * pow(10,5);
            bool operator <(const SHTagInfo &rhs) const {
                return this->weight < rhs.weight;
            }
        };
        
        struct sortNode {
            NodeID node_id;
            double weight;
            bool operator <(const sortNode &rhs) const {
                return this->weight < rhs.weight;
            }
        };

        std::vector<Shortcut> _new_shortcuts;
        std::vector<int> _edge_diffs;
        std::vector<NodeID> _remove;
        std::vector<bool> _to_remove;
        std::mutex _new_shortcuts_mutex;


        void _initVectors();
        void _contract(NodeID node);
        std::vector<Shortcut> _contract(NodeID node, ThreadData& td) const;
        void _quickContract(NodeID node);
        std::vector<Shortcut> _calcShortcuts(Shortcut const& start_edge, NodeID center_node,
                EdgeType direction, ThreadData& td) const;
        
        //tagging functions
        const std::vector<SHTagInfo> _calcTagInfos(NodeID node_id, EdgeType edge_type, uint lastRoundLvl) const;
        void _tagShortcutsOfNode(NodeID node_id, EdgeType edge_type, uint lastRoundLvl, ThreadData &td);       
        void _taggingLastRoundShortcuts(uint round);        
        bool _otherPathExist(ThreadData& td, NodeID start_node, NodeID end_node, uint radius) const;
                
        void _calcShortestDists(ThreadData& td, NodeID start_node, EdgeType direction,
                uint radius) const;
        
        Shortcut _createShortcut(Shortcut const& edge1, Shortcut const& edge2,
                EdgeType direction = EdgeType::OUT) const;
        

        
        void _chooseRemoveNodes(std::vector<NodeID> const& independent_set);
        void _chooseAllForRemove(std::vector<NodeID> const& independent_set);
        void _removeNodes(std::vector<NodeID>& nodes);
    public:
        void _markNeighbours(NodeID node, std::vector<bool>& marked) const;
        
        CHConstructor(CHGraphT& base_graph, uint num_threads = 1);

        /* functions for contraction */
        void quickContract(std::vector<NodeID>& nodes, uint max_degree,
                uint max_rounds);
        void contract(std::vector<NodeID>& nodes);
        void contract(std::vector<NodeID>& nodes, Prioritizer& prioritizer, bool taggingVisuallyUnpleasantSh);
        void rebuildCompleteGraph();

        /* const functions that use algorithms from the CHConstructor */
        std::vector<NodeID> calcIndependentSet(std::vector<sortNode> const& sorted_nodes,
                uint max_degree = MAX_UINT) const;
        std::vector<NodeID> calcIndependentSet(std::vector<NodeID> const& nodes,
                uint max_degree = MAX_UINT) const;
        int calcEdgeDiff(NodeID node) const;            
        std::vector<int> calcEdgeDiffs(std::vector<NodeID> const& nodes) const;
        std::vector<int> calcMaxStreetTypes(std::vector<NodeID> const& nodes) const;
        
        void calcNodeWeights(std::vector<NodeID> const& nodes, std::vector<double> &node_weights) const;
        std::vector<sortNode> calcNodeWeights(std::vector<NodeID> const& nodes) const;
        std::vector<int> calcWeightedEdgeDiffs(std::vector<NodeID> const& nodes) const;
        std::vector<int> calcGeoImportance(std::vector<NodeID> const& nodes) const;
        std::vector<NodeWeight> calcGeoImportance2(std::vector<NodeID> const& nodes, ls::ErrorMeasureType errorMeasure_type) const;
        std::vector<Shortcut> getShortcutsOfContracting(NodeID node) const;
        std::vector<std::vector<Shortcut>> getShortcutsOfContracting(std::vector<NodeID> const& nodes) const;
        std::vector<Shortcut> getShortcutsOfQuickContracting(NodeID node) const;
        std::vector<std::vector<Shortcut>> getShortcutsOfQuickContracting(std::vector<NodeID> const& nodes) const;

        double chainNodeShortcutWeight(NodeID node_id) const;                
        
        friend void unit_tests::testCHConstructor();   
        
        
    };

    /*
     * private
     */

    template <typename NodeT, typename EdgeT>
    struct CHConstructor<NodeT, EdgeT>::CompInOutProduct {
        CHGraphT const& g;

        CompInOutProduct(CHGraphT const& g)
        : g(g) {
        }

        bool operator()(NodeID node1, NodeID node2) const {
            uint edge_product1(g.getNrOfEdges(node1, EdgeType::IN)
                    * g.getNrOfEdges(node1, EdgeType::OUT));
            uint edge_product2(g.getNrOfEdges(node2, EdgeType::IN)
                    * g.getNrOfEdges(node2, EdgeType::OUT));

            return edge_product1 < edge_product2;
        }
    };

    template <typename NodeT, typename EdgeT>
    struct CHConstructor<NodeT, EdgeT>::PQElement {
        NodeID node;
        uint _dist;

        PQElement(NodeID node, uint dist)
        : node(node), _dist(dist) {
        }

        bool operator>(PQElement const& other) const {
            return _dist > other._dist;
        }

        /* make interface look similar to an edge */
        uint distance() const {
            return _dist;
        }
    };

    template <typename NodeT, typename EdgeT>
    auto CHConstructor<NodeT, EdgeT>::_myThreadData() -> ThreadData& {
        return _thread_data[omp_get_thread_num()];
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_initVectors() {
        _new_shortcuts.clear();
        _remove.clear();
        _to_remove.assign(_base_graph.getNrOfNodes(), false);
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_contract(NodeID node) {
        ThreadData & td(_myThreadData());
        auto shortcuts(_contract(node, td));

        _edge_diffs[node] = (int) shortcuts.size() - (int) _base_graph.getNrOfEdges(node);

        std::unique_lock<std::mutex> lock(_new_shortcuts_mutex);
        _new_shortcuts.insert(_new_shortcuts.end(), shortcuts.begin(), shortcuts.end());
    }

    template <typename NodeT, typename EdgeT>
    auto CHConstructor<NodeT, EdgeT>::_contract(NodeID node, ThreadData& td) const -> std::vector<Shortcut> {
        EdgeType search_direction;

        if (_base_graph.getNrOfEdges(node, EdgeType::IN) <= _base_graph.getNrOfEdges(node, EdgeType::OUT)) {
            search_direction = EdgeType::OUT;
        } else {
            search_direction = EdgeType::IN;
        }

        std::vector<Shortcut> shortcuts;
        for (auto const& edge : _base_graph.nodeEdges(node, !search_direction)) {
            if (edge.tgt == edge.src) continue; /* skip loops */
            auto new_shortcuts(_calcShortcuts(edge, node, search_direction, td));
            shortcuts.insert(shortcuts.end(), new_shortcuts.begin(), new_shortcuts.end());
        }

        return shortcuts;
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_quickContract(NodeID node) {
        auto shortcuts(getShortcutsOfQuickContracting(node));

        std::unique_lock<std::mutex> lock(_new_shortcuts_mutex);
        _new_shortcuts.insert(_new_shortcuts.end(), shortcuts.begin(), shortcuts.end());
    }

    template <typename NodeT, typename EdgeT>
    auto CHConstructor<NodeT, EdgeT>::_calcShortcuts(Shortcut const& start_edge, NodeID center_node,
            EdgeType direction, ThreadData& td) const -> std::vector<Shortcut> {

                struct Target {
                    NodeID end_node;
                    Shortcut const& end_edge;
                };
                std::vector<Target> targets;
                targets.reserve(_base_graph.getNrOfEdges(center_node, direction));

                NodeID start_node(otherNode(start_edge, !direction));
                uint radius = 0;

                for (auto const& edge : _base_graph.nodeEdges(center_node, direction)) {
                    if (edge.tgt == edge.src) continue; /* skip loops */
                    auto const end_node = otherNode(edge, direction);
                    if (start_node == end_node) continue; /* don't create loops */

                    radius = std::max(radius, edge.distance());
                    targets.emplace_back(Target{end_node, edge});
                }
                radius += start_edge.distance();

                _calcShortestDists(td, start_node, direction, radius);

                /* abort if start_edge wasn't a shortest path from start_node to center_node */
                if (td.dists[center_node] != start_edge.distance()) return std::vector<Shortcut>();

                std::vector<Shortcut> shortcuts;
                for (auto const& target : targets) {
                    /* we know a path within radius - so _calcShortestDists must have found one */
                    assert(c::NO_DIST != td.dists[target.end_node]);

                    uint center_node_dist(start_edge.distance() + target.end_edge.distance());
                    if (td.dists[target.end_node] == center_node_dist) {
                        shortcuts.push_back(_createShortcut(start_edge, target.end_edge, direction));
                    }
                }

                return shortcuts;
            }      
    
    template <typename NodeT, typename EdgeT>
    const std::vector<typename CHConstructor<NodeT, EdgeT>::SHTagInfo> CHConstructor<NodeT, EdgeT>::_calcTagInfos(NodeID node_id, EdgeType edge_type, uint lastRoundLvl) const {
        std::vector<CHConstructor<NodeT, EdgeT>::SHTagInfo> sh_tag_infos;
        for (const Shortcut &shortcut: _base_graph.nodeEdges(node_id, EdgeType::OUT)) {
            const Shortcut * shortcut_p = &shortcut;
            if(_base_graph.isShortcutOfRound(shortcut, lastRoundLvl)) {
                sh_tag_infos.push_back(SHTagInfo(shortcut_p, _base_graph));
            }
        }
        std::sort(sh_tag_infos.begin(), sh_tag_infos.end());            
        return sh_tag_infos;
    }
    
    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_tagShortcutsOfNode(NodeID node_id, EdgeType edge_type, uint lastRoundLvl, ThreadData &td) {        
        std::vector<SHTagInfo> sh_tag_infos = _calcTagInfos(node_id, edge_type, lastRoundLvl);                
        for (CHConstructor<NodeT, EdgeT>::SHTagInfo sh_tag_info: sh_tag_infos) {                        
            const Shortcut &shortcut = *sh_tag_info.shortcut_p;                        
            if (!_otherPathExist(td, shortcut.src, shortcut.tgt, (1+sh_tag_info.weight)*shortcut.dist)) {
                _base_graph.setEdgeFlag(shortcut.id, true);
            }
        }  
    }
    
    
    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_taggingLastRoundShortcuts(uint round) {
        if(round <= 1) return; //nothing to do for round 1        
        uint lastRound = round -1;        
        uint lastRoundLvl = lastRound-1;
        
        //initialize thread data
        ThreadData td;
        uint nr_of_nodes = _base_graph.getNrOfNodes();
        td.dists.assign(nr_of_nodes, c::NO_DIST);
        td.reset_dists.reserve(nr_of_nodes);
        
        //initialize all last shortcuts with unpleasing flag        
        for (const EdgeT &edge: _base_graph.getAllEdges()) {            
            if(_base_graph.isShortcutOfRound(edge.id, lastRoundLvl)) {
                _base_graph.setEdgeFlag(edge.id, false);
            }
        }
        
        for (int node_id = 0; node_id < _base_graph.getNrOfNodes(); node_id++) {             
            _tagShortcutsOfNode(node_id, EdgeType::IN, lastRoundLvl, td);
            _tagShortcutsOfNode(node_id, EdgeType::OUT, lastRoundLvl, td );
        }                        
    }
    
    template <typename NodeT, typename EdgeT>
    bool CHConstructor<NodeT, EdgeT>::_otherPathExist(ThreadData& td, NodeID start_node, NodeID end_node,
            /*EdgeType direction,*/ uint radius) const {        
        
        /* clear thread data first */
        td.pq = PQ();        
        for (auto node_id : td.reset_dists) {
            td.dists[node_id] = c::NO_DIST;
        }
        td.reset_dists.clear();
         

        /* now initialize with start node */
        td.pq.push(PQElement(start_node, 0));
        td.dists[start_node] = 0;
        td.reset_dists.push_back(start_node);

        while (!td.pq.empty() && td.pq.top().distance() <= radius && td.dists[end_node]== c::NO_DIST) {
            auto top = td.pq.top();
            td.pq.pop();
            if (td.dists[top.node] != top.distance()) continue;

            for (auto const& edge : _base_graph.nodeEdges(top.node, EdgeType::OUT)) {
                if (edge.speed == -1000) break; //visually unpleasing edges are not taken

                NodeID tgt_node(otherNode(edge, EdgeType::OUT));
                uint new_dist(top.distance() + edge.distance());

                if (new_dist < td.dists[tgt_node]) {
                    if (td.dists[tgt_node] == c::NO_DIST) {
                        td.reset_dists.push_back(tgt_node);
                    }
                    td.dists[tgt_node] = new_dist;
                    td.pq.push(PQElement(tgt_node, new_dist));
                }
            }
        }
        return td.dists[end_node]!= c::NO_DIST;
    }
    
    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_calcShortestDists(ThreadData& td, NodeID start_node,
            EdgeType direction, uint radius) const {
        /* calculates all shortest paths within radius distance from start_node */

        /* clear thread data first */
        td.pq = PQ();
        for (auto node_id : td.reset_dists) {
            td.dists[node_id] = c::NO_DIST;
        }
        td.reset_dists.clear();

        /* now initialize with start node */
        td.pq.push(PQElement(start_node, 0));
        td.dists[start_node] = 0;
        td.reset_dists.push_back(start_node);

        while (!td.pq.empty() && td.pq.top().distance() <= radius) {
            auto top = td.pq.top();
            td.pq.pop();
            if (td.dists[top.node] != top.distance()) continue;

            for (auto const& edge : _base_graph.nodeEdges(top.node, direction)) {
                NodeID tgt_node(otherNode(edge, direction));
                uint new_dist(top.distance() + edge.distance());

                if (new_dist < td.dists[tgt_node]) {
                    if (td.dists[tgt_node] == c::NO_DIST) {
                        td.reset_dists.push_back(tgt_node);
                    }
                    td.dists[tgt_node] = new_dist;
                    td.pq.push(PQElement(tgt_node, new_dist));
                }
            }
        }
    }

    template <typename NodeT, typename EdgeT>
    auto CHConstructor<NodeT, EdgeT>::_createShortcut(Shortcut const& edge1, Shortcut const& edge2,
            EdgeType direction) const -> Shortcut {
                /* make sure no "loop" edges are used */
                assert(edge1.src != edge1.tgt && edge2.src != edge2.tgt);

                if (direction == EdgeType::OUT) {
                    /* make sure no "loop" edges are created */
                    assert(edge1.src != edge2.tgt);
                    return make_shortcut(edge1, edge2);
                } else {
                    /* make sure no "loop" edges are created */
                    assert(edge2.src != edge1.tgt);
                    return make_shortcut(edge2, edge1);
                }
            }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_markNeighbours(NodeID node, std::vector<bool>& marked) const {
        for (uint i(0); i < 2; i++) {
            for (auto const& edge : _base_graph.nodeEdges(node, (EdgeType) i)) {
                marked[otherNode(edge, (EdgeType) i)] = true;
            }
        }
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_chooseRemoveNodes(std::vector<NodeID> const& independent_set) {
        double edge_diff_mean(0);
        for (NodeID node : independent_set) {
            edge_diff_mean += _edge_diffs[node];
        }
        edge_diff_mean /= independent_set.size();
        Print("The average edge difference is " << edge_diff_mean << ".");

        assert(_remove.empty());
        for (NodeID node : independent_set) {
            if (_edge_diffs[node] <= edge_diff_mean) {
                _remove.push_back(node);
                _to_remove[node] = true;
            }
        }
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_chooseAllForRemove(std::vector<NodeID> const& independent_set) {
        assert(_remove.empty());
        for (NodeID node : independent_set) {
            _remove.push_back(node);
            _to_remove[node] = true;
        }
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::_removeNodes(std::vector<NodeID>& nodes) {
        size_t remaining_nodes(nodes.size());
        size_t i(0);
        while (i < remaining_nodes) {
            NodeID node(nodes[i]);
            if (_to_remove[node]) {
                remaining_nodes--;
                nodes[i] = nodes[remaining_nodes];
                nodes[remaining_nodes] = node;
            } else {
                i++;
            }
        }

        nodes.resize(remaining_nodes);
    }

    /*
     * public
     */

    template <typename NodeT, typename EdgeT>
    CHConstructor<NodeT, EdgeT>::CHConstructor(CHGraphT& base_graph, uint num_threads)
    : _base_graph(base_graph), _num_threads(num_threads) {
        if (!_num_threads) {
            _num_threads = 1;
        }

        uint nr_of_nodes(_base_graph.getNrOfNodes());

        _thread_data.resize(_num_threads);
        _edge_diffs.resize(nr_of_nodes);
        _to_remove.resize(nr_of_nodes);

        for (auto& td : _thread_data) {
            td.dists.assign(nr_of_nodes, c::NO_DIST);
            td.reset_dists.reserve(nr_of_nodes);
        }
        _new_shortcuts.reserve(_base_graph.getNrOfEdges());
        _remove.reserve(nr_of_nodes);
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::quickContract(std::vector<NodeID>& nodes, uint max_degree, uint max_rounds) {
        using namespace std::chrono;

        Print("\nStarting the quick_contraction of nodes with degree smaller than " << max_degree << ".\n");

        for (uint round(1); round <= max_rounds; ++round) {
            steady_clock::time_point t1 = steady_clock::now();
            Print("Starting round " << round);
            Debug("Initializing the vectors for a new round.");
            _initVectors();

            Print("Sorting the remaining " << nodes.size() << " nodes.");
            std::sort(nodes.begin(), nodes.end(), CompInOutProduct(_base_graph));

            Debug("Constructing the independent set.");
            auto independent_set = calcIndependentSet(nodes, max_degree);
            Print("The independent set has size " << independent_set.size() << ".");

            if (independent_set.empty()) break;

            Debug("Quick-contracting all the nodes in the independent set.");
            uint size(independent_set.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
            for (uint i = 0; i < size; i++) {
                uint node(independent_set[i]);
                _quickContract(node);
            }
            Print("Number of possible new Shortcuts: " << _new_shortcuts.size());

            Debug("Remove the nodes with low edge difference.");
            _chooseAllForRemove(independent_set);
            _removeNodes(nodes);
            Print("Removed " << _remove.size() << " nodes with low edge difference.");

            Debug("Restructuring the graph.");
            _base_graph.restructure(_remove, _to_remove, _new_shortcuts);

            Print("Graph info:");
            _base_graph.printInfo(nodes);

            duration<double> time_span = duration_cast<duration<double>>(steady_clock::now() - t1);
            Print("Round took " << time_span.count() << " seconds.\n");
            Unused(time_span);
        }
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::contract(std::vector<NodeID>& nodes) {
        using namespace std::chrono;

        Print("\nStarting the contraction of " << nodes.size() << " nodes.\n");

        for (uint round(1); !nodes.empty(); ++round) {
            steady_clock::time_point t1 = steady_clock::now();
            Print("Starting round " << round);            
            Debug("Initializing the vectors for a new round.");
            _initVectors();

            Print("Sorting the remaining " << nodes.size() << " nodes.");
            std::sort(nodes.begin(), nodes.end(), CompInOutProduct(_base_graph));

            Debug("Constructing the independent set.");
            auto independent_set = calcIndependentSet(nodes);
            Print("The independent set has size " << independent_set.size() << ".");

            Debug("Contracting all the nodes in the independent set.");
            uint size(independent_set.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
            for (uint i = 0; i < size; i++) {
                uint node(independent_set[i]);
                _contract(node);
            }
            Print("Number of possible new Shortcuts: " << _new_shortcuts.size());

            Debug("Remove the nodes with low edge difference.");
            _chooseRemoveNodes(independent_set);
            _removeNodes(nodes);
            Print("Removed " << _remove.size() << " nodes with low edge difference.");

            Debug("Restructuring the graph.");
            _base_graph.restructure(_remove, _to_remove, _new_shortcuts);

            Print("Graph info:");
            _base_graph.printInfo(nodes);

            duration<double> time_span = duration_cast<duration<double>>(steady_clock::now() - t1);
            Print("Round took " << time_span.count() << " seconds.\n");
            Unused(time_span);
        }
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::contract(std::vector<NodeID>& nodes, Prioritizer& prioritizer, bool taggingVisuallyUnpleasantSh) {
        using namespace std::chrono;

        Print("\nStarting the contraction of " << nodes.size() << " nodes.\n");

        prioritizer.init(nodes);

        uint round(1);
        while (prioritizer.hasNodesLeft()) {
            steady_clock::time_point t1 = steady_clock::now();
            Print("Starting round " << round);
            
            if(taggingVisuallyUnpleasantSh) {
                Print("Tagging visually unpleasant for shortcuts of last round");
                _taggingLastRoundShortcuts(round);
            }
            
            Debug("Initializing the vectors for a new round.");
            _initVectors();
            
            Debug("Calculating list of nodes to be contracted next.");
            auto next_nodes(prioritizer.extractNextNodes());
            Print("There are " << next_nodes.size() << " nodes to be contracted in this round.");

            Debug("Contracting all the nodes in the independent set.");
            uint size(next_nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
            for (uint i = 0; i < size; i++) {
                uint node(next_nodes[i]);
                _contract(node);
            }
            Print("Number of new Shortcuts: " << _new_shortcuts.size());

            Debug("Mark nodes for removal from graph.");
            _chooseAllForRemove(next_nodes);
            Print("Marked " << _remove.size() << " nodes.");

            Debug("Restructuring the graph.");
            _base_graph.restructure(_remove, _to_remove, _new_shortcuts);

            Print("Graph info:");
            _base_graph.printInfo();

            duration<double> time_span = duration_cast<duration<double>>(steady_clock::now() - t1);
            Print("Round took " << time_span.count() << " seconds.\n");
            Unused(time_span);

            round++;
        }
    }

    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::rebuildCompleteGraph() {
        Print("Restoring edges from contracted nodes.");

        _base_graph.rebuildCompleteGraph();
    }

    template <typename NodeT, typename EdgeT>
    std::vector<NodeID> CHConstructor<NodeT, EdgeT>::calcIndependentSet(std::vector<NodeID> const& nodes,
            uint max_degree) const {
        std::vector<NodeID> independent_set;
        std::vector<bool> marked(_base_graph.getNrOfNodes(), false);
        independent_set.reserve(nodes.size());

        for (NodeID node : nodes) {
            if (!marked[node] && max_degree >= _base_graph.getNrOfEdges(node)) {
                marked[node] = true;
                _markNeighbours(node, marked);
                independent_set.push_back(node);
            }
        }
        return independent_set;
    }
    
    template <typename NodeT, typename EdgeT>
    std::vector<NodeID> CHConstructor<NodeT, EdgeT>::calcIndependentSet(std::vector<sortNode> const& sorted_nodes,
            uint max_degree) const {
        std::vector<NodeID> independent_set;
        std::vector<bool> marked(_base_graph.getNrOfNodes(), false);
        independent_set.reserve(sorted_nodes.size());

        for (const sortNode &sort_node : sorted_nodes) {
            NodeID node = sort_node.node_id;
            if (!marked[node] && max_degree >= _base_graph.getNrOfEdges(node)) {
                marked[node] = true;
                _markNeighbours(node, marked);
                independent_set.push_back(node);
            }
        }
        return independent_set;
    }

    template <typename NodeT, typename EdgeT>
    int CHConstructor<NodeT, EdgeT>::calcEdgeDiff(NodeID node) const {
        auto shortcuts(getShortcutsOfContracting(node));
        return (int) shortcuts.size() - (int) _base_graph.getNrOfEdges(node);
    }
       
    template <typename NodeT, typename EdgeT>
    std::vector<int> CHConstructor<NodeT, EdgeT>::calcEdgeDiffs(std::vector<NodeID> const& nodes) const {
        std::vector<int> edge_diffs(nodes.size());
        auto shortcuts(getShortcutsOfContracting(nodes));

        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            edge_diffs[i] = (int) shortcuts[i].size() - (int) _base_graph.getNrOfEdges(nodes[i]);
        }

        return edge_diffs;
    }
    
    template <typename NodeT, typename EdgeT>
    std::vector<int> CHConstructor<NodeT, EdgeT>::calcMaxStreetTypes(std::vector<NodeID> const& nodes) const {
        std::vector<int> max_street_types(nodes.size());
        //auto shortcuts(getShortcutsOfContracting(nodes));

        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            max_street_types[i] = (int) _base_graph.getMaxStreetType(nodes[i]);
        }

        return max_street_types;
    }
    
    template <typename NodeT, typename EdgeT>
    void CHConstructor<NodeT, EdgeT>::calcNodeWeights(std::vector<NodeID> const& nodes, std::vector<double> &node_weights) const {                        
        auto shortcuts(getShortcutsOfContracting(nodes));
        
        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            NodeID node_id = nodes[i];
            node_weights[node_id] = (int) shortcuts[i].size() - (int) _base_graph.getNrOfEdges(node_id);            
        }        
    } 
    
    
    template <typename NodeT, typename EdgeT>
    std::vector<typename CHConstructor<NodeT, EdgeT>::sortNode> CHConstructor<NodeT, EdgeT>::calcNodeWeights(std::vector<NodeID> const& nodes) const {
        std::vector<sortNode> node_weights(nodes.size());        
        auto shortcuts(getShortcutsOfContracting(nodes));

        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            node_weights[i].node_id = i;
            node_weights[i].weight = (int) shortcuts[i].size() - (int) _base_graph.getNrOfEdges(nodes[i]);
        }

        return node_weights;
    }
    
    template <typename NodeT, typename EdgeT>
    std::vector<int> CHConstructor<NodeT, EdgeT>::calcWeightedEdgeDiffs(std::vector<NodeID> const& nodes) const {
        std::vector<int> edge_diffs(nodes.size());
        auto shortcuts(getShortcutsOfContracting(nodes));

        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            /*
            auto edges = _base_graph. nodeEdges(nodes[i]);
            int old_weight = 0;            
            for (Shortcut shortcut: shortcuts[i]) {
                old_weight += pow(10, shortcut.type);               
            }*/
            
            int sh_weight = 0;            
            for (Shortcut shortcut: shortcuts[i]) {
                sh_weight += pow(0.1, shortcut.type);               
            }
            edge_diffs[i] = sh_weight;// - (int) _base_graph.getNrOfEdges(nodes[i]);
        }

        return edge_diffs;
    }
    
    template <typename NodeT, typename EdgeT>
    std::vector<NodeWeight> CHConstructor<NodeT, EdgeT>::calcGeoImportance2(std::vector<NodeID> const& nodes,
            ls::ErrorMeasureType errorMeasure_type) const {
        std::vector<NodeWeight> node_weights(nodes.size());
        auto shortcuts(getShortcutsOfContracting(nodes));

        std::unique_ptr<ls::ErrorMeasure> errorMeasure = ls::createErrorMeasure(errorMeasure_type);
        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            //NodeID node_id = nodes[i];            
            
            NodeWeight node_weight(0, 0, 0);
            
            for (const Shortcut& shortcut : shortcuts[i]) {                
                const auto& src_node = _base_graph.getNode(shortcut.src);
                const auto& tgt_node = _base_graph.getNode(shortcut.tgt);
                const auto& center_node = _base_graph.getNode(shortcut.center_node);
                
                //double geo_dist = geo::geoDist(src_node, tgt_node);
                assert(shortcut.dist> 0);
                //geo_importance += geo_dist * pow(geo_dist/edge.dist, 2);                            
                //node_weight.geo_measure = std::max(node_weight.geo_measure, geo_dist * pow(geo_dist/shortcut.dist, 1));                            
                node_weight.geo_measure = std::max(node_weight.geo_measure, ew::calcWeight(shortcut.dist, src_node, tgt_node));
                
                node_weight.error = std::max(node_weight.error, errorMeasure->calcError(src_node, tgt_node, center_node));                            
            }            
                        
            
            /*
            double divisor = shortcuts[i].size();
            if(divisor != 0) {
                //node_weight.geo_measure = geo_importance/_base_graph.getNrOfEdges(nodes[i]);
                node_weight.geo_measure = geo_importance/divisor;
            } else {
                node_weight.geo_measure = 0;
            }            */
            //node_weight.geo_measure = geo_importance;
            node_weight.edge_diff = (int) shortcuts[i].size(); //- (int) _base_graph.getNrOfEdges(nodes[i]);
            node_weights[i] = node_weight;
        }

        return node_weights;
    }
    
    template <typename NodeT, typename EdgeT>
    std::vector<int> CHConstructor<NodeT, EdgeT>::calcGeoImportance(std::vector<NodeID> const& nodes) const {
        std::vector<int> edge_diffs(nodes.size());
        auto shortcuts(getShortcutsOfContracting(nodes));

        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            NodeID node_id = nodes[i];
            /*
            auto edges = _base_graph. nodeEdges(nodes[i]);
            int old_weight = 0;            
            for (Shortcut shortcut: shortcuts[i]) {
                old_weight += pow(10, shortcut.type);               
            }*/
            
            double weight = 0;
            uint assert_counter = 0;            
            for (auto const& edge : _base_graph.nodeEdges(node_id, EdgeType::IN)) {
                auto other_node = otherNode(edge, EdgeType::IN);
                double geo_dist = geo::geoDist(_base_graph.getNode(node_id), _base_graph.getNode(other_node));
                weight += geo_dist * geo_dist/edge.dist;
                //double speed = geo_dist/edge.dist;
                //double fraction = edge.speed/speed;
                assert_counter++;
                //30000000 Umrechnungsfaktor
            }
            for (auto const& edge : _base_graph.nodeEdges(node_id, EdgeType::OUT)) {
                auto other_node = otherNode(edge, EdgeType::OUT);
                double geo_dist = geo::geoDist(_base_graph.getNode(node_id), _base_graph.getNode(other_node));
                weight += geo_dist * geo_dist/edge.dist;
                assert_counter++;
            }
            assert(assert_counter == _base_graph.getNrOfEdges(nodes[i]));
            
            
            edge_diffs[i] = weight;
        }

        return edge_diffs;
    }
    
       
    template <typename NodeT, typename EdgeT>
    auto CHConstructor<NodeT, EdgeT>::getShortcutsOfContracting(NodeID node) const -> std::vector<Shortcut> {
        ThreadData td;
        return _contract(node, td);
    }

    template <typename NodeT, typename EdgeT>
    auto CHConstructor<NodeT, EdgeT>::getShortcutsOfContracting(std::vector<NodeID> const& nodes) const
    -> std::vector<std::vector<Shortcut>>
    {
        std::vector<std::vector < Shortcut >> shortcuts(nodes.size());

        /* init thread data */
        std::vector<ThreadData> thread_data;
        thread_data.resize(_num_threads);
        auto nr_of_nodes(_base_graph.getNrOfNodes());
        for (auto& td : thread_data) {
            td.dists.assign(nr_of_nodes, c::NO_DIST);
            td.reset_dists.reserve(nr_of_nodes);
        }

        /* calc shortcuts */
        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            uint node(nodes[i]);
            
            auto td(thread_data[omp_get_thread_num()]);
            shortcuts[i] = _contract(node, td);
        }

        return shortcuts;
    }

    template <typename NodeT, typename EdgeT>
    auto CHConstructor<NodeT, EdgeT>::getShortcutsOfQuickContracting(NodeID node) const -> std::vector<Shortcut> {
        std::vector<Shortcut> shortcuts;
        for (auto const& in_edge : _base_graph.nodeEdges(node, EdgeType::IN)) {
            if (in_edge.tgt == in_edge.src) continue; /* skip loops */
            for (auto const& out_edge : _base_graph.nodeEdges(node, EdgeType::OUT)) {
                if (out_edge.tgt == out_edge.src) continue; /* skip loops */
                if (in_edge.src != out_edge.tgt) { /* don't create loops */
                    shortcuts.push_back(_createShortcut(in_edge, out_edge));
                }
            }
        }

        return shortcuts;
    }

    template <typename NodeT, typename EdgeT>
    auto CHConstructor<NodeT, EdgeT>::getShortcutsOfQuickContracting(std::vector<NodeID> const& nodes) const -> std::vector<std::vector<Shortcut>>
    {
        std::vector<std::vector < Shortcut >> shortcuts(nodes.size());

        /* calc shortcuts */
        uint size(nodes.size());
#pragma omp parallel for num_threads(_num_threads) schedule(dynamic)
        for (uint i = 0; i < size; i++) {
            shortcuts[i] = getShortcutsOfQuickContracting(nodes[i]);
        }

        return shortcuts;
    }
    
    //get weight for shortcut of a chainnode
    template <typename NodeT, typename EdgeT>
    double CHConstructor<NodeT, EdgeT>::chainNodeShortcutWeight(NodeID node_id) const
    {    
        assert(_base_graph.nodeNeighbours(node_id).size() <= 2);

        EdgeT in_edge;
        bool in_edge_found= false;
        NodeID in_edge_other_node = c::NO_NID;
        EdgeT out_edge;
        bool out_edge_found= false;
        NodeID out_edge_other_node = c::NO_NID;
        for (const EdgeT& edge: _base_graph.nodeEdges(node_id, EdgeType::IN)) {
            in_edge_other_node = otherNode(edge, EdgeType::IN);
            in_edge = edge;
            in_edge_found = true;        
            break;
        }
        if (in_edge_found) {
            for (const EdgeT& edge: _base_graph.nodeEdges(node_id, EdgeType::OUT)) {
                out_edge_other_node = otherNode(edge, EdgeType::OUT);
                if (out_edge_other_node != in_edge_other_node) {
                    out_edge = edge;
                    out_edge_found = true;            
                    break;
                }            
            }
            if (out_edge_found) {
                double shortcut_dist = in_edge.dist + out_edge.dist;
                //double geo_dist = geo::geoDist(getNode(in_edge_other_node), getNode(out_edge_other_node));
                return ew::calcWeight(shortcut_dist, _base_graph.getNode(in_edge_other_node), _base_graph.getNode(out_edge_other_node));
            }
        }
        return 0;
    }

}