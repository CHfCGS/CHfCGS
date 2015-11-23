#include <vector>
#include <limits>
#include <queue>

#include "chgraph.h"
#include "nodes_and_edges.h"

namespace chm {

template <class GraphT>
class CHDijkstra {
    GraphT myGraph;
    std::vector<double> distFwd;
    std::vector<double> distBwd;
    std::vector<bool> settledFwd;
    std::vector<bool> settledBwd;
       
    typedef std::priority_queue<BDPQElement, std::vector<BDPQElement>, std::greater<BDPQElement> > BDPQ;
    BDPQ myQueue;
    std::vector<NodeID> touchedNodes;    
    int nofTouchedNodes;

    CHDijkstra(GraphT _myGraph):myGraph(_mygraph) {        
        distFwd.resize(myGraph.getNrOfNodes);
        distBwd.resize(myGraph.getNrOfNodes);
        settledFwd.resize(myGraph.getNrOfNodes);;
        settledBwd.resize(myGraph.getNrOfNodes);;
        touchedNodes.resize(myGraph.getNrOfNodes);;
        for (int i = 0; i < myGraph.resize(myGraph.getNrOfNodes); i++) {
            distFwd[i] = numeric_limits<double>::max();
            distBwd[i] = numeric_limits<double>::max();
            settledFwd[i] = settledBwd[i] = false;
        }
        myQueue = new PriorityQueue<BDPQElement>();
        nofTouchedNodes = 0;
    }

    void labelFwd(int v, int d) {
        if ((distFwd[v] == numeric_limits<double>::max())&&(distBwd[v] == numeric_limits<double>::max())) {
            touchedNodes[nofTouchedNodes++] = v;
        }
        distFwd[v] = d;
        myQueue.add(BDPQElement(d, v, 0));
    }

    void labelBwd(int v, int d) {
        if ((distFwd[v] == numeric_limits<double>::max())&&(distBwd[v] == numeric_limits<double>::max())) {
            touchedNodes[nofTouchedNodes++] = v;
        }
        distBwd[v] = d;
        myQueue.add(BDPQElement(d, v, 1));
    }   

    int runDijkstra(int src, int trg) // returns distance from src to trg
    {

        // clean up previously touched nodes
        for (int i = 0; i < nofTouchedNodes; i++) {
            distFwd[touchedNodes[i]] = distBwd[touchedNodes[i]] = Integer.MAX_VALUE;
            settledFwd[touchedNodes[i]] = settledBwd[touchedNodes[i]] = false;
        }
        nofTouchedNodes = 0;
        myQueue.clear();
        // start with src
        labelFwd(src, 0);
        labelBwd(trg, 0);

        double bestDist = numeric_limits<double>::max();

        int edgeCount = 0;
        // otherwise we have to process pq until settling trg
        bool phaseFinished = false;
        while ((!myQueue.isEmpty())&&(phaseFinished == false)) {
            //System.out.print(".");
            BDPQElement cur = myQueue.top();
            myQueue.remove();
            double cur_dist = cur.key;
            NodeID cur_node = cur.value;
            int cur_side = cur.queue;

            //if (cur_dist>bestDist)
            //	phaseFinished=true;

            if (cur_side == 0) // we are in forward search
            {
                if (cur_dist == distFwd[cur_node]) { //suspicious
                    settledFwd[cur_node] = true;

                    bool stalled = false;

                    // check for stalling (if there is a node tmp_node (ABOVE) and an edge (tmp_node,cur_node)
                    // which sum to a smaller distance (!)
                    //for(int j=0; j<myGraph.nofInEdges(cur_node); j++)
                    for (auto &tmp_edge: myGraph.nodeEdges(cur_node, EdgeType::IN)) {                        
                        int tmp_wgt = tmp_edge.
                        NodeID tmp_node = tmp_edge.src;
                        if (distFwd[cur_node] - tmp_wgt > distFwd[tmp_node]) {
                            stalled = true;
                            break;
                        }
                        if (myGraph.level(tmp_node) < myGraph.level(cur_node))
                            break;

                    }
                    /*
                    for (int j = myGraph.getNrOFNodes(cur_node) - 1; j >= 0; j--) {
                        int tmp_edge = myGraph.inEdgeID(cur_node, j);
                        int tmp_wgt = myGraph.edgeWeight(tmp_edge);
                        int tmp_node = myGraph.edgeSource(tmp_edge);
                        if (distFwd[cur_node] - tmp_wgt > distFwd[tmp_node]) {
                            stalled = true;
                            break;
                        }
                        if (myGraph.level(tmp_node) < myGraph.level(cur_node))
                            break;
                    }
                     * */
                    if ((settledBwd[cur_node]) &&(distFwd[cur_node] + distBwd[cur_node] < bestDist))
                        bestDist = distFwd[cur_node] + distBwd[cur_node];

                    if (stalled == false)
                        for (int i = myGraph.nofOutEdges(cur_node) - 1; i >= 0; i--) {
                            int cur_edge = myGraph.outEdgeID(cur_node, i);
                            int cur_trg = myGraph.edgeTarget(cur_edge);
                            int cur_weight = myGraph.edgeWeight(cur_edge);
                            if (myGraph.level(cur_trg) >= myGraph.level(cur_node))
                                edgeCount++;
                            else
                                break;
                            if ((myGraph.level(cur_trg) >= myGraph.level(cur_node))&&(distFwd[cur_trg] > cur_dist + cur_weight)) {
                                labelFwd(cur_trg, cur_dist + cur_weight);
                            }
                        }
                }
            } else // we are in backward search
            {
                if (cur_dist == distBwd[cur_node]) {
                    settledBwd[cur_node] = true;
                    boolean stalled = false;

                    // check for stalling: if there is a node ABOVE cur_node ...
                    //for(int j=0; j<myGraph.nofOutEdges(cur_node); j++)
                    for (int j = myGraph.nofOutEdges(cur_node) - 1; j >= 0; j--) {
                        int tmp_edge = myGraph.outEdgeID(cur_node, j);
                        int tmp_wgt = myGraph.edgeWeight(tmp_edge);
                        int tmp_node = myGraph.edgeTarget(tmp_edge);
                        if (distBwd[cur_node] - tmp_wgt > distBwd[tmp_node]) {
                            stalled = true;
                            break;
                        }
                        if (myGraph.level(cur_node) > myGraph.level(tmp_node))
                            break;
                    }


                    if ((settledFwd[cur_node]) &&(distFwd[cur_node] + distBwd[cur_node] < bestDist))
                        bestDist = distFwd[cur_node] + distBwd[cur_node];

                    if (stalled == false)
                        for (int i = myGraph.nofInEdges(cur_node) - 1; i >= 0; i--) {
                            int cur_edge = myGraph.inEdgeID(cur_node, i);
                            int cur_trg = myGraph.edgeSource(cur_edge);
                            int cur_weight = myGraph.edgeWeight(cur_edge);
                            if (myGraph.level(cur_trg) >= myGraph.level(cur_node))
                                edgeCount++;
                            else
                                break;
                            if ((myGraph.level(cur_trg) >= myGraph.level(cur_node))&&(distBwd[cur_trg] > cur_dist + cur_weight)) {
                                labelBwd(cur_trg, cur_dist + cur_weight);
                            }
                        }
                }
            }
        }
        System.out.println("CH-BD has touched " + nofTouchedNodes + " and looked at " + edgeCount + " edges");

        if (1 == 0) {
            // ONLY FOR DEBUGGING
            // now check how many nodes bear correct distances
            int bestTmpDist = Integer.MAX_VALUE;
            Dijkstra myTmpDijkstraF = new Dijkstra(myGraph);
            Dijkstra myTmpDijkstraB = new BwDijkstra(myGraph);
            int fwdOK = 0, fwdTotal = 0, bwdOK = 0, bwdTotal = 0, bothOK = 0;

            for (int i = 0; i < nofTouchedNodes; i++) {
                //if ((i%100)==0)
                //	System.out.println(i+" ");
                int curNode = touchedNodes[i];
                boolean fwdUse = false, bwdUse = false;

                int fwdDist = distFwd[curNode];
                if (fwdDist != Integer.MAX_VALUE) {
                    fwdTotal++;
                    myTmpDijkstraF.runDijkstra(src, curNode);
                    assert(myTmpDijkstraF.dist[curNode] != Integer.MAX_VALUE);
                    if (myTmpDijkstraF.dist[curNode] == fwdDist) {
                        if (myGraph.level(curNode) > 80)
                            System.out.print(curNode + "(" + myGraph.level(curNode) + ") ");
                        fwdOK++;
                        fwdUse = true;
                    }
                }
                int bwdDist = distBwd[curNode];
                if (bwdDist != Integer.MAX_VALUE) {
                    bwdTotal++;
                    myTmpDijkstraB.runDijkstra(trg, curNode);
                    assert(myTmpDijkstraB.dist[curNode] != Integer.MAX_VALUE);
                    if (myTmpDijkstraB.dist[curNode] == bwdDist) {
                        bwdOK++;
                        bwdUse = true;
                    }
                }
                if ((bwdDist != Integer.MAX_VALUE)&&(fwdDist != Integer.MAX_VALUE)) {
                    if (myTmpDijkstraF.dist[curNode] + myTmpDijkstraB.dist[curNode] < bestTmpDist)
                        bestTmpDist = myTmpDijkstraF.dist[curNode] + myTmpDijkstraB.dist[curNode];
                }
                if (bwdUse || fwdUse)
                    bothOK++;
            }
            System.out.println("\n Forward search: " + fwdOK + "/" + fwdTotal);
            System.out.println("Backward search: " + bwdOK + "/" + bwdTotal);
            System.out.println("Best Distance:" + bestTmpDist + " and total nodes: " + bothOK);
        }
        return bestDist;

    }

    void checkCHreqs() {
        for (int i = 0; i < myGraph.nofNodes(); i++) {
            for (int j = 0; j < myGraph.nofOutEdges(i) - 1; j++) {
                int cur_edge = myGraph.outEdgeID(i, j);
                int next_edge = myGraph.outEdgeID(i, j + 1);
                int trg_node = myGraph.edgeTarget(cur_edge);
                int next_trg_node = myGraph.edgeTarget(next_edge);
                assert(myGraph.level(trg_node) <= myGraph.level(next_trg_node));
            }
            for (int j = 0; j < myGraph.nofInEdges(i) - 1; j++) {
                int cur_edge = myGraph.inEdgeID(i, j);
                int next_edge = myGraph.inEdgeID(i, j + 1);
                int src_node = myGraph.edgeSource(cur_edge);
                int next_src_node = myGraph.edgeSource(next_edge);
                assert(myGraph.level(src_node) <= myGraph.level(next_src_node));
            }
        }
        System.out.println("CH Reqs ok!");
    }
   
    
    struct BDPQElement {
        double key;
        NodeID value;
        int queue;

        BDPQElement(double a, NodeID b, int c) {
            key = a;
            value = b;
            queue = c;
        }
        
        bool operator>(BDPQElement const& other) const
            {               
            if (key > other.key) return true;           
            else return false;
        }
    }
}