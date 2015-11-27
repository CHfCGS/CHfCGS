#pragma once

#include "geoFunctions.h"
#include "nodes_and_edges.h"
#include "chgraph.h"

struct ILP_data {
    
    CHGraph<CHNode, CHEdge> &graph;
    const Chain &chain1;
    const Chain &chain2;
        
    ILP_Chain ilp_chain1;
    ILP_Chain ilp_chain2;
    LineID nextLineID = 0;
    std::vector<Line> potEdges1;
    std::vector<Line> potEdges2;
    std::vector<Line> allPotEdges;
    std::vector<Line> followerLines1;
    std::vector<Line> followerLines2;
    std::vector<Line> allFollowerLines;
    std::vector<Intersection> edgeIntersections;
    std::vector<Intersection> followerLinesUnorderings;
    
    bool followerIsPartial = true; //if this is false, the structure is invalid

    
    static bool testIntersection(const Line line1, const Line line2) {
        
        //4 Orientation tests
        double line2StartArea = geo::calcArea(line1.start.ch_node, line1.end.ch_node, line2.start.ch_node);
        double line2EndArea = geo::calcArea(line1.start.ch_node, line1.end.ch_node, line2.end.ch_node);
        bool line1BetweenLine2 = geo::differentSign(line2StartArea, line2EndArea);

        double line1StartArea = geo::calcArea(line2.start.ch_node, line2.end.ch_node, line1.start.ch_node);
        double line1EndArea = geo::calcArea(line2.start.ch_node, line2.end.ch_node, line1.end.ch_node);
        bool line2BetweenLine1 = geo::differentSign(line1StartArea, line1EndArea);

        return line1BetweenLine2 && line2BetweenLine1;
    }

    static std::vector<Intersection> calculateIntersections(const std::vector<Line> &lines) {
        //possible upgrade: Bentley-Ottmann
        //possible upgrade: could be adapted for 2 input vectors
        std::vector<Intersection> intersections;  
        if (lines.size() > 1) {
            for (uint i = 0; i < lines.size() - 1; i++) {
                for (uint j = i + 1; j < lines.size(); j++) {
                    if (testIntersection(lines.at(i), lines.at(j))) {
                        Intersection intersection(lines.at(i), lines.at(j));
                        intersections.push_back(intersection);
                    }
                }
            }
        }
        return intersections;
    }
    

    static bool testUnordering(const Line &line1, const Line &line2) {
        uint line1Chain1Pos, line1Chain2Pos, line2Chain1Pos, line2Chain2Pos;
        assert(line1.start.side != line1.end.side);
        if (line1.start.side) {            
            line1Chain1Pos = line1.start.posInChain;
            line1Chain2Pos = line1.end.posInChain;
        } else {
            line1Chain1Pos = line1.end.posInChain;
            line1Chain2Pos = line1.start.posInChain;
        }
        
        assert(line2.start.side != line2.end.side);
        if (line2.start.side) {            
            line2Chain1Pos = line2.start.posInChain;
            line2Chain2Pos = line2.end.posInChain;
        } else {
            line2Chain1Pos = line2.end.posInChain;
            line2Chain2Pos = line2.start.posInChain;
        }
        
        
        if ((line1Chain1Pos > line2Chain1Pos && line1Chain2Pos < line2Chain2Pos)
                || (line1Chain1Pos < line2Chain1Pos && line1Chain2Pos > line2Chain2Pos)) {
            return true;
        }
        return false;
    }
    
    static std::vector<Intersection> calculate_p_followerUnorderings(const std::vector<Line> &allFollowerLines) {
        std::vector<Intersection> followerLinesUnorderings;       
        assert(allFollowerLines.size() > 0);
        for (uint i = 0; i < allFollowerLines.size() - 1; i++) {
            for (uint j = i + 1; j < allFollowerLines.size(); j++) {
                const Line line1 = allFollowerLines.at(i);
                const Line line2 = allFollowerLines.at(j);                
                if (testUnordering(line1, line2)) {                    
                    Intersection unordering = Intersection(line1, line2);
                    followerLinesUnorderings.push_back(unordering);
                }                
            }
        }
        return followerLinesUnorderings;
    }

    

    static double computeLineError(const Line &line, const ILP_Chain &ilp_chain, ILP_Chain::iterator srcIt, ILP_Chain::iterator tgtIt) {
        assert(srcIt != tgtIt);
        assert(line.start.posInChain < line.end.posInChain);
        assert(ilp_chain.size() > 2);
        double maxError = 0;
        
        //calculate error for all nodes between src and tgt
        for (auto it = ++srcIt; it != tgtIt; it++) {
            double error = geo::calcPerpendicularLength(line.start.ch_node, line.end.ch_node, it->ch_node);
            if (error > maxError) {
                maxError = error;
            }
        }
        return maxError;
    }

    std::vector<Line> createPotentialEdges(ILP_Chain &ilp_chain, double epsilon) {

        std::vector<Line> lines;

        assert(ilp_chain.size() >= 2);
        for (auto srcIt = ilp_chain.begin(); srcIt != --ilp_chain.end(); srcIt++) {            
            auto incrSrcIt = srcIt;
            for (auto tgtIt = ++incrSrcIt; tgtIt != ilp_chain.end(); tgtIt++) {
                Line line(*srcIt, *tgtIt, nextLineID);
                if (computeLineError(line, ilp_chain, srcIt, tgtIt) < epsilon) {
                    lines.push_back(line);
                    nextLineID++;                    
                }
            }
            
        }
        return lines;
    }

    std::vector<Line> createFollowerLines(ILP_Chain &ilp_chainSrc, ILP_Chain& ilp_chainTgt, double eta) {
        std::vector<Line> followerLines;
        for (ILP_Node src : ilp_chainSrc) {
            bool isFollowing = false; //every node must have at least one other node which it could follow
            for (ILP_Node tgt : ilp_chainTgt) {                                
                if (geo::geoDist(src.ch_node, tgt.ch_node) < eta) {
                    Line line(src, tgt, nextLineID);
                    followerLines.push_back(line);
                    nextLineID++;
                    isFollowing = true;
                }
            }            
            assert(isFollowing);
            //followerIsPartial = !isFollowing;            
        }
        return followerLines;
    }    
    
    ILP_Chain ChainToILP_Chain(const Chain &chain, bool sameDirection) const {
        ILP_Chain ilp_chain;
        uint posInChain = 0;
        if (sameDirection) {            
            for (Chain::const_iterator it = chain.begin(); it != chain.end(); it++) {
                NodeID ch_node_id = *it;
                ILP_Node ilp_node = ILP_Node(graph.getNode(ch_node_id), ch_node_id, posInChain, true);
                ilp_chain.push_back(ilp_node);
                posInChain++;
            }
        } else {            
            for (Chain::const_reverse_iterator it = chain.rbegin(); it != chain.rend(); it++) {
                NodeID ch_node_id = *it;
                ILP_Node ilp_node = ILP_Node(graph.getNode(ch_node_id), ch_node_id, posInChain, false);
                ilp_chain.push_back(ilp_node);
                posInChain++;
            }
        }        
        return ilp_chain;

    }
    
    static std::vector<Line> concatLines(const std::vector<Line> &line1, const std::vector<Line> &line2)  {
        std::vector<Line> concatenatedLines;        
        for (Line line: line1) {
            concatenatedLines.push_back(line);
        }
        for (Line line: line2) {
            concatenatedLines.push_back(line);
        }
        return concatenatedLines;
    }
    
    ILP_data(CHGraph<CHNode, CHEdge> &graph, const Chain &chain1, const Chain &chain2, double epsilon, double eta, bool parallel):
        graph(graph), chain1(chain1), chain2(chain2) {
        
        ilp_chain1 = ChainToILP_Chain(chain1, true);                
        potEdges1 = createPotentialEdges(ilp_chain1, epsilon);
        
        if (parallel) {
            ilp_chain2 = ChainToILP_Chain(chain2, false);
            potEdges2 = createPotentialEdges(ilp_chain2, epsilon);
            followerLines1 = createFollowerLines(ilp_chain1, ilp_chain2, eta);
            followerLines2 = createFollowerLines(ilp_chain2, ilp_chain1, eta);            
            allFollowerLines = concatLines(followerLines1, followerLines2); 
            followerLinesUnorderings = calculate_p_followerUnorderings(allFollowerLines);            
        }
        allPotEdges = concatLines(potEdges1, potEdges2);
        edgeIntersections = calculateIntersections(allPotEdges);        
    }
};
