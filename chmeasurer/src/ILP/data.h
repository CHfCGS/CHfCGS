/*
 * File:   data.h
 * Author: tobias
 *
 * Created on July 30, 2017, 1:20 PM
 */

#ifndef DATA_H
#define	DATA_H

struct Data
{
    const CHGraph<CHNode, CHEdge> &graph;
    const Chain &chain1;
    const Chain &chain2;
    ILP_Chain ilp_chain1;
    ILP_Chain ilp_chain2;

    static bool testIntersection(const Line line1, const Line line2)
    {

        //4 Orientation tests
        double line2StartArea = geo::calcSignedArea(line1.start.ch_node, line1.end.ch_node, line2.start.ch_node);
        double line2EndArea = geo::calcSignedArea(line1.start.ch_node, line1.end.ch_node, line2.end.ch_node);
        bool line1BetweenLine2 = geo::differentSign(line2StartArea, line2EndArea);

        double line1StartArea = geo::calcSignedArea(line2.start.ch_node, line2.end.ch_node, line1.start.ch_node);
        double line1EndArea = geo::calcSignedArea(line2.start.ch_node, line2.end.ch_node, line1.end.ch_node);
        bool line2BetweenLine1 = geo::differentSign(line1StartArea, line1EndArea);

        return line1BetweenLine2 && line2BetweenLine1;
    }

    static double computeLineError(const Line &line, const ILP_Chain &ilp_chain, ILP_Chain::const_iterator srcIt, ILP_Chain::const_iterator tgtIt)
    {
        assert(srcIt != tgtIt);
        assert(line.start.posInChain < line.end.posInChain);
        assert(ilp_chain.size() >= 2);
        double maxError = 0;

        //calculate error for all nodes between src and tgt
        for (auto it = ++srcIt; it != tgtIt; it++)
        {
            double error = geo::calcPerpendicularLength(line.start.ch_node, line.end.ch_node, it->ch_node);
            if (error > maxError)
            {
                maxError = error;
            }
        }
        return maxError;
    }

    ILP_Chain ChainToILP_Chain(const Chain &chain, bool sameDirection) const
    {
        ILP_Chain ilp_chain;
        uint posInChain = 0;
        if (sameDirection)
        {
            for (Chain::const_iterator it = chain.begin(); it != chain.end(); it++)
            {
                NodeID ch_node_id = *it;
                ILP_Node ilp_node = ILP_Node(graph.getNode(ch_node_id), ch_node_id, posInChain, true);
                ilp_chain.push_back(ilp_node);
                posInChain++;
            }
        } else
        {
            for (Chain::const_reverse_iterator it = chain.rbegin(); it != chain.rend(); it++)
            {
                NodeID ch_node_id = *it;
                ILP_Node ilp_node = ILP_Node(graph.getNode(ch_node_id), ch_node_id, posInChain, false);
                ilp_chain.push_back(ilp_node);
                posInChain++;
            }
        }
        return ilp_chain;
    }

    static std::vector<Line> concatLines(const std::vector<Line> &line1, const std::vector<Line> &line2)
    {
        std::vector<Line> concatenatedLines;
        for (Line line : line1)
        {
            concatenatedLines.push_back(line);
        }
        for (Line line : line2)
        {
            concatenatedLines.push_back(line);
        }
        return concatenatedLines;
    }

    Data(const CHGraph<CHNode, CHEdge> &graph, const Chain &chain1, const Chain &chain2)
    : graph(graph), chain1(chain1), chain2(chain2) { }

    virtual ~Data() { }

};


#endif	/* DATA_H */

