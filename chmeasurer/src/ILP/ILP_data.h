#pragma once

#include "data.h"
#include "../geoFunctions.h"
#include "../nodes_and_edges.h"
#include "../chgraph.h"

struct ILP_data : Data
{
    std::vector<Line> followerLines1;
    std::vector<Line> followerLines2;
    std::vector<Line> allFollowerLines;

    std::vector<Intersection> followerLinesUnorderings;

    //bool followerIsPartial = true; //if this is false, the structure is invalid

    static std::vector<Intersection> calculateIntersections(const std::vector<Line> &lines)
    {
        //possible upgrade: Bentley-Ottmann
        //possible upgrade: could be adapted for 2 input vectors
        std::vector<Intersection> intersections;
        if (lines.size() > 1)
        {
            for (uint i = 0; i < lines.size() - 1; i++)
            {
                for (uint j = i + 1; j < lines.size(); j++)
                {
                    if (testIntersection(lines.at(i), lines.at(j)))
                    {
                        Intersection intersection(lines.at(i), lines.at(j));
                        intersections.push_back(intersection);
                    }
                }
            }
        }
        return intersections;
    }

    static bool testUnordering(const Line &line1, const Line &line2)
    {
        uint line1Chain1Pos, line1Chain2Pos, line2Chain1Pos, line2Chain2Pos;
        assert(line1.start.side != line1.end.side);
        if (line1.start.side)
        {
            line1Chain1Pos = line1.start.posInChain;
            line1Chain2Pos = line1.end.posInChain;
        } else
        {
            line1Chain1Pos = line1.end.posInChain;
            line1Chain2Pos = line1.start.posInChain;
        }

        assert(line2.start.side != line2.end.side);
        if (line2.start.side)
        {
            line2Chain1Pos = line2.start.posInChain;
            line2Chain2Pos = line2.end.posInChain;
        } else
        {
            line2Chain1Pos = line2.end.posInChain;
            line2Chain2Pos = line2.start.posInChain;
        }


        if ((line1Chain1Pos > line2Chain1Pos && line1Chain2Pos < line2Chain2Pos)
                || (line1Chain1Pos < line2Chain1Pos && line1Chain2Pos > line2Chain2Pos))
        {
            return true;
        }
        return false;
    }

    static std::vector<Intersection> calculate_p_followerUnorderings(const std::vector<Line> &allFollowerLines)
    {
        std::vector<Intersection> followerLinesUnorderings;

        assert(allFollowerLines.size() > 0);
        for (uint i = 0; i < allFollowerLines.size() - 1; i++)
        {
            for (uint j = i + 1; j < allFollowerLines.size(); j++)
            {
                const Line line1 = allFollowerLines.at(i);
                const Line line2 = allFollowerLines.at(j);
                if (testUnordering(line1, line2))
                {
                    Intersection unordering = Intersection(line1, line2);
                    followerLinesUnorderings.push_back(unordering);
                }
            }
        }
        return followerLinesUnorderings;
    }

    std::vector<Line> createFollowerLines(const ILP_Chain &ilp_chainSrc, const ILP_Chain& ilp_chainTgt, double eta)
    {
        std::vector<Line> followerLines;
        for (ILP_Node src : ilp_chainSrc)
        {
            //bool isFollowing = false; //every node must have at least one other node which it could follow
            // but we already have an eta
            for (ILP_Node tgt : ilp_chainTgt)
            {
                if (geo::geoDist(src.ch_node, tgt.ch_node) < eta)
                {
                    Line line(src, tgt, nextLineID);
                    followerLines.push_back(line);
                    nextLineID++;
                    //                    isFollowing = true;
                }
            }
            //assert(isFollowing);
            //followerIsPartial = !isFollowing;
        }
        return followerLines;
    }

    ILP_data(const CHGraph<CHNode, CHEdge> &graph, const Chain &chain1, const Chain &chain2, double epsilon, double eta, bool parallel)
    : Data(graph, chain1, chain2)
    {

        ilp_chain1 = ChainToILP_Chain(chain1, true);
        potEdges1 = createPotentialEdges(ilp_chain1, epsilon);

        if (parallel)
        {
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
