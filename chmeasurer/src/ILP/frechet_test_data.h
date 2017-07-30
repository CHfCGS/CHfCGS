#pragma once

#include "data.h"
#include "../geoFunctions.h"
#include "../nodes_and_edges.h"
#include "../chgraph.h"

#include <limits.h>

struct FrechetTest_data : Data
{
    const double LENGTH_FACTOR = 10e5;

    LineID nextLineID = 0;
    std::vector<Line> originalEdges1;
    std::vector<Line> originalEdges2;
    std::vector<Line> allOriginalEdges;
    std::vector<CrossLink> crossLinks1;
    std::vector<CrossLink> crossLinks2;
    std::vector<CrossLink> allCrossLinks;
    //std::vector<Intersection> edgeIntersections;
    std::vector<CrossLinkUnordering> crossLinksUnorderings;

    double max_link_length = 0;

    bool followerIsPartial = true; //if this is false, the structure is invalid

    void setMaxLinkLength()
    {
        for (const CrossLink& link : allCrossLinks)
        {
            max_link_length = std::max(max_link_length, link.length);
        }
    }

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

    static bool liesBetween(const uint pos, const Line& line)
    {
        return (pos > line.start.posInChain && pos < line.end.posInChain);
    }

    static bool liesPartiallyIn(const Line& line1, const Line& line2)
    {
        return liesBetween(line1.start.posInChain, line2)
                || liesBetween(line1.end.posInChain, line2);
    }

    static bool areOverlapping(const Line& line1, const Line& line2)
    {
        return liesPartiallyIn(line1, line2) || liesPartiallyIn(line2, line1);
    }

    static bool canCoexist(const CrossLink &link1, const CrossLink &link2)
    {
        //sources on same side
        if (link1.src.side == link2.src.side)
        {
            return !areOverlapping(link1.line, link2.line);
        } else
        { //different sides
            if (liesBetween(link1.src.posInChain, link2.line))
            {
                return false;
            } else if (liesBetween(link2.src.posInChain, link1.line))
            {
                return false;
            } else
            {
                return true;
            }
        }
    }

    static const uint NOPOS = UINT_MAX;

    struct continuousChainPos
    {
        uint nextpos;
        double signed_fraction;

        uint otherPos; //next from Fraction, needed because of unexact operations

        void init(uint pos)
        {
            this->nextpos = pos;
            signed_fraction = 0;
            otherPos = NOPOS;
        };

        void init(const CrossLink &link)
        {
            if (link.fraction < 0.5)
            {
                this->nextpos = link.line.start.posInChain;
                this->otherPos = link.line.end.posInChain;
                this->signed_fraction = link.fraction;
            } else
            {
                this->nextpos = link.line.end.posInChain;
                this->otherPos = link.line.start.posInChain;
                this->signed_fraction = link.fraction - 1;
            }
            assert(fabs(this->signed_fraction) <= 0.5);
        };

        bool operator<(const continuousChainPos &rhs) const
        {
            const double tolerance = 10000 * std::numeric_limits<double>::epsilon();
            //case: between the lines
            if (this->nextpos < rhs.nextpos)
            {
                if ((this->nextpos == rhs.otherPos) && (this->otherPos == rhs.nextpos))
                {
                    assert(this->signed_fraction >= 0);
                    assert(rhs.signed_fraction <= 0);
                    if (this->signed_fraction - rhs.signed_fraction >= 1 - tolerance)
                    {
                        return false;
                    } else
                    {
                        return true;
                    }
                } else
                {
                    return true;
                }
                //case: the positions lie around one node
            } else if (this->nextpos == rhs.nextpos)
            {
                if (this->signed_fraction < rhs.signed_fraction)
                {
                    //exception
                    if (fabs(this->signed_fraction) + fabs(rhs.signed_fraction) <= tolerance)
                    {
                        return false;
                    } else
                    {
                        return true;
                    }
                } else
                {
                    return false;
                }
            } else
            {
                return false;
            }
        }


        /*
        if (this->pos < rhs.pos) {
            return true;
        } else if (this->pos == rhs.pos) {
            return this->fraction < rhs.fraction;
        }
        else {
            return false;
        }*/

        /*
        bool operator >(const continuousChainPos &rhs) const {
            if (this->pos > rhs.pos) {
                return true;
            } else if (this->pos == rhs.pos) {
                return this->fraction > rhs.fraction;
            }
            else {
                return false;
            }
        }*/

        bool operator>(const continuousChainPos &rhs) const
        {
            return (rhs < *this);
        }
    };

    struct ChainPosPerLink
    {
        continuousChainPos chain1pos;
        continuousChainPos chain2pos;

        ChainPosPerLink(const CrossLink &link)
        {
            if (link.src.side)
            {
                chain1pos.init(link.src.posInChain);
                chain2pos.init(link);
            } else
            {
                chain1pos.init(link);
                chain2pos.init(link.src.posInChain);
            }
        }
        /*
        bool crosses(const ChainPosPerLink &rhs) const {
            return ((this->chain1pos > rhs.chain1pos && this->chain2pos < rhs.chain2pos)
                    || (this->chain1pos < rhs.chain1pos && this->chain2pos > rhs.chain2pos));
        }*/
    };

    static bool crosses(const ChainPosPerLink &lhs, const ChainPosPerLink &rhs)
    {
        if (lhs.chain1pos > rhs.chain1pos && lhs.chain2pos < rhs.chain2pos)
        {
            /*
            Print("lhs.chain1pos.pos" << lhs.chain1pos.nextpos);
            Print("lhs.chain1pos.fraction" << lhs.chain1pos.signed_fraction);
            Print("rhs.chain1pos.pos" << rhs.chain1pos.nextpos);
            Print("rhs.chain1pos.fraction" << rhs.chain1pos.signed_fraction);

            Print("lhs.chain2pos.pos" << lhs.chain2pos.nextpos);
            Print("lhs.chain2pos.fraction" << lhs.chain2pos.signed_fraction);
            Print("rhs.chain2pos.pos" << rhs.chain2pos.nextpos);
            Print("rhs.chain2pos.fraction" << rhs.chain2pos.signed_fraction);

            Print(">");
             */
            return true;
        } else if (lhs.chain1pos < rhs.chain1pos && lhs.chain2pos > rhs.chain2pos)
        {
            /*x
            Print("lhs.chain2pos.pos" << lhs.chain2pos.nextpos);
            Print("lhs.chain2pos.fraction" << lhs.chain2pos.signed_fraction);
            Print("rhs.chain2pos.pos" << rhs.chain2pos.nextpos);
            Print("rhs.chain2pos.fraction" << rhs.chain2pos.signed_fraction);
            Print("<");
             * */
            return true;
        } else
        {
            //Print("asdf");
            return false;
        }

        //return ((lhs.chain1pos > rhs.chain1pos && lhs.chain2pos < rhs.chain2pos)
        //      || (lhs.chain1pos < rhs.chain1pos && lhs.chain2pos > rhs.chain2pos));
    }

    static bool testUnordering(const CrossLink &link1, const CrossLink &link2)
    {
        assert(link1.line.start.side == link1.line.end.side);
        assert(link1.src.side != link1.line.start.side);
        assert(link2.line.start.side == link2.line.end.side);
        assert(link2.src.side != link2.line.start.side);


        if (canCoexist(link1, link2))
        {
            ChainPosPerLink cppl1(link1);
            ChainPosPerLink cppl2(link2);
            //return cppl1.crosses(cppl2);
            return crosses(cppl1, cppl2);
        } else
        {
            assert(false);
            return false;
        }
    }

    static std::vector<CrossLinkUnordering> calculate_crossLinkUnorderings(const std::vector<CrossLink> &allCrossLinksLines)
    {
        std::vector<CrossLinkUnordering> crossLinkUnorderings;
        assert(allCrossLinksLines.size() > 0);
        for (uint i = 0; i < allCrossLinksLines.size() - 1; i++)
        {
            for (uint j = i + 1; j < allCrossLinksLines.size(); j++)
            {
                const CrossLink& cross_link1 = allCrossLinksLines.at(i);
                const CrossLink& cross_link2 = allCrossLinksLines.at(j);
                if (testUnordering(cross_link1, cross_link2))
                {
                    CrossLinkUnordering unordering = CrossLinkUnordering(cross_link1, cross_link2);
                    crossLinkUnorderings.push_back(unordering);
                }
            }
        }
        return crossLinkUnorderings;
    }

    std::vector<Line> createOriginalEdges(const ILP_Chain &ilp_chain)
    {
        std::vector<Line> lines;

        assert(ilp_chain.size() >= 2);
        for (auto srcIt = ilp_chain.begin(); srcIt != --ilp_chain.end(); srcIt++)
        {
            auto nextIt = srcIt;
            nextIt++;
            Line line(*srcIt, *nextIt, nextLineID);
            lines.push_back(line);
            nextLineID++;
        }
        return lines;
    }

    std::vector<Line> createPotentialEdges(const ILP_Chain &ilp_chain, double epsilon)
    {
        std::vector<Line> lines;

        assert(ilp_chain.size() >= 2);
        for (auto srcIt = ilp_chain.begin(); srcIt != --ilp_chain.end(); srcIt++)
        {
            auto incrSrcIt = srcIt;
            for (auto tgtIt = ++incrSrcIt; tgtIt != ilp_chain.end(); tgtIt++)
            {
                Line line(*srcIt, *tgtIt, nextLineID);
                if (computeLineError(line, ilp_chain, srcIt, tgtIt) < epsilon)
                {
                    lines.push_back(line);
                    nextLineID++;
                }
            }

        }
        return lines;
    }

    std::vector<CrossLink> createCrossLinks(const ILP_Chain &ilp_chainSrc, const std::vector<Line>& edges_tgt, double eta)
    {
        std::vector<CrossLink> cross_links;
        for (ILP_Node src : ilp_chainSrc)
        {
            bool isFollowing = false; //every node must have at least one other edge to be linked with
            for (const Line& tgt_line : edges_tgt)
            {
                double geo_dist_start = geo::geoDist(src.ch_node, tgt_line.start.ch_node);
                double geo_dist_end = geo::geoDist(src.ch_node, tgt_line.end.ch_node);
                double perpendicular = geo::calcPerpendicularLength(tgt_line.start.ch_node, tgt_line.end.ch_node, src.ch_node);

                double fraction;
                double link_length;
                if (geo::isBetween(CHLine(tgt_line), graph.getNode(src.ch_node_id)))
                {
                    link_length = perpendicular;
                    double start_intersection_length = sqrt(pow(geo_dist_start, 2) - pow(perpendicular, 2));
                    double line_length = geo::geoDist(tgt_line.start.ch_node, tgt_line.end.ch_node);
                    assert(line_length != 0);
                    fraction = start_intersection_length / line_length;
                } else
                {
                    if (geo_dist_start < geo_dist_end)
                    {
                        link_length = geo_dist_start;
                        fraction = 0;
                    } else
                    {
                        link_length = geo_dist_end;
                        fraction = 1;
                    }
                }

                if (link_length <= eta)
                { //TODO remove
                    CrossLink link(src, tgt_line, fraction, LENGTH_FACTOR * link_length, nextLineID);
                    cross_links.push_back(link);
                    nextLineID++;
                    isFollowing = true;
                }
            }
            assert(isFollowing);
            //followerIsPartial = !isFollowing;
        }
        return cross_links;
    }

    std::vector<Line> createFollowerLines(const ILP_Chain &ilp_chainSrc, const ILP_Chain& ilp_chainTgt, double eta)
    {
        std::vector<Line> followerLines;
        for (ILP_Node src : ilp_chainSrc)
        {
            bool isFollowing = false; //every node must have at least one other node which it could follow
            for (ILP_Node tgt : ilp_chainTgt)
            {
                if (geo::geoDist(src.ch_node, tgt.ch_node) < eta)
                {
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

    static std::vector<CrossLink> concatLinks(const std::vector<CrossLink> &links1, const std::vector<CrossLink> &links2)
    {
        std::vector<CrossLink> concatenatedLinks;
        for (CrossLink link : links1)
        {
            concatenatedLinks.push_back(link);
        }
        for (CrossLink link : links2)
        {
            concatenatedLinks.push_back(link);
        }
        return concatenatedLinks;
    }

    FrechetTest_data(const CHGraph<CHNode, CHEdge> &graph, const Chain &chain1, const Chain &chain2, double epsilon, double eta, bool parallel)
    : Data(graph, chain1, chain2)
    {

        ilp_chain1 = ChainToILP_Chain(chain1, true);
        ilp_chain2 = ChainToILP_Chain(chain2, false);

        originalEdges1 = createOriginalEdges(ilp_chain1);
        originalEdges2 = createOriginalEdges(ilp_chain2);

        crossLinks1 = createCrossLinks(ilp_chain1, originalEdges2, eta);
        crossLinks2 = createCrossLinks(ilp_chain2, originalEdges1, eta);

        allCrossLinks = concatLinks(crossLinks1, crossLinks2);
        /*
        if (parallel) {

            potEdges2 = createPotentialEdges(ilp_chain2, epsilon);
            followerLines1 = createFollowerLines(ilp_chain1, ilp_chain2, eta);
            followerLines2 = createFollowerLines(ilp_chain2, ilp_chain1, eta);
            allFollowerLines = concatLines(followerLines1, followerLines2);
            followerLinesUnorderings = calculate_p_followerUnorderings(allFollowerLines);
        }*/
        crossLinksUnorderings = calculate_crossLinkUnorderings(allCrossLinks);

        allOriginalEdges = concatLines(originalEdges1, originalEdges2);
        setMaxLinkLength();
        //edgeIntersections = calculateIntersections(allPotEdges);
    }
};
