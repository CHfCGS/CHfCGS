#pragma once

#include <glpk.h>
#include <vector>             /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <sstream>

#include "frechet_test_data.h"

class PrllelLineSimplificFrechet : FrechetLP
{

    void setAllDegreeCoefficients(const FrechetTest_data &ilp_data)
    {
        size_t nofNodes = ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size();
        size_t chainPosOffset = ilp_data.ilp_chain1.size();

        for (const Line &line : ilp_data.potEdges1)
        {
            assert(line.start.ch_node_id != line.end.ch_node_id);
            setDegreeCoefficient(line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(line.end.posInChain + 1, line.id + 1, 1.0);
            setDegreeCoefficient(nofNodes + line.start.posInChain + 1, line.id + 1, -1.0);
        }

        for (const Line &line : ilp_data.potEdges2)
        {
            assert(line.start.ch_node_id != line.end.ch_node_id);
            setDegreeCoefficient(chainPosOffset + line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(chainPosOffset + line.end.posInChain + 1, line.id + 1, 1.0);
            setDegreeCoefficient(nofNodes + chainPosOffset + line.start.posInChain + 1, line.id + 1, -1.0);
        }

        for (const CrossLink &crossLink : ilp_data.crossLinks1)
        {
            //out
            setDegreeCoefficient(nofNodes + crossLink.src.posInChain + 1, crossLink.line.id + 1, 1.0);
        }

        for (const CrossLink &crossLink : ilp_data.crossLinks2)
        {
            //out
            setDegreeCoefficient(nofNodes + chainPosOffset + crossLink.src.posInChain + 1, crossLink.line.id + 1, 1.0);
        }
    }

    void setCrossLinksExistenceCoefficients(const FrechetTest_data &ilp_data)
    {
        size_t nofNodes = ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size();

        for (const CrossLink &crossLink : ilp_data.allCrossLinks)
        {
            setDegreeCoefficient(nofNodes + crossLink.id, crossLink.line.id, -1);
            setDegreeCoefficient(nofNodes + crossLink.id, crossLink.id, 1);
        }
    }

    void addInOutDegreeRows(const ILP_Chain& ilp_chain)
    {
        assert(ilp_chain.size() > 1);
        //nodes 0 and n-1 should have degree 1
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.front().ch_node_id, 1);
        //glp_set_row_name(lp, glp_get_num_rows(lp), "node: 0");
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, -1.0, 0);
        //setInOutDegreeCoefficients(lines, ilp_chain.front().ch_node_id);

        //all other nodes should have indegree = outdegree
        //for (uint node_id = 1; node_id < graph.nodes.size()-1; node_id++) {
        for (auto it = ++ilp_chain.begin(); it != --ilp_chain.end(); it++)
        {
            glp_add_rows(lp, 1);
            setRowName(it->ch_node_id, 1);
            //glp_set_row_name(lp, glp_get_num_rows(lp), row_name);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);
            //setInOutDegreeCoefficients(lines, node_id);
        }


        glp_add_rows(lp, 1);
        setRowName(ilp_chain.back().ch_node_id, 1);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);
        //setInOutDegreeCoefficients(lines, ilp_chain.back().ch_node_id);
    }


    //outdegrees has too be equal except for ends

    void addOutDegreeRows(const ILP_Chain& ilp_chain)
    {
        assert(ilp_chain.size() > 1);
        //all other nodes should have edge_outdegree = follower_outdegree
        for (auto it = ilp_chain.begin(); it != --ilp_chain.end(); it++)
        {
            glp_add_rows(lp, 1);
            setRowName(it->ch_node_id, 2);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 0.0, 0.0);
        }

        //except the last one
        glp_add_rows(lp, 1);
        setRowName(ilp_chain.back().ch_node_id, 2);
        glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0);
    }

    void addCrossLinksExistenceRows(const FrechetTest_data &ilp_data)
    {
        for (auto it = ++ilp_data.allPotEdges.begin(); it != ilp_data.allPotEdges.end(); it++)
        {
            glp_add_rows(lp, 1);
            setRowName(it->id, 3);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_LO, 0.0, 0.0);
        }
    }

    void addAllDegreeRows(FrechetTest_data &ilp_data) override
    {
        addInOutDegreeRows(ilp_data.ilp_chain1);
        addInOutDegreeRows(ilp_data.ilp_chain2);
        addOutDegreeRows(ilp_data.ilp_chain1);
        addOutDegreeRows(ilp_data.ilp_chain2);
    }

    void addAllIntersectionRows(const FrechetTest_data &ilp_data) override
    {
        addIntersectionRows(ilp_data.edgeIntersections);
        addLinkUnorderingRows(ilp_data.crossLinksUnorderings);
    }

    void addColumns(const std::vector<Line> &lines, double obj_coef)
    {
        debug_assert(lines.size() > 0);
        glp_add_cols(lp, lines.size());
        for (uint i = 0; i < lines.size(); i++)
        {
            const Line &line = lines.at(i);
            uint col_id = line.id + 1;
            std::string colType = obj_coef == 1 ? "Edge" : "followerLine";
            std::stringstream ss("");
            ss << colType << ": (" << line.start.ch_node_id << "," << line.end.ch_node_id << ")";
            std::string s = ss.str();
            char const *col_name = s.c_str();
            glp_set_col_name(lp, col_id, col_name);
            glp_set_col_kind(lp, col_id, GLP_BV);
            glp_set_obj_coef(lp, col_id, obj_coef); //objective is number of used edges
        }
    }

    void addCrossLinkColumns(const std::vector<CrossLink> &crossLinks, double obj_coef)
    {
        debug_assert(crossLinks.size() > 0);
        glp_add_cols(lp, crossLinks.size());
        for (uint i = 0; i < crossLinks.size(); i++)
        {
            const CrossLink &crossLink = crossLinks.at(i);
            uint col_id = crossLink.id + 1;
            std::string colType = "CrossLink";
            std::stringstream ss("");
            ss << colType << ": (" << crossLink.src.ch_node_id << "->("
                    << crossLink.line.end.ch_node_id << "," << crossLink.line.end.ch_node_id << "))";
            std::string s = ss.str();
            char const *col_name = s.c_str();
            glp_set_col_name(lp, col_id, col_name);
            glp_set_col_kind(lp, col_id, GLP_BV);
            glp_set_obj_coef(lp, col_id, obj_coef); //objective is number of used edges
        }
    }

    void setColumns(const FrechetTest_data &ilp_data)
    {
        addColumns(ilp_data.potEdges1, 1.0);
        addColumns(ilp_data.potEdges2, 1.0);
        addCrossLinkColumns(ilp_data.crossLinks1, 0);
        addCrossLinkColumns(ilp_data.crossLinks1, 0);
    }

    void prepareArray(const FrechetTest_data &ilp_data)
    {
        assert(ilp_data.ilp_chain1.size() > 1);
        assert(ilp_data.ilp_chain2.size() > 1);

        preNofCols = ilp_data.potEdges1.size() + ilp_data.potEdges2.size()
                + ilp_data.crossLinks1.size() + ilp_data.crossLinks2.size();
        preNofRows = 2 * (ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size()) //2 degree equations for each node
                + ilp_data.allCrossLinks.size()
                + ilp_data.edgeIntersections.size()
                + ilp_data.crossLinksUnorderings.size();
        preNofNonZeros = 3 * (ilp_data.potEdges1.size() + ilp_data.potEdges2.size())
                + 3 * ilp_data.allCrossLinks.size()
                + 2 * ilp_data.edgeIntersections.size()
                + 2 * ilp_data.crossLinksUnorderings.size();
        nofNonZeros = 0;
        enoughSpace = preNofNonZeros < (int) size;
        assert(enoughSpace);
    }


public:

    PrllelLineSimplificFrechet(const CHGraph<CHNode, CHEdge> &graph) : FrechetLP(graph) { }

    ~PrllelLineSimplificFrechet() { }

    double simplify(const Chain& chain1, const Chain& chain2, double epsilon, double eta)
    {

        FrechetTest_data ilp_data = FrechetTest_data(graph, chain1, chain2, 0, eta, false);
        prepareArray(ilp_data);
        lp = glp_create_prob();
        glp_term_out(GLP_OFF); //make glp not verbose

        glp_set_prob_name(lp, "parallel line generalization");

        glp_set_obj_dir(lp, GLP_MIN); //minimize number of used edges

        setColumns(ilp_data);

        addCrossLinksExistenceRows(ilp_data);
        setCrossLinksExistenceCoefficients(ilp_data);

        addAllDegreeRows(ilp_data);
        setAllDegreeCoefficients(ilp_data);

        addAllIntersectionRows(ilp_data);

        assert(preNofRows == glp_get_num_rows(lp));
        assert(preNofNonZeros == (int) nofNonZeros);

        return solve();
    }
};