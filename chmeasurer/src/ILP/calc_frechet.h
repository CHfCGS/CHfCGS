#pragma once

#include <glpk.h>
#include <vector>             /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <sstream>

#include "frechetLP.h"

class CalcFrechetILP : public FrechetLP
{
protected:

    void setAllDegreeCoefficients(const FrechetTest_data& ilp_data)
    {

        //size_t nofNodes = ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size();
        size_t chainPosOffset = ilp_data.ilp_chain1.size();
        //const double INDEGREEMULTIPLIER = 1000.0;

        /*
        for (Line &line : ilp_data.potEdges1) {
            assert(line.start.ch_node_id != line.end.ch_node_id);
            setDegreeCoefficient(line.start.posInChain + 1, line.id + 1 , -1.0);
            setDegreeCoefficient(line.end.posInChain + 1, line.id + 1 , 1.0);
            setDegreeCoefficient(nofNodes + line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(2*nofNodes + line.end.posInChain + 1, line.id + 1, INDEGREEMULTIPLIER);
        }

        for (Line &line : ilp_data.potEdges2) {
            assert(line.start.ch_node_id != line.end.ch_node_id);
            setDegreeCoefficient(chainPosOffset + line.start.posInChain + 1, line.id + 1 , -1.0);
            setDegreeCoefficient(chainPosOffset + line.end.posInChain + 1, line.id + 1 , 1.0);
            setDegreeCoefficient(nofNodes + chainPosOffset + line.start.posInChain + 1, line.id + 1, -1.0);
            setDegreeCoefficient(2*nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1, INDEGREEMULTIPLIER);
        }*/

        for (const CrossLink &link : ilp_data.crossLinks1)
        {
            //assert(link.start.ch_node_id != link.end.ch_node_id);
            //out
            setDegreeCoefficient(link.src.posInChain + 1, link.id + 1, 1.0);
            //in
            //setDegreeCoefficient(2*nofNodes + chainPosOffset + line.end.posInChain + 1, line.id + 1, -1.0);
        }

        for (const CrossLink &link : ilp_data.crossLinks2)
        {
            //assert(line.start.ch_node_id != line.end.ch_node_id);
            //out
            setDegreeCoefficient(chainPosOffset + link.src.posInChain + 1, link.id + 1, 1.0);
            //in
            //setDegreeCoefficient(2*nofNodes + line.end.posInChain + 1, line.id + 1, -1.0);
        }

    }

    void addOutDegreeRows(ILP_Chain ilp_chain)
    {
        assert(ilp_chain.size() > 1);
        for (auto it = ilp_chain.begin(); it != ilp_chain.end(); it++)
        {
            glp_add_rows(lp, 1);
            setRowName(it->ch_node_id, 2);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_FX, 1.0, 0.0);
        }
    }

    void addAllDegreeRows(FrechetTest_data &ilp_data) override
    {
        //addInOutDegreeRows(ilp_data.ilp_chain1);
        //addInOutDegreeRows(ilp_data.ilp_chain2);
        addOutDegreeRows(ilp_data.ilp_chain1);
        addOutDegreeRows(ilp_data.ilp_chain2);
        //addInDegreeRows(ilp_data.ilp_chain1);
        //addInDegreeRows(ilp_data.ilp_chain2);
    }

    void setLimitEntries(const CrossLink& link)
    {
        //entry for link
        nofNonZeros++;
        size_t index = nofNonZeros;
        ia[index] = glp_get_num_rows(lp);
        ja[index] = link.id + 1;
        ar[index] = link.length;

        //entry for limit
        nofNonZeros++;
        index = nofNonZeros;
        ia[index] = glp_get_num_rows(lp);
        ja[index] = glp_get_num_cols(lp);
        ar[index] = -1.0;
    }

    void addAllIntersectionRows(const FrechetTest_data &ilp_data) override
    {
        //addIntersectionRows(ilp_data.edgeIntersections);
        addLinkUnorderingRows(ilp_data.crossLinksUnorderings);
    }

    void addLimitRows(const FrechetTest_data &ilp_data)
    {
        for (const CrossLink &link : ilp_data.allCrossLinks)
        {
            glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, glp_get_num_rows(lp), GLP_UP, 0.0, 0.0);

            setLimitEntries(link);

            std::stringstream ss("");
            ss << "Limit: (" << link.src.ch_node_id << "->" << link.line.start.ch_node_id << "," << link.line.end.ch_node_id
                    << "length: " << link.length << "),";
            std::string s = ss.str();
            char const *row_name = s.c_str();
            glp_set_row_name(lp, glp_get_num_rows(lp), row_name);
            //assert(line1.start.ch_node_id != line2.end.ch_node_id);
        }
    }

    //is the last column

    void addLimitColumn(uint nextNodeID, double max_link_length)
    {
        glp_add_cols(lp, 1);
        uint col_id = glp_get_num_cols(lp);
        assert(col_id == nextNodeID + 1);
        //CrossLink &link = links.at(i);
        //uint col_id = line.id + 1;
        std::string colType = "Limit"; // obj_coef == 1 ? "Edge" : "followerLine";
        std::stringstream ss("");
        //ss << colType << ": (" << link.src.ch_node_id << "," << link.line.start.ch_node_id << "," << link.line.end.ch_node_id << ")";
        ss << colType;
        std::string s = ss.str();
        char const *col_name = s.c_str();
        glp_set_col_name(lp, col_id, col_name);
        glp_set_col_kind(lp, col_id, GLP_CV);
        glp_set_col_bnds(lp, col_id, GLP_DB, 0.0, max_link_length);
        glp_set_obj_coef(lp, col_id, 1); //objective is max length
    }

    void setColumns(const FrechetTest_data &ilp_data) override
    {
        addColumns(ilp_data.potEdges1, 0);
        addColumns(ilp_data.potEdges2, 0);
        addLinkColumns(ilp_data.crossLinks1);
        addLinkColumns(ilp_data.crossLinks2);
        addLimitColumn(ilp_data.nextLineID, ilp_data.max_link_length);
    }

    void prepareArray(const FrechetTest_data &ilp_data)
    {
        assert(ilp_data.ilp_chain1.size() > 1);
        assert(ilp_data.ilp_chain2.size() > 1);

        preNofCols = ilp_data.allPotEdges.size()
                + ilp_data.crossLinks1.size() + ilp_data.crossLinks2.size()
                + 1; //limit
        preNofRows = ilp_data.ilp_chain1.size() + ilp_data.ilp_chain2.size() //out degree
                + ilp_data.crossLinksUnorderings.size()
                + ilp_data.allCrossLinks.size(); //limits
        preNofNonZeros = 3 * ilp_data.allCrossLinks.size() //one for outgoing, two for limitation row
                + 2 * ilp_data.crossLinksUnorderings.size();
        nofNonZeros = 0;
        enoughSpace = preNofNonZeros < (int) size;
        assert(enoughSpace);
    }


public:

    CalcFrechetILP(const CHGraph<CHNode, CHEdge> &graph) : FrechetLP(graph) { }

    ~CalcFrechetILP() { }

    double simplify(const Chain& chain1, const Chain& chain2, double eta)
    {

        FrechetTest_data ilp_data = FrechetTest_data(graph, chain1, chain2, 0, eta, false);
        prepareArray(ilp_data);
        lp = glp_create_prob();
        glp_term_out(GLP_ON); //make glp not verbose

        glp_set_prob_name(lp, "calc frechet");

        glp_set_obj_dir(lp, GLP_MIN); //minimize length of longest cross link

        setColumns(ilp_data);

        addAllDegreeRows(ilp_data);
        setAllDegreeCoefficients(ilp_data);

        addAllIntersectionRows(ilp_data);
        addLimitRows(ilp_data);

        assert(preNofCols == glp_get_num_cols(lp));
        assert(preNofRows == glp_get_num_rows(lp));
        assert(preNofNonZeros == (int) nofNonZeros);

        return solve();
    }
};