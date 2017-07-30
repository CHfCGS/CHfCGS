/*
 * File:  calculation.h
 * Author: tobias
 *
 * Created on July 30, 2017, 1:23 PM
 */

#ifndef CALCULATION_H
#define	CALCULATION_H

class Calculation
{



protected:
    const CHGraph<CHNode, CHEdge> &graph;

    glp_prob *lp;

    const static size_t size = 4000000;
    size_t nofNonZeros = 0;

    //numbers of needed sizes beforehand
    int preNofRows;
    int preNofCols;
    int preNofNonZeros = 0;

    //whole ILP is invalid if size is not big enough
    bool enoughSpace;

    std::vector<int> ia;
    std::vector<int> ja;
    std::vector<double> ar;

    double objective_value = std::numeric_limits<double>::max();

    bool isInRange(size_t index)
    {
        return (1 <= index && index <= size); //0 is not valid for ILP
    }

    void setRowName(NodeID node_id, uint nr)
    {
        std::stringstream ss("");
        ss << "node: " << nr << ". " << node_id;
        std::string s = ss.str();
        char const *row_name = s.c_str();
        glp_set_row_name(lp, glp_get_num_rows(lp), row_name);
    }

    Calculation(const CHGraph<CHNode, CHEdge> &graph) : graph(graph)
    {
        ia.resize(size + 1);
        ja.resize(size + 1);
        ar.resize(size + 1);
    }

    virtual ~Calculation() { }
public:

};


#endif	/* CALCULATION_H */

