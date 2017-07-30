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

    Data(const CHGraph<CHNode, CHEdge> &graph, const Chain &chain1, const Chain &chain2)
    : graph(graph), chain1(chain1), chain2(chain2) { }

    virtual ~Data() { }

};


#endif	/* DATA_H */

