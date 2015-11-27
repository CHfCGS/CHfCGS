#pragma once

#include "../nodes_and_edges.h"

struct Kink {
    const chc::NodeID src;
    const chc::NodeID peak;
    const chc::NodeID tgt;
    
};

class DiscreteCurveEvolution {
    double turnAngle(const Kink kink);
};