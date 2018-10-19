#include "deltaPhi.h"
#include <cmath>

float deltaR2(float e1, float p1, float e2, float p2) {

    float de = e1 - e2;
    float dp = deltaPhi(p1,p2);
    return de*de+dp*dp;
}

float deltaR(float e1, float p1, float e2, float p2) {
    return sqrt(deltaR2(e1,p1,e2,p2));
}

