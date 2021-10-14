#ifndef CRR_H
#define CRR_H

#include <iostream>
#include <cmath>
#include <vector>
class CRR {
    
    private:

    float S_0;
    float u;
    float d;
    float r;

    public:
    CRR(){
        S_0 = 1;
        u = 0.5;
        d = - 0.5;
        r = 0.0;

    }

    CRR(float& price, float& up, float& down, float& freerate);

    int CheckInput() const ;

    float GetRiskNeutral() const;

    float S(int & n, int & i) const;

    float Payoff(int& T, float& K, const char& O) const;

    float Getr() const;

};

#endif