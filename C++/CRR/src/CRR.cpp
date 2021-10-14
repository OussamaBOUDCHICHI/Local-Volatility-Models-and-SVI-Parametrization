#include "../include/CRR.hpp"

using namespace std;

CRR::CRR(float& price, float& up, float& down, float& freerate){
    S_0 = price;
    u = up;
    d = down;
    r = freerate;
}

int CRR::CheckInput() const {
    if (S_0 <= 0 || u <= -1.0 || d <= -1.0 || u <= d || r <= -1.0)
    {
        cout << "Illegal data ranges" << endl;
        cout << "Terminating program" << endl;
        return 1;
    }

    // Checking for arbitrage

    if (r >= u || r <= d)
    {
        cout << "Arbitrage Exists" << endl;
        cout << "Terminating program" << endl;
        return 1;
    }
    cout << "Input data checked" << endl;
    cout << "There's no arbitrage, all clear !" << endl;
    return 0;
}

float CRR::GetRiskNeutral() const {

    return (r-d)/(u-d);
}

float CRR::S(int& n, int& i) const {
    return S_0*pow(1+u, i)*pow(1+d, n-i);

}

float CRR::Payoff(int& T, float& K, const char& O) const {
    float q = this->GetRiskNeutral();
    vector <double> ArbitPrice(T);
    for(int i = 0; i <= T; i++){
        if(O == 'C') {
        ArbitPrice.push_back(max(this->S(T,i) - K, float(0)));
    }
    if(O == 'P') {
        ArbitPrice.push_back(max(K - this->S(T,i), float(0)));
    }
    }


    for (int i=T-1; i >=0; i--){

        for (int j=0; j <= i; j++){
            ArbitPrice[j] = (q*ArbitPrice[j+1] + (1-q)*ArbitPrice[j]) / (1 + r);
        }

    }
    
    return ArbitPrice[0];

}

float CRR::Getr() const {
    return r;
}

