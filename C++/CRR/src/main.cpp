#include "../include/CRR.hpp"
#include <iostream>

using namespace std;

int main(){
    float K;
    int T;
    float u, d, r, S_0;

    cout << "Enter Strike K : " << endl;
    cin >> K;
    cout << "Enter Maturity T :" << endl;
    cin >> T;
    cout << "Enter  u : " << endl;
    cin >> u;
    cout << "Enter d : " << endl;
    cin >> d;
    cout << "Enter r : " << endl;
    cin >> r;
    cout << "Enter S_0 : " << endl;
    cin >> S_0;

    CRR mdl(S_0, u, d, r);

    cout << mdl.CheckInput() << endl;
    cout << "Risk-Neutral probability : " << mdl.GetRiskNeutral() << endl;

    cout << "Call Payoff : " << mdl.Payoff(T, K, 'C') << endl;
    cout << "Put Payoff : " << mdl.Payoff(T, K, 'P') << endl;


    return 0;

}