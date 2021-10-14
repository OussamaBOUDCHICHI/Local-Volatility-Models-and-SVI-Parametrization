#include "../include/BlackScholesModel.hpp"
#include <iostream>

using namespace std;



int main() {

    float S_0 ;
    float r ;
    float T ;
    float t_0 ;
    float K ;
    float sigma ;
    char flag ;
    double P_market;
    std::string method;
    double sig_0(0.2);
    int maxIter(50);
    

    cout << "Enter S_0 : " << endl;
    cin >> S_0;

    cout << "Enter r : " << endl;
    cin >> r;

    cout << "Enter T : " << endl;
    cin >> T;

    cout << "Enter t_0 : " << endl;
    cin >> t_0;

    cout << "Enter K : " << endl;
    cin >> K;

    cout << "Enter sigma : " << endl;
    cin >> sigma;

    cout << "Enter Vanilla type (C / P) : " << endl;
    cin >> flag;

    cout << "Enter Method : " << endl;
    cin >> method;

    cout << "Enter Market : " << endl;
    cin >> P_market;

    
    BlackScholesModel BS(S_0, t_0, T, r, K, sigma, flag);
    double C = BS.price();
    
    double Delta = BlackScholesModel::OptionGreeks().compDelta(S_0, t_0, T, r, K, sigma, flag);
    double vol = BS.impVolatility(sig_0, P_market, maxIter, method);
    cout << "C : " << C << endl;
    cout << "sigma : " << vol << endl;
    cout << "Delta : " << Delta << endl;
    

    return 0;
}