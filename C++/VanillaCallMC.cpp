#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>

using namespace std;

float stoploss(double x, double  K){
    return max(x - K, 0.0);
}

int main() {

        mt19937 G(time(NULL));
        normal_distribution <double> N(0,1);


        float T ;
        int n = 10000;
        float K ;
        float r ;
        float S ;
        float sigma ;

        cout << "Enter T : " << endl;
        cin >> T;

        cout << "Enter K : " << endl;
        cin >> K;

        cout << "Enter r : " << endl;
        cin >> r;

        cout << "Enter S_0 : " << endl;
        cin >> S;

        cout << "Enter sigma : " << endl;
        cin >> sigma;



        float S1 = S * exp((r - 0.5 * (sigma*sigma)) * T);
        float sum = 0.0;


        for (int i=0; i < n; i++){
            float S2 = S1*exp(sigma *  sqrt(T) * N(G));
            sum += stoploss(S2, K);
        }

        float premium = exp(-r*T) * (sum /(float) n);

cout << "The price of the vanilla is :  " << premium << endl;

}