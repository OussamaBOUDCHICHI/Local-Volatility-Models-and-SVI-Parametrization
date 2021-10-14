#ifndef BLACKSCHOLES_H
#define BLACKSCHOLES_H

#include <cmath>
#include <iostream>
#include <string>
/*************************
BS Model Class
**************************/

class BlackScholesModel {
    private:
        double r, T, S_0, sigma, q, K, t_0;
        char type ;

    public:
        BlackScholesModel();

        BlackScholesModel(double S_0, double t_0, double T, double r, double K, 
               double sigma, char type, double q = 0.0):
        S_0(S_0), t_0(t_0), T(T), r(r), K(K), sigma(sigma), q(q), type(type) {}

        virtual ~BlackScholesModel() {}
        friend class OptionGreeks;
        friend class Normal;

        double price() const;
        double impVolatility(double sigma_0, double P_market, int maxIter, std::string method = "N-R", double a = 0.0001, double b = 2.0, double tolerance=1e-5) const;

    // utils class

    class Normal {
        private:
            double x;
        
        public:
        Normal() {}
        virtual ~Normal() {}
        Normal(double x) : x(x) {}

        double CDF() {
            return 0.5 * erfc(-x * M_SQRT1_2);
        }

        double PDF() {
            return 0.5 * M_2_SQRTPI * M_SQRT1_2 * exp(- 0.5 * x * x);
        }
    };
    
    //Option Greeks::

    class OptionGreeks {

        private:
            double delta, gamma, theta, vega, rho;
        
        public:
            OptionGreeks() {}
            double compVega(double S_0, double K, double r, double sigma,
            double T, double t_0, char type, double q = 0.0);

            double compDelta(double S_0, double K, double r, double sigma,
            double T, double t_0, char type, double q = 0.0);
            
            double compGamma(double S_0, double K, double r, double sigma,
            double T, double t_0, char type , double q = 0.0);
            
            double compRho(double S_0, double K, double r, double sigma,
            double T, double t_0, char type,  double q = 0.0);
            
            double compTheta(double S_0, double K, double r, double sigma,
            double T, double t_0, char type, double q = 0.0);

    };

};


#endif