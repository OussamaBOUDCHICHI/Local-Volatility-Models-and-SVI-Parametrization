#include "../include/BlackScholesModel.hpp"

using namespace std;

BlackScholesModel::BlackScholesModel() : S_0(100.0), t_0(0.0), T(0.5), r(0.02), K(80.0), sigma(0.2), q(0.0), type('C') {}


// Delta
double BlackScholesModel::OptionGreeks::compDelta(double S_0, double K, double r, double sigma,
            double T, double t_0, char type, double q) {

                double d1 = (log(S_0 / K) + (r + 0.5 * sigma * sigma) * (T - t_0)) / (sigma * sqrt(T - t_0));

                return Normal(d1).CDF();
            }

// Vega
double BlackScholesModel::OptionGreeks::compVega(double S_0, double K, double r, double sigma,
            double T, double t_0, char type, double q) {

                double d1 = (log(S_0 / K) + (r + 0.5 * sigma * sigma) * (T - t_0)) / (sigma * sqrt(T - t_0));

                return S_0 * sqrt(T - t_0) * Normal(d1).PDF();
            }


// Gamma

double BlackScholesModel::OptionGreeks::compGamma(double S_0, double K, double r, double sigma,
            double T, double t_0, char type, double q) {

                double d1 = (log(S_0 / K) + (r + 0.5 * sigma * sigma) * (T - t_0)) / (sigma * sqrt(T - t_0));

                return Normal(d1).PDF() / (S_0 * sigma * sqrt(T - t_0));

            }

// Theta

double BlackScholesModel::OptionGreeks::compTheta(double S_0, double K, double r, double sigma,
            double T, double t_0, char type, double q) {

                double d1 = (log(S_0 / K) + (r + 0.5 * sigma * sigma) * (T - t_0)) / (sigma * sqrt(T - t_0));
                double d2 = d1 - sigma * sqrt(T - t_0);
                return (S_0 * sigma * Normal(d1).PDF()) / (2 * sqrt(T - t_0)) + r * K * exp(-r * (T - t_0)) * Normal(d2).CDF() ;
            }

// Rho 

double BlackScholesModel::OptionGreeks::compRho(double S_0, double K, double r, double sigma,
            double T, double t_0, char type, double q) {

                double d1 = (log(S_0 / K) + (r + 0.5 * sigma * sigma) * (T - t_0)) / (sigma * sqrt(T - t_0));
                double d2 = d1 - sigma * sqrt(T - t_0);
                return K * (T - t_0) * exp(-r * (T - t_0)) * Normal(d2).CDF();

            }



double BlackScholesModel::price() const {

        double d1 = (log(S_0 / K) + (r + 0.5 * sigma * sigma) * (T - t_0)) / (sigma * sqrt(T - t_0));
        double d2 = d1 - sigma * sqrt(T - t_0);

        double C =  S_0 * Normal(d1).CDF() - K * exp(-r * (T - t_0)) * Normal(d2).CDF();

        if (this -> type == 'P') C = C - S_0 + K *  exp(-r * (T - t_0));

        return C;
}

double BlackScholesModel::impVolatility(double sigma_0, double P_market, int maxIter, std::string method, double a , double b, double tolerance ) const {

        double sig = sigma_0;
        

        int iter = 0;
        std::string message;
        cout << "Using method : " + method + " ... \n *****************************************" << endl;;

       if (method == "N-R") {
           
           double Price = this ->price();
           double stoppingCriterion = abs(Price - P_market);
           while((stoppingCriterion > tolerance) && (iter < maxIter)) {
               iter += 1;

               double vega = OptionGreeks().compVega(this->S_0, this->K, this->r, this->q, this->sigma, this->T, this->t_0, this->type);
               if (vega == 0.0){

               message = "Vega equals to 0 at iter : " + std::to_string(iter) + ".\nI suggest to use another methode(Dichotomy).\nTerminating Program with 0.0 ..";
               sig =  0.0;
               break;

                }

               sig = sig - (Price - P_market) / vega;
               if (sig < 0.0 || isnan(sig)) {
                   message = "Solution out of bounds at iter : " + std::to_string(iter) + ".\n I suggest to use another methode(Dichotomy). \n Terminating Program with 0.0 ..";
                   sig =  0.0;
                   break;

               }

               BlackScholesModel BS(this->S_0, this->K, this->r, this->q, sig, this->T, this->t_0, this->type);
               Price = BS.price();
               
               stoppingCriterion = abs(Price - P_market);

               message = "Algorithm converged in : " + std::to_string(iter) + " iterations.";
            }

                    
        }

        if (method == "Dichotomy") {

                BlackScholesModel BS_min(this->S_0, this->K, this->r, this->q, a, this->T, this->t_0, this->type);
                BlackScholesModel BS_max(this->S_0, this->K, this->r, this->q, b, this->T, this->t_0, this->type);
                
                double P_min = BS_min.price(); 
                double P_max = BS_max.price();
                

                if (P_min > P_market || P_max < P_market){
                    message = "Illegal bounds. Try to change bounds.\n Terminating program with 0.0 ..." ;
                    return 0.0;
                    
                }

                double sig_min = a;
                double sig_max = b;

                sig = 0.5 * (sig_min + sig_max);
                BlackScholesModel BS(this->S_0, this->K, this->r, this->q, sig, this->T, this->t_0, this->type);
                double Price = BS.price();
                
            
                double stoppingCriterion = abs(Price - P_market);

                while((stoppingCriterion > tolerance) && (iter < maxIter)) {
                        iter += 1;
                        if (Price - P_market > 0.0) {
                            sig_max= sig;
                            sig = 0.5 * (sig_min + sig_max);
                        } 

                        else {
                            sig_min = sig;
                            sig = 0.5 * (sig_min + sig_max);
                        }

                        BlackScholesModel BS(this->S_0, this->K, this->r, this->q, sig, this->T, this->t_0, this->type);
                        Price = BS.price();
                        stoppingCriterion = abs(Price - P_market);

                        message = "Algorithm converged in : " + std::to_string(iter) + " iterations.";
                }

        }

        cout << message << endl;
        return sig;
}