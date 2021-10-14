#ifndef DIFFUSIONS_H
#define DIFFUSIONS_H
#include <cmath>

/***************************************************
This class describes a stochastic process governed by :
dX_t = b(t, X_t)dt + sigma(t, X_t)dW_t ; 
st (W_t)_t is a standard brownian motion
***************************************************/

class Diffusions {

    private:
        double x_0;
    
    public:
        Diffusions(double x0) : x_0(x0) {}
        virtual ~Diffusions() {}

        double getx_0() const {return x_0;}
        virtual double drift(double t, double x) const = 0;
        virtual double diffusion(double t, double x) const = 0;
        
        // Euler approx of expect
        //E(X_{t_0 + delta} | X_{t_0} = x_0) = x_0 + mu(t_0, x_0) delta

        virtual double expectation(double t_0, double x_0, double dt) {
            return x_0 + drift(t_0, x_0) * dt;
        }

        //Euler approx of Var:
        //V(X_{t_0 + delta } | X_{t_0} = x_0) = sigma(t_0, x_0)^2 * dt

        virtual double variance(double t_0, double x_0, double dt) {
            double sig = diffusion(t_0, x_0);
            return sig * sig * dt;
        }
};

/***************************************************
This class describes a stochastic process governed by :
dX_t = (r - 0.5 * sigma^2)dt + sigma * dW_t ; 
st (W_t)_t is a standard brownian motion
***************************************************/

class BlackScholesProcess : public Diffusions {
    
    private:
        double r, sigma;
    
    public:
        BlackScholesProcess(double r, double volatility, double s_0 = 0.01):
        Diffusions(s_0), r(r), sigma(volatility) {}

        double drift(double t, double x) const {
            return r - 0.5 * sigma * sigma;
        }

        double diffusion(double t, double x) const {
            return sigma;
        }
};

/***************************************************
This class describes a stochastic process governed by :
dX_t = -a * X_t * dt + sigma * dW_t ; 
st (W_t)_t is a standard brownian motion
***************************************************/

class OrsteinUhlenbeckProcess : public Diffusions {
    private:
        double a, sigma;
    
    public:
        OrsteinUhlenbeckProcess(double a, double sigma, double x_0 = 0.01)
        : Diffusions(x_0), a(a), sigma(sigma) {}


        double drift(double t, double x) const {
            return - a * x ;
        }

        double diffusion(double t, double x) const {

            return sigma;
        }

        double expectation(double t_0, double x_0, double dt) const {
            return x_0 * exp(- a * dt);
        }

        double variance(double t_0, double x_0, double dt) const {
            return 0.5 * sigma * sigma / a * (1.0 - exp(-2.0 * a * dt));
        }

};

/***************************************************
This class describes a stochastic process governed by :
dX_t = a *(b - X_t) * dt + sigma * sqrt(X_t) * dW_t ; 
st (W_t)_t is a standard brownian motion
***************************************************/

class SquareRootProcess : public Diffusions {
    private:
        double b, a, sigma;
    
    public:
        SquareRootProcess(double b, double a, double sigma, double x_0 = 0.01):
        Diffusions(x_0), b(b), a(a), sigma(sigma) {}

        double drift(double t, double x) const {

            return a * (b - x);
        }

        double diffusion(double t, double x) const {
            return sigma * sqrt(x);
        }
};


#endif