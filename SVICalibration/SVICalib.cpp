#include <nlopt.hpp>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include "gnuplot-iostream.h"

/***********************************************
*                                              *
*          @brief   SVI CALIBRATION            *
*          @author  BOUDCHICHI Oussama         *
*                                              *
* *********************************************/         



using namespace std;

// % Begin

// % Global variables to store data and communicate with functions
std::vector <double> lMoneyness;
std::vector <double> totalVar;

// % SVI Function:
vector<double> SVI(const vector<double>& theta, const vector<double>& x) 
{
    double a, b, rho, m, sigma;
    a = theta[0]; b = theta[1]; rho = theta[2];
    m = theta[3]; sigma = theta[4];

    vector<double> w(x.size());
    for(int i=0; i < x.size(); i++) {
        double y = a + b * (rho * (x[i] - m) + std::sqrt((x[i] - m) * (x[i] - m) + sigma * sigma));
        w.push_back(y);

    }

    return w;
}

// % Least-Squares objective function
double LSQ(const vector<double>& theta,  std::vector<double> &grad, void *my_func_data) {

    double a, b, rho, m, sigma;
    a = theta[0]; b = theta[1]; rho = theta[2];
    m = theta[3]; sigma = theta[4];

    double* data = static_cast <double*>(my_func_data);
    //vector <double> W = SVI(theta, lMoneyness);
    double sum = 0;
    int n  = sizeof(&data) / sizeof(&data[0]);
    for (int i=0; i < n; i++) {

        double W = a + b * (rho * (data[i] - m) + sqrt(pow(data[i] - m, 2.) + pow(sigma, 2.))); 
        sum += (W - totalVar[i]) * (W - totalVar[i]);

    }

    return sqrt(sum);

}


// % Convexity Test g function : 
vector<double> g_v(const vector<double>& theta) {

    vector <double> W = SVI(theta, lMoneyness);
    double a, b, rho, m, sigma;
    a = theta[0]; b = theta[1]; rho = theta[2];
    m = theta[3]; sigma = theta[4];
    vector<double> ret(lMoneyness.size());

    for(int i=0; i < W.size(); i++) {
        double wp = b * rho + (b * (lMoneyness[i] - m)) / sqrt((lMoneyness[i] - m) * (lMoneyness[i] - m) + sigma * sigma);
        double wpp = (b * sigma * sigma) / sqrt(pow((lMoneyness[i] - m) + sigma * sigma, 3.0));
        double y = 0.5 * wpp + pow((1 - (lMoneyness[i] * wp) / (2 * W[i])), 2.0) - pow((0.5 * wp), 2.0) * (1 / (W[i]) + .25);
        ret.push_back(y);

    }

    return ret;

}

// % Multi Constraints function. NLopt convention : c(x) <= 0
void multi_constraint(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)
{
    double *data = static_cast <double*>(f_data);
    // int s = sizeof(data);
    // vector<double> V(data, data + s);

    for (int i = 0; i < m; i++) {
        double a, b, rho, m, sigma;
        a = x[0]; b = x[1]; rho = x[2];
        m = x[3]; sigma = x[4];
        double w = a  + b * (rho * (data[i] - m ) + sqrt(data[i] * data[i] + sigma * sigma));
        double wp = b * rho + (b * (data[i] - m)) / sqrt((data[i] - m) * (data[i] - m) + sigma * sigma);
        double wpp = (b * sigma * sigma) / sqrt(pow((data[i] - m) + sigma * sigma, 3.0));
        result[i] = 5e-3 - (0.5 * wpp + pow((1 - (data[i] * wp) / (2 * w)), 2.0) - pow((0.5 * wp), 2.0) * (1 / (w) + .25));
    }
}


// % Utils : 
vector<double> getMinMax(vector<double>& V) {
    auto [min, max] = std::minmax_element(begin(V), end(V));
    vector<double> result {*min, *max};
    return result;

}

// % 
vector<double> getSlice(vector<double>& V, vector<double>& TT, const double& T) {
    vector<double> res;
    for(int i=0; i < V.size(); i++) {
        if (TT[i] == T) res.push_back(V[i]);
    }

    return res;
}




int main() {

    // %  Read data :
    ifstream indata; indata.open("data.txt");
    if(!indata) {
        cerr << "Error: file coukd not be opened" << endl;
        exit(1);
    }

    vector<double> x, w,  TT;
    double tempx, tempw, tempt;
    while(!indata.eof()) {
        indata >> tempx >> tempw >> tempt;
        x.push_back(tempx); w.push_back(tempw); TT.push_back(tempt);
    }

    indata.close();
    // % Finish Reading.
    

     
    // % Get desired slice, T = T_0
    double Maturity = TT[0];
    lMoneyness = getSlice(x, TT, Maturity);
    totalVar   = getSlice(w, TT, Maturity);
    
    // Initialize Optimizer, algo : COBYLA
    nlopt::opt opt(nlopt::LN_COBYLA, 5);

    // % Get min_x(w), max_x(w) 
    vector<double> m_M_x, m_M_w;
    m_M_x = getMinMax(lMoneyness); m_M_w = getMinMax(totalVar);

    // % Define Bounds
    vector <double> lb {1e-5, 1e-2,  -0.9999999999, 2. * m_M_x[0], 1e-2};
    vector <double> ub {m_M_w[1], 1., 0.9999999999, 2. * m_M_x[1], 1.};

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    // % Constraints tolerance
    vector<double> tol(lMoneyness.size(), 0.);

    // % log-Moneyness Data feeded to L-SQ
    double* data = &lMoneyness[0];
    opt.set_min_objective(LSQ, data);

    // % Add inequality constraints
    opt.add_inequality_mconstraint(multi_constraint, data, tol);
    
    // % Set solution relative tolerance
    opt.set_xtol_rel(1e-6);
    
    // % Initial value
    vector<double> theta {0.5 * m_M_w[0] + 1e-5, 0.1, -0.5, m_M_x[0] + m_M_x[1], 0.1};

    // % Variable to store the minimimun attained by LSQ
    double minf;

    auto t1 = chrono::system_clock::now();
    // % Optimize 
    nlopt::result result = opt.optimize(theta, minf);

    auto t2 = chrono::system_clock::now();
    chrono::duration<double> diff = t2-t1;

    // % Compute fitted values, and g function
    vector<double> W = SVI(theta, lMoneyness);
    vector<double> G = g_v(theta);

    // % Export data
    ofstream out("svi.dat");
    for(int i=0; i< W.size(); i++)
    {
        
        out << lMoneyness[i]<< "\t" << W[i] << "\t" << totalVar[i] <<  "\t" << G[i] << endl;
    }
    out.close();

    cout << "\n****************************************************\n";
    cout << "Solution, for maturity " << Maturity << " : \n";
    for (double x : theta) cout << x << " " << setprecision(5);
    

    cout << "\n****************************************************" << endl;
    cout << "\nObj min : " << minf << endl;
        
    cout << "Obj evals : " << opt.get_numevals() << endl;
        
    cout << "Method : " << opt.get_algorithm_name() << endl;

    cout << "Esplaped time : " << diff.count() << " .s" << endl;

    // % Plot results :

    Gnuplot gp;
    vector<pair<double, double> > xsvi_pnts;
    vector<pair<double, double> > xmark_pnts;
    
    for(int i=0; i < W.size(); i++) {
        xsvi_pnts.push_back(make_pair(lMoneyness[i], W[i]));
        xmark_pnts.push_back(make_pair(lMoneyness[i], totalVar[i]));
        
    }
    gp << "set style line 1 linecolor rgb 'black' linetype 1 dt 2  linewidth 2\n";
    gp << "set style line 2 lw 2 lc rgb 'blue' ps 1 pt 6 pi 1\n";

    gp << "plot"<< gp.file1d(xsvi_pnts) << "w l ls 1 title 'SVI' ,"
    << gp.file1d(xmark_pnts) << "ls 2 title 'Market'" << endl;   



    return 0;

    // % End. 
}