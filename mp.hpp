#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

const int ARRAY_NUM = 30;    // number of members in array.
const int ARRAY_NUM_BIG = 60;

// Classes
//// Class Parameter: store parameters in chemical equations.
class Parameter {
public:
    // S + E <b==f> [SE] <k2-=k1> S* + E
    double *f = new double [ARRAY_NUM];    // forward reaction coefficient for reaction 1.
    double *b = new double [ARRAY_NUM];    // backward reaction coefficient for reaction 1.
    double *kf = new double [ARRAY_NUM];    // forward reaction coefficient for reaction 2.
    double *kb = new double [ARRAY_NUM];    // backward reaction coefficient for reaction 2, usually very small.
    Parameter();
    Parameter(size_t dimension, double input_f[], double input_b[], double input_kf[], double input_kb[]);    // constructor.
    
    // Member Function declaration.
    void setParameter(size_t dimension, double input_f[], double input_b[], double input_kf[], double input_kb[]);    // set the values.
    void printParameter(size_t dimension);    // print the values.
};    // remember the ;

//// Class Concentrate: store Concentrate values with time.
class Concentrate {
public:
    double time;    // store the time.
    double *reactant = double [ARRAY_NUM_BIG];    // store reactant concentrate at time.
    Concentrate();
Concentrate(size_t dimension, double now_t, double y[]);    // constructor.

    // Member Function declaration.
    void setTime(double now_t);    // set the time of Concentrate.
    void setConcentrate(size_t dimension, double y[]);    // set the concentrate values.
    void printConcentrate(size_t dimension);    // print the values.
};


// Functions
//// Function odeFunction1s1p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
int odeFunction1s1p(double t, const double y[], double f[], void *para);
//// Function odeJacobian1s1p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
int odeJacobian1s1p(double t, const double y[], double *dfdy, double dydt[], void *para);
//// Function odeFunction1s2p:  needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
int odeFunction1s2p(double t, const double y[], double f[], void *para);
//// Function odeJacobian1s2p:  needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
int odeJacobian1s2p(double t, const double y[], double *dfdy, double dydt[], void *para);
//// Function odeRun: using gsl to solve ode.
std::vector<Concentrate> odeRun(int func, double start_t, double end_t, Parameter *params, double y[], size_t dimension, double pacelen);

//// Function initiateParameter: set the parameter array to a certain velue.
int initiateParameter(double par[], int len, double value);


