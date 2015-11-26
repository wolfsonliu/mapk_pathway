#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


// Classes

//// Class
//// Parameter: store parameters in chemical equations.
class Parameter {
public:
    size_t function_num;
    // S + E <b==f> [SE] <k2-=k1> S* + E
    double *f;    // forward reaction coefficient for reaction 1.
    double *b;    // backward reaction coefficient for reaction 1.
    double *kf;    // forward reaction coefficient for reaction 2.
    double *kb;    // backward reaction coefficient for reaction 2, usually very small.
    Parameter();
    Parameter(size_t dimension, double input_f[], double input_b[], double input_kf[], double input_kb[]);    // constructor.
    
    // Member Function declaration.
    void setParameter(double input_f[], double input_b[], double input_kf[], double input_kb[]);    // set the values.
    void printParameter();    // print the values.
};    // remember the ;

//// Class
//// Concentration: store Concentration values with time.
class Concentration {
private:
    bool is_new;
public:
    size_t reactant_num;
    double time;    // store the time.
    double *reactant;    // store reactant concentrate at time.
    Concentration();
    Concentration(size_t dimension, double now_t, double y[]);    // constructor.

    // Member Function Declaration.
    void setTime(double now_t);    // set the time of Concentration.
    void setConcentration(size_t dimension, double y[]);    // set the concentrate values.
    void printConcentration();    // print the values.
    void delConcentration();    // set free the memory.
};

//// Class
//// ReactantConcentration: store the list of reactant concentration changes with time.
class ReactantConcentration {
public:
    int reaction_type;
    std::vector<Concentration> list;
    std::vector<std::string> name;
    size_t reactant_num;
    ReactantConcentration(int type, std::vector<std::string> react_name, size_t react_num);
    // Member-Function
    //// Member Function
    //// odeFunction1s1p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
    static int odeFunction1s1p(double t, const double y[], double f[], void *para);
    //// Member Function
    //// odeJacobian1s1p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
    static int odeJacobian1s1p(double t, const double y[], double *dfdy, double dydt[], void *para);
    //// Member Function
    //// odeFunction1s2p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
    static int odeFunction1s2p(double t, const double y[], double f[], void *para);
    //// Member Function
    //// odeJacobian1s2p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
    static int odeJacobian1s2p(double t, const double y[], double *dfdy, double dydt[], void *para);   

    // setReactantName: set the name of reactant in certain reactions.
    void setReactantName(const std::vector<std::string> &namelist);
    //// Member-Function
    //// odeRun: using gsl to solve ode.
    void odeRun(double start_t, double end_t, Parameter *params, double reactant[], double pacelen);

};


// Functions

//// Function
//// initiateParameter: set the parameter array to a certain velue.
int initiateParameter(double par[], int len, double value);

//// Function
//// reactantName: return the vector of strings in a kind of reaction.
std::vector<std::string> reactantName(int functype);

