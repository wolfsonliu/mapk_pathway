////////////////////////////////////////////////////////////////////////////////
// This file is the head file for MAPK pathway ordinary differential equation.
// Author: Wolfson
// Algorithm: using GNU Scientific Library.
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

////////////////////
// Classes
////////////////////


class Parameter
//// Class: Parameter
//// Store parameters in chemical equations.
{
public:
    size_t function_num;
    // S + E <b==f> [SE] <k2-=k1> S* + E
    double *f;    // forward reaction coefficient for reaction 1.
    double *b;    // backward reaction coefficient for reaction 1.
    double *kf;    // forward reaction coefficient for reaction 2.
    double *kb;    // backward reaction coefficient for reaction 2, usually very small.
    Parameter():
        function_num(0),
        f(0),
        b(0),
        kf(0),
        kb(0){}//
    Parameter(size_t dimension);     // constructor.
    Parameter(size_t dimension,
              double input_f[],
              double input_b[],
              double input_kf[],
              double input_kb[]);    // constructor.
    Parameter(size_t dimension,
              const double &input_f,
              const double &input_b,
              const double &input_kf,
              const double &input_kb);    // constructor.
    
    // Member Function declaration.
    //// initParameter:
    //// Initiate parameter not in declaration.
    void initParameter(size_t dimension);
    void initParameter(size_t dimension,
                       const double &input_f,
                       const double &input_b,
                       const double &input_kf,
                       const double &input_kb);
    void initParameter(size_t dimension,
                       double input_f[],
                       double input_b[],
                       double input_kf[],
                       double input_kb[]);
    //// setParameter:
    //// Set the values.
    void setParameter(double input_f[],
                      double input_b[],
                      double input_kf[],
                      double input_kb[]);
    void setParameter(const double &input_f,
                      const double &input_b,
                      const double &input_kf,
                      const double &input_kb);
    void setParameter(double input[]);
    void setParameter(std::vector<double> &input);
    //// printParameter:
    //// Print the values.
    void printParameter();
    //// outputParameter:
    //// Output parameter to file.
    void outputParameter(std::ofstream & outfile);
};    // remember the ;


class Concentration
//// Class: Concentration
//// store Concentration values with time.
{
private:
    bool is_new;
public:
    size_t reactant_num;
    double time;    // store the time.
    double *reactant;    // store reactant concentrate at time.
    Concentration();
    Concentration(size_t dimension,
                  double now_t,
                  double y[]);    // constructor.


    // Member Function Declaration.
    //// setTime: 
    //// Set the time of Concentration.
    void setTime(double now_t);
    //// setConcentration:
    //// Set the concentrate values.
    void setConcentration(size_t dimension, double y[]);
    //// printConcentration:
    //// Print the values.
    void printConcentration();    
    //// delConcentration:
    //// Set free the memory.
    void delConcentration();
};


class ReactantConcentration
//// Class: ReactantConcentration
//// store the list of reactant concentration changes with time.
{
    bool reach_stable;
public:
    int reaction_type;
    int function_num;
    std::vector<Concentration> list;
    std::vector<double> dissipation;    // energy dissipation.
    std::vector<std::string> name;
    size_t reactant_num;
    int final_out;    // the order of final output.
    ReactantConcentration(int type);

    //// setReactantName:
    //// Set the name of reactant in certain reactions.
    void setReactantName(const std::vector<std::string> &namelist);
    //// odeRun:
    //// Using gsl to solve ode.
    void odeRun(double start_t,
                double end_t,
                Parameter *params,
                double reactant[],
                double pacelen);
    //// outputList:
    //// Output result to file.
    void outputList(std::ofstream & outfile);
    //// outputStable:
    //// Output the stable situation to file.
    void outputStable(std::ofstream & outfile);
    
private:
    //// reactantName:
    //// Return the vector of strings in a kind of reaction.
    void reactantName();
    //// funcNumber:
    //// Set the number of chemical reaction according to reaction type.
    void funcNumber();
    //// finalOrder:
    //// Set the number of final output order according to reaction type.
    void finalOrder();

    //// odeFunction1s1p:
    //// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
    static int odeFunction1s1p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian1s1p:
    //// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
    static int odeJacobian1s1p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction1s2p:
    //// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
    static int odeFunction1s2p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian1s2p:
    //// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
    static int odeJacobian1s2p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction2s1p:
    //// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 1 phosphorylation.
    static int odeFunction2s1p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian2s1p:
    //// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 1 phosphorylation.
    static int odeJacobian2s1p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction2s2p:
    //// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 2 phosphorylations.
    static int odeFunction2s2p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian2s1p:
    //// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 2 phosphorylations.
    static int odeJacobian2s2p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction3s1p:
    //// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 1 phosphorylation.
    static int odeFunction3s1p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian3s1p:
    //// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 1 phosphorylation.
    static int odeJacobian3s1p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction3s2p:
    //// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
    static int odeFunction3s2p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian3s2p:
    //// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
    static int odeJacobian3s2p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// judgeStable:
    //// Whether the ordinary differential equations reach stable.
    int isStable();
    //// calculateDissipation:
    //// Energy dissipation of each time.
    double calculateDissipation(const double yarray[],
                                void * param);
    //// propensity1s1p:
    //// Needed by calculateDissipation, 1 stage, 1 phosphorylation.
    static void propensityFunction1s1p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity1s2p:
    //// Needed by calculateDissipation, 1 stage, 2 phosphorylations.
    static void propensityFunction1s2p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity2s1p:
    //// Needed by calculateDissipation, 2 stages, 1 phosphorylation.
    static void propensityFunction2s1p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity2s2p:
    //// Needed by calculateDissipation, 2 stages, 2 phosphorylations.
    static void propensityFunction2s2p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity3s1p:
    //// Needed by calculateDissipation, 3 stages, 1 phosphorylation.
    static void propensityFunction3s1p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity3s2p:
    //// Needed by calculateDissipation, 3 stages, 2 phosphorylations.
    static void propensityFunction3s2p(const double y[],
                                       double prob[],
                                       void *para);
};


////////////////////
// Functions
////////////////////


