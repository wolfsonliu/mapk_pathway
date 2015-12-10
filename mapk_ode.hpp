////////////////////////////////////////////////////////////////////////////////
// This file is the head file for MAPK pathway ordinary differential equation.
// Author: Wolfson
// Algorithm: using GNU Scientific Library.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>


////////////////////
// Classes
////////////////////


class Parameter
//// Class: Parameter
//// Store parameters in chemical equations.
{
    size_t function_num;
public:
    // S + E <b==f> [SE] <k2-=k1> S* + E
    double *f;    // forward reaction coefficient for reaction 1.
    double *b;    // backward reaction coefficient for reaction 1.
    double *kf;    // forward reaction coefficient for reaction 2.
    double *kb;    // backward reaction coefficient for reaction 2, usually very small.
    // Constructor
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
    
    // Member-Function
    //// initParameter:
    //////// Initiate parameter not in declaration.
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
    //////// Set the values.
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
    //////// Print the values.
    void printParameter();
    //// outputParameter:
    //////// Output parameter to file.
    void outputParameter(std::ofstream & outfile);
};    // remember the ;


class Concentration
//// Class: Concentration
//// store Concentration values with time.
{
    bool is_new;
    size_t reactant_num;
    double *reactant;    // store reactant concentrate at time.
public:
    // Constructor
    Concentration():
        is_new(true),
        reactant_num(0),
        reactant(0){}
    Concentration(size_t dimension);
    Concentration(size_t dimension,
                  double y[]);


    // Member-Function
    //// setConcentration:
    //////// Set the concentrate values.
    void setConcentration(double y[]);
    void setConcentration(size_t dimension,
                          double y[]);
    //// printConcentration:
    //////// Print the values.
    void printConcentration();    
    //// reset:
    //////// Set free the memory.
    void reset();
    //// size:
    //////// Return the size of reactant.
    inline size_t size();
    //// operator[]:
    //////// Used for subscript operator overloading.
    inline double &operator[](int i);
};


class ReactantConcentration
//// Class: ReactantConcentration
//// Store the list of reactant concentration changes with time.
{
    bool                       reach_stable;
public:
    int                        reaction_type;
    int                        function_num;
    size_t                     reactant_num;
    std::vector<double>        time;
    std::vector<Concentration> list;
    std::vector<double>        dissipation;
    // energy dissipation.
    std::vector<std::string>   name;
    int                        final_out;
    // the order of final output.

    // Constructor
    ReactantConcentration(int type);

    // Member-Function
    //// setReactantName:
    //////// Set the name of reactant in certain reactions.
    void setReactantName(const std::vector<std::string> &namelist);
    //// outputList:
    //////// Output result to file.
    void outputList(std::ofstream & outfile);
    //// outputStable:
    /////// Output the stable situation to file.
    void outputFinal(std::ofstream & outfile);
    //// getFinalOutput:
    //////// Return the output reactant concentration values at stable stage.
    double getFinalOutput();
    //// getFinalTime:
    //////// Return the time at stable stage.
    double getFinalTime();
    //// getFinalDissipation:
    //////// Return the energy dissipation at stable stage.
    double getFinalDissipation();
    //// odeRun:
    //////// Using gsl to solve ode.
    void odeRun(double     start_t,
                double     end_t,
                Parameter *params,
                double     reactant[],
                double     pacelen);
private:
    // Member-Function
    //// reactantName:
    //////// Return the vector of strings in a kind of reaction.
    void reactantName();
    //// funcNumber:
    //////// Set the number of chemical reaction according to reaction type.
    void funcNumber();
    //// finalOrder:
    //////// Set the number of final output order according to reaction type.
    void finalOrder();
    //// judgeStable:
    //////// Whether the ordinary differential equations reach stable.
    bool isStable();
    //// isMonotone:
    //////// Whether the values are monotonous.
    bool isMonotone(int order);

    //// odeFunction1s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
    static int odeFunction1s1p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian1s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
    static int odeJacobian1s1p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction1s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
    static int odeFunction1s2p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian1s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
    static int odeJacobian1s2p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction2s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 1 phosphorylation.
    static int odeFunction2s1p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian2s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 1 phosphorylation.
    static int odeJacobian2s1p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction2s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 2 phosphorylations.
    static int odeFunction2s2p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian2s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 2 phosphorylations.
    static int odeJacobian2s2p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction3s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 1 phosphorylation.
    static int odeFunction3s1p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian3s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 1 phosphorylation.
    static int odeJacobian3s1p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);
    //// odeFunction3s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
    static int odeFunction3s2p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian3s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
    static int odeJacobian3s2p(double t, const double y[],
                               double *dfdy, double dydt[], void *para);

    //// calculateDissipation:
    //////// Energy dissipation of each time.
    double calculateDissipation(const double yarray[],
                                void * param);
    //// propensity1s1p:
    //////// Needed by calculateDissipation, 1 stage, 1 phosphorylation.
    static void propensityFunction1s1p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity1s2p:
    //////// Needed by calculateDissipation, 1 stage, 2 phosphorylations.
    static void propensityFunction1s2p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity2s1p:
    //////// Needed by calculateDissipation, 2 stages, 1 phosphorylation.
    static void propensityFunction2s1p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity2s2p:
    //////// Needed by calculateDissipation, 2 stages, 2 phosphorylations.
    static void propensityFunction2s2p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity3s1p:
    //////// Needed by calculateDissipation, 3 stages, 1 phosphorylation.
    static void propensityFunction3s1p(const double y[],
                                       double prob[],
                                       void *para);
    //// propensity3s2p:
    //////// Needed by calculateDissipation, 3 stages, 2 phosphorylations.
    static void propensityFunction3s2p(const double y[],
                                       double prob[],
                                       void *para);
};


class Stable
//// Class: Stable
//// Use in StableVector.
{

public:
    double input;
    double output;
    double time;
    double dissipation;

    // Constructor
    Stable():
        input(0.0),
        output(0.0),
        time(0.0),
        dissipation(0.0) {}
    Stable(double i, double o, double t, double d);

    // Member-Function
    //// setValue
    //////// Set the value of instant.
    void setValue(double i, double o, double t, double d);
    //// getValue
    //////// Return value of according.
    double getValue(int according);
    //// lt
    //////// lt is less than <.
    static bool lt(Stable& a, Stable& b, int according);
    //// gt
    //////// gt is greater than >.
    static bool gt(Stable& a, Stable& b, int according);
    //// eq
    //////// eq is equal =.
    static bool eq(Stable& a, Stable& b, int according);
};

    
class StableList
//// Class: StableList
//// Store stable situation of different input.
{
    int               reaction_type;
    int               function_num;
    size_t            reactant_num;
    static bool compa(Stable& a, Stable& b);
    static bool compb(Stable& a, Stable& b);
    static bool compc(Stable& a, Stable& b);
    static bool compd(Stable& a, Stable& b);
public:
    std::vector<Stable> situation;


    // Constructor
    StableList():
        reaction_type(0),
        function_num(0),
        reactant_num(0) {}
    StableList(int type);
    StableList(int         type,
               int         functionnum,
               int         reactantnum);

    typedef std::vector<Stable>::iterator iterator;
    // Member-Function
    //// addBack:
    //////// Add values to the end.
    void addBack(double iput,
                 double oput,
                 double t,
                 double e);
    //// addInter:
    //////// Insert values according to iterator.
    void addInter(iterator it,
                  double                      iput,
                  double                      oput,
                  double                      t,
                  double                      e);
    void addInter(iterator   it,
                  std::vector<Stable>::value_type value);
    //// notSufficient:
    //////// Return the first pointer if the gap between output is too large.
    std::vector<Stable>::iterator unSufficient();
    //// maxOutput:
    //////// Return the max output value.
    std::vector<Stable>::iterator max(int according = 2);
    //// sortList:
    //////// Sort situation.
    void sortList(int according = 1);

    //// outputList:
    //////// Output result to file.
    void outputList(std::ofstream & outfile);
    //// getPercent
    //////// Return the iterator pointed to situation with value near specific percentage.
    std::vector<Stable>::iterator getPercent(double percent, int according = 2);
    //// getArray
    //////// Return double[] of values using according to specific.
    double* getArray(int according);
    //// begin:
    //////// Return the begin of situation.
    inline std::vector<Stable>::iterator begin() {return situation.begin();}
    //// end:
    //////// Return the end of situation.
    inline std::vector<Stable>::iterator end() {return situation.end();}
    //// size:
    //////// Return the size of situation.
    inline size_t size();
};


class DataInterpolation
//// Class: StableList
//// Store stable situation of different input.
{
public:
    std::vector<double> x;
    std::vector<double> y;

    // Member-Function
    //// interpolate:
    //////// Interpolate data.
    void interpolate(const double* xarray,
                     const double* yarray,
                     long          length,
                     long          number);
    //// size:
    //////// Return the size of data.
    inline size_t size();
};


////////////////////
// Member-Functions
////////////////////


size_t Concentration::size()
//// Member-Function: Concentration::size
//// Return the size of reactant.
{
    return reactant_num;
}


double& Concentration::operator[](int i)
//// Member-Function: Concentration::operator[]
//// Used for subscript operator overloading.
{
    if (i >= reactant_num) {
        std::cout << "Concentration: subscript operator out of bounds."
                  << "\n";
        exit(1);
    }
    return reactant[i];
}



size_t StableList::size()
//// Member-Function: StableList::size
//// Return the size of situation.
{
    return situation.size();
}


size_t DataInterpolation::size()
//// Member-Function: DataInterpolation::size
//// Return the size of DataInterpolation.
{
    return x.size();
}
