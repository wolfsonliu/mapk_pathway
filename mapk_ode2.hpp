////////////////////////////////////////////////////////////////////////////////
// This file is the head file for MAPK pathway ordinary differential equation.
// Author: Wolfson
// Algorithm: using GNU Scientific Library.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
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


class Parameter;
class Concentration;
class ReactantConcentration;
class Stable;
class StableList;
class DataInterpolation;

////////////////////
// Functions:
////////////////////

std::ostream& operator<<(std::ostream&          out,
                         Parameter&             para);
std::ostream& operator<<(std::ostream&          out,
                         Concentration&         conc);
std::ostream& operator<<(std::ostream&          out,
                         ReactantConcentration& rc);
std::ostream& operator<<(std::ostream&          out,
                         StableList&            sl);
std::ostream& operator<<(std::ostream&          out,
                         DataInterpolation&     di);
std::ostream& operator<<(std::ostream&          out,
                         Stable&                stb);
bool isStableUseEndSome(std::unique_ptr<double[]>&  x,
                        std::unique_ptr<double[]>&  y,
                        size_t                      size,
                        double                      threshold = 0.0001,
                        int                         endnum    = 5);
bool isMonotoneUseTwoVct(std::vector<double>& x,
                         std::vector<double>& y);
bool isOscillationUseOneVct(//std::vector<double>& x,
                            std::vector<double>& y);
void reactantName(std::vector<std::string>& name, int& reaction_type);
//// isMonotoneUseTwoVct:
//////// Whether the values are monotonous.

//// Function: isStableUseEndSome
//// Whether reach stable,
//// return true when stable, return false when not stable.
////////////////////


////////////////////
// Classes
////////////////////


class Parameter: public std::enable_shared_from_this<Parameter>
//// Class: Parameter
//// Store parameters in chemical equations.
{
    // shared_from_this();
    size_t function_num;
    double calculateATP();
public:
    // S + E <b==f> [SE] <k2-=k1> S* + E
    std::unique_ptr<double[]> f;
    // forward reaction coefficient for reaction 1.
    std::unique_ptr<double[]> b;
    // backward reaction coefficient for reaction 1.
    std::unique_ptr<double[]> kf;
    // forward reaction coefficient for reaction 2.
    std::unique_ptr<double[]> kb;
    double                    atp;
    // backward reaction coefficient for reaction 2, usually very small.
    // Constructor
    Parameter():
        function_num(0),
        atp(0.0) {}//
    Parameter(const int dimension);
    Parameter(const int dimension,
              const std::unique_ptr<double[]>& input_f,
              const std::unique_ptr<double[]>& input_b,
              const std::unique_ptr<double[]>& input_kf,
              const std::unique_ptr<double[]>& input_kb);
    Parameter(const int dimension,
              const double& input_f,
              const double& input_b,
              const double& input_kf,
              const double& input_kb);
    Parameter(const Parameter& para);
    
    // Destructor
    virtual ~Parameter() {}
    // Friend
    friend std::ostream& operator<<(std::ostream&    out,
                                    Parameter&       para);
    // Operator
    //// operator=
    //////// Overload operator=.
    Parameter& operator=(const Parameter& para);
    
    // Member-Function
    //// initParameter:
    //////// Initiate parameter not in declaration.
    void initParameter(const int& dimension);
    void initParameter(const int& dimension,
                       const double& input_f,
                       const double& input_b,
                       const double& input_kf,
                       const double& input_kb);
    void initParameter(const int& dimension,
                       const std::unique_ptr<double[]>& input_f,
                       const std::unique_ptr<double[]>& input_b,
                       const std::unique_ptr<double[]>& input_kf,
                       const std::unique_ptr<double[]>& input_kb);
    //// setParameter:
    //////// Set the values.
    void setParameter(const std::unique_ptr<double[]>& input_f,
                      const std::unique_ptr<double[]>& input_b,
                      const std::unique_ptr<double[]>& input_kf,
                      const std::unique_ptr<double[]>& input_kb);
    void setParameter(const double& input_f,
                      const double& input_b,
                      const double& input_kf,
                      const double& input_kb);
    void setParameter(const std::unique_ptr<double[]>& input);
    void setParameter(const std::vector<double>& input);
    //// printParameter:
    //////// Print the values.
    void printParameter();

    
};    // remember the ;


class Concentration: public std::enable_shared_from_this<Concentration>
//// Class: Concentration
//// store Concentration values with time.
{
    bool is_new;
    size_t reactant_num;
    std::unique_ptr<double[]> reactant;
    // store reactant concentrate at time.
public:
    // Constructor
    Concentration():
        is_new(true),
        reactant_num(0) {}//
    Concentration(const size_t& dimension);
    Concentration(const size_t& dimension,
                  const std::unique_ptr<double[]>& y);
    Concentration(const Concentration& conc);
    // Destructor
    virtual ~Concentration() {}
    // Friend
    friend std::ostream& operator<<(std::ostream&  out,
                                    Concentration& conc);

    // Operator
    //// operator=
    //////// Overload operator=.
    Concentration& operator=(const Concentration& conc);
    // Member-Function
    //// setConcentration:
    //////// Set the concentrate values.
    void setConcentration(const std::unique_ptr<double[]>& y);
    void setConcentration(double y[]);
    void setConcentration(const size_t& dimension,
                          const std::unique_ptr<double[]>& y);
    //// printConcentration:
    //////// Print the values.
    void printConcentration();    
    //// reset:
    //////// Set free the memory.
    void reset();
    //// size:
    //////// Return the size of reactant.
    size_t size();
    //// get():
    //////// Used for subscript operator overloading.
    double get(const int ord);
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
    //std::vector<std::string>   name;
    int                        final_out;
    // the order of final output.

    // Constructor
    ReactantConcentration():
        reach_stable(true),
        reaction_type(0),
        function_num(0),
        final_out(0) {}
    ReactantConcentration(const int& type);

    // Destructor
    virtual ~ReactantConcentration() {}
    // Friend
    friend std::ostream& operator<<(std::ostream&          out,
                                    ReactantConcentration& rc);
    friend bool isStableUseEndSome(std::unique_ptr<double[]>& thevector,
                                   size_t size,
                                   double threshold,
                                   int    endnum);
    friend bool isMonotoneUseTwoVct(std::vector<double>& x,
                                    std::vector<double>& y);
    friend bool isOscillationUseOneVct(//std::vector<double>& x,
                                       std::vector<double>& y);
    // Member-Function
    //// setReactantName:
    //////// Set the name of reactant in certain reactions.
    //void setReactantName(const std::vector<std::string> &namelist);
    //// outputStable:
    /////// Output the stable situation to file.
    void outputFinal(std::ofstream& outfile);
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
                Parameter* params,
                double     reactant[],
                double     pacelen);
    //// reachStable:
    //////// Check whether reach stable.
    bool reachStable();
    //// isMonotone:
    //////// Check whether result monotone.
    bool isMonotone();
    //// isOscillation:
    //////// Check whether exist oscillaiton.
    bool isOscillation();
    
private:
    // Member-Function
    //// reactantName:
    //////// Return the vector of strings in a kind of reaction.
    //void reactantName();
    //// funcNumber:
    //////// Set the number of chemical reaction according to reaction type.
    void funcNumber();
    //// finalOrder:
    //////// Set the number of final output order according to reaction type.p
    void finalOrder();
    //// judgeStable:
    //////// Whether the ordinary differential equations reach stable.
    bool isStable();
    //// odeFunction1s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
    static int odeFunction1s1p(double t, const double y[],
                               double f[], void *para);
    //// odeJacobian1s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
    static int odeJacobian1s1p(double t, const double y[],
                               double dfdy[], double dydt[],
                               void* para);
    //// odeFunction1s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
    static int odeFunction1s2p(double t, const double y[],
                               double f[], void* para);
    //// odeJacobian1s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
    static int odeJacobian1s2p(double t, const double y[],
                               double dfdy[], double dydt[],
                               void* para);
    //// odeFunction2s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 1 phosphorylation.
    static int odeFunction2s1p(double t, const double y[],
                               double f[], void* para);
    //// odeJacobian2s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 1 phosphorylation.
    static int odeJacobian2s1p(double t, const double y[],
                               double dfdy[], double dydt[],
                               void* para);
    //// odeFunction2s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 2 phosphorylations.
    static int odeFunction2s2p(double t, const double y[],
                               double f[], void* para);
    //// odeJacobian2s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 2 phosphorylations.
    static int odeJacobian2s2p(double t, const double y[],
                               double dfdy[], double dydt[],
                               void* para);
    //// odeFunction3s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 1 phosphorylation.
    static int odeFunction3s1p(double t, const double y[],
                               double f[], void* para);
    //// odeJacobian3s1p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 1 phosphorylation.
    static int odeJacobian3s1p(double t, const double y[],
                               double dfdy[], double dydt[],
                               void *para);
    //// odeFunction3s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
    static int odeFunction3s2p(double t, const double y[],
                               double f[], void* para);
    //// odeJacobian3s2p:
    //////// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
    static int odeJacobian3s2p(double t, const double y[],
                               double dfdy[], double dydt[],
                               void* para);

    //// calculateDissipation:
    //////// Energy dissipation of each time.
    double calculateDissipation(double yarray[],
                                void*  param);
    //// propensity1s1p:
    //////// Needed by calculateDissipation, 1 stage, 1 phosphorylation.
    static void propensityFunction1s1p(const double y[],
                                       double       prob[],
                                       void*        para);
    //// propensity1s2p:
    //////// Needed by calculateDissipation, 1 stage, 2 phosphorylations.
    static void propensityFunction1s2p(const double y[],
                                       double       prob[],
                                       void*        para);
    //// propensity2s1p:
    //////// Needed by calculateDissipation, 2 stages, 1 phosphorylation.
    static void propensityFunction2s1p(const double y[],
                                       double       prob[],
                                       void*        para);
    //// propensity2s2p:
    //////// Needed by calculateDissipation, 2 stages, 2 phosphorylations.
    static void propensityFunction2s2p(const double y[],
                                       double       prob[],
                                       void*        para);
    //// propensity3s1p:
    //////// Needed by calculateDissipation, 3 stages, 1 phosphorylation.
    static void propensityFunction3s1p(const double y[],
                                       double       prob[],
                                       void*        para);
    //// propensity3s2p:
    //////// Needed by calculateDissipation, 3 stages, 2 phosphorylations.
    static void propensityFunction3s2p(const double y[],
                                       double       prob[],
                                       void*        para);
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
    Stable(const double& i,
           const double& o,
           const double& t,
           const double& d):
        input(i),
        output(o),
        time(t),
        dissipation(d) {}
    Stable(const Stable& stb):
        input(stb.input),
        output(stb.output),
        time(stb.time),
        dissipation(stb.dissipation) {}
    
    // Destructor
    virtual ~Stable() {}

    // Friend
    friend std::ostream& operator<<(std::ostream& out,
                                    Stable&       stb);
    // Member-Function
    //// setValue
    //////// Set the value of instant.
    void setValue(const double& i,
                  const double& o,
                  const double& t,
                  const double& d);
    //// getValue
    //////// Return value of according.
    double getValue(const int& according);
    //// lt
    //////// lt is less than <.
    static bool lt(const Stable& a,
                   const Stable& b,
                   const int&    according);
    //// gt
    //////// gt is greater than >.
    static bool gt(const Stable& a,
                   const Stable& b,
                   const int&    according);
    //// eq
    //////// eq is equal =.
    static bool eq(const Stable& a,
                   const Stable& b,
                   const int&    according);
};

    
class StableList
//// Class: StableList
//// Store stable situation of different input.
{
    int               reaction_type;
    int               function_num;
    size_t            reactant_num;
    static bool ltcompa(Stable& a, Stable& b);
    static bool ltcompb(Stable& a, Stable& b);
    static bool ltcompc(Stable& a, Stable& b);
    static bool ltcompd(Stable& a, Stable& b);
public:
    std::vector<Stable> situation;


    // Constructor
    StableList():
        reaction_type(0),
        function_num(0),
        reactant_num(0) {}
    StableList(const int& type);
    StableList(const int& type,
               const int& functionnum,
               const int& reactantnum):
        reaction_type(type),
        function_num(functionnum),
        reactant_num(reactantnum) {}
    // Destructor
    virtual ~StableList() {}

    // Friend
    friend std::ostream& operator<<(std::ostream& out,
                                    StableList&   sl);
    friend bool isStableUseEndSome(std::unique_ptr<double[]>&  thevector,
                                   size_t size,
                                   double threshold,
                                   int    endnum);
    // Type Define
    typedef std::vector<Stable>::iterator iterator;
    typedef Stable value_type;
    // Operator
    //// operator[]
    Stable& operator[](const size_t index);
    
    // Member-Function
    //// addBack:
    //////// Add values to the end.
    void addBack(const double& iput,
                 const double& oput,
                 const double& t,
                 const double& e);
    //// addInter:
    //////// Insert values according to iterator.
    void addInter(iterator      it,
                  const double& iput,
                  const double& oput,
                  const double& t,
                  const double& e);
    void addInter(iterator   it,
                  std::vector<Stable>::value_type value);
    //// unSufficient:
    //////// Return the first pointer if the gap between output is too large.
    std::vector<Stable>::iterator unSufficient();
    //// isMonotone
    //////// Return true if exists monotone.
    bool isMonotone();
    //// max:
    //////// Return the max value.
    std::vector<Stable>::iterator max(int according = 2);
    //// min:
    //////// Return the min value.
    std::vector<Stable>::iterator min(int according = 2);
    //// isStable:
    //////// Judge whether reach stable.
    bool isStable(int xaccording = 1, int yaccording = 2);
    //// sortList:
    //////// Sort situation.
    void sortList(int according = 1);

    //// getPercent
    //////// Return the iterator pointed to situation with value near specific percentage.
    std::vector<Stable>::iterator getPercent(const double& percent,
                                             int according = 2);
    //// getArray
    //////// Return double[] of values using according to specific.
    std::unique_ptr<double[]> getArray(const int& according);
    //// begin:
    //////// Return the begin of situation.
    std::vector<Stable>::iterator begin() {return situation.begin();}
    //// end:
    //////// Return the end of situation.
    std::vector<Stable>::iterator end() {return situation.end();}
    //// size:
    //////// Return the size of situation.
    size_t size();
};


class DataInterpolation
//// Class: StableList
//// Store stable situation of different input.
{
public:
    std::vector<double> x;
    std::vector<double> y;

    // Destructor
    virtual ~DataInterpolation() {}
    // Friend
    friend std::ostream& operator<<(std::ostream&    out,
                                    DataInterpolation& di);
    // Member-Function
    //// interpolate:
    //////// Interpolate data.
    void interpolate(std::unique_ptr<double[]> xarray,
                     std::unique_ptr<double[]> yarray,
                     long                      length,
                     long                      number);
    void interpolateL(std::unique_ptr<double[]> xarray,
                      std::unique_ptr<double[]> yarray,
                      long                      length,
                      long                      number);
    //// size:
    //////// Return the size of data.
    size_t size();
    //// nearPercentNum
    //////// Return the order of value near specific percentage.
    int nearPercentNum(const double& percent,
                       const int     according = 2);
    //// nearValueNum
    //////// Return the order of value near specific percentage.
    int nearValueNum(const double& value,
                     const int     according = 1);    
};


////////////////////
// Member-Functions
////////////////////







