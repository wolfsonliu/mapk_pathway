//-*-coding:utf-8-*-
////////////////////////////////////////
// This file is the function file for MAPK pathway ordinary differential equation.
// Author: Wolfson
// Date: Nov.29, 2015
// Edit: Jan.5, 2015
// Algorithm: using GNU Scientific Library.
////////////////////////////////////////

#include "mapk_ode2.hpp"



////////////////////
// Functions
////////////////////  
// {{{ Function

std::ostream& operator<<(std::ostream& out, Parameter& para)
// {{{ operator<<
//// Output parameter to stream.
{
    
    for (size_t i = 0; i != para.function_num; ++i) {
        out << i << ","
            << para.f[i] << ","
            << para.b[i] << ","
            << para.kf[i] << ","
            << para.kb[i] << "\n";
    }

    return out;
}
// }}}


std::ostream& operator<<(std::ostream& out, Concentration& conc)
// {{{ operator<<
//// Output parameter to stream.
{
    out << conc.reactant[0];
    for (size_t i = 1; i != conc.reactant_num; ++i) {
        out << "," << conc.reactant[i];
    }

    return out;
}
// }}}


std::ostream& operator<<(std::ostream& out, ReactantConcentration& rc)
// {{{ operator<<
//// Output result to stream.
{
    for (std::vector<Concentration>::size_type vi = 0;
         vi != rc.list.size();
         ++vi) {
        out << rc.time[vi];
        for (size_t i = 0; i != rc.reactant_num; ++i) {
            out << "," << rc.list[vi].get(i);
        }//for

        out << "," << rc.dissipation[vi];
        out << "\n";
    }//for

    return out;
}
// }}}


std::ostream& operator<<(std::ostream& out, Stable& stb)
// {{{ operator<<

//// Output Stable to stream.
{
    out << stb.input << ","
        << stb.output << ","
        << stb.time << ","
        << stb.dissipation;

    return out;
}

// }}}


std::ostream& operator<<(std::ostream& out, StableList& sl)
// {{{ operator<<
//// Output result to stream.
{
    
    for (size_t i = 0;
         i != sl.situation.size();
         ++i) {
        out << sl.situation[i].input << ","
            << sl.situation[i].output << ","
            << sl.situation[i].time << ","
            << sl.situation[i].dissipation << "\n";
    } // for.
    return out;
}
// }}}


std::ostream& operator<<(std::ostream& out, DataInterpolation& di)
// {{{ operator<<
//// Output result to stream.
{
    for (size_t i = 0; i != di.size(); ++i) {
        out << di.x[i] << ","
            << di.y[i] << "\n";
    }// for.
    return out;
}
// }}}


bool isStableUseEndSome(std::unique_ptr<double[]>& x,
                        std::unique_ptr<double[]>& y,
                        size_t                     size,
                        double                     threshold,
                        size_t                     endnum)
// {{{ Function: isStableUseEndSome

//// Whether reach stable,
//// return true when stable, return false when not stable.
{

    if (size < endnum + 1) {
        return false;
    }
    bool   judge    = true;
    double dy       = 0.0;
    double dx       = 0.0;

    for (size_t i = 0; i != endnum; ++i) {

        dx = std::fabs(x[size - i - 1] -
                       x[size - i - 2]);
        dy = std::fabs(y[size - i - 1] -
                       y[size - i - 2]);
        if (dx != 0 && y[size -i - 1] != 0) {
            judge = (dy * x[size - i - 1]) / (dx * y[size - i - 1])
                < threshold ? true : false;
        }
        if (!judge) {
            return judge;
        }
    } // for.
    return judge;
}

// }}}


bool isMonotoneUseTwoVct(std::vector<double>& x,
                         std::vector<double>& y)
// {{{ Function: isMonotoneWithTwoVct

//// Whether monotonous.
{
    double dx        = 0.0;
    double dy        = 0.0;
    //double kvalue    = 0.0;
    //std::vector<double> k;
    long   increase  = 0;
    long   decrease  = 0;
    long   totalsize = x.size() - 1;
    //int    sign      = 0;
    
    for (size_t i = 0; i != x.size() - 2; ++i) {
        // calculate the k.
        dx = x[i + 2] - x[i];
        dy = y[i + 2] - y[i];
        //k.push_back(dy / dx);
        if ((dy / dx) >= 0) {
            ++increase;
        } else if ((dy / dx) < 0) {
            ++decrease;
        }
    }

    if (//increase > static_cast<int>(totalsize * 0.6) &&
        decrease < static_cast<int>(totalsize * 0.01)) {
        return true;
    } else if (//decrease > static_cast<int>(totalsize * 0.6) &&
               increase < static_cast<int>(totalsize * 0.01)) {
        return true;
    } else {
        return false;
    }
}

// }}}


bool isOscillationUseOneVct(//std::vector<double>& x,
                            std::vector<double>& y)
// {{{ Function: isOscillation
//// Return true if exists Oscillation.
{
    //bool hasvalley = false;
    //bool haspeak   = false;
    bool isoscillation = false;
    std::vector<double>::iterator firstmax_i =
        std::max_element(y.begin(), y.end());
    double first_max = *firstmax_i;
    
    //std::vector<double>::iterator firstminafter_fm_i =
    //    std::min_element(firstmax_i + 1, y.end()); 
    //double first_min_after_fm = *firstminafter_fm_i;
    
    //std::vector<double>::iterator secondmax_i =
    //    std::max_element(firstminafter_fm_i, y.end());
    //double second_max = *secondmax_i;
    
    
    if (first_max > 1.05 * y.back()) {
        isoscillation = true;
    }
                                            
    //for (std::vector<double>::iterator vi = firstmax_i;
    //     vi != y.end() - 1;
    //     ++vi) {
    //    if (*vi < *(vi - 1) && *vi < *(vi + 1)) {
    //        hasvalley = true;
    //    }
    //} // for.
    if (isoscillation) {
        return true;
    }
    
    return false;
}
// }}}


void reactantName(std::vector<std::string>& name, int& reaction_type)
// {{{ Function: reactantName
//// Return the vector of strings in a kind of reaction.
{
    if (reaction_type == 1) {
    // 1s1p
        name.push_back("map3k");
        name.push_back("map3kp");
        name.push_back("map3k_ras");
        name.push_back("map3kp_m3kp");
        name.push_back("ras");
        name.push_back("m3kp");

    } else if (reaction_type == 2) {
    // 1s2p
        name.push_back("map3k");
        name.push_back("map3kp");
        name.push_back("map3k2p");
        name.push_back("map3k_ras");
        name.push_back("map3kp_m3kp");
        name.push_back("map3kp_ras");
        name.push_back("map3k2p_m3kp");
        name.push_back("ras");
        name.push_back("m3kp");       
    } else if (reaction_type == 3) {
    // 2s1p
        name.push_back("map3k");
        name.push_back("map3kp");
        name.push_back("map2k");
        name.push_back("map2kp");
        name.push_back("map3k_ras");
        name.push_back("map3kp_m3kp");
        name.push_back("map2k_map3kp");
        name.push_back("map2kp_m2kp");
        name.push_back("ras");
        name.push_back("m3kp");
        name.push_back("m2kp");
    } else if (reaction_type == 4) {
    // 2s2p
        name.push_back("map3k");
        name.push_back("map3kp");
        name.push_back("map2k");
        name.push_back("map2kp");
        name.push_back("map2k2p");
        name.push_back("map3k_ras");
        name.push_back("map3kp_m3kp");
        name.push_back("map2k_map3kp");
        name.push_back("map2kp_m2kp");
        name.push_back("map2kp_map3kp");
        name.push_back("map2k2p_m2kp");
        name.push_back("ras");
        name.push_back("m3kp");
        name.push_back("m2kp");
    } else if (reaction_type == 5) {
    // 3s1p
        name.push_back("map3k");
        name.push_back("map3kp");
        name.push_back("map2k");
        name.push_back("map2kp");
        name.push_back("mapk");
        name.push_back("mapkp");
        name.push_back("map3k_ras");
        name.push_back("map3kp_m3kp");
        name.push_back("map2k_map3kp");
        name.push_back("map2kp_m2kp");
        name.push_back("mapk_map2kp");
        name.push_back("mapkp_mkp");
        name.push_back("ras");
        name.push_back("m3kp");
        name.push_back("m2kp");
        name.push_back("mkp");
    } else if (reaction_type == 6) {
        name.push_back("map3k");
        name.push_back("map3kp");
        name.push_back("map2k");
        name.push_back("map2kp");
        name.push_back("map2k2p");
        name.push_back("mapk");
        name.push_back("mapkp");
        name.push_back("mapk2p");
        name.push_back("map3k_ras");
        name.push_back("map3kp_m3kp");
        name.push_back("map2k_map3kp");
        name.push_back("map2kp_m2kp");
        name.push_back("map2kp_map3kp");
        name.push_back("map2k2p_m2kp");
        name.push_back("mapk_map2k2p");
        name.push_back("mapkp_mkp");
        name.push_back("mapkp_map2k2p");
        name.push_back("mapk2p_mkp");
        name.push_back("ras");
        name.push_back("m3kp");
        name.push_back("m2kp");
        name.push_back("mkp");
    }
}
// }}}

////////////////////

// }}}


////////////////////
// Class: Parameter
////////////////////
// {{{ Class: Parameter

Parameter::Parameter(const int dimension)
// {{{ Constructor: Parameter
{
    function_num = dimension;
    f.reset(new double[function_num]);
    b.reset(new double[function_num]);
    kf.reset(new double[function_num]);
    kb.reset(new double[function_num]);

}
// }}}


Parameter::Parameter(const int dimension,
                     const std::unique_ptr<double[]>& input_f,
                     const std::unique_ptr<double[]>& input_b,
                     const std::unique_ptr<double[]>& input_kf,
                     const std::unique_ptr<double[]>& input_kb)
// {{{ Constructor: Parameter
{
    function_num = dimension;
    f.reset(new double[function_num]);
    b.reset(new double[function_num]);
    kf.reset(new double[function_num]);
    kb.reset(new double[function_num]);
 
    for (size_t i = 0; i != function_num; ++i) {
        f[i] = input_f[i];
        b[i] = input_b[i];
        kf[i] = input_kf[i];
        kb[i] = input_kb[i];
    }
    atp = calculateATP();
}
// }}}


Parameter::Parameter(const int dimension,
                     const double& input_f,
                     const double& input_b,
                     const double& input_kf,
                     const double& input_kb)
// {{{ Constructor: Parameter
{
    function_num = dimension;
    f.reset(new double[function_num]);
    b.reset(new double[function_num]);
    kf.reset(new double[function_num]);
    kb.reset(new double[function_num]);
 
    for (size_t i = 0; i != function_num; ++i) {
        f[i] = input_f;
        b[i] = input_b;
        kf[i] = input_kf;
        kb[i] = input_kb;
    }
    atp = calculateATP();
}
// }}}


Parameter::Parameter(const Parameter& para)
// {{{ Constructor: Parameter
{
    function_num = para.function_num;
    f.reset(new double[function_num]);
    b.reset(new double[function_num]);
    kf.reset(new double[function_num]);
    kb.reset(new double[function_num]);
 
    for (size_t i = 0; i != function_num; ++i) {
        f[i] = para.f[i];
        b[i] = para.b[i];
        kf[i] = para.kf[i];
        kb[i] = para.kb[i];
    }
    atp = calculateATP();
}
// }}}


Parameter& Parameter::operator=(const Parameter& para)
// {{{ operator=
//// Overload operator=.
{
    function_num = para.function_num;
    f.reset(new double[function_num]);
    b.reset(new double[function_num]);
    kf.reset(new double[function_num]);
    kb.reset(new double[function_num]);
 
    for (size_t i = 0; i != function_num; ++i) {
        f[i] = para.f[i];
        b[i] = para.b[i];
        kf[i] = para.kf[i];
        kb[i] = para.kb[i];
    }
    atp = para.atp;
    return *shared_from_this();
}   
// }}}


void Parameter::initParameter(const int& dimension)
// {{{ Member-Function: initParameter
//// initiate parameter no in declaration.
{
    function_num = dimension;
    f.reset(new double[function_num]);
    b.reset(new double[function_num]);
    kf.reset(new double[function_num]);
    kb.reset(new double[function_num]);
}
// }}}


void Parameter::initParameter(const int&    dimension,
                              const double& input_f,
                              const double& input_b,
                              const double& input_kf,
                              const double& input_kb)
// {{{ Member-Function: initParameter
//// initiate parameter no in declaration.
{
    function_num = dimension;
    f.reset(new double[function_num]);
    b.reset(new double[function_num]);
    kf.reset(new double[function_num]);
    kb.reset(new double[function_num]);
    for (size_t i = 0; i != function_num; ++i) {
        f[i] = input_f;
        b[i] = input_b;
        kf[i] = input_kf;
        kb[i] = input_kb;
    }
    atp = calculateATP();
}
// }}}


void Parameter::initParameter(const int& dimension,
                              const std::unique_ptr<double[]>& input_f,
                              const std::unique_ptr<double[]>& input_b,
                              const std::unique_ptr<double[]>& input_kf,
                              const std::unique_ptr<double[]>& input_kb)
// {{{ Member-Function: initParameter
//// initiate parameter no in declaration.
{
    function_num = dimension;
    f.reset(new double[function_num]);
    b.reset(new double[function_num]);
    kf.reset(new double[function_num]);
    kb.reset(new double[function_num]);
    for (size_t i = 0; i != function_num; ++i) {
        f[i] = input_f[i];
        b[i] = input_b[i];
        kf[i] = input_kf[i];
        kb[i] = input_kb[i];
    }
    atp = calculateATP();
}
// }}}


void Parameter::setParameter(const std::unique_ptr<double[]>& input_f,
                             const std::unique_ptr<double[]>& input_b,
                             const std::unique_ptr<double[]>& input_kf,
                             const std::unique_ptr<double[]>& input_kb)
// {{{ Member-Function: setParameter

//// set the value.
{
    if (function_num == 0) {
        std::cout << "function number is not set.";
        exit(1);
    }
    for (size_t i = 0; i != function_num; ++i) {
        f[i] = input_f[i];
        b[i] = input_b[i];
        kf[i] = input_kf[i];
        kb[i] = input_kb[i];
    }
    atp = calculateATP();
}

// }}}



void Parameter::setParameter(const double& input_f,
                             const double& input_b,
                             const double& input_kf,
                             const double& input_kb)
// {{{ Member-Function: setParameter

//// set the parameter array to a certain velue.
{
    if (function_num == 0) {
        std::cout << "function number is not set.";
        exit(1);
    }

    for (size_t i = 0; i != function_num; ++i) {
        f[i] = input_f;
        b[i] = input_b;
        kf[i] = input_kf;
        kb[i] = input_kb;
    }
    atp = calculateATP();
}

// }}}


void Parameter::setParameter(const std::unique_ptr<double[]>& input)
// {{{ Member-Function: setParameter
//// Set the parameter array to a certain velue.
{
    if (function_num == 0) {
        std::cout << "function number is not set.";
        exit(1);
    }

    for (size_t i = 0; i != function_num; ++i) {
        f[i] = input[0];
        b[i] = input[1];
        kf[i] = input[2];
        kb[i] = input[3];
    }
    atp = calculateATP();
}
// }}}


void Parameter::setParameter(const std::vector<double>& input)
// {{{ Member-Function: setParameter

//// Set the parameter array to a certain velue.
{
    if (function_num == 0 || function_num != input.size()) {
        std::cout << "function number is not set.";
        exit(1);
    }
    for (size_t i = 0; i != function_num; ++i) {
        f[i] = input[0];
        b[i] = input[1];
        kf[i] = input[2];
        kb[i] = input[3];
    }
    atp = calculateATP();
}

// }}}


void Parameter::printParameter()
// {{{ Member-Function: printParameter
//// print the values.
{
    for (size_t i = 0; i != function_num; ++i) {
        std::cout << i << ": ";
        std::cout << f[i] << "\t"
                  << b[i] << "\t"
                  << kf[i] << "\t"
                  << kb[i] << "\n";
    }
}
// }}}


double Parameter::calculateATP()
// {{{ Member-Function: calculateATP
//// Calculate the free energy of ATP with these parameters.    
{
    double energy = 0.0;
    for(int i = 0; i != function_num / 2.0; ++i ) {
        energy = energy +
            std::log10((f[2*i] * f[2*i+1] * kf[2*i] * kf[2*i+1]) /
                       (b[2*i] * b[2*i+1] * kb[2*i] * kb[2*i+1]));
    }
    return energy;
}
// }}}

////////////////////

// }}}


////////////////////
// Class: Concentration
////////////////////
// {{{ Class: Concentration
Concentration::Concentration(const size_t& dimension)
// {{{ Constructor: Concentration
//// Constructor of Class Concentration.
{
    reactant_num = dimension;
    reactant.reset(new double[dimension]);
    for (size_t i = 0; i != reactant_num; ++i) {
        reactant[i] = 0.0;
    }
    is_new = true;
    //std::cout << "Concentration constructor:\n dimension " << dimension
    //        << "\nreactant_num " << reactant_num << "\n";
}
// }}}


Concentration::Concentration(const size_t& dimension,
                             const std::unique_ptr<double[]>& y)
// {{{ Constructor: Concentration
//// Constructor of Class Concentration.
{
    reactant_num = dimension;
    reactant.reset(new double [dimension]);
    for (size_t i = 0; i != reactant_num; ++i) {
        reactant[i] = y[i];
    }
    is_new = false;
}
// }}}


Concentration::Concentration(const Concentration& conc)
// {{{ Constructor: Concentration
//// Constructor of Class Concentration.
{
    is_new = false;
    reactant_num = conc.reactant_num;
    reactant.reset(new double [reactant_num]);
    for (size_t i = 0; i != reactant_num; ++i) {
        reactant[i] = conc.reactant[i];
    } // for
}
// }}}


Concentration& Concentration::operator=(const Concentration& conc)
// {{{ Operator: =
//// Operator = overload.
{
    is_new = false;
    reactant_num = conc.reactant_num;
    reactant.reset(new double [reactant_num]);
    for (size_t i = 0; i != reactant_num; ++i) {
        reactant[i] = conc.reactant[i];
    } // for
    return *shared_from_this();
}
// }}}


void Concentration::setConcentration(const std::unique_ptr<double[]>& y)
// {{{ Member-Function: setConcentration
//// Set the value of Concentrations of reactant.
{
    if (reactant_num == 0) {
        std::cout << "Concentration reactant_num is 0.\n";
        exit(1);
    }
    for (size_t i = 0; i != reactant_num; ++i) {
        reactant[i] = y[i];
    }
    is_new = false;
}
// }}}


void Concentration::setConcentration(double y[])
// {{{ Member-Function: setConcentration
//// Set the value of Concentrations of reactant.
{
    if (reactant_num == 0) {
        std::cout << "Concentration reactant_num is 0.\n";
        exit(1);
    }
    for (size_t i = 0; i != reactant_num; ++i) {
        reactant[i] = y[i];
    }
    is_new = false;
}
// }}}


void Concentration::setConcentration(const size_t& dimension,
                                     const std::unique_ptr<double[]>& y)
// {{{ Member-Function: setConcentration
//// Set the value of Concentrations of reactant.
{
    reactant_num = dimension;
    reactant.reset(new double [reactant_num]);
    for (size_t i = 0; i != reactant_num; ++i) {
        reactant[i] = y[i];
    }
    is_new = false;
}
// }}}


void Concentration::printConcentration()
// {{{ Member-Function: printConcentration
//// Print the value to screen.
{
    std::cout << "Value: ";
    for (size_t i = 0; i != reactant_num; ++i) {
        std::cout << reactant[i] << "\t";
    }
    std::cout << std::endl;

}
// }}}


void Concentration::reset()
// {{{ Member-Function: reset
//// Delete the array.
{
    if (!is_new) {
        reactant.release();
    }
    is_new = true;
}
// }}}


size_t Concentration::size()
// {{{ Member-Function: Concentration::size
//// Return the size of reactant.
{
    return reactant_num;
}
// }}}


double Concentration::get(const size_t& ord)
// {{{ Member-Function: Concentration::get
//// Get value.
{
    if (ord >= reactant_num) {
        std::cout << "Concentration: subscript operator out of bounds.\n"
                  << "reactant_num: " << reactant_num
                  << "\norder: " << ord
                  << "\n";
        exit(1);
    }
    return reactant[ord];
}
// }}}

////////////////////
// }}}


////////////////////
// Class: ReactantCocentration
////////////////////
// {{{ Class: ReactantConcentration

ReactantConcentration::ReactantConcentration(const int& type)
// {{{ Constructor: ReactantConcentration
{
    reaction_type = type;
    reach_stable = true;
    //reactantName();
    
    funcNumber();
    finalOrder();
}
// }}}


bool ReactantConcentration::isStable()
// {{{ Member-Function: isStable
//// Whether the ordinary differential equations reach stable,
//// return 1 when stable, return 0 when not stable.
{
    std::vector<Concentration>::size_type vct_len = list.size();
    std::unique_ptr<double[]> x_array(new double[vct_len]);
    std::unique_ptr<double[]> y_array(new double[vct_len]);
    for (size_t si = 0; si != vct_len; si++) {
        x_array[si] = time[si];
        y_array[si] = list[si].get(final_out);
    } // for.
    return isStableUseEndSome(x_array, y_array, vct_len, 0.02, 5);
}
// }}}


bool ReactantConcentration::isMonotone()
// {{{ Member-Function: isMonotone
//// Check whether result monotone.
{
    std::vector<double> output;
    for (size_t si = 0; si != list.size(); ++si) {
        output.push_back(list[si].get(final_out));
    }
    return isMonotoneUseTwoVct(time, output);

}
// }}}


bool ReactantConcentration::isOscillation()
// {{{ Member-Function: isOscillation
//// Check whether exist oscillation.
{
    std::vector<double> output;
    for (size_t si = 0; si != list.size(); ++si) {
        output.push_back(list[si].get(final_out));
    } // for.
    return isOscillationUseOneVct(output);
}
// }}}


bool ReactantConcentration::reachStable()
// {{{ Member-Function: reachStable
//// Check whether reach stable.
{
    return reach_stable;
}
// }}}


void ReactantConcentration::outputFinal(std::ofstream& outfile)
// {{{ Member-Function: outputStable
//// Output the stable situation to file.
{
    if (reach_stable == true && isMonotone()) {
        for (std::vector<Concentration>::size_type vi = 0;
             vi != list.size();
             ++vi) {
            outfile << time[vi];
            for (size_t i = 0; i != reactant_num; ++i) {
                outfile << "," << list[vi].get(i);
            }
            outfile << "," << dissipation[vi];
            outfile << "\n";
        }
    }
}
// }}}


double ReactantConcentration::getFinalOutput()
// {{{ Member-Function: getFinalOutput
//// Return the output reactant concentration values at stable stage.
{

    return list.back().get(final_out);
}
// }}}


double ReactantConcentration::getFinalTime()
// {{{ Member-Function: getFinalTime
//// Return the time at stable stage.
{

    return time.back();
}
// }}}


double ReactantConcentration::getFinalDissipation()
// {{{ Member-Function: getFinalDissipation
//// Return the energy dissipation at stable stage.
{

    return dissipation.back();
}
// }}}


void ReactantConcentration::odeRun(double     start_t,
                                   double     end_t,
                                   Parameter* params,
                                   double*    reactant,
                                   double     pacelen)
// {{{ Member-Function: odeRun
//// Using gsl to solve ode.
{
    // initiate.
    double* y = reactant;
    double time_span = end_t - start_t;    // used to reset start_t and end_t.
    bool continue_flag = true;
    typedef int (* Func)(double       t,
                         const double y[],
                         double       f[],
                         void*        para);
    typedef int (* Jac)(double        t,
                        const double  y[],
                        double        dfdy[],
                        double        dydt[],
                        void*         para);

    Func odeFunction = 0;    // function pointer.
    Jac odeJacobian = 0;    // function pointer.

    // choose function  by reaction_type.
    if (reaction_type == 1) {
        odeFunction = odeFunction1s1p;
        odeJacobian = odeJacobian1s1p;
    } else if (reaction_type == 2) {
        odeFunction = odeFunction1s2p;
        odeJacobian = odeJacobian1s2p;
    } else if (reaction_type == 3) {
        odeFunction = odeFunction2s1p;
        odeJacobian = odeJacobian2s1p;
    } else if (reaction_type == 4) {
        odeFunction = odeFunction2s2p;
        odeJacobian = odeJacobian2s2p;
    } else if (reaction_type == 5) {
        odeFunction = odeFunction3s1p;
        odeJacobian = odeJacobian3s1p;
    } else if (reaction_type == 6) {
        odeFunction = odeFunction3s2p;
        odeJacobian = odeJacobian3s2p;
    }// HERE!!!!!!!!!!!!!!!!!!!!!!signal handler.
    gsl_odeiv2_system sys = {odeFunction,
                             odeJacobian,
                             reactant_num,
                             &*params};
    gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys,
                                                         gsl_odeiv2_step_rk8pd,
                                                         1e-6,
                                                         1e-6,
                                                         0.0);
    

    int count_while = 0;
    // calculate.
    do {
        Concentration reaction(reactant_num);    // temporary store y.
        long times = static_cast<long>((end_t - start_t)/pacelen);
        for (long i = 0; i != times; ++i) {
            double ti = pacelen + start_t;
            int    status = gsl_odeiv2_driver_apply(d, &start_t, ti, y);
            // start_t and y will be change after one run.
            if (status != GSL_SUCCESS) {
                std::cout << "Error, return value = " << status << "\n";
                reach_stable = false;
                continue_flag = false;
                break;
            }
            reaction.setConcentration(y);
            // set the concentrate of reactant to y[].
            time.push_back(ti);
            list.push_back(reaction);    // record the situation at ti.
            dissipation.push_back(calculateDissipation(y, params));
            if (isStable()) {
                reach_stable = true;
                continue_flag = false;
                break;
            }
        }


        if (isStable()) {
            // if continue_flag is false, end while loop.
            break;
        }
        
        if (count_while > 100 * times) {
            reach_stable = false;
            break;
        }
       
        end_t = start_t + time_span;
        ++count_while;
    } while (continue_flag);

    gsl_odeiv2_driver_free(d);
}
// }}}


double ReactantConcentration::calculateDissipation(double yarray[],
                                                   void*  param)
// {{{ Member-Function: calculateDissipation
//// Energy dissipation of each time.
{
   double  energy = 0;
   double* probarray = new double[function_num * 4];

   typedef void (*propensity)(const double y[],
                              double       prob[],
                              void*        para);
   propensity profunc = 0;
   // choose function by reaction_type.
   if (reaction_type == 1) {
       profunc = propensityFunction1s1p;
   } else if (reaction_type == 2) {
       profunc = propensityFunction1s2p;
   } else if (reaction_type == 3) {
       profunc = propensityFunction2s1p;
   } else if (reaction_type == 4) {
       profunc = propensityFunction2s2p;
   } else if (reaction_type == 5) {
       profunc = propensityFunction3s1p;
   } else if (reaction_type == 6) {
       profunc = propensityFunction3s2p;
   }
   profunc(yarray, probarray, param);
   // calculate energy dissipation.
   for (int i = 0; i != static_cast<int>(function_num / 2); ++i) {
       for (int j = 0; j != 4; ++j) {
           if ((probarray[8 * i + 2 * j] > 0) &
               (probarray[8 * i + 2 * j + 1] > 0)) {
               energy = energy +
                   (
                    probarray[8 * i + 2 * j] - probarray[8 * i + 2 * j + 1]
                    ) *
                   log(probarray[8 * i + 2 * j] / probarray[8 * i + 2 * j + 1]);
           }
       }
   }
   delete [] probarray;
   return energy;
}  
// }}}


void ReactantConcentration::finalOrder()
// {{{ Member-Function: finalOrder
//// Set the number of final output order according to reaction type.
{
    switch (reaction_type) {
    case 1:
        final_out = 1;
        break;
    case 2:
        final_out = 2;
        break;
    case 3:
        final_out = 3;
        break;
    case 4:
        final_out = 4;
        break;
    case 5:
        final_out = 5;
        break;
    case 6:
        final_out = 7;
        break;
    default:
        std::cout << "Wrong input of reaction type in Number!"
                  << std::endl;
    }
}
// }}}


void ReactantConcentration::funcNumber()
// {{{ Member-Function: funcumber
//// Return the number of chemical reaction according to reaction type.
{
    switch (reaction_type) {
    case 1:
        reactant_num = 6;
        function_num = 2;
        break;
    case 2:
        reactant_num = 9;
        function_num = 4;
        break;
    case 3:
        reactant_num = 11;
        function_num = 4;
        break;
    case 4:
        reactant_num = 14;
        function_num = 6;
        break;
    case 5:
        reactant_num = 16;
        function_num = 6;
        break;
    case 6:
        reactant_num = 22;
        function_num = 10;
        break;
    default:
        std::cout << "Wrong input of reaction type in chemNumber!"
                  << std::endl;
    }
}
// }}}


int ReactantConcentration::odeFunction1s1p(double t,
                                           const double y[],
                                           double f[],
                                           void *para)
// {{{ Member-Function: odeFunction1s1p
//// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k       = y[0];    // MAPKKK
    double map3kp      = y[1];    // MAPKKK-phosphoryl
    double map3k_ras   = y[2];    // intermedia product
    double map3kp_m3kp = y[3];    // intermedia product
    double ras         = y[4];    // MAPKKK kinase
    double m3kp        = y[5];    // MAPKKK-p phosphatase 
    
    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map3k
    f[0] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // map3kp
    f[1] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras;
    // map3k_ras
    f[2] = params.f[0] * map3k * ras - params.b[0] * map3k_ras -
        params.kf[0] * map3k_ras + params.kb[0] * map3kp * ras;
    // map3kp_m3kp
    f[3] = params.f[1] * map3kp * m3kp - params.b[1] * map3kp_m3kp -
        params.kf[1] * map3kp_m3kp + params.kb[1] * map3k * m3kp;
    // ras: dras = -dmap3k_ras
    f[4] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras;
    // m3kp: dm3kp = -dmap3kp_m3kp
    f[5] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;    
    
    return GSL_SUCCESS;    // in gsl header.
}
// }}}


int ReactantConcentration::odeJacobian1s1p(double       t,
                                           const double y[],
                                           double       dfdy[],
                                           double       dfdt[],
                                           void*        para)
// {{{ Member-Function: odeJacobian1s1p

//// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    // initiate.
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k       = y[0];    // MAPKKK
    double map3kp      = y[1];    // MAPKKK-phosphoryl
    double map3k_ras   = y[2];    // intermedia product
    map3k_ras          = y[2];    // get rid of warning
    double map3kp_m3kp = y[3];    // intermedia product
    map3kp_m3kp        = y[3];    // get rid of warning
    double ras         = y[4];    // MAPKKK kinase
    double m3kp        = y[5];    // MAPKKK-p phosphatase 
    
    
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 6, 6);
    gsl_matrix*            m = &dfdy_mat.matrix;
    
    // calculate matrix.
    // map3k.
    gsl_matrix_set(m, 0, 0, -params.f[0] * ras -
                   params.kb[1] * m3kp);
    gsl_matrix_set(m, 0, 1, 0.0);
    gsl_matrix_set(m, 0, 2, params.b[0]);
    gsl_matrix_set(m, 0, 3, params.kf[1]);
    gsl_matrix_set(m, 0, 4, -params.f[0] * map3k);
    gsl_matrix_set(m, 0, 5, -params.kb[1] * map3k);
    // map3kp.
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, -params.f[1] * m3kp -
                   params.kb[0] * ras);
    gsl_matrix_set(m, 1, 2, params.kf[0]);
    gsl_matrix_set(m, 1, 3, params.b[1]);
    gsl_matrix_set(m, 1, 4, -params.kb[0] * map3kp);
    gsl_matrix_set(m, 1, 5, -params.f[1] * map3kp);
    // map3k_ras.
    gsl_matrix_set(m, 2, 0, params.f[0] * ras);
    gsl_matrix_set(m, 2, 1, params.kb[0] * ras);
    gsl_matrix_set(m, 2, 2, -params.b[0] - params.kf[0]);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 2, 4, params.f[0] * map3k +
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 2, 5, 0.0);
    // map3kp_m3kp.
    gsl_matrix_set(m, 3, 0, params.kb[1] * m3kp);
    gsl_matrix_set(m, 3, 1, params.f[1] * m3kp);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -params.b[1] - params.kf[1]);
    gsl_matrix_set(m, 3, 4, 0.0);
    gsl_matrix_set(m, 3, 5, params.f[1] * map3kp +
                   params.kb[1] * map3k);
    // ras.
    gsl_matrix_set(m, 4, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 4, 1, -params.kb[0] * ras);
    gsl_matrix_set(m, 4, 2, params.b[0] + params.kf[0]);
    gsl_matrix_set(m, 4, 3, 0.0);
    gsl_matrix_set(m, 4, 4, -params.f[0] * map3k -
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 4, 5, 0.0);
    // m3kp.
    gsl_matrix_set(m, 5, 0, -params.kb[1] * m3kp);
    gsl_matrix_set(m, 5, 1, -params.f[1] * m3kp);
    gsl_matrix_set(m, 5, 2, 0.0);
    gsl_matrix_set(m, 5, 3, params.b[1] + params.kf[1]);
    gsl_matrix_set(m, 5, 4, 0.0);
    gsl_matrix_set(m, 5, 5, -params.f[1] * map3kp -
                   params.kb[1] * map3k);

    for (int i = 0; i != 6; ++i) {
        dfdt[i] = 0.0;
    }

    return GSL_SUCCESS;
}

// }}}


int ReactantConcentration::odeFunction1s2p(double       t,
                                           const double y[],
                                           double       f[],
                                           void*        para)
// {{{ Member-Function: odeFunction1s2p
//// Needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k        = y[0];    // MAPKKK
    double map3kp       = y[1];    // MAPKKK-phosphoryl
    double map3k2p      = y[2];    // MAPKKK-phosphoryl-phosphoryl
    double map3k_ras    = y[3];    // intermedia product
    double map3kp_m3kp  = y[4];    // intermedia product
    double map3kp_ras   = y[5];    // intermedia product
    double map3k2p_m3kp = y[6];    // intermedia product
    double ras          = y[7];    // MAPKKK kinase
    double m3kp         = y[8];    // MAPKKK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map3kp + ras <b2==f2> map3kp_ras <kb2-=kf2> map3k2p + ras
    // map3k2p + m3kp <b3==f3> map3k2p_m3kp <kb3-=kf3> map3kp + m3kp
    
    // map3k
    f[0] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // map3kp
    f[1] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras -
        params.f[2] * map3kp * ras + params.b[2] * map3kp_ras +
        params.kf[3] * map3k2p_m3kp - params.kb[3] * map3kp * m3kp;
    // map3k2p
    f[2] = -params.f[3] * map3k2p * m3kp + params.b[3] * map3k2p_m3kp +
        params.kf[2] * map3kp_ras - params.kb[2] * map3k2p * ras;
    // map3k_ras
    f[3] = params.f[0] * map3k * ras - params.b[0] * map3k_ras -
        params.kf[0] * map3k_ras + params.kb[0] * map3kp * ras;
    // map3kp_m3kp
    f[4] = params.f[1] * map3kp * m3kp - params.b[1] * map3kp_m3kp -
        params.kf[1] * map3kp_m3kp + params.kb[1] * map3k * m3kp;
    // map3kp_ras
    f[5] = params.f[2] * map3kp * ras - params.b[2] * map3kp_ras -
        params.kf[2] * map3kp_ras + params.kb[2] * map3k2p * ras;
    // map3k2p_m3kp
    f[6] = params.f[3] * map3k2p * m3kp - params.b[3] * map3k2p_m3kp -
        params.kf[3] * map3k2p_m3kp + params.kb[3] * map3kp * m3kp;
    // ras: dras = - dmap3k_ras - dma3kp_ras
    f[7] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras -
        params.f[2] * map3kp * ras + params.b[2] * map3kp_ras +
        params.kf[2] * map3kp_ras - params.kb[2] * map3k2p * ras;
    // m3kp: dm3kp = - dmap3kp_m3kp - dmap3k2p_m3kp
    f[8] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp -
        params.f[3] * map3k2p * m3kp + params.b[3] * map3k2p_m3kp +
        params.kf[3] * map3k2p_m3kp - params.kb[3] * map3kp * m3kp;
    
    return GSL_SUCCESS;    // in gsl header.
}
// }}}


int ReactantConcentration::odeJacobian1s2p(double       t,
                                           const double y[],
                                           double       dfdy[],
                                           double       dfdt[],
                                           void*        para)
// {{{ Member-Function: odeJacobian1s2p

//// needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    // initiate.
    Parameter params = *static_cast<Parameter *>(para);
    
    double map3k        = y[0];    // MAPKKK
    double map3kp       = y[1];    // MAPKKK-phosphoryl
    double map3k2p      = y[2];    // MAPKKK-phosphoryl-phosphoryl
    double map3k_ras    = y[3];    // intermedia product
    map3k_ras           = y[3];    // get rid of warning
    double map3kp_m3kp  = y[4];    // intermedia product
    map3kp_m3kp         = y[4];    // get rid of warning
    double map3kp_ras   = y[5];    // intermedia product
    map3kp_ras          = y[5];    // get rid of warning
    double map3k2p_m3kp = y[6];    // intermedia product
    map3k2p_m3kp        = y[6];    // get rid of warning
    double ras          = y[7];    // MAPKKK kinase
    double m3kp         = y[8];    // MAPKKK-p phosphatase
    
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 9, 9);
    gsl_matrix *m = &dfdy_mat.matrix;

    // calculate matrix.
    // map3k.
    gsl_matrix_set(m, 0, 0, -params.f[0] * ras -
                   params.kb[1] * m3kp);
    gsl_matrix_set(m, 0, 1, 0.0);
    gsl_matrix_set(m, 0, 2, 0.0);
    gsl_matrix_set(m, 0, 3, params.b[0]);
    gsl_matrix_set(m, 0, 4, params.kf[1]);
    gsl_matrix_set(m, 0, 5, 0.0);
    gsl_matrix_set(m, 0, 6, 0.0);
    gsl_matrix_set(m, 0, 7, -params.f[0] * map3k);
    gsl_matrix_set(m, 0, 8, -params.kb[1] * map3k);
    // map3kp.
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, -params.f[1] * m3kp -
                   params.kb[0] * ras -
                   params.f[2] * ras -
                   params.kb[3] * m3kp);
    gsl_matrix_set(m, 1, 2, 0.0);
    gsl_matrix_set(m, 1, 3, params.kf[0]);
    gsl_matrix_set(m, 1, 4, params.b[1]);
    gsl_matrix_set(m, 1, 5, params.b[2]);
    gsl_matrix_set(m, 1, 6, params.kf[3]);
    gsl_matrix_set(m, 1, 7, -params.kb[0] * map3kp -
                   params.f[2] * map3kp);
    gsl_matrix_set(m, 1, 8, -params.f[1] * map3kp -
                   params.kb[3] * map3kp);
    // map3k2p.
    gsl_matrix_set(m, 2, 0, 0.0);
    gsl_matrix_set(m, 2, 1, 0.0);
    gsl_matrix_set(m, 2, 2, -params.f[3] * m3kp -
                   params.kb[2] * ras);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 2, 4, 0.0);
    gsl_matrix_set(m, 2, 5, params.kf[2]);
    gsl_matrix_set(m, 2, 6, params.b[3]);
    gsl_matrix_set(m, 2, 7, -params.kb[2] * map3k2p);
    gsl_matrix_set(m, 2, 8, -params.f[3] * map3k2p);
    // map3k_ras.
    gsl_matrix_set(m, 3, 0, params.f[0] * ras);
    gsl_matrix_set(m, 3, 1, params.kb[0] * ras);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -params.b[0] - params.kf[0]);
    gsl_matrix_set(m, 3, 4, 0.0);
    gsl_matrix_set(m, 3, 5, 0.0);
    gsl_matrix_set(m, 3, 6, 0.0);
    gsl_matrix_set(m, 3, 7, params.f[0] * map3k +
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 3, 8, 0.0);
    // map3kp_m3kp.
    gsl_matrix_set(m, 4, 0, params.kb[1] * m3kp);
    gsl_matrix_set(m, 4, 1, params.f[1] * m3kp);
    gsl_matrix_set(m, 4, 2, 0.0);
    gsl_matrix_set(m, 4, 3, 0.0);
    gsl_matrix_set(m, 4, 4, -params.b[1] - params.kf[1]);
    gsl_matrix_set(m, 4, 5, 0.0);
    gsl_matrix_set(m, 4, 6, 0.0);
    gsl_matrix_set(m, 4, 7, 0.0);
    gsl_matrix_set(m, 4, 8, params.f[1] * map3kp +
                   params.kb[1] * map3k);
    // map3kp_ras.
    gsl_matrix_set(m, 5, 0, 0.0);
    gsl_matrix_set(m, 5, 1, params.f[2] * ras);
    gsl_matrix_set(m, 5, 2, params.kb[2] * ras);
    gsl_matrix_set(m, 5, 3, 0.0);
    gsl_matrix_set(m, 5, 4, 0.0);
    gsl_matrix_set(m, 5, 5, -params.b[2] - params.kf[2]);
    gsl_matrix_set(m, 5, 6, 0.0);
    gsl_matrix_set(m, 5, 7, params.f[2] * map3kp +
                   params.kb[2] * map3k2p);
    gsl_matrix_set(m, 5, 8, 0.0);
    // map3k2p_m3kp.
    gsl_matrix_set(m, 6, 0, 0.0);
    gsl_matrix_set(m, 6, 1, params.kb[3] * m3kp);
    gsl_matrix_set(m, 6, 2, params.f[3] * m3kp);
    gsl_matrix_set(m, 6, 3, 0.0);
    gsl_matrix_set(m, 6, 4, 0.0);
    gsl_matrix_set(m, 6, 5, 0.0);
    gsl_matrix_set(m, 6, 6, -params.b[3] - params.kf[3]);
    gsl_matrix_set(m, 6, 7, 0.0);
    gsl_matrix_set(m, 6, 8, params.f[3] * map3k2p +
                   params.kb[3] * map3kp);
    // ras.
    gsl_matrix_set(m, 7, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 7, 1, -params.kb[0] * ras -
                   params.f[2] * ras);
    gsl_matrix_set(m, 7, 2, -params.kb[2] * ras);
    gsl_matrix_set(m, 7, 3, params.b[0] + params.kf[0]);
    gsl_matrix_set(m, 7, 4, 0.0);
    gsl_matrix_set(m, 7, 5, params.b[2] + params.kf[2]);
    gsl_matrix_set(m, 7, 6, 0.0);
    gsl_matrix_set(m, 7, 7, -params.f[0] * map3k -
                   params.kb[0] * map3kp -
                   params.f[2] * map3kp -
                   params.kb[2] * map3k2p);
    gsl_matrix_set(m, 7, 8, 0.0);
    // m3kp.
    gsl_matrix_set(m, 8, 0, -params.kb[1] * m3kp);
    gsl_matrix_set(m, 8, 1, -params.f[1] * m3kp -
                   params.kb[3] * m3kp);
    gsl_matrix_set(m, 8, 2, -params.f[3] * m3kp);
    gsl_matrix_set(m, 8, 3, 0.0);
    gsl_matrix_set(m, 8, 4, params.b[1] + params.kf[1]);
    gsl_matrix_set(m, 8, 5, 0.0);
    gsl_matrix_set(m, 8, 6, params.b[3] + params.kf[3]);
    gsl_matrix_set(m, 8, 7, 0.0);
    gsl_matrix_set(m, 8, 8, -params.f[1] * map3kp -
                   params.kb[1] * map3k -
                   params.f[3] * map3k2p -
                   params.kb[3] * map3kp);

    for (int i = 0; i != 9; ++i) {
        dfdt[i] = 0.0;
    }

    return GSL_SUCCESS;
}

// }}}


int ReactantConcentration::odeFunction2s1p(double       t,
                                           const double y[],
                                           double       f[],
                                           void*        para)
// {{{ Member-Function: odeFunction2s1p
//// needed by <gsl/gsl_odeiv2.h>, 2 stages, 1 phosphorylation.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k        = y[0];    // MAPKKK
    double map3kp       = y[1];    // MAPKKK-phosphoryl
    double map2k        = y[2];    // MAPKK
    double map2kp       = y[3];    // MAPKK-phosphoryl
    double map3k_ras    = y[4];    // intermedia product
    double map3kp_m3kp  = y[5];    // intermedia product
    double map2k_map3kp = y[6];    // intermedia product
    double map2kp_m2kp  = y[7];    // intermedia product
    double ras          = y[8];    // MAPKKK kinase
    double m3kp         = y[9];    // MAPKKK-p phosphatase
    double m2kp         = y[10];   // MAPKK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map2k + map3kp <b2==f2> map2k_map3kp <kb2-=kf2> map2kp + map3kp
    // map2kp + m2kp <b3==f3> map2kp_m2kp <kb3-=kf3> map2k + m2kp
    
    // map3k
    f[0] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // map3kp
    f[1] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras -
        params.f[2] * map2k * map3kp + params.b[2] * map2k_map3kp +
        params.kf[2] * map2k_map3kp - params.kb[2] * map2kp * map3kp;
    // map2k
    f[2] = -params.f[2] * map2k * map3kp + params.b[2] * map2k_map3kp +
        params.kf[3] * map2kp_m2kp - params.kb[3] * map2k * m2kp;
    // map2kp
    f[3] = -params.f[3] * map2kp * m2kp + params.b[3] * map2kp_m2kp +
        params.kf[2] * map2k_map3kp - params.kb[2] * map2kp * map3kp;
    // map3k_ras
    f[4] = params.f[0] * map3k * ras - params.b[0] * map3k_ras -
        params.kf[0] * map3k_ras + params.kb[0] * map3kp * ras;
    // map3kp_m3kp
    f[5] = params.f[1] * map3kp * m3kp - params.b[1] * map3kp_m3kp -
        params.kf[1] * map3kp_m3kp + params.kb[1] * map3k * m3kp;
    // map2k_map3kp
    f[6] = params.f[2] * map2k * map3kp - params.b[2] * map2k_map3kp -
        params.kf[2] * map2k_map3kp + params.kb[2] * map2kp * map3kp;
    // map2kp_m2kp
    f[7] = params.f[3] * map2kp * m2kp - params.b[3] * map2kp_m2kp -
        params.kf[3] * map2kp_m2kp + params.kb[3] * map2k * m2kp;
    // ras: dras = -dmap3k_ras
    f[8] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras;
    // m3kp: dm3kp = -dmap3kp_m3kp
    f[9] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // m2kp: dm2kp = -dmap2kp_m2kp
    f[10] = -params.f[3] * map2kp * m2kp + params.b[3] * map2kp_m2kp +
        params.kf[3] * map2kp_m2kp - params.kb[3] * map2k * m2kp;
    
    return GSL_SUCCESS;    // in gsl header.
}
// }}}


int ReactantConcentration::odeJacobian2s1p(double       t,
                                           const double y[],
                                           double       dfdy[],
                                           double       dfdt[],
                                           void*        para)
// {{{ Member-Function: odeJacobian2s1p

//// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 1 phosphorylation.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    // initiate.
    Parameter params = *static_cast<Parameter *>(para);
    
    double map3k        = y[0];    // MAPKKK
    double map3kp       = y[1];    // MAPKKK-phosphoryl
    double map2k        = y[2];    // MAPKK
    double map2kp       = y[3];    // MAPKK-phosphoryl
    double map3k_ras    = y[4];    // intermedia product
    map3k_ras           = y[4];    // get rid of warning
    double map3kp_m3kp  = y[5];    // intermedia product
    map3kp_m3kp         = y[5];    // get rid of warning
    double map2k_map3kp = y[6];    // intermedia product
    map2k_map3kp        = y[6];    // get rid of warning
    double map2kp_m2kp  = y[7];    // intermedia product
    map2kp_m2kp         = y[7];    // get rid of warning
    double ras          = y[8];    // MAPKKK kinase
    double m3kp         = y[9];    // MAPKKK-p phosphatase
    double m2kp         = y[10];   // MAPKK-p phosphatase
    
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 11, 11);
    gsl_matrix *m = &dfdy_mat.matrix;

    // calculate matrix.
    // map3k.
    gsl_matrix_set(m, 0, 0, -params.f[0] * ras -
                   params.kb[1] * m3kp);
    gsl_matrix_set(m, 0, 1, 0.0);
    gsl_matrix_set(m, 0, 2, 0.0);    
    gsl_matrix_set(m, 0, 3, 0.0);
    gsl_matrix_set(m, 0, 4, params.b[0]);
    gsl_matrix_set(m, 0, 5, params.kf[1]);
    gsl_matrix_set(m, 0, 6, 0.0);
    gsl_matrix_set(m, 0, 7, 0.0);
    gsl_matrix_set(m, 0, 8, -params.f[0] * map3k);
    gsl_matrix_set(m, 0, 9, -params.kb[1] * map3k);
    gsl_matrix_set(m, 0, 10, 0.0);
    // map3kp.
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, -params.f[1] * m3kp -
                   params.kb[0] * ras -
                   params.f[2] * map2k -
                   params.kb[2] * map2kp);
    gsl_matrix_set(m, 1, 2, -params.f[2] * map3kp);
    gsl_matrix_set(m, 1, 3, -params.kb[2] * map3kp);
    gsl_matrix_set(m, 1, 4, params.kf[0]);
    gsl_matrix_set(m, 1, 5, params.b[1]);
    gsl_matrix_set(m, 1, 6, params.b[2] + params.kf[2]);
    gsl_matrix_set(m, 1, 7, 0.0);
    gsl_matrix_set(m, 1, 8, -params.kb[0] * map3kp);
    gsl_matrix_set(m, 1, 9, -params.f[1] * map3kp);
    gsl_matrix_set(m, 1, 10, 0.0);
    // map2k.
    gsl_matrix_set(m, 2, 0, 0.0);
    gsl_matrix_set(m, 2, 1, -params.f[2] * map2k);
    gsl_matrix_set(m, 2, 2, -params.f[2] * map3kp -
                   params.kb[3] * m2kp);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 2, 4, 0.0);
    gsl_matrix_set(m, 2, 5, 0.0);
    gsl_matrix_set(m, 2, 6, params.b[2]);
    gsl_matrix_set(m, 2, 7, params.kf[3]);
    gsl_matrix_set(m, 2, 8, 0.0);
    gsl_matrix_set(m, 2, 9, 0.0);
    gsl_matrix_set(m, 2, 10, -params.kb[3] * map2k);
    // map2kp.
    gsl_matrix_set(m, 3, 0, 0.0);
    gsl_matrix_set(m, 3, 1, -params.kb[2] * map2kp);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -params.f[3] * m2kp -
                   params.kb[2] * map3kp);
    gsl_matrix_set(m, 3, 4, 0.0);
    gsl_matrix_set(m, 3, 5, 0.0);
    gsl_matrix_set(m, 3, 6, params.kf[2]);
    gsl_matrix_set(m, 3, 7, params.b[3]);
    gsl_matrix_set(m, 3, 8, 0.0);
    gsl_matrix_set(m, 3, 9, 0.0);
    gsl_matrix_set(m, 3, 10, -params.f[3] * map2kp);
    // map3k_ras.
    gsl_matrix_set(m, 4, 0, params.f[0] * ras);
    gsl_matrix_set(m, 4, 1, params.kb[0] * ras);
    gsl_matrix_set(m, 4, 2, 0.0);
    gsl_matrix_set(m, 4, 3, 0.0);
    gsl_matrix_set(m, 4, 4, -params.b[0] - params.kf[0]);
    gsl_matrix_set(m, 4, 5, 0.0);
    gsl_matrix_set(m, 4, 6, 0.0);
    gsl_matrix_set(m, 4, 7, 0.0);
    gsl_matrix_set(m, 4, 8, params.f[0] * map3k +
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 4, 9, 0.0);
    gsl_matrix_set(m, 4, 10, 0.0);           
    // map3kp_m3kp.
    gsl_matrix_set(m, 5, 0, params.kb[1] * m3kp);
    gsl_matrix_set(m, 5, 1, params.f[1] * m3kp);
    gsl_matrix_set(m, 5, 2, 0.0);
    gsl_matrix_set(m, 5, 3, 0.0);
    gsl_matrix_set(m, 5, 4, 0.0);
    gsl_matrix_set(m, 5, 5, -params.b[1] - params.kf[1]);
    gsl_matrix_set(m, 5, 6, 0.0);
    gsl_matrix_set(m, 5, 7, 0.0);
    gsl_matrix_set(m, 5, 8, 0.0);
    gsl_matrix_set(m, 5, 9, params.f[1] * map3kp +
                   params.kb[1] * map3k);
    gsl_matrix_set(m, 5, 10, 0.0);           
    // map2k_map3kp.
    gsl_matrix_set(m, 6, 0, 0.0);
    gsl_matrix_set(m, 6, 1, params.f[2] * map2k +
                   params.kb[2] * map2kp);
    gsl_matrix_set(m, 6, 2, params.f[2] * map3kp);
    gsl_matrix_set(m, 6, 3, params.kb[2] * map3kp);
    gsl_matrix_set(m, 6, 4, 0.0);
    gsl_matrix_set(m, 6, 5, 0.0);
    gsl_matrix_set(m, 6, 6, -params.b[2] - params.kf[2]);
    gsl_matrix_set(m, 6, 7, 0.0);
    gsl_matrix_set(m, 6, 8, 0.0);
    gsl_matrix_set(m, 6, 9, 0.0);
    gsl_matrix_set(m, 6, 10, 0.0);       
    // map2kp_m2kp.
    gsl_matrix_set(m, 7, 0, 0.0);
    gsl_matrix_set(m, 7, 1, 0.0);
    gsl_matrix_set(m, 7, 2, params.kb[3] * m2kp);
    gsl_matrix_set(m, 7, 3, params.f[3] * m2kp);
    gsl_matrix_set(m, 7, 4, 0.0);
    gsl_matrix_set(m, 7, 5, 0.0);
    gsl_matrix_set(m, 7, 6, 0.0);
    gsl_matrix_set(m, 7, 7, -params.b[3] - params.kf[3]);
    gsl_matrix_set(m, 7, 8, 0.0);
    gsl_matrix_set(m, 7, 9, 0.0);
    gsl_matrix_set(m, 7, 10, params.f[3] * map2kp +
                   params.kb[3] * map2k);
    // ras.
    gsl_matrix_set(m, 8, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 8, 1, -params.kb[0] * ras);
    gsl_matrix_set(m, 8, 2, 0.0);
    gsl_matrix_set(m, 8, 3, 0.0);
    gsl_matrix_set(m, 8, 4, params.b[0] + params.kf[0]);
    gsl_matrix_set(m, 8, 5, 0.0);
    gsl_matrix_set(m, 8, 6, 0.0);
    gsl_matrix_set(m, 8, 7, 0.0);
    gsl_matrix_set(m, 8, 8, -params.f[0] * map3k -
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 8, 9, 0.0);
    gsl_matrix_set(m, 8, 10, 0.0);
    // m3kp.
    gsl_matrix_set(m, 9, 0, -params.kb[1] * m3kp);
    gsl_matrix_set(m, 9, 1, -params.f[1] * m3kp);
    gsl_matrix_set(m, 9, 2, 0.0);
    gsl_matrix_set(m, 9, 3, 0.0);
    gsl_matrix_set(m, 9, 4, 0.0);
    gsl_matrix_set(m, 9, 5, params.b[1] + params.kf[1]);
    gsl_matrix_set(m, 9, 6, 0.0);
    gsl_matrix_set(m, 9, 7, 0.0);
    gsl_matrix_set(m, 9, 8, 0.0);
    gsl_matrix_set(m, 9, 9, -params.f[1] * map3kp -
                   params.kb[1] * map3k);
    gsl_matrix_set(m, 9, 10, 0.0);
    // m2kp.
    gsl_matrix_set(m, 10, 0, 0.0);
    gsl_matrix_set(m, 10, 1, 0.0);
    gsl_matrix_set(m, 10, 2, -params.kb[3] * m2kp);
    gsl_matrix_set(m, 10, 3, -params.f[3] * m2kp);
    gsl_matrix_set(m, 10, 4, 0.0);
    gsl_matrix_set(m, 10, 5, 0.0);
    gsl_matrix_set(m, 10, 6, 0.0);
    gsl_matrix_set(m, 10, 7, params.b[3] + params.kf[3]);
    gsl_matrix_set(m, 10, 8, 0.0);
    gsl_matrix_set(m, 10, 9, 0.0);
    gsl_matrix_set(m, 10, 10, -params.f[3] * map2kp -
                   params.kb[3] * map2k);
            
    for (int i = 0; i != 11; ++i) {
        dfdt[i] = 0.0;
    }

    return GSL_SUCCESS;
}

// }}}


int ReactantConcentration::odeFunction2s2p(double       t,
                                           const double y[],
                                           double       f[],
                                           void*        para)
// {{{ Member-Function: odeFunction2s2p

//// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 2 phosphorylations.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double map2k2p       = y[4];    // MAPKK-phosphoryl-phosphoryl
    double map3k_ras     = y[5];    // intermedia product
    double map3kp_m3kp   = y[6];    // intermedia product
    double map2k_map3kp  = y[7];    // intermedia product
    double map2kp_m2kp   = y[8];    // intermedia product
    double map2kp_map3kp = y[9];    // intermedia product
    double map2k2p_m2kp  = y[10];   // intermedia product
    double ras           = y[11];   // MAPKKK kinase
    double m3kp          = y[12];   // MAPKKK-p phosphatase
    double m2kp          = y[13];   // MAPKK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map2k + map3kp <b2==f2> map2k_map3kp <kb2-=kf2> map2kp + map3kp
    // map2kp + m2kp <b3==f3> map2kp_m2kp <kb3-=kf3> map2k + m2kp
    // map2kp + map3kp <b4==f4> map2kp_map3kp <kb4-=kf4> map2k2p + map3kp
    // map2k2p + m2kp <b5==f5> map2k2p_m2kp <kb5-=kf5> map2kp + m2kp
    
    // map3k
    f[0] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // map3kp
    f[1] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras -
        params.f[2] * map2k * map3kp + params.b[2] * map2k_map3kp +
        params.kf[2] * map2k_map3kp - params.kb[2] * map2kp * map3kp -
        params.f[4] * map2kp * map3kp + params.b[4] * map2kp_map3kp +
        params.kf[4] * map2kp_map3kp - params.kb[4] * map2k2p * map3kp;
    // map2k
    f[2] = -params.f[2] * map2k * map3kp + params.b[2] * map2k_map3kp +
        params.kf[3] * map2kp_m2kp - params.kb[3] * map2k * m2kp;
    // map2kp
    f[3] = -params.f[3] * map2kp * m2kp + params.b[3] * map2kp_m2kp +
        params.kf[2] * map2k_map3kp - params.kb[2] * map2kp * map3kp -
        params.f[4] * map2kp * map3kp + params.b[4] * map2kp_map3kp +
        params.kf[5] * map2k2p_m2kp - params.kb[5] * map2kp * m2kp;
    // map2k2p
    f[4] = -params.f[5] * map2k2p * m2kp + params.b[5] * map2k2p_m2kp +
        params.kf[4] * map2kp_map3kp - params.kb[4] * map2k2p * map3kp;
    // map3k_ras
    f[5] = params.f[0] * map3k * ras - params.b[0] * map3k_ras -
        params.kf[0] * map3k_ras + params.kb[0] * map3kp * ras;
    // map3kp_m3kp
    f[6] = params.f[1] * map3kp * m3kp - params.b[1] * map3kp_m3kp -
        params.kf[1] * map3kp_m3kp + params.kb[1] * map3k * m3kp;
    // map2k_map3kp
    f[7] = params.f[2] * map2k * map3kp - params.b[2] * map2k_map3kp -
        params.kf[2] * map2k_map3kp + params.kb[2] * map2kp * map3kp;
    // map2kp_m2kp
    f[8] = params.f[3] * map2kp * m2kp - params.b[3] * map2kp_m2kp -
        params.kf[3] * map2kp_m2kp + params.kb[3] * map2k * m2kp;
    // map2kp_map3kp
    f[9] = params.f[4] * map2kp * map3kp - params.b[4] * map2kp_map3kp -
        params.kf[4] * map2kp_map3kp + params.kb[4] * map2k2p * map3kp;
    // map2k2p_m2kp
    f[10] = params.f[5] * map2k2p * m2kp - params.b[5] * map2k2p_m2kp -
        params.kf[5] * map2k2p_m2kp + params.kb[5] * map2kp * m2kp;
    // ras: dras = -dmap3k_ras
    f[11] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras;
    // m3kp: dm3kp = -dmap3kp_m3kp
    f[12] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // m2kp: dm2kp = -dmap2kp_m2kp - dmap2k2p_m2kp
    f[13] = -params.f[3] * map2kp * m2kp + params.b[3] * map2kp_m2kp +
        params.kf[3] * map2kp_m2kp - params.kb[3] * map2k * m2kp -
        params.f[5] * map2k2p * m2kp + params.b[5] * map2k2p_m2kp +
        params.kf[5] * map2k2p_m2kp - params.kb[5] * map2kp * m2kp;
    
    return GSL_SUCCESS;    // in gsl header.
}

// }}}


int ReactantConcentration::odeJacobian2s2p(double       t,
                                           const double y[],
                                           double       dfdy[],
                                           double       dfdt[],
                                           void*        para)
// {{{ Member-Function: odeJacobian2s2p

//// Needed by <gsl/gsl_odeiv2.h>, 2 stages, 2 phosphorylations.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    // initiate.
    Parameter params = *static_cast<Parameter *>(para);
    
    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double map2k2p       = y[4];    // MAPKK-phosphoryl-phosphoryl
    double map3k_ras     = y[5];    // intermedia product
    map3k_ras            = y[5];    // get rid of warning
    double map3kp_m3kp   = y[6];    // intermedia product
    map3kp_m3kp          = y[6];    // get rid of warning
    double map2k_map3kp  = y[7];    // intermedia product
    map2k_map3kp         = y[7];    // get rid of warning
    double map2kp_m2kp   = y[8];    // intermedia product
    map2kp_m2kp          = y[8];    // get rid of warning
    double map2kp_map3kp = y[9];    // intermedia product
    map2kp_map3kp        = y[9];    // get rid of warning
    double map2k2p_m2kp  = y[10];   // intermedia product
    map2k2p_m2kp         = y[10];   // get rid of warning
    double ras           = y[11];   // MAPKKK kinase
    double m3kp          = y[12];   // MAPKKK-p phosphatase
    double m2kp          = y[13];   // MAPKK-p phosphatase
    
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 14, 14);
    gsl_matrix *m = &dfdy_mat.matrix;

    // calculate matrix.
    // map3k.
    gsl_matrix_set(m, 0, 0, -params.f[0] * ras -
                   params.kb[1] * m3kp);
    gsl_matrix_set(m, 0, 1, 0.0);
    gsl_matrix_set(m, 0, 2, 0.0);    
    gsl_matrix_set(m, 0, 3, 0.0);
    gsl_matrix_set(m, 0, 4, 0.0);
    gsl_matrix_set(m, 0, 5, params.b[0]);
    gsl_matrix_set(m, 0, 6, params.kf[1]);
    gsl_matrix_set(m, 0, 7, 0.0);
    gsl_matrix_set(m, 0, 8, 0.0);
    gsl_matrix_set(m, 0, 9, 0.0);
    gsl_matrix_set(m, 0, 10, 0.0);
    gsl_matrix_set(m, 0, 11, -params.f[0] * map3k);
    gsl_matrix_set(m, 0, 12, -params.kb[1] * map3k);
    gsl_matrix_set(m, 0, 13, 0.0);
    // map3kp.
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, -params.f[1] * m3kp -
                   params.kb[0] * ras -
                   params.f[2] * map2k -
                   params.kb[2] * map2kp -
                   params.f[4] * map2kp -
                   params.kb[4] * map2k2p);
    gsl_matrix_set(m, 1, 2, -params.f[2] * map3kp);
    gsl_matrix_set(m, 1, 3, -params.kb[2] * map3kp -
                   params.f[4] * map3kp);
    gsl_matrix_set(m, 1, 4, -params.kb[4] * map3kp);
    gsl_matrix_set(m, 1, 5, params.kf[0]);
    gsl_matrix_set(m, 1, 6, params.b[1]);
    gsl_matrix_set(m, 1, 7, params.b[2] + params.kf[2]);
    gsl_matrix_set(m, 1, 8, 0.0);
    gsl_matrix_set(m, 1, 9, params.b[4] + params.kf[4]);
    gsl_matrix_set(m, 1, 10, 0.0);
    gsl_matrix_set(m, 1, 11, -params.kb[0] * map3kp);
    gsl_matrix_set(m, 1, 12, -params.f[1] * map3kp);
    gsl_matrix_set(m, 1, 13, 0.0);
    // map2k.
    gsl_matrix_set(m, 2, 0, 0.0);
    gsl_matrix_set(m, 2, 1, -params.f[2] * map2k);
    gsl_matrix_set(m, 2, 2, -params.f[2] * map3kp -
                   params.kb[3] * m2kp);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 2, 4, 0.0);
    gsl_matrix_set(m, 2, 5, 0.0);
    gsl_matrix_set(m, 2, 6, 0.0);
    gsl_matrix_set(m, 2, 7, params.b[2]);
    gsl_matrix_set(m, 2, 8, params.kf[3]);
    gsl_matrix_set(m, 2, 9, 0.0);
    gsl_matrix_set(m, 2, 10, 0.0);
    gsl_matrix_set(m, 2, 11, 0.0);
    gsl_matrix_set(m, 2, 12, 0.0);
    gsl_matrix_set(m, 2, 13, -params.kb[3] * map2k);
    // map2kp.
    gsl_matrix_set(m, 3, 0, 0.0);
    gsl_matrix_set(m, 3, 1, -params.kb[2] * map2kp);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -params.f[3] * m2kp -
                   params.kb[2] * map3kp -
                   params.f[4] * map3kp -
                   params.kb[5] * m2kp);
    gsl_matrix_set(m, 3, 4, 0.0);
    gsl_matrix_set(m, 3, 5, 0.0);
    gsl_matrix_set(m, 3, 6, 0.0);
    gsl_matrix_set(m, 3, 7, params.kf[2]);
    gsl_matrix_set(m, 3, 8, params.b[3]);
    gsl_matrix_set(m, 3, 9, params.b[4]);
    gsl_matrix_set(m, 3, 10, params.kf[5]);
    gsl_matrix_set(m, 3, 11, 0.0);
    gsl_matrix_set(m, 3, 12, 0.0);
    gsl_matrix_set(m, 3, 13, -params.f[3] * map2kp -
                   params.kb[5] * map2kp);
    // map2k2p
    gsl_matrix_set(m, 4, 0, 0.0);
    gsl_matrix_set(m, 4, 1, -params.kb[4] * map2k2p);
    gsl_matrix_set(m, 4, 2, 0.0);
    gsl_matrix_set(m, 4, 3, 0.0);
    gsl_matrix_set(m, 4, 4, -params.f[5] * m2kp -
                   params.kb[4] * map3kp);
    gsl_matrix_set(m, 4, 5, 0.0);
    gsl_matrix_set(m, 4, 6, 0.0);
    gsl_matrix_set(m, 4, 7, 0.0);
    gsl_matrix_set(m, 4, 8, 0.0);
    gsl_matrix_set(m, 4, 9, params.kf[4]);
    gsl_matrix_set(m, 4, 10, params.b[5]);
    gsl_matrix_set(m, 4, 11, 0.0);
    gsl_matrix_set(m, 4, 12, 0.0);
    gsl_matrix_set(m, 4, 13, -params.f[5] * map2k2p);
    // map3k_ras.
    gsl_matrix_set(m, 5, 0, params.f[0] * ras);
    gsl_matrix_set(m, 5, 1, params.kb[0] * ras);
    gsl_matrix_set(m, 5, 2, 0.0);
    gsl_matrix_set(m, 5, 3, 0.0);
    gsl_matrix_set(m, 5, 4, 0.0);
    gsl_matrix_set(m, 5, 5, -params.b[0] - params.kf[0]);
    gsl_matrix_set(m, 5, 6, 0.0);
    gsl_matrix_set(m, 5, 7, 0.0);
    gsl_matrix_set(m, 5, 8, 0.0);
    gsl_matrix_set(m, 5, 9, 0.0);
    gsl_matrix_set(m, 5, 10, 0.0);
    gsl_matrix_set(m, 5, 11, params.f[0] * map3k +
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 5, 12, 0.0);
    gsl_matrix_set(m, 5, 13, 0.0);
    // map3kp_m3kp.
    gsl_matrix_set(m, 6, 0, params.kb[1] * m3kp);
    gsl_matrix_set(m, 6, 1, params.f[1] * m3kp);
    gsl_matrix_set(m, 6, 2, 0.0);
    gsl_matrix_set(m, 6, 3, 0.0);
    gsl_matrix_set(m, 6, 4, 0.0);
    gsl_matrix_set(m, 6, 5, 0.0);
    gsl_matrix_set(m, 6, 6, -params.b[1] - params.kf[1]);
    gsl_matrix_set(m, 6, 7, 0.0);
    gsl_matrix_set(m, 6, 8, 0.0);
    gsl_matrix_set(m, 6, 9, 0.0);
    gsl_matrix_set(m, 6, 10, 0.0);
    gsl_matrix_set(m, 6, 11, 0.0);
    gsl_matrix_set(m, 6, 12, params.f[1] * map3kp +
                   params.kb[1] * map3k);
    gsl_matrix_set(m, 6, 13, 0.0);
    // map2k_map3kp.
    gsl_matrix_set(m, 7, 0, 0.0);
    gsl_matrix_set(m, 7, 1, params.f[2] * map2k +
                   params.kb[2] * map2kp);
    gsl_matrix_set(m, 7, 2, params.f[2] * map3kp);
    gsl_matrix_set(m, 7, 3, params.kb[2] * map3kp);
    gsl_matrix_set(m, 7, 4, 0.0);
    gsl_matrix_set(m, 7, 5, 0.0);
    gsl_matrix_set(m, 7, 6, 0.0);
    gsl_matrix_set(m, 7, 7, -params.b[2] - params.kf[2]);
    gsl_matrix_set(m, 7, 8, 0.0);
    gsl_matrix_set(m, 7, 9, 0.0);
    gsl_matrix_set(m, 7, 10, 0.0);
    gsl_matrix_set(m, 7, 11, 0.0);
    gsl_matrix_set(m, 7, 12, 0.0);
    gsl_matrix_set(m, 7, 13, 0.0);
    // map2kp_m2kp.
    gsl_matrix_set(m, 8, 0, 0.0);
    gsl_matrix_set(m, 8, 1, 0.0);
    gsl_matrix_set(m, 8, 2, params.kb[3] * m2kp);
    gsl_matrix_set(m, 8, 3, params.f[3] * m2kp);
    gsl_matrix_set(m, 8, 4, 0.0);
    gsl_matrix_set(m, 8, 5, 0.0);
    gsl_matrix_set(m, 8, 6, 0.0);
    gsl_matrix_set(m, 8, 7, 0.0);
    gsl_matrix_set(m, 8, 8, -params.b[3] - params.kf[3]);
    gsl_matrix_set(m, 8, 9, 0.0);
    gsl_matrix_set(m, 8, 10, 0.0);
    gsl_matrix_set(m, 8, 11, 0.0);
    gsl_matrix_set(m, 8, 12, 0.0);
    gsl_matrix_set(m, 8, 13, params.f[3] * map2kp +
                   params.kb[3] * map2k);
    // map2kp_map3kp
    gsl_matrix_set(m, 9, 0, 0.0);
    gsl_matrix_set(m, 9, 1, params.f[4] * map2kp +
                   params.kb[4] * map2k2p);
    gsl_matrix_set(m, 9, 2, 0.0);
    gsl_matrix_set(m, 9, 3, params.f[4] * map3kp);
    gsl_matrix_set(m, 9, 4, params.kb[4] * map3kp);
    gsl_matrix_set(m, 9, 5, 0.0);
    gsl_matrix_set(m, 9, 6, 0.0);
    gsl_matrix_set(m, 9, 7, 0.0);
    gsl_matrix_set(m, 9, 8, 0.0);
    gsl_matrix_set(m, 9, 9, -params.b[4] - params.kf[4]);
    gsl_matrix_set(m, 9, 10, 0.0);
    gsl_matrix_set(m, 9, 11, 0.0);
    gsl_matrix_set(m, 9, 12, 0.0);
    gsl_matrix_set(m, 9, 13, 0.0);
    // map2k2p_m2kp
    gsl_matrix_set(m, 10, 0, 0.0);
    gsl_matrix_set(m, 10, 1, 0.0);
    gsl_matrix_set(m, 10, 2, 0.0);
    gsl_matrix_set(m, 10, 3, params.kb[5] * m2kp);
    gsl_matrix_set(m, 10, 4, params.f[5] * m2kp);
    gsl_matrix_set(m, 10, 5, 0.0);
    gsl_matrix_set(m, 10, 6, 0.0);
    gsl_matrix_set(m, 10, 7, 0.0);
    gsl_matrix_set(m, 10, 8, 0.0);
    gsl_matrix_set(m, 10, 9, 0.0);
    gsl_matrix_set(m, 10, 10, -params.b[5] - params.kf[5]);
    gsl_matrix_set(m, 10, 11, 0.0);
    gsl_matrix_set(m, 10, 12, 0.0);
    gsl_matrix_set(m, 10, 13, params.f[5] * map2k2p +
                   params.kb[5] * map2kp);
    // ras.
    gsl_matrix_set(m, 11, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 11, 1, -params.kb[0] * ras);
    gsl_matrix_set(m, 11, 2, 0.0);
    gsl_matrix_set(m, 11, 3, 0.0);
    gsl_matrix_set(m, 11, 4, 0.0);
    gsl_matrix_set(m, 11, 5, params.b[0] + params.kf[0]);
    gsl_matrix_set(m, 11, 6, 0.0);
    gsl_matrix_set(m, 11, 7, 0.0);
    gsl_matrix_set(m, 11, 8, 0.0);
    gsl_matrix_set(m, 11, 9, 0.0);
    gsl_matrix_set(m, 11, 10, 0.0);
    gsl_matrix_set(m, 11, 11, -params.f[0] * map3k -
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 11, 12, 0.0);
    gsl_matrix_set(m, 11, 13, 0.0);
    // m3kp.
    gsl_matrix_set(m, 12, 0, -params.kb[1] * m3kp);
    gsl_matrix_set(m, 12, 1, -params.f[1] * m3kp);
    gsl_matrix_set(m, 12, 2, 0.0);
    gsl_matrix_set(m, 12, 3, 0.0);
    gsl_matrix_set(m, 12, 4, 0.0);
    gsl_matrix_set(m, 12, 5, 0.0);
    gsl_matrix_set(m, 12, 6, params.b[1] + params.kf[1]);
    gsl_matrix_set(m, 12, 7, 0.0);
    gsl_matrix_set(m, 12, 8, 0.0);
    gsl_matrix_set(m, 12, 9, 0.0);
    gsl_matrix_set(m, 12, 10, 0.0);
    gsl_matrix_set(m, 12, 11, 0.0);
    gsl_matrix_set(m, 12, 12, -params.f[1] * map3kp -
                   params.kb[1] * map3k);
    gsl_matrix_set(m, 12, 13, 0.0);
    // m2kp.
    gsl_matrix_set(m, 13, 0, 0.0);
    gsl_matrix_set(m, 13, 1, 0.0);
    gsl_matrix_set(m, 13, 2, -params.kb[3] * m2kp);
    gsl_matrix_set(m, 13, 3, -params.f[3] * m2kp -
                   params.kb[5] * m2kp);
    gsl_matrix_set(m, 13, 4, -params.f[5] * m2kp);
    gsl_matrix_set(m, 13, 5, 0.0);
    gsl_matrix_set(m, 13, 6, 0.0);
    gsl_matrix_set(m, 13, 7, 0.0);
    gsl_matrix_set(m, 13, 8, params.b[3]  + params.kf[3]);
    gsl_matrix_set(m, 13, 9, 0.0);
    gsl_matrix_set(m, 13, 10, params.b[5] + params.kb[5]);
    gsl_matrix_set(m, 13, 11, 0.0);
    gsl_matrix_set(m, 13, 12, 0.0);
    gsl_matrix_set(m, 13, 13, -params.f[3] * map2kp -
                   params.kb[3] * map2k -
                   params.f[5] * map2k2p -
                   params.kb[5] * map2kp);
            
    for (int i = 0; i != 14; ++i) {
        dfdt[i] = 0.0;
    }

    return GSL_SUCCESS;
}

// }}}


int ReactantConcentration::odeFunction3s1p(double       t,
                                           const double y[],
                                           double       f[],
                                           void*        para)
// {{{ Member-Function: odeFunction3s1p
//// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 1 phosphorylation.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double mapk          = y[4];    // MAPK
    double mapkp         = y[5];    // MAPK-phosphoryl
    double map3k_ras     = y[6];    // intermedia product
    double map3kp_m3kp   = y[7];    // intermedia product
    double map2k_map3kp  = y[8];    // intermedia product
    double map2kp_m2kp   = y[9];    // intermedia product
    double mapk_map2kp   = y[10];   // intermedia product
    double mapkp_mkp     = y[11];   // intermedia product
    double ras           = y[12];   // MAPKKK kinase
    double m3kp          = y[13];   // MAPKKK-p phosphatase
    double m2kp          = y[14];   // MAPKK-p phosphatase
    double mkp           = y[15];   // MAPK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map2k + map3kp <b2==f2> map2k_map3kp <kb2-=kf2> map2kp + map3kp
    // map2kp + m2kp <b3==f3> map2kp_m2kp <kb3-=kf3> map2k + m2kp
    // mapk + map2kp <b4==f4> mapk_map2kp <kb4-=kf4> mapkp + map2kp
    // mapkp + mkp <b5==f5> mapkp_mkp <kb5-=kf5> mapk + mkp

    
    // map3k
    f[0] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // map3kp
    f[1] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras -
        params.f[2] * map2k * map3kp + params.b[2] * map2k_map3kp +
        params.kf[2] * map2k_map3kp - params.kb[2] * map2kp * map3kp;
    // map2k
    f[2] = -params.f[2] * map2k * map3kp + params.b[2] * map2k_map3kp +
        params.kf[3] * map2kp_m2kp - params.kb[3] * map2k * m2kp;
    // map2kp
    f[3] = -params.f[3] * map2kp * m2kp + params.b[3] * map2kp_m2kp +
        params.kf[2] * map2k_map3kp - params.kb[2] * map2kp * map3kp -
        params.f[4] * mapk * map2kp + params.b[4] * mapk_map2kp +
        params.kf[4] * mapk_map2kp - params.kb[4] * mapkp * map2kp;
    // mapk
    f[4] = -params.f[4] * mapk * map2kp + params.b[4] * mapk_map2kp +
        params.kf[5] * mapkp_mkp - params.kb[5] * mapk * mkp;
    // mapkp
    f[5] = -params.f[5] * mapkp * mkp + params.b[5] * mapkp_mkp +
        params.kf[4] * mapk_map2kp - params.kb[4] * mapkp * map2kp;
    // map3k_ras
    f[6] = params.f[0] * map3k * ras - params.b[0] * map3k_ras -
        params.kf[0] * map3k_ras + params.kb[0] * map3kp * ras;
    // map3kp_m3kp
    f[7] = params.f[1] * map3kp * m3kp - params.b[1] * map3kp_m3kp -
        params.kf[1] * map3kp_m3kp + params.kb[1] * map3k * m3kp;
    // map2k_map3kp
    f[8] = params.f[2] * map2k * map3kp - params.b[2] * map2k_map3kp -
        params.kf[2] * map2k_map3kp + params.kb[2] * map2kp * map3kp;
    // map2kp_m2kp
    f[9] = params.f[3] * map2kp * m2kp - params.b[3] * map2kp_m2kp -
        params.kf[3] * map2kp_m2kp + params.kb[3] * map2k * m2kp;
    // mapk_map2kp
    f[10] = params.f[4] * mapk * map2kp - params.b[4] * mapk_map2kp -
        params.kf[4] * mapk_map2kp + params.kb[4] * mapkp * map2kp;
    // mapkp_mkp
    f[11] = params.f[5] * mapkp * mkp - params.b[5] * mapkp_mkp -
        params.kf[5] * mapkp_mkp + params.kb[5] * mapk * mkp;
    // ras: dras = -dmap3k_ras
    f[12] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras;
    // m3kp: dm3kp = -dmap3kp_m3kp
    f[13] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // m2kp: dm2kp = -dmap2kp_m2kp
    f[14] = -params.f[3] * map2kp * m2kp + params.b[3] * map2kp_m2kp +
        params.kf[3] * map2kp_m2kp - params.kb[3] * map2k * m2kp;
    // mkp: dmkp = -dmapkp_mkp
    f[15] = -params.f[5] * mapkp * mkp + params.b[5] * mapkp_mkp +
        params.kf[5] * mapkp_mkp - params.kb[5] * mapk * mkp;
    
    return GSL_SUCCESS;    // in gsl header.
}
// }}}


int ReactantConcentration::odeJacobian3s1p(double       t,
                                           const double y[],
                                           double       dfdy[],
                                           double       dfdt[],
                                           void*        para)
// {{{ Member-Function: odeJacobian3s1p
//// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 1 phosphorylation.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    // initiate.
    Parameter params = *static_cast<Parameter *>(para);

    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double mapk          = y[4];    // MAPK
    double mapkp         = y[5];    // MAPK-phosphoryl
    double map3k_ras     = y[6];    // intermedia product
    map3k_ras            = y[6];    // get rid of warning
    double map3kp_m3kp   = y[7];    // intermedia product
    map3kp_m3kp          = y[7];    // get rid of warning
    double map2k_map3kp  = y[8];    // intermedia product
    map2k_map3kp         = y[8];    // get rid of warning
    double map2kp_m2kp   = y[9];    // intermedia product
    map2kp_m2kp          = y[9];    // get rid of warning
    double mapk_map2kp   = y[10];   // intermedia product
    mapk_map2kp          = y[10];   // get rid of warning
    double mapkp_mkp     = y[11];   // intermedia product
    mapkp_mkp            = y[11];   // get rid of warning
    double ras           = y[12];   // MAPKKK kinase
    double m3kp          = y[13];   // MAPKKK-p phosphatase
    double m2kp          = y[14];   // MAPKK-p phosphatase
    double mkp           = y[15];   // MAPK-p phosphatase

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 16, 16);
    gsl_matrix *m = &dfdy_mat.matrix;

    // calculate matrix.
    // map3k.
    gsl_matrix_set(m, 0, 0, -params.f[0] * ras -
                   params.kb[1] * m3kp);
    gsl_matrix_set(m, 0, 1, 0.0);
    gsl_matrix_set(m, 0, 2, 0.0);    
    gsl_matrix_set(m, 0, 3, 0.0);
    gsl_matrix_set(m, 0, 4, 0.0);
    gsl_matrix_set(m, 0, 5, 0.0);
    gsl_matrix_set(m, 0, 6, params.b[0]);
    gsl_matrix_set(m, 0, 7, params.kf[1]);
    gsl_matrix_set(m, 0, 8, 0.0);
    gsl_matrix_set(m, 0, 9, 0.0);
    gsl_matrix_set(m, 0, 10, 0.0);
    gsl_matrix_set(m, 0, 11, 0.0);
    gsl_matrix_set(m, 0, 12, -params.f[0] * map3k);
    gsl_matrix_set(m, 0, 13, -params.kb[1] * map3k);
    gsl_matrix_set(m, 0, 14, 0.0);
    gsl_matrix_set(m, 0, 15, 0.0);
    // map3kp.
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, -params.f[1] * m3kp -
                   params.kb[0] * ras -
                   params.f[2] * map2k -
                   params.kb[2] * map2kp);
    gsl_matrix_set(m, 1, 2, -params.f[2] * map3kp);
    gsl_matrix_set(m, 1, 3, -params.kb[2] * map3kp);
    gsl_matrix_set(m, 1, 4, 0.0);
    gsl_matrix_set(m, 1, 5, 0.0);
    gsl_matrix_set(m, 1, 6, params.kf[0]);
    gsl_matrix_set(m, 1, 7, params.b[1]);
    gsl_matrix_set(m, 1, 8, params.b[2] + params.kf[2]);
    gsl_matrix_set(m, 1, 9, 0.0);
    gsl_matrix_set(m, 1, 10, 0.0);
    gsl_matrix_set(m, 1, 11, 0.0);
    gsl_matrix_set(m, 1, 12, -params.kb[0] * map3kp);
    gsl_matrix_set(m, 1, 13, -params.f[1] * map3kp);
    gsl_matrix_set(m, 1, 14, 0.0);
    gsl_matrix_set(m, 1, 15, 0.0);
    // map2k.
    gsl_matrix_set(m, 2, 0, 0.0);
    gsl_matrix_set(m, 2, 1, -params.f[2] * map2k);
    gsl_matrix_set(m, 2, 2, -params.f[2] * map3kp -
                   params.kb[3] * m2kp);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 2, 4, 0.0);
    gsl_matrix_set(m, 2, 5, 0.0);
    gsl_matrix_set(m, 2, 6, 0.0);
    gsl_matrix_set(m, 2, 7, 0.0);
    gsl_matrix_set(m, 2, 8, params.b[2]);
    gsl_matrix_set(m, 2, 9, params.kf[3]);
    gsl_matrix_set(m, 2, 10, 0.0);
    gsl_matrix_set(m, 2, 11, 0.0);
    gsl_matrix_set(m, 2, 12, 0.0);
    gsl_matrix_set(m, 2, 13, 0.0);
    gsl_matrix_set(m, 2, 14, -params.kb[3] * map2k);
    gsl_matrix_set(m, 2, 15, 0.0);
    // map2kp.
    gsl_matrix_set(m, 3, 0, 0.0);
    gsl_matrix_set(m, 3, 1, -params.kb[2] * map2kp);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -params.f[3] * m2kp -
                   params.kb[2] * map3kp -
                   params.f[4] * mapk -
                   params.kb[4] * mapkp);
    gsl_matrix_set(m, 3, 4, -params.f[4] * map2kp);
    gsl_matrix_set(m, 3, 5, -params.kb[4] * map2kp);
    gsl_matrix_set(m, 3, 6, 0.0);
    gsl_matrix_set(m, 3, 7, 0.0);
    gsl_matrix_set(m, 3, 8, params.kf[2]);
    gsl_matrix_set(m, 3, 9, params.b[3]);
    gsl_matrix_set(m, 3, 10, params.b[4] + params.kf[4]);
    gsl_matrix_set(m, 3, 11, 0.0);
    gsl_matrix_set(m, 3, 12, 0.0);
    gsl_matrix_set(m, 3, 13, 0.0);
    gsl_matrix_set(m, 3, 14, -params.f[3] * map2kp);
    gsl_matrix_set(m, 3, 15, 0.0);
    // mapk
    gsl_matrix_set(m, 4, 0, 0.0);
    gsl_matrix_set(m, 4, 1, 0.0);
    gsl_matrix_set(m, 4, 2, 0.0);
    gsl_matrix_set(m, 4, 3, -params.f[4] * mapk);
    gsl_matrix_set(m, 4, 4, -params.f[4] * map2kp -
                   params.kb[5] * mkp);
    gsl_matrix_set(m, 4, 5, 0.0);
    gsl_matrix_set(m, 4, 6, 0.0);
    gsl_matrix_set(m, 4, 7, 0.0);
    gsl_matrix_set(m, 4, 8, 0.0);
    gsl_matrix_set(m, 4, 9, 0.0);
    gsl_matrix_set(m, 4, 10, params.b[4]);
    gsl_matrix_set(m, 4, 11, params.kf[5]);
    gsl_matrix_set(m, 4, 12, 0.0);
    gsl_matrix_set(m, 4, 13, 0.0);
    gsl_matrix_set(m, 4, 14, 0.0);
    gsl_matrix_set(m, 4, 15, -params.kb[5] * mapk);
    // mapkp
    gsl_matrix_set(m, 5, 0, 0.0);
    gsl_matrix_set(m, 5, 1, 0.0);
    gsl_matrix_set(m, 5, 2, 0.0);
    gsl_matrix_set(m, 5, 3, -params.kb[4] * mapkp);
    gsl_matrix_set(m, 5, 4, 0.0);
    gsl_matrix_set(m, 5, 5, -params.f[5] * mkp -
                   params.kb[4] * map2kp);
    gsl_matrix_set(m, 5, 6, 0.0);
    gsl_matrix_set(m, 5, 7, 0.0);
    gsl_matrix_set(m, 5, 8, 0.0);
    gsl_matrix_set(m, 5, 9, 0.0);
    gsl_matrix_set(m, 5, 10, params.kf[4]);
    gsl_matrix_set(m, 5, 11, params.b[5]);
    gsl_matrix_set(m, 5, 12, 0.0);
    gsl_matrix_set(m, 5, 13, 0.0);
    gsl_matrix_set(m, 5, 14, 0.0);
    gsl_matrix_set(m, 5, 15, -params.f[5] * mapkp);
    // map3k_ras.
    gsl_matrix_set(m, 6, 0, params.f[0] * ras);
    gsl_matrix_set(m, 6, 1, params.kb[0] * ras);
    gsl_matrix_set(m, 6, 2, 0.0);
    gsl_matrix_set(m, 6, 3, 0.0);
    gsl_matrix_set(m, 6, 4, 0.0);
    gsl_matrix_set(m, 6, 5, 0.0);
    gsl_matrix_set(m, 6, 6, -params.b[0] - params.kf[0] );
    gsl_matrix_set(m, 6, 7, 0.0);
    gsl_matrix_set(m, 6, 8, 0.0);
    gsl_matrix_set(m, 6, 9, 0.0);
    gsl_matrix_set(m, 6, 10, 0.0);
    gsl_matrix_set(m, 6, 11, 0.0);
    gsl_matrix_set(m, 6, 12, params.f[0] * map3k +
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 6, 13, 0.0);
    gsl_matrix_set(m, 6, 14, 0.0);
    gsl_matrix_set(m, 6, 15, 0.0);
    // map3kp_m3kp.
    gsl_matrix_set(m, 7, 0, params.kb[1] * m3kp);
    gsl_matrix_set(m, 7, 1, params.f[1] * m3kp);
    gsl_matrix_set(m, 7, 2, 0.0);
    gsl_matrix_set(m, 7, 3, 0.0);
    gsl_matrix_set(m, 7, 4, 0.0);
    gsl_matrix_set(m, 7, 5, 0.0);
    gsl_matrix_set(m, 7, 6, 0.0);
    gsl_matrix_set(m, 7, 7, -params.b[1] - params.kf[1]);
    gsl_matrix_set(m, 7, 8, 0.0);
    gsl_matrix_set(m, 7, 9, 0.0);
    gsl_matrix_set(m, 7, 10, 0.0);
    gsl_matrix_set(m, 7, 11, 0.0);
    gsl_matrix_set(m, 7, 12, 0.0);
    gsl_matrix_set(m, 7, 13, params.f[1] * map3kp +
                   params.kb[1] * map3k);
    gsl_matrix_set(m, 7, 14, 0.0);
    gsl_matrix_set(m, 7, 15, 0.0);
    // map2k_map3kp.
    gsl_matrix_set(m, 8, 0, 0.0);
    gsl_matrix_set(m, 8, 1, params.f[2] * map2k +
                   params.kb[2] * map2kp);
    gsl_matrix_set(m, 8, 2, params.f[2] * map3kp);
    gsl_matrix_set(m, 8, 3, params.kb[2] * map3kp);
    gsl_matrix_set(m, 8, 4, 0.0);
    gsl_matrix_set(m, 8, 5, 0.0);
    gsl_matrix_set(m, 8, 6, 0.0);
    gsl_matrix_set(m, 8, 7, 0.0);
    gsl_matrix_set(m, 8, 8, -params.b[2] - params.kf[2]);
    gsl_matrix_set(m, 8, 9, 0.0);
    gsl_matrix_set(m, 8, 10, 0.0);
    gsl_matrix_set(m, 8, 11, 0.0);
    gsl_matrix_set(m, 8, 12, 0.0);
    gsl_matrix_set(m, 8, 13, 0.0);
    gsl_matrix_set(m, 8, 14, 0.0);
    gsl_matrix_set(m, 8, 15, 0.0);
    // map2kp_m2kp.
    gsl_matrix_set(m, 9, 0, 0.0);
    gsl_matrix_set(m, 9, 1, 0.0);
    gsl_matrix_set(m, 9, 2, params.kb[3] * m2kp);
    gsl_matrix_set(m, 9, 3, params.f[3] * m2kp);
    gsl_matrix_set(m, 9, 4, 0.0);
    gsl_matrix_set(m, 9, 5, 0.0);
    gsl_matrix_set(m, 9, 6, 0.0);
    gsl_matrix_set(m, 9, 7, 0.0);
    gsl_matrix_set(m, 9, 8, 0.0);
    gsl_matrix_set(m, 9, 9, -params.b[3] - params.kf[3]);
    gsl_matrix_set(m, 9, 10, 0.0);
    gsl_matrix_set(m, 9, 11, 0.0);
    gsl_matrix_set(m, 9, 12, 0.0);
    gsl_matrix_set(m, 9, 13, 0.0);
    gsl_matrix_set(m, 9, 14, params.f[3] * map2kp +
                   params.kb[3] * map2k);
    gsl_matrix_set(m, 9, 15, 0.0);
    // mapk_map2kp
    gsl_matrix_set(m, 10, 0, 0.0);
    gsl_matrix_set(m, 10, 1, 0.0);
    gsl_matrix_set(m, 10, 2, 0.0);
    gsl_matrix_set(m, 10, 3, params.f[4] * mapk +
                   params.kb[4] * mapkp);
    gsl_matrix_set(m, 10, 4, params.f[4] * map2kp);
    gsl_matrix_set(m, 10, 5, params.kb[4] * map2kp);
    gsl_matrix_set(m, 10, 6, 0.0);
    gsl_matrix_set(m, 10, 7, 0.0);
    gsl_matrix_set(m, 10, 8, 0.0);
    gsl_matrix_set(m, 10, 9, 0.0);
    gsl_matrix_set(m, 10, 10, -params.b[4] - params.kf[4]);
    gsl_matrix_set(m, 10, 11, 0.0);
    gsl_matrix_set(m, 10, 12, 0.0);
    gsl_matrix_set(m, 10, 13, 0.0);
    gsl_matrix_set(m, 10, 14, 0.0);
    gsl_matrix_set(m, 10, 15, 0.0);
    // mapkp_mkp
    gsl_matrix_set(m, 11, 0, 0.0);
    gsl_matrix_set(m, 11, 1, 0.0);
    gsl_matrix_set(m, 11, 2, 0.0);
    gsl_matrix_set(m, 11, 3, 0.0);
    gsl_matrix_set(m, 11, 4, params.kb[5] * mkp);
    gsl_matrix_set(m, 11, 5, params.f[5] * mkp);
    gsl_matrix_set(m, 11, 6, 0.0);
    gsl_matrix_set(m, 11, 7, 0.0);
    gsl_matrix_set(m, 11, 8, 0.0);
    gsl_matrix_set(m, 11, 9, 0.0);
    gsl_matrix_set(m, 11, 10, 0.0);
    gsl_matrix_set(m, 11, 11, -params.b[5] - params.kf[5]);
    gsl_matrix_set(m, 11, 12, 0.0);
    gsl_matrix_set(m, 11, 13, 0.0);
    gsl_matrix_set(m, 11, 14, 0.0);
    gsl_matrix_set(m, 11, 15, params.f[5] * mapkp +
                   params.kb[5] * mapk);
    // ras.
    gsl_matrix_set(m, 12, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 12, 1, -params.kb[0] * ras);
    gsl_matrix_set(m, 12, 2, 0.0);
    gsl_matrix_set(m, 12, 3, 0.0);
    gsl_matrix_set(m, 12, 4, 0.0);
    gsl_matrix_set(m, 12, 5, 0.0);
    gsl_matrix_set(m, 12, 6, params.b[0] + params.kf[0]);
    gsl_matrix_set(m, 12, 7, 0.0);
    gsl_matrix_set(m, 12, 8, 0.0);
    gsl_matrix_set(m, 12, 9, 0.0);
    gsl_matrix_set(m, 12, 10, 0.0);
    gsl_matrix_set(m, 12, 11, 0.0);
    gsl_matrix_set(m, 12, 12, -params.f[0] * map3k -
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 12, 13, 0.0);
    gsl_matrix_set(m, 12, 14, 0.0);
    gsl_matrix_set(m, 12, 15, 0.0);
    // m3kp.
    gsl_matrix_set(m, 13, 0, -params.kb[1] * m3kp);
    gsl_matrix_set(m, 13, 1, -params.f[1] * m3kp);
    gsl_matrix_set(m, 13, 2, 0.0);
    gsl_matrix_set(m, 13, 3, 0.0);
    gsl_matrix_set(m, 13, 4, 0.0);
    gsl_matrix_set(m, 13, 5, 0.0);
    gsl_matrix_set(m, 13, 6, 0.0);
    gsl_matrix_set(m, 13, 7, params.b[1] + params.kf[1]);
    gsl_matrix_set(m, 13, 8, 0.0);
    gsl_matrix_set(m, 13, 9, 0.0);
    gsl_matrix_set(m, 13, 10, 0.0);
    gsl_matrix_set(m, 13, 11, 0.0);
    gsl_matrix_set(m, 13, 12, 0.0);
    gsl_matrix_set(m, 13, 13, -params.f[1] * map3kp -
                   params.kb[1] * map3k);
    gsl_matrix_set(m, 13, 14, 0.0);
    gsl_matrix_set(m, 13, 15, 0.0);
    // m2kp.
    gsl_matrix_set(m, 14, 0, 0.0);
    gsl_matrix_set(m, 14, 1, 0.0);
    gsl_matrix_set(m, 14, 2, -params.kb[3] * m2kp);
    gsl_matrix_set(m, 14, 3, -params.f[3] * m2kp);
    gsl_matrix_set(m, 14, 4, 0.0);
    gsl_matrix_set(m, 14, 5, 0.0);
    gsl_matrix_set(m, 14, 6, 0.0);
    gsl_matrix_set(m, 14, 7, 0.0);
    gsl_matrix_set(m, 14, 8, 0.0);
    gsl_matrix_set(m, 14, 9, params.b[3] + params.kf[3]);
    gsl_matrix_set(m, 14, 10, 0.0);
    gsl_matrix_set(m, 14, 11, 0.0);
    gsl_matrix_set(m, 14, 12, 0.0);
    gsl_matrix_set(m, 14, 13, 0.0);
    gsl_matrix_set(m, 14, 14, -params.f[3] * map2kp -
                   params.kb[3] * map2k);
    gsl_matrix_set(m, 14, 15, 0.0);
    // mkp.
    gsl_matrix_set(m, 15, 0, 0.0);
    gsl_matrix_set(m, 15, 1, 0.0);
    gsl_matrix_set(m, 15, 2, 0.0);
    gsl_matrix_set(m, 15, 3, 0.0);
    gsl_matrix_set(m, 15, 4, -params.kb[5] * mkp);
    gsl_matrix_set(m, 15, 5, -params.f[5] * mkp);
    gsl_matrix_set(m, 15, 6, 0.0);
    gsl_matrix_set(m, 15, 7, 0.0);
    gsl_matrix_set(m, 15, 8, 0.0);
    gsl_matrix_set(m, 15, 9, 0.0);
    gsl_matrix_set(m, 15, 10, 0.0);
    gsl_matrix_set(m, 15, 11, params.b[5] + params.kf[5]);
    gsl_matrix_set(m, 15, 12, 0.0);
    gsl_matrix_set(m, 15, 13, 0.0);
    gsl_matrix_set(m, 15, 14, 0.0);
    gsl_matrix_set(m, 15, 15, -params.f[5] * mapkp -
                   params.kb[5] * mapk);

     
    for (int i = 0; i != 16; ++i) {
    dfdt[i] = 0.0;
    }

    return GSL_SUCCESS;
}
// }}}


int ReactantConcentration::odeFunction3s2p(double t,
                                           const double y[],
                                           double f[],
                                           void *para)
// {{{ Member-Function: odeFunction3s2p

//// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double map2k2p       = y[4];    // MAPKK-phosphoryl-phosphoryl
    double mapk          = y[5];    // MAPK
    double mapkp         = y[6];    // MAPK-phosphoryl
    double mapk2p        = y[7];    // MAPK-phosphoryl-phosphoryl
    double map3k_ras     = y[8];    // intermedia product
    map3k_ras            = y[8];    // get rid of warning
    double map3kp_m3kp   = y[9];    // intermedia product
    map3kp_m3kp          = y[9];    // get rid of warning
    double map2k_map3kp  = y[10];   // intermedia product
    map2k_map3kp         = y[10];   // get rid of warning
    double map2kp_m2kp   = y[11];   // intermedia product
    map2kp_m2kp          = y[11];   // get rid of warning
    double map2kp_map3kp = y[12];   // intermedia product
    map2kp_map3kp        = y[12];   // get rid of warning
    double map2k2p_m2kp  = y[13];   // intermedia product
    map2k2p_m2kp         = y[13];   // get rid of warning
    double mapk_map2k2p  = y[14];   // intermedia product
    mapk_map2k2p         = y[14];   // get rid of warning
    double mapkp_mkp     = y[15];   // intermedia product
    mapkp_mkp            = y[15];   // get rid of warning
    double mapkp_map2k2p = y[16];   // intermedia product
    mapkp_map2k2p        = y[16];   // get rid of warning
    double mapk2p_mkp    = y[17];   // intermedia product
    mapk2p_mkp           = y[17];   // get rid of warning
    double ras           = y[18];   // MAPKKK kinase
    double m3kp          = y[19];   // MAPKKK-p phosphatase
    double m2kp          = y[20];   // MAPKK-p phosphatase
    double mkp           = y[21];   // MAPK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map2k + map3kp <b2==f2> map2k_map3kp <kb2-=kf2> map2kp + map3kp
    // map2kp + m2kp <b3==f3> map2kp_m2kp <kb3-=kf3> map2k + m2kp
    // map2kp + map3kp <b4==f4> map2kp_map3kp <kb4-=kf4> map2k2p + map3kp
    // map2k2p + m2kp <b5==f5> map2k2p_m2kp <kb5-=kf5> map2kp + m2kp
    // mapk + map2k2p <b6==f6> mapk_map2k2p <kb6-=kf6> mapkp + map2k2p
    // mapkp + mkp <b7==f7> mapkp_mkp <kb7-=kf7> mapk + mkp
    // mapkp + map2k2p <b8==f8> mapkp_map2k2p <kb8-=kf8> mapk2p + map2k2p
    // mapk2p + mkp <b9==f9> mapkp_mkp <kb9-=kf9> mapkp + mkp

    
    // map3k
    f[0] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // map3kp
    f[1] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras -
        params.f[2] * map2k * map3kp + params.b[2] * map2k_map3kp +
        params.kf[2] * map2k_map3kp - params.kb[2] * map2kp * map3kp -
        params.f[4] * map2kp * map3kp + params.b[4] * map2kp_map3kp +
        params.kf[4] * map2kp_map3kp - params.kb[4] * map2k2p * map3kp;
    // map2k
    f[2] = -params.f[2] * map2k * map3kp + params.b[2] * map2k_map3kp +
        params.kf[3] * map2kp_m2kp - params.kb[3] * map2k * m2kp;
    // map2kp
    f[3] = -params.f[3] * map2kp * m2kp + params.b[3] * map2kp_m2kp +
        params.kf[2] * map2k_map3kp - params.kb[2] * map2kp * map3kp -
        params.f[4] * map2kp * map3kp + params.b[4] * map2kp_map3kp +
        params.kf[5] * map2k2p_m2kp - params.kb[5] * map2kp * m2kp;
    // map2k2p
    f[4] = -params.f[5] * map2k2p * m2kp + params.b[5] * map2k2p_m2kp +
        params.kf[4] * map2kp_map3kp - params.kb[4] * map2k2p * map3kp -
        params.f[6] * mapk * map2k2p + params.b[6] * mapk_map2k2p +
        params.kf[6] * mapk_map2k2p - params.kb[6] * mapkp * map2k2p -
        params.f[8] * mapkp * map2k2p + params.b[8] * mapkp_map2k2p +
        params.kf[8] * mapkp_map2k2p - params.kb[8] * mapk2p * map2k2p;
    // mapk
    f[5] = -params.f[6] * mapk * map2k2p + params.b[6] * mapk_map2k2p +
        params.kf[7] * mapkp_mkp - params.kb[7] * mapk * mkp;
    // mapkp
    f[6] = -params.f[7] * mapkp * mkp + params.b[7] * mapkp_mkp +
        params.kf[6] * mapk_map2k2p - params.kb[6] * mapkp * map2k2p -
        params.f[8] * mapkp * map2k2p + params.b[8] * mapkp_map2k2p +
        params.kf[9] * mapk2p_mkp - params.kb[9] * mapkp * mkp;
    // mapk2p
    f[7] = -params.f[9] * mapk2p * mkp + params.b[9] * mapk2p_mkp +
        params.kf[8] * mapkp_map2k2p - params.kb[8] * mapk2p * map2k2p;
    
    // map3k_ras
    f[8] = params.f[0] * map3k * ras - params.b[0] * map3k_ras -
        params.kf[0] * map3k_ras + params.kb[0] * map3kp * ras;
    // map3kp_m3kp
    f[9] = params.f[1] * map3kp * m3kp - params.b[1] * map3kp_m3kp -
        params.kf[1] * map3kp_m3kp + params.kb[1] * map3k * m3kp;
    // map2k_map3kp
    f[10] = params.f[2] * map2k * map3kp - params.b[2] * map2k_map3kp -
        params.kf[2] * map2k_map3kp + params.kb[2] * map2kp * map3kp;
    // map2kp_m2kp
    f[11] = params.f[3] * map2kp * m2kp - params.b[3] * map2kp_m2kp -
        params.kf[3] * map2kp_m2kp + params.kb[3] * map2k * m2kp;
    // map2kp_map3kp
    f[12] = params.f[4] * map2kp * map3kp - params.b[4] * map2kp_map3kp -
        params.kf[4] * map2kp_map3kp + params.kb[4] * map2k2p * map3kp;
    // map2k2p_m2kp
    f[13] = params.f[5] * map2k2p * m2kp - params.b[5] * map2k2p_m2kp -
        params.kf[5] * map2k2p_m2kp + params.kb[5] * map2kp * m2kp;
    // mapk_map2k2p
    f[14] = params.f[6] * mapk * map2k2p - params.b[6] * mapk_map2k2p -
        params.kf[6] * mapk_map2k2p + params.kb[6] * mapkp * map2k2p;
    // mapkp_mkp
    f[15] = params.f[7] * mapkp * mkp - params.b[7] * mapkp_mkp -
        params.kf[7] * mapkp_mkp + params.kb[7] * mapk * mkp;
    // mapkp_map2k2p
    f[16] = params.f[8] * mapkp * map2k2p - params.b[8] * mapkp_map2k2p -
        params.kf[8] * mapkp_map2k2p + params.kb[8] * mapk2p * map2k2p;
    // mapk2p_mkp
    f[17] = params.f[9] * mapk2p * mkp - params.b[9] * mapk2p_mkp -
        params.kf[9] * mapk2p_mkp + params.kb[9] * mapkp * mkp;
    // ras: dras = -dmap3k_ras
    f[18] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras +
        params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras;
    // m3kp: dm3kp = -dmap3kp_m3kp
    f[19] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp +
        params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // m2kp: dm2kp = -dmap2kp_m2kp - dmap2k2p_m2kp
    f[20] = -params.f[3] * map2kp * m2kp + params.b[3] * map2kp_m2kp +
        params.kf[3] * map2kp_m2kp - params.kb[3] * map2k * m2kp -
        params.f[5] * map2k2p * m2kp + params.b[5] * map2k2p_m2kp +
        params.kf[5] * map2k2p_m2kp - params.kb[5] * map2kp * m2kp;
    // mkp: dmkp = -dmapkp_mkp - dmapk2p_mkp
    f[21] = -params.f[7] * mapkp * mkp + params.b[7] * mapkp_mkp +
        params.kf[7] * mapkp_mkp - params.kb[7] * mapk * mkp -
        params.f[9] * mapk2p * mkp + params.b[9] * mapk2p_mkp +
        params.kf[9] * mapk2p_mkp - params.kb[9] * mapkp * mkp;
    
    return GSL_SUCCESS;    // in gsl header.
}

// }}}


int ReactantConcentration::odeJacobian3s2p(double t,
                                           const double y[],
                                           double *dfdy,
                                           double dfdt[],
                                           void *para)
// {{{ Member-Function: odeJacobian3s2p
//// Needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    // initiate.
    Parameter params = *static_cast<Parameter *>(para);

    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double map2k2p       = y[4];    // MAPKK-phosphoryl-phosphoryl
    double mapk          = y[5];    // MAPK
    double mapkp         = y[6];    // MAPK-phosphoryl
    double mapk2p        = y[7];    // MAPK-phosphoryl-phosphoryl
    double map3k_ras     = y[8];    // intermedia product
    map3k_ras            = y[8];    // get rid of warning
    double map3kp_m3kp   = y[9];    // intermedia product
    map3kp_m3kp          = y[9];    // get rid of warning
    double map2k_map3kp  = y[10];   // intermedia product
    map2k_map3kp         = y[10];   // get rid of warning
    double map2kp_m2kp   = y[11];   // intermedia product
    map2kp_m2kp          = y[11];   // get rid of warning
    double map2kp_map3kp = y[12];   // intermedia product
    map2kp_map3kp        = y[12];   // get rid of warning
    double map2k2p_m2kp  = y[13];   // intermedia product
    map2k2p_m2kp         = y[13];   // get rid of warning
    double mapk_map2k2p  = y[14];   // intermedia product
    mapk_map2k2p         = y[14];   // get rid of warning
    double mapkp_mkp     = y[15];   // intermedia product
    mapkp_mkp            = y[15];   // get rid of warning
    double mapkp_map2k2p = y[16];   // intermedia product
    mapkp_map2k2p        = y[16];   // get rid of warning
    double mapk2p_mkp    = y[17];   // intermedia product
    mapk2p_mkp           = y[17];   // get rid of warning
    double ras           = y[18];   // MAPKKK kinase
    double m3kp          = y[19];   // MAPKKK-p phosphatase
    double m2kp          = y[20];   // MAPKK-p phosphatase
    double mkp           = y[21];   // MAPK-p phosphatase

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 21, 21);
    gsl_matrix *m = &dfdy_mat.matrix;

    // calculate matrix.
    // map3k.
    gsl_matrix_set(m, 0, 0, -params.f[0] * ras -
                   params.kb[1] * m3kp);
    gsl_matrix_set(m, 0, 1, 0.0);
    gsl_matrix_set(m, 0, 2, 0.0);    
    gsl_matrix_set(m, 0, 3, 0.0);
    gsl_matrix_set(m, 0, 4, 0.0);
    gsl_matrix_set(m, 0, 5, 0.0);
    gsl_matrix_set(m, 0, 6, 0.0);
    gsl_matrix_set(m, 0, 7, 0.0);
    gsl_matrix_set(m, 0, 8, params.b[0]);
    gsl_matrix_set(m, 0, 9, params.kf[1]);
    gsl_matrix_set(m, 0, 10, 0.0);
    gsl_matrix_set(m, 0, 11, 0.0);
    gsl_matrix_set(m, 0, 12, 0.0);
    gsl_matrix_set(m, 0, 13, 0.0);
    gsl_matrix_set(m, 0, 14, 0.0);
    gsl_matrix_set(m, 0, 15, 0.0);
    gsl_matrix_set(m, 0, 16, 0.0);
    gsl_matrix_set(m, 0, 17, 0.0);
    gsl_matrix_set(m, 0, 18, -params.f[0] * map3k);
    gsl_matrix_set(m, 0, 19, -params.kb[1] * map3k);
    gsl_matrix_set(m, 0, 20, 0.0);
    gsl_matrix_set(m, 0, 21, 0.0);
    // map3kp.
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, -params.f[1] * m3kp -
                   params.kb[0] * ras -
                   params.f[2] * map2k -
                   params.kb[2] * map2kp -
                   params.f[4] * map2kp -
                   params.kb[4] * map2k2p);
    gsl_matrix_set(m, 1, 2, -params.f[2] * map3kp);
    gsl_matrix_set(m, 1, 3, -params.kb[2] * map3kp -
                   params.f[4] * map3kp);
    gsl_matrix_set(m, 1, 4, -params.kb[4]);
    gsl_matrix_set(m, 1, 5, 0.0);
    gsl_matrix_set(m, 1, 6, 0.0);
    gsl_matrix_set(m, 1, 7, 0.0);
    gsl_matrix_set(m, 1, 8, params.kf[0]);
    gsl_matrix_set(m, 1, 9, params.b[1]);
    gsl_matrix_set(m, 1, 10, params.b[2] + params.kf[2]);
    gsl_matrix_set(m, 1, 11, 0.0);
    gsl_matrix_set(m, 1, 12, params.b[4] + params.kf[4]);
    gsl_matrix_set(m, 1, 13, 0.0);
    gsl_matrix_set(m, 1, 14, 0.0);
    gsl_matrix_set(m, 1, 15, 0.0);
    gsl_matrix_set(m, 1, 16, 0.0);
    gsl_matrix_set(m, 1, 17, 0.0);
    gsl_matrix_set(m, 1, 18, -params.kb[0] * map3kp);
    gsl_matrix_set(m, 1, 19, -params.f[1] * map3kp);
    gsl_matrix_set(m, 1, 20, 0.0);
    gsl_matrix_set(m, 1, 21, 0.0);
    // map2k.
    gsl_matrix_set(m, 2, 0, 0.0);
    gsl_matrix_set(m, 2, 1, -params.f[2] * map2k);
    gsl_matrix_set(m, 2, 2, -params.f[2] * map3kp -
                   params.kb[3] * m2kp);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 2, 4, 0.0);
    gsl_matrix_set(m, 2, 5, 0.0);
    gsl_matrix_set(m, 2, 6, 0.0);
    gsl_matrix_set(m, 2, 7, 0.0);
    gsl_matrix_set(m, 2, 8, 0.0);
    gsl_matrix_set(m, 2, 9, 0.0);
    gsl_matrix_set(m, 2, 10, params.b[2]);
    gsl_matrix_set(m, 2, 11, params.kf[3]);
    gsl_matrix_set(m, 2, 12, 0.0);
    gsl_matrix_set(m, 2, 13, 0.0);
    gsl_matrix_set(m, 2, 14, 0.0);
    gsl_matrix_set(m, 2, 15, 0.0);
    gsl_matrix_set(m, 2, 16, 0.0);
    gsl_matrix_set(m, 2, 17, 0.0);
    gsl_matrix_set(m, 2, 18, 0.0);
    gsl_matrix_set(m, 2, 19, 0.0);
    gsl_matrix_set(m, 2, 20, -params.kb[3] * map2k );
    gsl_matrix_set(m, 2, 21, 0.0);
    // map2kp.
    gsl_matrix_set(m, 3, 0, 0.0);
    gsl_matrix_set(m, 3, 1, -params.kb[2] * map2kp -
                   params.f[4] * map2kp);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -params.f[3] * m2kp -
                   params.kb[2] * map3kp -
                   params.f[4] * map3kp -
                   params.kb[5] * m2kp);
    gsl_matrix_set(m, 3, 4, 0.0);
    gsl_matrix_set(m, 3, 5, 0.0);
    gsl_matrix_set(m, 3, 6, 0.0);
    gsl_matrix_set(m, 3, 7, 0.0);
    gsl_matrix_set(m, 3, 8, 0.0);
    gsl_matrix_set(m, 3, 9, 0.0);
    gsl_matrix_set(m, 3, 10, params.kf[2]);
    gsl_matrix_set(m, 3, 11, params.b[3]);
    gsl_matrix_set(m, 3, 12, params.b[4]);
    gsl_matrix_set(m, 3, 13, params.kf[5]);
    gsl_matrix_set(m, 3, 14, 0.0);
    gsl_matrix_set(m, 3, 15, 0.0);
    gsl_matrix_set(m, 3, 16, 0.0);
    gsl_matrix_set(m, 3, 17, 0.0);
    gsl_matrix_set(m, 3, 18, 0.0);
    gsl_matrix_set(m, 3, 19, 0.0);
    gsl_matrix_set(m, 3, 20, -params.f[3] * map2kp -
                   params.kb[5] * map2kp);
    gsl_matrix_set(m, 3, 21, 0.0);
    // map2k2p
    gsl_matrix_set(m, 4, 0, 0.0);
    gsl_matrix_set(m, 4, 1, -params.kb[4] * map2k2p);
    gsl_matrix_set(m, 4, 2, 0.0);
    gsl_matrix_set(m, 4, 3, 0.0);
    gsl_matrix_set(m, 4, 4, -params.f[5] * m2kp -
                   params.kb[4] * map3kp -
                   params.f[6] * mapk -
                   params.kb[6] * mapkp -
                   params.f[8] * mapkp -
                   params.kb[8] * mapk2p);
    gsl_matrix_set(m, 4, 5, -params.f[6] * map2k2p);
    gsl_matrix_set(m, 4, 6, -params.kb[6] * map2k2p -
                   params.f[8] * map2k2p);
    gsl_matrix_set(m, 4, 7, -params.kb[8] * map2k2p);
    gsl_matrix_set(m, 4, 8, 0.0);
    gsl_matrix_set(m, 4, 9, 0.0);
    gsl_matrix_set(m, 4, 10, 0.0);
    gsl_matrix_set(m, 4, 11, 0.0);
    gsl_matrix_set(m, 4, 12, params.kf[4]);
    gsl_matrix_set(m, 4, 13, params.b[5]);
    gsl_matrix_set(m, 4, 14, params.b[6] + params.kf[6]);
    gsl_matrix_set(m, 4, 15, 0.0);
    gsl_matrix_set(m, 4, 16, params.b[8] + params.kf[8]);
    gsl_matrix_set(m, 4, 17, 0.0);
    gsl_matrix_set(m, 4, 18, 0.0);
    gsl_matrix_set(m, 4, 19, 0.0);
    gsl_matrix_set(m, 4, 20, -params.f[5] * map2k2p);
    gsl_matrix_set(m, 4, 21, 0.0);
    // mapk
    gsl_matrix_set(m, 5, 0, 0.0);
    gsl_matrix_set(m, 5, 1, 0.0);
    gsl_matrix_set(m, 5, 2, 0.0);
    gsl_matrix_set(m, 5, 3, 0.0);
    gsl_matrix_set(m, 5, 4, -params.f[6] * mapk);
    gsl_matrix_set(m, 5, 5, -params.f[6] * map2k2p -
                   params.kb[7] * mkp);
    gsl_matrix_set(m, 5, 6, 0.0);
    gsl_matrix_set(m, 5, 7, 0.0);
    gsl_matrix_set(m, 5, 8, 0.0);
    gsl_matrix_set(m, 5, 9, 0.0);
    gsl_matrix_set(m, 5, 10, 0.0);
    gsl_matrix_set(m, 5, 11, 0.0);
    gsl_matrix_set(m, 5, 12, 0.0);
    gsl_matrix_set(m, 5, 13, 0.0);
    gsl_matrix_set(m, 5, 14, params.b[6]);
    gsl_matrix_set(m, 5, 15, params.kf[7]);
    gsl_matrix_set(m, 5, 16, 0.0);
    gsl_matrix_set(m, 5, 17, 0.0);
    gsl_matrix_set(m, 5, 18, 0.0);
    gsl_matrix_set(m, 5, 19, 0.0);
    gsl_matrix_set(m, 5, 20, 0.0);
    gsl_matrix_set(m, 5, 21, -params.kb[7] * mapk);
    // mapkp
    gsl_matrix_set(m, 6, 0, 0.0);
    gsl_matrix_set(m, 6, 1, 0.0);
    gsl_matrix_set(m, 6, 2, 0.0);
    gsl_matrix_set(m, 6, 3, 0.0);
    gsl_matrix_set(m, 6, 4, -params.kb[6] * mapkp -
                   params.f[8] * mapkp);
    gsl_matrix_set(m, 6, 5, 0.0);
    gsl_matrix_set(m, 6, 6, -params.f[7] * mkp -
                   params.kb[6] * map2k2p -
                   params.f[8] * map2k2p -
                   params.kb[9] * mkp);
    gsl_matrix_set(m, 6, 7, 0.0);
    gsl_matrix_set(m, 6, 8, 0.0);
    gsl_matrix_set(m, 6, 9, 0.0);
    gsl_matrix_set(m, 6, 10, 0.0);
    gsl_matrix_set(m, 6, 11, 0.0);
    gsl_matrix_set(m, 6, 12, 0.0);
    gsl_matrix_set(m, 6, 13, 0.0);
    gsl_matrix_set(m, 6, 14, params.kf[6]);
    gsl_matrix_set(m, 6, 15, params.b[7]);
    gsl_matrix_set(m, 6, 16, params.b[8]);
    gsl_matrix_set(m, 6, 17, params.kf[9]);
    gsl_matrix_set(m, 6, 18, 0.0);
    gsl_matrix_set(m, 6, 19, 0.0);
    gsl_matrix_set(m, 6, 20, 0.0);
    gsl_matrix_set(m, 6, 21, -params.f[7] * mapkp -
                   params.kb[9] * mapkp);
    // mapk2p
    gsl_matrix_set(m, 7, 0, 0.0);
    gsl_matrix_set(m, 7, 1, 0.0);
    gsl_matrix_set(m, 7, 2, 0.0);
    gsl_matrix_set(m, 7, 3, 0.0);
    gsl_matrix_set(m, 7, 4, -params.kb[8] * mapk2p);
    gsl_matrix_set(m, 7, 5, 0.0);
    gsl_matrix_set(m, 7, 6, 0.0);
    gsl_matrix_set(m, 7, 7, -params.f[9] * mkp -
                   params.kb[8] * map2k2p);
    gsl_matrix_set(m, 7, 8, 0.0);
    gsl_matrix_set(m, 7, 9, 0.0);
    gsl_matrix_set(m, 7, 10, 0.0);
    gsl_matrix_set(m, 7, 11, 0.0);
    gsl_matrix_set(m, 7, 12, 0.0);
    gsl_matrix_set(m, 7, 13, 0.0);
    gsl_matrix_set(m, 7, 14, 0.0);
    gsl_matrix_set(m, 7, 15, 0.0);
    gsl_matrix_set(m, 7, 16, params.kf[8]);
    gsl_matrix_set(m, 7, 17, params.b[9]);
    gsl_matrix_set(m, 7, 18, 0.0);
    gsl_matrix_set(m, 7, 19, 0.0);
    gsl_matrix_set(m, 7, 20, 0.0);
    gsl_matrix_set(m, 7, 21, -params.f[9] * mapk2p);    
    // map3k_ras.
    gsl_matrix_set(m, 8, 0, params.f[0] * ras);
    gsl_matrix_set(m, 8, 1, params.kb[0] * ras);
    gsl_matrix_set(m, 8, 2, 0.0);
    gsl_matrix_set(m, 8, 3, 0.0);
    gsl_matrix_set(m, 8, 4, 0.0);
    gsl_matrix_set(m, 8, 5, 0.0);
    gsl_matrix_set(m, 8, 6, 0.0);
    gsl_matrix_set(m, 8, 7, 0.0);
    gsl_matrix_set(m, 8, 8, -params.b[0] - params.kf[0]);
    gsl_matrix_set(m, 8, 9, 0.0);
    gsl_matrix_set(m, 8, 10, 0.0);
    gsl_matrix_set(m, 8, 11, 0.0);
    gsl_matrix_set(m, 8, 12, 0.0);
    gsl_matrix_set(m, 8, 13, 0.0);
    gsl_matrix_set(m, 8, 14, 0.0);
    gsl_matrix_set(m, 8, 15, 0.0);
    gsl_matrix_set(m, 8, 16, 0.0);
    gsl_matrix_set(m, 8, 17, 0.0);
    gsl_matrix_set(m, 8, 18, params.f[0] * map3k +
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 8, 19, 0.0);
    gsl_matrix_set(m, 8, 20, 0.0);
    gsl_matrix_set(m, 8, 21, 0.0);
    // map3kp_m3kp.
    gsl_matrix_set(m, 9, 0, params.kb[1] * m3kp);
    gsl_matrix_set(m, 9, 1, params.f[1] * m3kp);
    gsl_matrix_set(m, 9, 2, 0.0);
    gsl_matrix_set(m, 9, 3, 0.0);
    gsl_matrix_set(m, 9, 4, 0.0);
    gsl_matrix_set(m, 9, 5, 0.0);
    gsl_matrix_set(m, 9, 6, 0.0);
    gsl_matrix_set(m, 9, 7, 0.0);
    gsl_matrix_set(m, 9, 8, 0.0);
    gsl_matrix_set(m, 9, 9, -params.b[1] - params.kf[1]);
    gsl_matrix_set(m, 9, 10, 0.0);
    gsl_matrix_set(m, 9, 11, 0.0);
    gsl_matrix_set(m, 9, 12, 0.0);
    gsl_matrix_set(m, 9, 13, 0.0);
    gsl_matrix_set(m, 9, 14, 0.0);
    gsl_matrix_set(m, 9, 15, 0.0);
    gsl_matrix_set(m, 9, 16, 0.0);
    gsl_matrix_set(m, 9, 17, 0.0);
    gsl_matrix_set(m, 9, 18, 0.0);
    gsl_matrix_set(m, 9, 19, params.f[1] * map3kp +
                   params.kb[1] * map3k);
    gsl_matrix_set(m, 9, 20, 0.0);
    gsl_matrix_set(m, 9, 21, 0.0);
    // map2k_map3kp.
    gsl_matrix_set(m, 10, 0, 0.0);
    gsl_matrix_set(m, 10, 1, params.f[2] * map2k +
                   params.kb[2] * map2kp);
    gsl_matrix_set(m, 10, 2, params.f[2] * map3kp);
    gsl_matrix_set(m, 10, 3, params.kb[2] * map3kp);
    gsl_matrix_set(m, 10, 4, 0.0);
    gsl_matrix_set(m, 10, 5, 0.0);
    gsl_matrix_set(m, 10, 6, 0.0);
    gsl_matrix_set(m, 10, 7, 0.0);
    gsl_matrix_set(m, 10, 8, 0.0);
    gsl_matrix_set(m, 10, 9, 0.0);
    gsl_matrix_set(m, 10, 10, -params.b[2] - params.kf[2]);
    gsl_matrix_set(m, 10, 11, 0.0);
    gsl_matrix_set(m, 10, 12, 0.0);
    gsl_matrix_set(m, 10, 13, 0.0);
    gsl_matrix_set(m, 10, 14, 0.0);
    gsl_matrix_set(m, 10, 15, 0.0);
    gsl_matrix_set(m, 10, 16, 0.0);
    gsl_matrix_set(m, 10, 17, 0.0);
    gsl_matrix_set(m, 10, 18, 0.0);
    gsl_matrix_set(m, 10, 19, 0.0);
    gsl_matrix_set(m, 10, 20, 0.0);
    gsl_matrix_set(m, 10, 21, 0.0);
    // map2kp_m2kp.
    gsl_matrix_set(m, 11, 0, 0.0);
    gsl_matrix_set(m, 11, 1, 0.0);
    gsl_matrix_set(m, 11, 2, params.kb[3] * m2kp);
    gsl_matrix_set(m, 11, 3, params.f[3] * m2kp);
    gsl_matrix_set(m, 11, 4, 0.0);
    gsl_matrix_set(m, 11, 5, 0.0);
    gsl_matrix_set(m, 11, 6, 0.0);
    gsl_matrix_set(m, 11, 7, 0.0);
    gsl_matrix_set(m, 11, 8, 0.0);
    gsl_matrix_set(m, 11, 9, 0.0);
    gsl_matrix_set(m, 11, 10, 0.0);
    gsl_matrix_set(m, 11, 11, -params.b[3] - params.kf[3]);
    gsl_matrix_set(m, 11, 12, 0.0);
    gsl_matrix_set(m, 11, 13, 0.0);
    gsl_matrix_set(m, 11, 14, 0.0);
    gsl_matrix_set(m, 11, 15, 0.0);
    gsl_matrix_set(m, 11, 16, 0.0);
    gsl_matrix_set(m, 11, 17, 0.0);
    gsl_matrix_set(m, 11, 18, 0.0);
    gsl_matrix_set(m, 11, 19, 0.0);
    gsl_matrix_set(m, 11, 20, params.f[3] * map2kp +
                   params.kb[3] * map2k);
    gsl_matrix_set(m, 11, 21, 0.0);
    // map2kp_map3kp
    gsl_matrix_set(m, 12, 0, 0.0);
    gsl_matrix_set(m, 12, 1, params.f[4] * map2kp +
                   params.kb[4] * map2k2p);
    gsl_matrix_set(m, 12, 2, 0.0);
    gsl_matrix_set(m, 12, 3, params.f[4] * map3kp);
    gsl_matrix_set(m, 12, 4, params.kb[4] * map3kp);
    gsl_matrix_set(m, 12, 5, 0.0);
    gsl_matrix_set(m, 12, 6, 0.0);
    gsl_matrix_set(m, 12, 7, 0.0);
    gsl_matrix_set(m, 12, 8, 0.0);
    gsl_matrix_set(m, 12, 9, 0.0);
    gsl_matrix_set(m, 12, 10, 0.0);
    gsl_matrix_set(m, 12, 11, 0.0);
    gsl_matrix_set(m, 12, 12, -params.b[4] - params.kf[4]);
    gsl_matrix_set(m, 12, 13, 0.0);
    gsl_matrix_set(m, 12, 14, 0.0);
    gsl_matrix_set(m, 12, 15, 0.0);
    gsl_matrix_set(m, 12, 16, 0.0);
    gsl_matrix_set(m, 12, 17, 0.0);
    gsl_matrix_set(m, 12, 18, 0.0);
    gsl_matrix_set(m, 12, 19, 0.0);
    gsl_matrix_set(m, 12, 20, 0.0);
    gsl_matrix_set(m, 12, 21, 0.0);
    // map2k2p_m2kp
    gsl_matrix_set(m, 13, 0, 0.0);
    gsl_matrix_set(m, 13, 1, 0.0);
    gsl_matrix_set(m, 13, 2, 0.0);
    gsl_matrix_set(m, 13, 3, params.kb[5] * m2kp);
    gsl_matrix_set(m, 13, 4, params.f[5] * m2kp);
    gsl_matrix_set(m, 13, 5, 0.0);
    gsl_matrix_set(m, 13, 6, 0.0);
    gsl_matrix_set(m, 13, 7, 0.0);
    gsl_matrix_set(m, 13, 8, 0.0);
    gsl_matrix_set(m, 13, 9, 0.0);
    gsl_matrix_set(m, 13, 10, 0.0);
    gsl_matrix_set(m, 13, 11, 0.0);
    gsl_matrix_set(m, 13, 12, 0.0);
    gsl_matrix_set(m, 13, 13, -params.b[5] - params.kf[5]);
    gsl_matrix_set(m, 13, 14, 0.0);
    gsl_matrix_set(m, 13, 15, 0.0);
    gsl_matrix_set(m, 13, 16, 0.0);
    gsl_matrix_set(m, 13, 17, 0.0);
    gsl_matrix_set(m, 13, 18, 0.0);
    gsl_matrix_set(m, 13, 19, 0.0);
    gsl_matrix_set(m, 13, 20, params.f[5] * map2k2p +
                   params.kb[5] * map2kp);
    gsl_matrix_set(m, 13, 21, 0.0);
    // mapk_map2k2p
    gsl_matrix_set(m, 14, 0, 0.0);
    gsl_matrix_set(m, 14, 1, 0.0);
    gsl_matrix_set(m, 14, 2, 0.0);
    gsl_matrix_set(m, 14, 3, 0.0);
    gsl_matrix_set(m, 14, 4, params.f[6] * mapk +
                  params.kb[6] * mapkp);
    gsl_matrix_set(m, 14, 5, params.f[6] * map2k2p);
    gsl_matrix_set(m, 14, 6, params.kb[6] * map2k2p);
    gsl_matrix_set(m, 14, 7, 0.0);
    gsl_matrix_set(m, 14, 8, 0.0);
    gsl_matrix_set(m, 14, 9, 0.0);
    gsl_matrix_set(m, 14, 10, 0.0);
    gsl_matrix_set(m, 14, 11, 0.0);
    gsl_matrix_set(m, 14, 12, 0.0);
    gsl_matrix_set(m, 14, 13, 0.0);
    gsl_matrix_set(m, 14, 14, -params.b[6] - params.kf[6]);
    gsl_matrix_set(m, 14, 15, 0.0);
    gsl_matrix_set(m, 14, 16, 0.0);
    gsl_matrix_set(m, 14, 17, 0.0);
    gsl_matrix_set(m, 14, 18, 0.0);
    gsl_matrix_set(m, 14, 19, 0.0);
    gsl_matrix_set(m, 14, 20, 0.0);
    gsl_matrix_set(m, 14, 21, 0.0);
    // mapkp_mkp
    gsl_matrix_set(m, 15, 0, 0.0);
    gsl_matrix_set(m, 15, 1, 0.0);
    gsl_matrix_set(m, 15, 2, 0.0);
    gsl_matrix_set(m, 15, 3, 0.0);
    gsl_matrix_set(m, 15, 4, 0.0);
    gsl_matrix_set(m, 15, 5, params.kb[7] * mkp);
    gsl_matrix_set(m, 15, 6, params.f[7] * mkp);
    gsl_matrix_set(m, 15, 7, 0.0);
    gsl_matrix_set(m, 15, 8, 0.0);
    gsl_matrix_set(m, 15, 9, 0.0);
    gsl_matrix_set(m, 15, 10, 0.0);
    gsl_matrix_set(m, 15, 11, 0.0);
    gsl_matrix_set(m, 15, 12, 0.0);
    gsl_matrix_set(m, 15, 13, 0.0);
    gsl_matrix_set(m, 15, 14, 0.0);
    gsl_matrix_set(m, 15, 15, -params.b[7] - params.kf[7]);
    gsl_matrix_set(m, 15, 16, 0.0);
    gsl_matrix_set(m, 15, 17, 0.0);
    gsl_matrix_set(m, 15, 18, 0.0);
    gsl_matrix_set(m, 15, 19, 0.0);
    gsl_matrix_set(m, 15, 20, 0.0);
    gsl_matrix_set(m, 15, 21, params.f[7] * mapkp +
                   params.kb[7] * mapk);
    // mapkp_map2k2p
    gsl_matrix_set(m, 16, 0, 0.0);
    gsl_matrix_set(m, 16, 1, 0.0);
    gsl_matrix_set(m, 16, 2, 0.0);
    gsl_matrix_set(m, 16, 3, 0.0);
    gsl_matrix_set(m, 16, 4, params.f[8] * mapkp +
                   params.kb[8] * mapk2p);
    gsl_matrix_set(m, 16, 5, 0.0);
    gsl_matrix_set(m, 16, 6, params.f[8] * map2k2p);
    gsl_matrix_set(m, 16, 7, params.kb[8] * map2k2p);
    gsl_matrix_set(m, 16, 8, 0.0);
    gsl_matrix_set(m, 16, 9, 0.0);
    gsl_matrix_set(m, 16, 10, 0.0);
    gsl_matrix_set(m, 16, 11, 0.0);
    gsl_matrix_set(m, 16, 12, 0.0);
    gsl_matrix_set(m, 16, 13, 0.0);
    gsl_matrix_set(m, 16, 14, 0.0);
    gsl_matrix_set(m, 16, 15, 0.0);
    gsl_matrix_set(m, 16, 16, -params.b[8] - params.kf[8]);
    gsl_matrix_set(m, 16, 17, 0.0);
    gsl_matrix_set(m, 16, 18, 0.0);
    gsl_matrix_set(m, 16, 19, 0.0);
    gsl_matrix_set(m, 16, 20, 0.0);
    gsl_matrix_set(m, 16, 21, 0.0);
    // mapk2p_mkp
    gsl_matrix_set(m, 17, 0, 0.0);
    gsl_matrix_set(m, 17, 1, 0.0);
    gsl_matrix_set(m, 17, 2, 0.0);
    gsl_matrix_set(m, 17, 3, 0.0);
    gsl_matrix_set(m, 17, 4, 0.0);
    gsl_matrix_set(m, 17, 5, 0.0);
    gsl_matrix_set(m, 17, 6, params.kb[9] * mkp);
    gsl_matrix_set(m, 17, 7, params.f[9] * mkp);
    gsl_matrix_set(m, 17, 8, 0.0);
    gsl_matrix_set(m, 17, 9, 0.0);
    gsl_matrix_set(m, 17, 10, 0.0);
    gsl_matrix_set(m, 17, 11, 0.0);
    gsl_matrix_set(m, 17, 12, 0.0);
    gsl_matrix_set(m, 17, 13, 0.0);
    gsl_matrix_set(m, 17, 14, 0.0);
    gsl_matrix_set(m, 17, 15, 0.0);
    gsl_matrix_set(m, 17, 16, 0.0);
    gsl_matrix_set(m, 17, 17, -params.b[9] - params.kf[9]);
    gsl_matrix_set(m, 17, 18, 0.0);
    gsl_matrix_set(m, 17, 19, 0.0);
    gsl_matrix_set(m, 17, 20, 0.0);
    gsl_matrix_set(m, 17, 21, params.f[9] * mapk2p +
                   params.kb[9] * mapkp);
    // ras.
    gsl_matrix_set(m, 18, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 18, 1, -params.kb[0] * ras);
    gsl_matrix_set(m, 18, 2, 0.0);
    gsl_matrix_set(m, 18, 3, 0.0);
    gsl_matrix_set(m, 18, 4, 0.0);
    gsl_matrix_set(m, 18, 5, 0.0);
    gsl_matrix_set(m, 18, 6, 0.0);
    gsl_matrix_set(m, 18, 7, 0.0);
    gsl_matrix_set(m, 18, 8, params.b[0] + params.kf[0]);
    gsl_matrix_set(m, 18, 9, 0.0);
    gsl_matrix_set(m, 18, 10, 0.0);
    gsl_matrix_set(m, 18, 11, 0.0);
    gsl_matrix_set(m, 18, 12, 0.0);
    gsl_matrix_set(m, 18, 13, 0.0);
    gsl_matrix_set(m, 18, 14, 0.0);
    gsl_matrix_set(m, 18, 15, 0.0);
    gsl_matrix_set(m, 18, 16, 0.0);
    gsl_matrix_set(m, 18, 17, 0.0);
    gsl_matrix_set(m, 18, 18, -params.f[0] * map3k -
                   params.kb[0] * map3kp);
    gsl_matrix_set(m, 18, 19, 0.0);
    gsl_matrix_set(m, 18, 20, 0.0);
    gsl_matrix_set(m, 18, 21, 0.0);
    // m3kp.
    gsl_matrix_set(m, 19, 0, -params.kb[1] * m3kp);
    gsl_matrix_set(m, 19, 1, -params.f[1] * m3kp);
    gsl_matrix_set(m, 19, 2, 0.0);
    gsl_matrix_set(m, 19, 3, 0.0);
    gsl_matrix_set(m, 19, 4, 0.0);
    gsl_matrix_set(m, 19, 5, 0.0);
    gsl_matrix_set(m, 19, 6, 0.0);
    gsl_matrix_set(m, 19, 7, 0.0);
    gsl_matrix_set(m, 19, 8, 0.0);
    gsl_matrix_set(m, 19, 9, params.b[1] + params.kf[1]);
    gsl_matrix_set(m, 19, 10, 0.0);
    gsl_matrix_set(m, 19, 11, 0.0);
    gsl_matrix_set(m, 19, 12, 0.0);
    gsl_matrix_set(m, 19, 13, 0.0);
    gsl_matrix_set(m, 19, 14, 0.0);
    gsl_matrix_set(m, 19, 15, 0.0);
    gsl_matrix_set(m, 19, 16, 0.0);
    gsl_matrix_set(m, 19, 17, 0.0);
    gsl_matrix_set(m, 19, 18, 0.0);
    gsl_matrix_set(m, 19, 19, -params.f[1] * map3kp -
                   params.kb[1] * map3k);
    gsl_matrix_set(m, 19, 20, 0.0);
    gsl_matrix_set(m, 19, 21, 0.0);
    // m2kp.
    gsl_matrix_set(m, 20, 0, 0.0);
    gsl_matrix_set(m, 20, 1, 0.0);
    gsl_matrix_set(m, 20, 2, -params.kb[3] * m2kp);
    gsl_matrix_set(m, 20, 3, -params.f[3] * m2kp - params.kb[5] * m2kp);
    gsl_matrix_set(m, 20, 4, -params.f[5] * m2kp);
    gsl_matrix_set(m, 20, 5, 0.0);
    gsl_matrix_set(m, 20, 6, 0.0);
    gsl_matrix_set(m, 20, 7, 0.0);
    gsl_matrix_set(m, 20, 8, 0.0);
    gsl_matrix_set(m, 20, 9, 0.0);
    gsl_matrix_set(m, 20, 10, 0.0);
    gsl_matrix_set(m, 20, 11, params.b[3] + params.kf[3]);
    gsl_matrix_set(m, 20, 12, 0.0);
    gsl_matrix_set(m, 20, 13, params.b[5] + params.kf[5]);
    gsl_matrix_set(m, 20, 14, 0.0);
    gsl_matrix_set(m, 20, 15, 0.0);
    gsl_matrix_set(m, 20, 16, 0.0);
    gsl_matrix_set(m, 20, 17, 0.0);
    gsl_matrix_set(m, 20, 18, 0.0);
    gsl_matrix_set(m, 20, 19, 0.0);
    gsl_matrix_set(m, 20, 20, -params.f[3] * map2kp -
                   params.kb[3] * map2k -
                   params.f[5] * map2k2p -
                   params.kb[5] * map2kp);
    gsl_matrix_set(m, 20, 21, 0.0);
    // mkp.
    gsl_matrix_set(m, 21, 0, 0.0);
    gsl_matrix_set(m, 21, 1, 0.0);
    gsl_matrix_set(m, 21, 2, 0.0);
    gsl_matrix_set(m, 21, 3, 0.0);
    gsl_matrix_set(m, 21, 4, 0.0);
    gsl_matrix_set(m, 21, 5, -params.kb[7] * mkp);
    gsl_matrix_set(m, 21, 6, -params.f[7] * mkp -
                   params.kb[9] * mkp);
    gsl_matrix_set(m, 21, 7, -params.f[9] * mkp);
    gsl_matrix_set(m, 21, 8, 0.0);
    gsl_matrix_set(m, 21, 9, 0.0);
    gsl_matrix_set(m, 21, 10, 0.0);
    gsl_matrix_set(m, 21, 11, 0.0);
    gsl_matrix_set(m, 21, 12, 0.0);
    gsl_matrix_set(m, 21, 13, 0.0);
    gsl_matrix_set(m, 21, 14, 0.0);
    gsl_matrix_set(m, 21, 15, params.b[7] + params.kf[7]);
    gsl_matrix_set(m, 21, 16, 0.0);
    gsl_matrix_set(m, 21, 17, params.b[9] + params.kf[9]);
    gsl_matrix_set(m, 21, 18, 0.0);
    gsl_matrix_set(m, 21, 19, 0.0);
    gsl_matrix_set(m, 21, 20, 0.0);
    gsl_matrix_set(m, 21, 21, -params.f[7] * mapkp -
                   params.kb[7] * mapk -
                   params.f[9] * mapk2p -
                   params.kb[9] * mapkp);

     
    for (int i = 0; i != 22; ++i) {
    dfdt[i] = 0.0;
    }

    return GSL_SUCCESS;
}
// }}}


void ReactantConcentration::propensityFunction1s1p(const double y[],
                                                   double       prob[],
                                                   void*        para)
// {{{ Member-Function: propensity1s1p
//// needed by <gsl/gsl_odeiv2.h>, 3 stages, 2 phosphorylations.
{
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k       = y[0];    // MAPKKK
    double map3kp      = y[1];    // MAPKKK-phosphoryl
    double map3k_ras   = y[2];    // intermedia product
    double map3kp_m3kp = y[3];    // intermedia product
    double ras         = y[4];    // MAPKKK kinase
    double m3kp        = y[5];    // MAPKKK-p phosphatase 
    
    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    
    prob[0] = params.f[0] * map3k * ras;
    prob[1] = params.b[0] * map3k_ras;
    prob[2] = params.kf[0] * map3k_ras;
    prob[3] = params.kb[0] * map3kp * ras;
    prob[4] = params.f[1] * map3kp * m3kp;
    prob[5] = params.b[1] * map3kp_m3kp;
    prob[6] = params.kf[1] * map3kp_m3kp;
    prob[7] = params.kb[1] * map3k * m3kp;
    
}
// }}}


void ReactantConcentration::propensityFunction1s2p(const double y[],
                                                   double       prob[],
                                                   void*        para)
// {{{ Member-Function: propensity1s2p
//// Needed by calculateDissipation, 1 stage, 2 phosphorylations.
{
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k        = y[0];    // MAPKKK
    double map3kp       = y[1];    // MAPKKK-phosphoryl
    double map3k2p      = y[2];    // MAPKKK-phosphoryl-phosphoryl
    double map3k_ras    = y[3];    // intermedia product
    double map3kp_m3kp  = y[4];    // intermedia product
    double map3kp_ras   = y[5];    // intermedia product
    double map3k2p_m3kp = y[6];    // intermedia product
    double ras          = y[7];    // MAPKKK kinase
    double m3kp         = y[8];    // MAPKKK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map3kp + ras <b2==f2> map3kp_ras <kb2-=kf2> map3k2p + ras
    // map3k2p + m3kp <b3==f3> map3k2p_m3kp <kb3-=kf3> map3kp + m3kp
    
    prob[0] = params.f[0] * map3k * ras;
    prob[1] = params.b[0] * map3k_ras;
    prob[2] = params.kf[0] * map3k_ras;
    prob[3] = params.kb[0] * map3kp * ras;
    prob[4] = params.f[1] * map3kp * m3kp;
    prob[5] = params.b[1] * map3kp_m3kp;
    prob[6] = params.kf[1] * map3kp_m3kp;
    prob[7] = params.kb[1] * map3k * m3kp;
    prob[8] = params.f[2] * map3kp * ras;
    prob[9] = params.b[2] * map3kp_ras;
    prob[10] = params.kf[2] * map3kp_ras;
    prob[11] = params.kb[2] * map3k2p * ras;
    prob[12] = params.f[3] * map3k2p;
    prob[13] = params.b[3] * map3k2p_m3kp;
    prob[14] = params.kf[3] * map3k2p_m3kp;
    prob[15] = params.kb[3] * map3kp * m3kp;
    
}
// }}}


void ReactantConcentration::propensityFunction2s1p(const double y[],
                                                   double       prob[],
                                                   void*        para)
// {{{ Member-Function: propensity2s1p
//// Nneeded by calculateDissipation, 2 stages, 1 phosphorylation.
{
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k        = y[0];    // MAPKKK
    double map3kp       = y[1];    // MAPKKK-phosphoryl
    double map2k        = y[2];    // MAPKK
    double map2kp       = y[3];    // MAPKK-phosphoryl
    double map3k_ras    = y[4];    // intermedia product
    double map3kp_m3kp  = y[5];    // intermedia product
    double map2k_map3kp = y[6];    // intermedia product
    double map2kp_m2kp  = y[7];    // intermedia product
    double ras          = y[8];    // MAPKKK kinase
    double m3kp         = y[9];    // MAPKKK-p phosphatase
    double m2kp         = y[10];   // MAPKK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map2k + map3kp <b2==f2> map2k_map3kp <kb2-=kf2> map2kp + map3kp
    // map2kp + m2kp <b3==f3> map2kp_m2kp <kb3-=kf3> map2k + m2kp
    
    prob[0] = params.f[0] * map3k * ras;
    prob[1] = params.b[0] * map3k_ras;
    prob[2] = params.kf[0] * map3k_ras;
    prob[3] = params.kb[0] * map3kp * ras;
    prob[4] = params.f[1] * map3kp * m3kp;
    prob[5] = params.b[1] * map3kp_m3kp;
    prob[6] = params.kf[1] * map3kp_m3kp;
    prob[7] = params.kb[1] * map3k * m3kp;
    prob[8] = params.f[2] * map2k * map3kp;
    prob[9] = params.b[2] * map2k_map3kp;
    prob[10] = params.kf[2] * map2k_map3kp;
    prob[11] = params.kb[2] * map2kp * map3kp;
    prob[12] = params.f[3] * map2kp * m2kp;
    prob[13] = params.b[3] * map2kp_m2kp;
    prob[14] = params.kf[3] * map2kp_m2kp;
    prob[15] = params.kb[3] * map2k * m2kp;
    
}
// }}}


void ReactantConcentration::propensityFunction2s2p(const double y[],
                                                   double       prob[],
                                                   void*        para)
// {{{ Member-Function: propensity2s2p
//// Needed by calculateDissipation, 2 stages, 2 phosphorylations.
{
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    
    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double map2k2p       = y[4];    // MAPKK-phosphoryl-phosphoryl
    double map3k_ras     = y[5];    // intermedia product
    double map3kp_m3kp   = y[6];    // intermedia product
    double map2k_map3kp  = y[7];    // intermedia product
    double map2kp_m2kp   = y[8];    // intermedia product
    double map2kp_map3kp = y[9];    // intermedia product
    double map2k2p_m2kp  = y[10];   // intermedia product
    double ras           = y[11];   // MAPKKK kinase
    double m3kp          = y[12];   // MAPKKK-p phosphatase
    double m2kp          = y[13];   // MAPKK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map2k + map3kp <b2==f2> map2k_map3kp <kb2-=kf2> map2kp + map3kp
    // map2kp + m2kp <b3==f3> map2kp_m2kp <kb3-=kf3> map2k + m2kp
    // map2kp + map3kp <b4==f4> map2kp_map3kp <kb4-=kf4> map2k2p + map3kp
    // map2k2p + m2kp <b5==f5> map2k2p_m2kp <kb5-=kf5> map2kp + m2kp
    
    prob[0] = params.f[0] * map3k * ras;
    prob[1] = params.b[0] * map3k_ras;
    prob[2] = params.kf[0] * map3k_ras;
    prob[3] = params.kb[0] * map3kp * ras;
    prob[4] = params.f[1] * map3kp * m3kp;
    prob[5] = params.b[1] * map3kp_m3kp;
    prob[6] = params.kf[1] * map3kp_m3kp;
    prob[7] = params.kb[1] * map3k * m3kp;
    prob[8] = params.f[2] * map2k * map3kp;
    prob[9] = params.b[2] * map2k_map3kp;
    prob[10] = params.kf[2] * map2k_map3kp;
    prob[11] = params.kb[2] * map2kp * map3kp;
    prob[12] = params.f[3] * map2kp * m2kp;
    prob[13] = params.b[3] * map2kp_m2kp;
    prob[14] = params.kf[3] * map2kp_m2kp;
    prob[15] = params.kb[3] * map2k * m2kp;
    prob[16] = params.f[4] * map2kp * map3kp;
    prob[17] = params.b[4] * map2kp_map3kp;
    prob[18] = params.kf[4] * map2kp_map3kp;
    prob[19] = params.kb[4] * map2k2p * map3kp;
    prob[20] = params.f[5] * map2k2p * m2kp;
    prob[21] = params.b[5] * map2k2p_m2kp;
    prob[22] = params.kf[5] * map2k2p_m2kp;
    prob[23] = params.kb[5] * map2kp * m2kp;
    
}
// }}}


void ReactantConcentration::propensityFunction3s1p(const double y[],
                                                   double       prob[],
                                                   void*        para)
// {{{ Member-Function: propensity3s1p
//// Needed by calculateDissipation, 3 stages, 1 phosphorylation.
{
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double mapk          = y[4];    // MAPK
    double mapkp         = y[5];    // MAPK-phosphoryl
    double map3k_ras     = y[6];    // intermedia product
    double map3kp_m3kp   = y[7];    // intermedia product
    double map2k_map3kp  = y[8];    // intermedia product
    double map2kp_m2kp   = y[9];    // intermedia product
    double mapk_map2kp   = y[10];   // intermedia product
    double mapkp_mkp     = y[11];   // intermedia product
    double ras           = y[12];   // MAPKKK kinase
    double m3kp          = y[13];   // MAPKKK-p phosphatase
    double m2kp          = y[14];   // MAPKK-p phosphatase
    double mkp           = y[15];   // MAPK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map2k + map3kp <b2==f2> map2k_map3kp <kb2-=kf2> map2kp + map3kp
    // map2kp + m2kp <b3==f3> map2kp_m2kp <kb3-=kf3> map2k + m2kp
    // mapk + map2kp <b4==f4> mapk_map2kp <kb4-=kf4> mapkp + map2kp
    // mapkp + mkp <b5==f5> mapkp_mkp <kb5-=kf5> mapk + mkp
    
    prob[0] = params.f[0] * map3k * ras;
    prob[1] = params.b[0] * map3k_ras;
    prob[2] = params.kf[0] * map3k_ras;
    prob[3] = params.kb[0] * map3kp * ras;
    prob[4] = params.f[1] * map3kp * m3kp;
    prob[5] = params.b[1] * map3kp_m3kp;
    prob[6] = params.kf[1] * map3kp_m3kp;
    prob[7] = params.kb[1] * map3k * m3kp;
    prob[8] = params.f[2] * map2k * map3kp;
    prob[9] = params.b[2] * map2k_map3kp;
    prob[10] = params.kf[2] * map2k_map3kp;
    prob[11] = params.kb[2] * map2kp * map3kp;
    prob[12] = params.f[3] * map2kp * m2kp;
    prob[13] = params.b[3] * map2kp_m2kp;
    prob[14] = params.kf[3] * map2kp_m2kp;
    prob[15] = params.kb[3] * map2k * m2kp;
    prob[16] = params.f[4] * mapk * map2kp;
    prob[17] = params.b[4] * mapk_map2kp;
    prob[18] = params.kf[4] * mapk_map2kp;
    prob[19] = params.kb[4] * mapkp * map2kp;
    prob[20] = params.f[5] * mapkp * mkp;
    prob[21] = params.b[5] * mapkp_mkp;
    prob[22] = params.kf[5] * mapkp_mkp;
    prob[23] = params.kb[5] * mapk * mkp;
    
}
// }}}


void ReactantConcentration::propensityFunction3s2p(const double y[],
                                                   double       prob[],
                                                   void*        para)
// {{{ Member-Function: propensity3s2p

//// Needed by calculateDissipation, 3 stages, 2 phosphorylations.
{
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k         = y[0];    // MAPKKK
    double map3kp        = y[1];    // MAPKKK-phosphoryl
    double map2k         = y[2];    // MAPKK
    double map2kp        = y[3];    // MAPKK-phosphoryl
    double map2k2p       = y[4];    // MAPKK-phosphoryl-phosphoryl
    double mapk          = y[5];    // MAPK
    double mapkp         = y[6];    // MAPK-phosphoryl
    double mapk2p        = y[7];    // MAPK-phosphoryl-phosphoryl
    double map3k_ras     = y[8];    // intermedia product
    double map3kp_m3kp   = y[9];    // intermedia product
    double map2k_map3kp  = y[10];   // intermedia product
    double map2kp_m2kp   = y[11];   // intermedia product
    double map2kp_map3kp = y[12];   // intermedia product
    double map2k2p_m2kp  = y[13];   // intermedia product
    double mapk_map2k2p  = y[14];   // intermedia product
    double mapkp_mkp     = y[15];   // intermedia product
    double mapkp_map2k2p = y[16];   // intermedia product
    double mapk2p_mkp    = y[17];   // intermedia product
    double ras           = y[18];   // MAPKKK kinase
    double m3kp          = y[19];   // MAPKKK-p phosphatase
    double m2kp          = y[20];   // MAPKK-p phosphatase
    double mkp           = y[21];   // MAPK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map2k + map3kp <b2==f2> map2k_map3kp <kb2-=kf2> map2kp + map3kp
    // map2kp + m2kp <b3==f3> map2kp_m2kp <kb3-=kf3> map2k + m2kp
    // map2kp + map3kp <b4==f4> map2kp_map3kp <kb4-=kf4> map2k2p + map3kp
    // map2k2p + m2kp <b5==f5> map2k2p_m2kp <kb5-=kf5> map2kp + m2kp
    // mapk + map2k2p <b6==f6> mapk_map2k2p <kb6-=kf6> mapkp + map2k2p
    // mapkp + mkp <b7==f7> mapkp_mkp <kb7-=kf7> mapk + mkp
    // mapkp + map2k2p <b8==f8> mapkp_map2k2p <kb8-=kf8> mapk2p + map2k2p
    // mapk2p + mkp <b9==f9> mapkp_mkp <kb9-=kf9> mapkp + mkp

    prob[0] = params.f[0] * map3k * ras;
    prob[1] = params.b[0] * map3k_ras;
    prob[2] = params.kf[0] * map3k_ras;
    prob[3] = params.kb[0] * map3kp * ras;
    prob[4] = params.f[1] * map3kp * m3kp;
    prob[5] = params.b[1] * map3kp_m3kp;
    prob[6] = params.kf[1] * map3kp_m3kp;
    prob[7] = params.kb[1] * map3k * m3kp;
    prob[8] = params.f[2] * map2k * map3kp;
    prob[9] = params.b[2] * map2k_map3kp;
    prob[10] = params.kf[2] * map2k_map3kp; 
    prob[11] = params.kb[2] * map2kp * map3kp;
    prob[12] = params.f[3] * map2kp * m2kp;
    prob[13] = params.b[3] * map2kp_m2kp;
    prob[14] = params.kf[3] * map2kp_m2kp;
    prob[15] = params.kb[3] * map2k * m2kp;
    prob[16] = params.f[4] * map2kp * map3kp;
    prob[17] = params.b[4] * map2kp_map3kp;
    prob[18] = params.kf[4] * map2kp_map3kp;
    prob[19] = params.kb[4] * map2k2p * map3kp;
    prob[20] = params.f[5] * map2k2p * m2kp;
    prob[21] = params.b[5] * map2k2p_m2kp;
    prob[22] = params.kf[5] * map2k2p_m2kp;
    prob[23] = params.kb[5] * map2kp * m2kp;
    prob[24] = params.f[6] * mapk * map2k2p;
    prob[25] = params.b[6] * mapk_map2k2p;
    prob[26] = params.kf[6] * mapk_map2k2p;
    prob[27] = params.kb[6] * mapkp * map2k2p;
    prob[28] = params.f[7] * mapkp * mkp;
    prob[29] = params.b[7] * mapkp_mkp;
    prob[30] = params.kf[7] * mapkp_mkp;
    prob[31] = params.kb[7] * mapk * mkp;
    prob[32] = params.f[8] * mapkp * map2k2p;
    prob[33] = params.b[8] * mapkp_map2k2p;
    prob[34] = params.kf[8] * mapkp_map2k2p;
    prob[35] = params.kb[8] * mapk2p * map2k2p;
    prob[36] = params.f[9] * mapk2p * mkp;
    prob[37] = params.b[9] * mapk2p_mkp;
    prob[38] = params.kf[9] * mapk2p_mkp;
    prob[39] = params.kb[9] * mapkp * mkp;
}

// }}}
////////////////////

// }}}


////////////////////
// Class: Stable
////////////////////
// {{{ Class: Stable

void Stable::setValue(const double& i,
                      const double& o,
                      const double& t,
                      const double& d)
// {{{ Member-Function: setValue
//// Return value of according.
{
    input       = i;
    output      = o;
    time        = t;
    dissipation = d;
}
// }}}


double Stable::getValue(const int& according)
// {{{ Member-Function: getValue
//// Return value of according.
{
    switch(according) {
    case 1:
        return input;
    case 2:
        return output;
    case 3:
        return time;
    case 4:
        return dissipation;
    default:
        return input;
    }    
}
// }}}


bool Stable::lt(const Stable& a,
                const Stable& b,
                const int&    according)
// {{{ Member-Function: compare
//// Compare the value of Stable according to one of member.
{
    switch(according) {
    case 1:
        return (a.input < b.input);
    case 2:
        return (a.output < b.output);
    case 3:
        return (a.time < b.time);
    case 4:
        return (a.dissipation < b.dissipation);
    default:
        return false;
    }
}
// }}}


bool Stable::gt(const Stable& a,
                const Stable& b,
                const int&    according)
// {{{ Member-Function: compare
//// Compare the value of Stable according to one of member.
{
    switch(according) {
    case 1:
        return (a.input > b.input);
    case 2:
        return (a.output > b.output);
    case 3:
        return (a.time > b.time);
    case 4:
        return (a.dissipation > b.dissipation);
    default:
        return false;
    }
}
// }}}


bool Stable::eq(const Stable& a,
                const Stable& b,
                const int&    according)
// {{{ Member-Function: compare
//// Compare the value of Stable according to one of member.
{
    switch(according) {
    case 1:
        return (a.input == b.input);
    case 2:
        return (a.output == b.output);
    case 3:
        return (a.time == b.time);
    case 4:
        return (a.dissipation == b.dissipation);
    default:
        return false;
    }
}
// }}}
////////////////////

// }}}


////////////////////
// Class: StableList
////////////////////
// {{{ Class: StableList

StableList::StableList(const int& type)
// {{{ Constructor: StableList
////
{
    reaction_type = type;
    switch (reaction_type) {
    case 1:
        function_num = 2;
        reactant_num = 6;
        break;
    case 2:
        function_num = 4;
        reactant_num = 9;
        break;
    case 3:
        function_num = 4;
        reactant_num = 11;
        break;
    case 4:
        function_num = 6;
        reactant_num = 14;
        break;
    case 5:
        function_num = 6;
        reactant_num = 16;
        break;
    case 6:
        function_num = 10;
        reactant_num = 22;
        break;
    default:
        std::cout << "Wrong input of reaction type in RCStable!"
                  << std::endl;
    }
}
// }}}


Stable& StableList::operator[](const size_t index)
// {{{ Operator: []
//// Return the Stable value.
{
    if (index >= size()) {
        std::cout << "StableList []: out of range.\n";
    }
    return situation[index];
}
// }}}


void StableList::addBack(const double& iput,
                         const double& oput,
                         const double& t,
                         const double& e)
// {{{ Member-Function: addBack
//// Add values to the end.
{
    Stable tmp(iput, oput, t, e);
    situation.push_back(tmp);
}
// }}}


void StableList::addInter(std::vector<Stable>::iterator it,
                          const double&                 iput,
                          const double&                 oput,
                          const double&                 t,
                          const double&                 e)
// {{{ Member-Function: addInter
//// Insert values according to iterator.
{
    Stable tmp(iput, oput, t, e);
    situation.insert(it, tmp);
}
// }}}


void StableList::addInter(StableList::iterator   it,
                          std::vector<Stable>::value_type value)
// {{{ Member-Function: addInter
//// Insert values according to iterator.
{
    situation.insert(it, value);
}
// }}}


bool StableList::ltcompa(Stable& a, Stable& b)
{
    return (a.input < b.input);
}
bool StableList::ltcompb(Stable& a, Stable& b)
{
    return (a.output < b.output);
}
bool StableList::ltcompc(Stable& a, Stable& b)
{
    return (a.time < b.time);
}
bool StableList::ltcompd(Stable& a, Stable& b)
{
    return (a.dissipation < b.dissipation);
}


void StableList::sortList(int according)
// {{{ Member-Function: sort
//// Sort the vector.
{
    typedef bool (* CompFunc)(Stable& a, Stable& b);
    CompFunc mycomp = 0;
    switch(according) {
    case 1:
        mycomp = ltcompa;
        break;
    case 2:
        mycomp = ltcompb;
        break;
    case 3:
        mycomp = ltcompc;
        break;
    case 4:
        mycomp = ltcompd;
        break;
    }
    std::sort(situation.begin(),
              situation.end(),
              mycomp);
}
// }}}


std::vector<Stable>::iterator StableList::max(int according)
// {{{ Member-Function: maxOutput
//// Return the max output value.
{
    typedef bool (* CompFunc)(Stable& a, Stable& b);
    CompFunc mycomp = 0;
    switch(according) {
    case 1:
        mycomp = ltcompa;
        break;
    case 2:
        mycomp = ltcompb;
        break;
    case 3:
        mycomp = ltcompc;
        break;
    case 4:
        mycomp = ltcompd;
        break;
    }
    return std::max_element(situation.begin(),
                            situation.end(),
                            mycomp);
}
// }}}


std::vector<Stable>::iterator StableList::min(int according)
// {{{ Member-Function: maxOutput
//// Return the max output value.
{
    typedef bool (* CompFunc)(Stable& a, Stable& b);
    CompFunc mycomp = 0;
    switch(according) {
    case 1:
        mycomp = ltcompa;
        break;
    case 2:
        mycomp = ltcompb;
        break;
    case 3:
        mycomp = ltcompc;
        break;
    case 4:
        mycomp = ltcompd;
        break;
    }
    return std::min_element(situation.begin(),
                            situation.end(),
                            mycomp);
}
// }}}


std::vector<Stable>::iterator StableList::unSufficient()
// {{{ Member-Function: notSufficient
//// Return the first pointer if the gap between output is too large.
{
    if (situation.size() == 1) {
        return situation.end();
    } else if (situation.size() == 2) {
        std::vector<Stable>::iterator iter = situation.begin();
        ++iter;
        return iter;
    }
    
    //double xthreshold = 0.05 * (max(1)->input - min(1)->input);
    double ythreshold = 0.1 * (max(2)->output - min(2)->output);

    for (std::vector<Stable>::iterator ii = situation.begin() + 1;
         ii != situation.end();
         ++ii) {
        if (//(std::fabs(ii->output - (ii - 1)->output) > ythreshold)||
            (std::fabs(ii->output - (ii - 1)->output) > ythreshold)) {
            if (std::fabs((ii->input - (ii - 1)->input) / ii->input) > 0.01) {
                // This if is used to check
                // whether the inserted input value is equal.
                // If the inputs have equal values,
                // the cubic spline interpolation by gsl would fail.
                return ii;
            } else {
                continue;
            }
        }
    }
    return situation.end();
}
// }}}


bool StableList::isStable(int xaccording, int yaccording)
// {{{ Member-Function: isStable
//// Whether reach stable,
//// return true when stable, return false when not stable.
{
    std::vector<Stable>::size_type vct_len = situation.size();
    std::unique_ptr<double[]> x_array(new double[vct_len]);
    std::unique_ptr<double[]> y_array(new double[vct_len]);
    for (size_t si = 0; si != vct_len; ++si) {
        x_array[si] = situation[si].getValue(xaccording);
        y_array[si] = situation[si].getValue(yaccording);
    } // for.
    return isStableUseEndSome(x_array, y_array, vct_len, 0.02, 5);
}
// }}}


std::unique_ptr<double[]> StableList::getArray(const int& according)
// {{{ Member-Function: getArray
//// Return double[] of values using according to specific.
{
    std::unique_ptr<double[]> tmp(new double[situation.size()]);
    for (size_t i = 0; i != situation.size(); ++i) {
        tmp[i] = situation[i].getValue(according);
    }
    return tmp;
}
// }}}


size_t StableList::size()
// {{{ Member-Function: size
//// Return the size of situation.
{
    return situation.size();
}
// }}}


bool StableList::isMonotone()
// {{{ Member-Function: isMonotone
//// Return true if exists Monotone.
{
    double tmp       = situation.back().getValue(2) -
        situation.front().getValue(2) ;
    //double maxoutput = (*max(2)).getValue(2);
    int    sign      = tmp/std::fabs(tmp);
    for (std::vector<Stable>::size_type si = 1;
         si != situation.size();
         ++si) {
        if (sign == 1 &&
            situation[si].getValue(2) <
            (situation[si - 1].getValue(2) * 0.98)) {
            return false;
        } else if (sign == -1 &&
                   situation[si].getValue(2) >
                   (situation[si - 1].getValue(2) * 1.02)) {
            return false;
        }// if.
    } // for.
    return true;
}
// }}}
////////////////////
// }}}


////////////////////
// Class: DataInterpolation
////////////////////
// {{{
void DataInterpolation::interpolate(std::unique_ptr<double[]> xarray,
                                    std::unique_ptr<double[]> yarray,
                                    long                      length,
                                    long                      number)
// {{{ Member-Function: interpolate
//// Interpolate data.
{
    double* xary = new double[length];
    double* yary = new double[length];
    for (int i = 0; i != length; ++i) {
        xary[i] = xarray[i];
        yary[i] = yarray[i];
    }
    double max_x = *std::max_element(xary, xary + length);
    double min_x = *std::min_element(xary, xary + length);
    double pace = (max_x - min_x) / number;
    double xi, yi;
    gsl_interp_accel* acc = gsl_interp_accel_alloc ();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, length);
    gsl_spline_init(spline, xary, yary, length);
    for (xi = min_x; xi < max_x; xi += pace) {
        yi = gsl_spline_eval(spline, xi, acc);
        x.push_back(xi);
        y.push_back(yi);
    }// for.
    // free memory.
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    delete [] xary;
    delete [] yary;
}
// }}}


void DataInterpolation::interpolateL(std::unique_ptr<double[]> xarray,
                                     std::unique_ptr<double[]> yarray,
                                     long                      length,
                                     long                      number)
// {{{ Member-Function: interpolate
//// Interpolate data.
{
    double* xaryl = new double[length];
    double* yaryl = new double[length];
    for (int i = 0; i != length; ++i) {
        xaryl[i] = std::log10(xarray[i]);
        yaryl[i] = std::log10(yarray[i]);
    }
    double max_x = *std::max_element(xaryl, xaryl + length);
    double min_x = *std::min_element(xaryl, xaryl + length);
    double pace = (max_x - min_x) / number;
    double xi, yi;
    gsl_interp_accel* acc = gsl_interp_accel_alloc ();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, length);
    gsl_spline_init(spline, xaryl, yaryl, length);
    for (xi = min_x; xi < max_x; xi += pace) {
        yi = gsl_spline_eval(spline, xi, acc);
        x.push_back(std::pow(10,xi));
        y.push_back(std::pow(10,yi));
    }// for.
    // free memory.
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    delete [] xaryl;
    delete [] yaryl;
}
// }}}


size_t DataInterpolation::size()
// {{{ Member-Function: DataInterpolation::size
//// Return the size of DataInterpolation.
{
    return x.size();
}
// }}}


int DataInterpolation::nearPercentNum(const double& percent,
                                      const int     according)
// {{{ Member-Function: nearPercentNum
//// Return the iterator pointed to situation with value near specific percentage.
{
    double value = 0.0;
    if (according == 1) {
        value = percent * ((*max_element(x.begin(), x.end())) -
                           (*min_element(x.begin(), x.end()))) +
            (*min_element(x.begin(), x.end()));
    } else if (according == 2) {
        value = percent * ((*max_element(y.begin(), y.end())) -
                           (*min_element(y.begin(), y.end()))) +
            (*min_element(y.begin(), y.end()));
    }
    double smallest = 0.0;
    std::vector<double> p;
    for (std::vector<double>::size_type i = 0; i != x.size(); ++i) {
        if (according == 1) {
            p.push_back(std::fabs(x[i] - value));
        } else if (according == 2) {
            p.push_back(std::fabs(y[i] - value));
        }
    } // for.
    smallest = *std::min_element(p.begin(), p.end());
    for (std::vector<double>::size_type i = 0; i != p.size(); ++i) {
        if (p[i] == smallest) {
            return i;
        } // if.
    } // for.
    return 0;
}
// }}}


int DataInterpolation::nearValueNum(const double& value,
                                    const int     according)
// {{{ Member-Function: nearValueNum
//// Return the iterator pointed to situation with value near specific percentage.
{

    double smallest = 0.0;
    std::vector<double> p;
    for (std::vector<double>::size_type i = 0; i != x.size(); ++i) {
        // Get the vector of difference between the number and the value.
        if (according == 1) {
            p.push_back(std::fabs(x[i] - value));
        } else if (according == 2) {
            p.push_back(std::fabs(y[i] - value));
        }
    } // for.
    smallest = *std::min_element(p.begin(), p.end());
    for (std::vector<double>::size_type i = 0; i != p.size(); ++i) {
        if (p[i] == smallest) {
            return i;
        } // if.
    } // for.
    return 0;
}
// }}}
////////////////////
// }}}


/////////
// EOF //
/////////
