#include "mp.hpp"
using namespace std;

// Member-Function
//// Member-Function Parameter::Parameter: constructor of Class Parameter.
Parameter::Parameter()
{
    function_num = 0;
    f[0] = 0.0;
    b[0] = 0.0;
    kf[0] = 0.0;
    kb[0] = 0.0;
}
//// Member-Function Parameter::Parameter
Parameter::Parameter(size_t dimension, double input_f[], double input_b[], double input_kf[], double input_kb[])
{
    function_num = dimension;
    f = new double[function_num];
    b = new double[function_num];
    kf = new double[function_num];
    kb = new double[function_num];
    for (size_t i = 0; i != function_num; ++i) {
	f[i] = input_f[i];
	b[i] = input_b[i];
	kf[i] = input_kf[i];
	kb[i] = input_kb[i];
    }
}
//// Member-Function Parameter::setParameter: set the value.
void Parameter::setParameter(double input_f[], double input_b[], double input_kf[], double input_kb[])
{
    if (function_num == 0) {
	throw std::runtime_error("function number is not set.");
    }
    for (int i = 0; i != function_num; ++i) {
	f[i] = input_f[i];
	b[i] = input_b[i];
	kf[i] = input_kf[i];
	kb[i] = input_kb[i];
    }
}
//// Member-Function Parameter::printParameter: print the values.
void Parameter::printParameter()
{
    for (size_t i = 0; i != function_num; ++i) {
	std::cout << i << ": ";
	std::cout << f[i] << "\t"
	     << b[i] << "\t"
	     << kf[i] << "\t"
	     << kb[i] << "\n";
    }
}

//// Member-Function
//// Concentration::Concentration: constructor of Class Concentration.
Concentration::Concentration()
{
    is_new = true;
    reactant_num = 0;
    time = 0;
    reactant = new double[1];
    reactant[0] = 0;
}

//// Member-Function
//// Concentration::Concentration: constructor of Class Concentration.
Concentration::Concentration(size_t dimension, double now_t, double y[])
{
    reactant_num = dimension;
    if (reactant_num == 0) {
	throw std::runtime_error("reactant number is not set.");
    }
    reactant = new double [reactant_num];
    time = now_t;
    for (int i = 0; i != reactant_num; ++i) {
	reactant[i] = y[i];
    }
}
//// Member-Function
//// Concentration::setTime: set the value of time
void Concentration::setTime(double now_t)
{
    time = now_t;
}
//// Member-Function
//// Concentration::setConcentration: set the value of Concentrations of reactant.
void Concentration::setConcentration(size_t dimension, double y[])
{
    reactant_num = dimension;
    reactant = new double [reactant_num];
    for (size_t i = 0; i != reactant_num; ++i) {
        reactant[i] = y[i];
    }
    is_new = false;
}
//// Member-Function
//// Concentration::printConcentration: print the value to screen.
void Concentration::printConcentration()
{
    std::cout << "Time: " << time << "\n";
    std::cout << "Value: ";
    for (size_t i = 0; i != reactant_num; ++i) {
	std::cout << reactant[i] << "\t";
    }
    std::cout << std::endl;

}
//// Member-Function Concentration::delConcentration: delete the array.
void Concentration::delConcentration()
// delete the array.
{
    if (!is_new) {
	delete [] reactant;
    }
    is_new = true;
}


//// Member-Function
//// ReactantConcentration::ReactantConcentration: initiation.
ReactantConcentration::ReactantConcentration(int type, std::vector<std::string> react_name, size_t react_num)
{
    reaction_type = type;
    name = react_name;
    reactant_num = react_num;
}


//// Member-Function
//// ReactantConcentration::setReactantName: store the names of reactants in the list.
void ReactantConcentration::setReactantName(const std::vector<std::string> &namelist)
{
    name = namelist;
    reactant_num = name.size();
}

//// Member-Function
//// ReactantConcentration::odeFunction1s1p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
int ReactantConcentration::odeFunction1s1p(double t, const double y[], double f[], void *para)
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k = y[0];    // MAPKKK
    double map3kp = y[1];    // MAPKKK-phosphoryl
    double map3k_ras = y[2];    // intermedia product
    double map3kp_m3kp = y[3];    // intermedia product
    double ras = y[4];    // MAPKKK kinase
    double m3kp = y[5];    // MAPKKK-p phosphatase 
    
    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map3k
    f[0] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras + params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // map3kp
    f[1] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp + params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras;
    // map3k_ras
    f[2] = params.f[0] * map3k * ras - params.b[0] * map3k_ras - params.kf[0] * map3k_ras + params.kb[0] * map3kp * ras;
    // map3kp_m3kp
    f[3] = params.f[1] * map3kp * m3kp - params.b[1] * map3kp_m3kp - params.kf[1] * map3kp_m3kp + params.kb[1] * map3k * m3kp;
    // ras: dras = -dmap3k_ras
    f[4] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras + params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras;
    // m3kp: dm3kp = -dmap3kp_m3kp
    f[5] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp + params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;    
    
    return GSL_SUCCESS;    // in gsl header.
}

//// Member-Function
//// ReactantConcentration::odeJacobian1s1p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
int ReactantConcentration::odeJacobian1s1p(double t, const double y[], double *dfdy, double dfdt[], void *para)
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    // initiate.
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, m3kp as the phosphatase for map3kp.
    double map3k = y[0];    // MAPKKK
    double map3kp = y[1];    // MAPKKK-phosphoryl
    double map3k_ras = y[2];    // intermedia product
    double map3kp_m3kp = y[3];    // intermedia product
    double ras = y[4];    // MAPKKK kinase
    double m3kp = y[5];    // MAPKKK-p phosphatase 

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 6, 6);
    gsl_matrix *m = &dfdy_mat.matrix;
    
    // calculate matrix.
    // map3k.
    gsl_matrix_set(m, 0, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 0, 1, params.b[0]);
    gsl_matrix_set(m, 0, 2, -params.f[0] * map3k);
    gsl_matrix_set(m, 0, 3, 0.0);
    gsl_matrix_set(m, 0, 4, params.kf[1]);
    gsl_matrix_set(m, 0, 5, 0.0);
    // map3kp.
    gsl_matrix_set(m, 1, 0, 0.0);
    gsl_matrix_set(m, 1, 1, params.kf[0]);
    gsl_matrix_set(m, 1, 2, 0.0);
    gsl_matrix_set(m, 1, 3, -params.f[1] * m3kp);
    gsl_matrix_set(m, 1, 4, params.b[1]);
    gsl_matrix_set(m, 1, 5, -params.f[1] * map3kp);
    // map3k_ras.
    gsl_matrix_set(m, 2, 0, params.f[0] * ras);
    gsl_matrix_set(m, 2, 1, -params.b[0] - params.f[0]);
    gsl_matrix_set(m, 2, 2, params.f[0] * map3k);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 2, 4, 0.0);
    gsl_matrix_set(m, 2, 5, 0.0);
    // map3kp_m3kp.
    gsl_matrix_set(m, 3, 0, 0.0);
    gsl_matrix_set(m, 3, 1, 0.0);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, params.f[1] * m3kp);
    gsl_matrix_set(m, 3, 4, -params.b[1] - params.kf[1]);
    gsl_matrix_set(m, 3, 5, params.f[1] * map3kp);
    // ras.
    gsl_matrix_set(m, 4, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 4, 1, params.b[0] + params.f[0]);
    gsl_matrix_set(m, 4, 2, -params.f[0] * map3k);
    gsl_matrix_set(m, 4, 3, 0.0);
    gsl_matrix_set(m, 4, 4, 0.0);
    gsl_matrix_set(m, 4, 5, 0.0);
    // m3kp.
    gsl_matrix_set(m, 5, 0, 0.0);
    gsl_matrix_set(m, 5, 1, 0.0);
    gsl_matrix_set(m, 5, 2, 0.0);
    gsl_matrix_set(m, 5, 3, -params.f[1] * m3kp);
    gsl_matrix_set(m, 5, 4, params.b[1] + params.kf[1]);
    gsl_matrix_set(m, 5, 5, -params.f[1] * map3kp);

    for (int i = 0; i != 6; ++i) {
	dfdt[i] = 0.0;
    }

    return GSL_SUCCESS;
}


//// Member-Function
//// ReactantConcentration::odeFunction1s2p:  needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
int ReactantConcentration::odeFunction1s2p(double t, const double y[], double f[], void *para)
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    Parameter params = *static_cast<Parameter *>(para);
    // use ras as the kinase for map3k, pase as the phosphatase for map3kp.
    double map3k = y[0];    // MAPKKK
    double map3kp = y[1];    // MAPKKK-phosphoryl
    double map3k2p = y[2];    // MAPKKK-phosphoryl-phosphoryl
    double map3k_ras = y[3];    // intermedia product
    double map3kp_m3kp = y[4];    // intermedia product
    double map3kp_ras = y[5];    // intermedia product
    double map3k2p_m3kp = y[6];    // intermedia product
    double ras = y[7];    // MAPKKK kinase
    double m3kp = y[8];    // MAPKKK-p phosphatase

    // map3k + ras <b0==f0> map3k_ras <kb0-=kf0> map3kp + ras
    // map3kp + m3kp <b1==f1> map3kp_m3kp <kb1-=kf1> map3k + m3kp
    // map3kp + ras <b2==f2> map3kp_ras <kb2-=kf2> map3k2p + ras
    // map3kp + m3kp <b3==f3> map3k2p_m3kp <kb3-=kf3> map3kp + m3kp
    // map3k
    f[0] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras + params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp;
    // map3kp
    f[1] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp + params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras -
	params.f[2] * map3kp * ras + params.b[2] * map3kp_ras + params.kf[3] * map3k2p_m3kp - params.kb[3] * map3kp * m3kp;
    // map3k2p
    f[2] = -params.f[3] * map3k2p * m3kp + params.b[3] * map3k2p_m3kp + params.kf[2] * map3kp_ras - params.kb[2] * map3k2p * ras;
    // map3k_ras
    f[3] = params.f[0] * map3k * ras - params.b[0] * map3k_ras - params.kf[0] * map3k_ras + params.kb[0] * map3kp * ras;
    // map3kp_m3kp
    f[4] = params.f[1] * map3kp * m3kp - params.b[1] * map3kp_m3kp - params.kf[1] * map3kp_m3kp + params.kb[1] * map3k * m3kp;
    // map3kp_ras
    f[5] = params.f[2] * map3kp * ras - params.b[2] * map3kp_ras - params.kf[2] * map3kp_ras + params.kb[2] * map3kp * ras;
    // map3k2p_m3kp
    f[6] = params.f[3] * map3k2p * m3kp - params.b[3] * map3k2p_m3kp - params.kf[3] * map3k2p_m3kp + params.kb[3] * map3kp * m3kp;
    // ras: dras = - dmap3k_ras - dma3kp_ras
    f[7] = -params.f[0] * map3k * ras + params.b[0] * map3k_ras + params.kf[0] * map3k_ras - params.kb[0] * map3kp * ras -
	params.f[2] * map3kp * ras + params.b[2] * map3kp_ras + params.kf[2] * map3kp_ras - params.kb[2] * map3k2p * ras;
    // m3kp: dm3kp = - dmap3kp_m3kp - dmap3k2p_m3kp
    f[8] = -params.f[1] * map3kp * m3kp + params.b[1] * map3kp_m3kp + params.kf[1] * map3kp_m3kp - params.kb[1] * map3k * m3kp -
	params.f[3] * map3k2p * m3kp + params.b[3] * map3k2p_m3kp + params.kf[3] * map3k2p_m3kp - params.kb[3] * map3kp * m3kp;
    
    return GSL_SUCCESS;    // in gsl header.
}


//// Member-Function
//// ReactantConcentration::odeJacobian1s2p:  needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
int ReactantConcentration::odeJacobian1s2p(double t, const double y[], double *dfdy, double dfdt[], void *para)
{
    static_cast<void>(t);    // convert type to avoid unused parameter warning
    // initiate.
    Parameter params = *static_cast<Parameter *>(para);
    double map3k = y[0];    // MAPKKK
    double map3kp = y[1];    // MAPKKK-phosphoryl
    double map3k2p = y[2];    // MAPKKK-phosphoryl-phosphoryl
    double map3k_ras = y[3];    // intermedia product
    double map3kp_m3kp = y[4];    // intermedia product
    double map3kp_ras = y[5];    // intermedia product
    double map3k2p_m3kp = y[6];    // intermedia product
    double ras = y[7];    // MAPKKK kinase
    double m3kp = y[8];    // MAPKKK-p phosphatase
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 6, 6);
    gsl_matrix *m = &dfdy_mat.matrix;

    // calculate matrix.
    // map3k.
    gsl_matrix_set(m, 0, 0, -params.f[0] * ras - params.kb[1] * m3kp);
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
    gsl_matrix_set(m, 1, 1, -params.f[1] * m3kp - params.kb[0] * ras - params.f[2] * ras - params.kb[3] * m3kp);
    gsl_matrix_set(m, 1, 2, 0.0);
    gsl_matrix_set(m, 1, 3, params.kf[0]);
    gsl_matrix_set(m, 1, 4, params.b[1]);
    gsl_matrix_set(m, 1, 5, params.b[2]);
    gsl_matrix_set(m, 1, 6, params.kf[3]);
    gsl_matrix_set(m, 1, 7, -params.kb[0] * map3k - params.f[2] * map3kp);
    gsl_matrix_set(m, 1, 8, -params.f[1] * map3kp - params.kb[3] * map3kp);
    // map3k2p.
    gsl_matrix_set(m, 2, 0, 0.0);
    gsl_matrix_set(m, 2, 1, 0.0);
    gsl_matrix_set(m, 2, 2, -params.f[3] * m3kp - params.kb[2] * ras);
    gsl_matrix_set(m, 2, 3, 0.0);
    gsl_matrix_set(m, 2, 4, 0.0);
    gsl_matrix_set(m, 2, 5, params.kf[2]);
    gsl_matrix_set(m, 2, 6, params.b[3]);
    gsl_matrix_set(m, 2, 7, -params.kb[2] * map3k2p);
    gsl_matrix_set(m, 2, 8, -params.f[3] * map3k2p);
    // map3k_ras.
    gsl_matrix_set(m, 3, 0, params.f[0] * ras + params.kb[1] * m3kp);
    gsl_matrix_set(m, 3, 1, params.kb[0] * ras);
    gsl_matrix_set(m, 3, 2, 0.0);
    gsl_matrix_set(m, 3, 3, -params.b[0] - params.kf[0]);
    gsl_matrix_set(m, 3, 4, 0.0);
    gsl_matrix_set(m, 3, 5, 0.0);
    gsl_matrix_set(m, 3, 6, 0.0);
    gsl_matrix_set(m, 3, 7, params.f[0] * map3k + params.kb[0] * map3kp);
    gsl_matrix_set(m, 3, 8, 0.0);
    // map3kp_m3kp.
    gsl_matrix_set(m, 4, 0, 0.0);
    gsl_matrix_set(m, 4, 1, params.f[1] * m3kp);
    gsl_matrix_set(m, 4, 2, 0.0);
    gsl_matrix_set(m, 4, 3, 0.0);
    gsl_matrix_set(m, 4, 4, -params.b[1] - params.kf[1]);
    gsl_matrix_set(m, 4, 5, 0.0);
    gsl_matrix_set(m, 4, 6, 0.0);
    gsl_matrix_set(m, 4, 7, 0.0);
    gsl_matrix_set(m, 4, 8, params.f[1] * map3kp + params.kb[1] * map3k);
    // map3kp_ras.
    gsl_matrix_set(m, 5, 0, 0.0);
    gsl_matrix_set(m, 5, 1, params.f[2] * ras + params.kb[2] * ras);
    gsl_matrix_set(m, 5, 2, 0.0);
    gsl_matrix_set(m, 5, 3, 0.0);
    gsl_matrix_set(m, 5, 4, 0.0);
    gsl_matrix_set(m, 5, 5, -params.b[2] - params.kf[2]);
    gsl_matrix_set(m, 5, 6, 0.0);
    gsl_matrix_set(m, 5, 7, params.f[2] * map3kp + params.kb[2] * map3kp);
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
    gsl_matrix_set(m, 6, 8, params.f[3] * map3k2p + params.kb[3] * map3kp);
    // ras.
    gsl_matrix_set(m, 7, 0, -params.f[0] * ras);
    gsl_matrix_set(m, 7, 1, -params.kb[0] * ras - params.f[2] * ras);
    gsl_matrix_set(m, 7, 2, -params.kb[2] * ras);
    gsl_matrix_set(m, 7, 3, params.b[0] + params.kf[0]);
    gsl_matrix_set(m, 7, 4, 0.0);
    gsl_matrix_set(m, 7, 5, params.b[2] + params.kf[2]);
    gsl_matrix_set(m, 7, 6, 0.0);
    gsl_matrix_set(m, 7, 7, -params.f[0] * map3k - params.kb[0] * map3kp - params.f[2] * map3kp - params.kb[2] * map3k2p);
    gsl_matrix_set(m, 7, 8, 0.0);
    // m3kp.
    gsl_matrix_set(m, 8, 0, -params.kb[1] * m3kp);
    gsl_matrix_set(m, 8, 1, -params.f[1] * m3kp - params.kb[3] * m3kp);
    gsl_matrix_set(m, 8, 2, -params.f[3] * m3kp);
    gsl_matrix_set(m, 8, 3, 0.0);
    gsl_matrix_set(m, 8, 4, params.b[1] + params.kf[1]);
    gsl_matrix_set(m, 8, 5, 0.0);
    gsl_matrix_set(m, 8, 6, params.b[3] + params.kf[3]);
    gsl_matrix_set(m, 8, 7, 0.0);
    gsl_matrix_set(m, 8, 8, -params.f[1] * map3kp - params.kb[1] * map3k - params.f[3] * map3k2p - params.kb[3] * map3kp);

    for (int i = 0; i != 9; ++i) {
	dfdt[i] = 0.0;
    }

    return GSL_SUCCESS;
}


//// Member-Function
//// ReactantConcentration::odeRun: using gsl to solve ode.
void ReactantConcentration::odeRun(double start_t, double end_t, Parameter *params, double reactant[], double pacelen)
{
    // initiate.
    double *y = reactant;
    typedef int (* Func)(double t, const double y[], double f[], void *para);
    typedef int (* Jac)(double t, const double y[], double *dfdy, double dydt[], void *para);

    Func odeFunction;    // function pointer.
    Jac odeJacobian;    // function pointer.
    
    if (reaction_type == 1) {
	odeFunction = odeFunction1s1p;
	odeJacobian = odeJacobian1s1p;
     } else if (reaction_type == 2) {
	odeFunction = odeFunction1s2p;
	odeJacobian = odeJacobian1s2p;
    }    // HERE!!!!!!!!!!!!!!!!!!!!!!signal handler.
    gsl_odeiv2_system sys = {odeFunction, odeJacobian, reactant_num, &*params};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    Concentration reaction;    // temporary store y.

    //std::ofstream outfile("infuncresult.csv", ofstream::out | ofstream::app);
    //outfile << "Time";
    //for (size_t i = 0; i != reactant_num; ++i) {
    //	outfile << "," << name[i];
    //}
    //outfile << "\n";
    
    // calculate.
    do {
	long times = static_cast<long>((end_t - start_t)/pacelen);
	for (long i = 0; i != times; ++i) {
	    double ti = i * pacelen + start_t;
	    int status = gsl_odeiv2_driver_apply(d, &start_t, ti, y);
	    if (status != GSL_SUCCESS) {
		printf("Error, return value = %d\n", status);
		break;
	    }
	    reaction.setTime(ti);    // set the time of reaction to ti.
	    reaction.setConcentration(reactant_num, y);    // set the concentrate of reactant to y[]
	    //outfile << reaction.time;
	    //for (int i = 0; i != reactant_num; ++i) {
	    //	outfile << "," << reaction.reactant[i];
	    //}
	    //outfile << "\n";

	    //reaction.printConcentration();
	    list.push_back(reaction);    // record the situation at ti.
	}
    } while (0);// *********************** add here

    //outfile.close();
    //outfile.clear();
    gsl_odeiv2_driver_free (d);
}


// Functions

//// Function
//// initiateParameter: set the parameter array to a certain velue.
int initiateParameter(double par[], int len, double value)
{
    for (int i = 0; i != len; ++i) {
	par[i] = value;
    }
    return 0;
}

//// Function
//// reactantName: return the vector of strings in a kind of reaction.
std::vector<std::string> reactantName(int functype)
{
    std::vector<std::string> name_vector;
    if (functype == 1) {
	name_vector.push_back("map3k");
	name_vector.push_back("map3kp");
	name_vector.push_back("map3k_ras");
	name_vector.push_back("map3kp_m3kp");
	name_vector.push_back("ras");
	name_vector.push_back("m3kp");

    } else if (functype == 2) {
	name_vector.push_back("map3k");
	name_vector.push_back("map3kp");
	name_vector.push_back("map3k2p");
	name_vector.push_back("map3k_ras");
	name_vector.push_back("map3kp_m3kp");
	name_vector.push_back("map3kp_ras");
	name_vector.push_back("map3k2p_m3kp");
	name_vector.push_back("ras");
	name_vector.push_back("m3kp");       
    }

    return name_vector;
}
