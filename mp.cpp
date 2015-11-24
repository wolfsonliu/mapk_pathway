#include "mp.hpp"

using namespace std;

// Member-Function
//// Member-Function Parameter::Parameter: constructor of Class Parameter.
Parameter::Parameter()
{
    f[0] = 0.0;
    b[0] = 0.0;
    kf[0] = 0.0;
    kb[0] = 0.0;
}

Parameter::Parameter(size_t dimension, double input_f[], double input_b[], double input_kf[], double input_kb[])
{
    if (dimension > ARRAY_NUM) {
	cout << "dimension too large";
    }
    for (size_t i = 0; i != dimension; ++i) {
	f[i] = input_f[i];
	b[i] = input_b[i];
	kf[i] = input_kf[i];
	kb[i] = input_kb[i];
    }
}
 
void Parameter::setParameter(size_t dimension, double input_f[], double input_b[], double input_kf[], double input_kb[])
{
    if (dimension > ARRAY_NUM) {
	cout << "dimension too large";
    }
    for (int i = 0; i != dimension; ++i) {
	f[i] = input_f[i];
	b[i] = input_b[i];
	kf[i] = input_kf[i];
	kb[i] = input_kb[i];
    }
}

void Parameter::printParameter(size_t dimension)
{
    for (size_t i = 0; i != dimension; ++i) {
	cout << i << ": ";
	cout << f[i] << "\t"
	     << b[i] << "\t"
	     << kf[i] << "\t"
	     << kb[i] << "\n";
    }
}
//// Member-Function Concentrate::Concentrate: constructor of Class Concentrate

Concentrate::Concentrate()
{
    time = 0;
    reactant[0] = 0.0;
}

Concentrate::Concentrate(size_t dimension, double now_t, double y[])
{
    if (dimension > ARRAY_NUM_BIG) {
	cout << "dimension too large";
    }
    time = now_t;
    for (int i = 0; i != dimension; ++i) {
	reactant[i] = y[i];
    }
}

void Concentrate::setTime(double now_t)
{
    time = now_t;
}

void Concentrate::setConcentrate(size_t dimension, double y[])
{
    if (dimension > ARRAY_NUM_BIG) {
	cout << "dimension too large";
    }
    for (size_t i = 0; i != dimension; ++i) {
        reactant[i] = y[i];
    }
}

void Concentrate::printConcentrate(size_t dimension)
{
    cout << "Time: " << time << "\n";
    cout << "Value: "
    for (size_t i = 0; i != dimension; ++i) {
	cout << reactant[i] << "\t";
    }
    cout << endl;

}

// Function
//// Function odeFunction1s1p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
int odeFunction1s1p(double t, const double y[], double f[], void *para)
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

//// Function odeJacobian1s1p: needed by <gsl/gsl_odeiv2.h>, 1 stage, 1 phosphorylation.
int odeJacobian1s1p(double t, const double y[], double *dfdy, double dydt[], void *para)
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

    for (int i = 0; i != NUM_OF_CHEM; ++i) {
	dydt[i] = 0.0;
    }

    return GSL_SUCCESS;
}


//// Function odeFunction1s2p:  needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
int odeFunction1s2p(double t, const double y[], double f[], void *para)
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


//// Function odeJacobian1s2p:  needed by <gsl/gsl_odeiv2.h>, 1 stage, 2 phosphorylations.
int odeJacobian1s2p(double t, const double y[], double *dfdy, double dydt[], void *para)
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

    for (int i = 0; i != NUM_OF_CHEM; ++i) {
	dydt[i] = 0.0;
    }

    return GSL_SUCCESS;
}



//// Function odeRun: using gsl to solve ode.
vector<Concentrate> odeRun(int func, double start_t, double end_t, Parameter *params, double y[], size_t dimension, double pacelen)
{
    // initiate.
    if (func == 1) {
	gsl_odeiv2_system sys = {odeFunction1s1p, odeJacobian1s1p, dimension, &*params};
    } else if (func == 2) {
	gsl_odeive2_system sys = {odeFunction1s2p, odeJacobian1s2p, dimension, &*params};
    }
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    vector<Concentrate>  return_record;
    Concentrate reaction;

    // calculate.
    long times = static_cast<long>((end_t - start_t)/pacelen);
    for (long i = 0; i != times; ++i) {
	double ti = i * pacelen + start_t;
	int status = gsl_odeiv2_driver_apply(d, &start_t, ti, y);

	if (status != GSL_SUCCESS) {
	    printf("Error, return value = %d\n", status);
	    break;
	}
	reaction.setTime(ti);    // set the time of reaction to ti.
	reaction.setConcentrate(dimension, y);    // set the concentrate of reactant to y[].
	return_record.push_back(reaction);    // record the situation at ti.
    }
    gsl_odeiv2_driver_free (d);
    return return_record;
}

	
//// Function initiateParameter: set the parameter array to a certain velue.
int initiateParameter(double par[], int len, double value)
{
    for (int i = 0; i != len; ++i) {
	par[i] = value;
    }
    return 0;
}

