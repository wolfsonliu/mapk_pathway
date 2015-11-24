#include "mp.hpp"

using namespace std;


int main(int argc, char *argv[])
{
    // using inside.
    double tstart = 0.0;    // start of time.
    double tend = 0.3;    // end of time.
    double lenpace = 0.01;    // length of pace.
    int num_reactant_tmp = 6;
    int num_chem_tmp = 2;
    int reaction_type = 1;
    // defination.

    //// variables
    vector<Concentrate> result;    // record result.
    double start_time = tstart;    // start of time.
    double end_time = tend;    // end of time.
    size_t num_reactant = num_reactant_tmp;    // dimension of Jacobian.
    size_t num_chem = num_chem_tmp;
    double pace_len = lenpace;    // length of pace.
    double f[num_chem];    // forward reaction coefficient for reaction 1.
    double b[num_chem];    // backward reaction coefficient for reaction 1.
    double kf[num_chem];    // forward reaction coefficient for reaction 2.
    double kb[num_chem];    // backward reaction coefficient for reaction 2, usually very small.
    initiateParameter(f, num_chem, 10.00);
    initiateParameter(b, num_chem, 1.50);
    initiateParameter(kf, num_chem, 1.50);
    initiateParameter(kb, num_chem, 0.03);
    cout << f[0] << b[0] << kf[0] << kb[0] << endl;
    Parameter para(num_chem, f, b, kf, kb);
    // para.printParameter(num_chem);
    Parameter *params = &para;
    double reactant[num_chem] = {3e-3, 0.0, 1, 0.0, 0.0, 3e-4};

    

    // calculate.
    result = odeRun(reaction_type, start_time, end_time, params, reactant, num_reactant, pace_len);
    ofstream outfile("result.csv", ofstream::out | ofstream::app);
    outfile << "Time,map3k,map3kp,map3k_ras,map3kp_m3kp,ras,m3kp\n";
    for (vector<Concentrate>::size_type vi = 0;
	 vi != result.size();
	 ++vi) {
	outfile << result[vi].time;
	for (int i = 0; i != NUM_OF_REACTANT; ++i) {
	    outfile << "," << result[vi].reactant[i];
	}
	outfile << "\n";
    }
    outfile.close();
    outfile.clear();
    return 0;
    
}
