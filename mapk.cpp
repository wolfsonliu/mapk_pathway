//-*-coding:utf-8-*-
//Author: Wolfson
//Date: Nov. 29, 2015
#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/config.hpp>
#include <boost/program_options/environment_iterator.hpp>
#include <boost/program_options/eof_iterator.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/version.hpp>
#include "mapk_ode.hpp"


using namespace std;




int main(int argc, char *argv[])
{
    // defination.
    int reaction_type;
    double start_time;    // start of time.
    double end_time;    // end of time.
    size_t num_chem;    // number of chemical equation, dimension of f,b,kf,kb.
    size_t num_reactant;    // number of reactions, dimension of Jacobian.
    vector<string> reactant_name;
    double *reactant;    // initiation of reactant.
    double pace_len;
    double *f;    // forward reaction coefficient for reaction 1.
    double *b;    // backward reaction coefficient for reaction 1.
    double *kf;    // forward reaction coefficient for reaction 2.
    double *kb;    // backward reaction coefficient for reaction 2, usually very small.
    int product_ord;    // the final out put number in reactant, start from 0.
    string infile_name;
    ifstream infile("mapk.cfg");
    string outfile_name;    // output file name.
    ofstream outfile;
    
    ////////////////////
    // boost::program_options
    ////////////////////

    //// used to parse arguments inputed by command line or config file.
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");    // set option description.
    desc.add_options()
	("help,h", "Produce help message.")
	("reaction-type,r", po::value<int>(&reaction_type), "Choose reaction type:\n"
                                                       "\t1: one step reaction with one phosphorylation;\n"
	                                               "\t2: one step reaction with two phosphorylation.")
	("time,t", po::value< vector<double> >()->multitoken(), "Input start time and end time with space as seperator.")
	("pace-length,p", po::value<double>(&pace_len), "Set the pace lenght of reactions.")
	("parameter-f,A", po::value< vector<double> >()->multitoken(),
	 "Forward reaction coefficient for the first <b==f>, must be inputed and seperate by space.")
	("parameter-b,B", po::value< vector<double> >()->multitoken(),
	 "Backward reaction coefficient for the first <b==f>, must be inputed and seperate by space.")
	("parameter-kf,C", po::value< vector<double> >()->multitoken(),
	 "Forward reaction coefficient for the second <kb-=kf>, must be inputed and seperate by space.")
	("parameter-kb,D", po::value< vector<double> >()->multitoken(),
	 "Backward reaction coefficient for the second <kb-=kf>, must be inputed and seperate by space.")
	("reactant,R", po::value< vector<double> >()->multitoken(),
	 "Initiation of reactants, must be inputed and seperate by space.")
	("product-order,P", po::value<int>(&product_ord), "The order of final out put in reactants, from 0.")
	//("input,i", po::value<string>(&infile_name), "The input file name.")
	("output,f", po::value<string>(&outfile_name), "The output file name.");
    
    po::variables_map vm;    // make variables map.
    po::store(po::parse_command_line(argc, argv, desc), vm);    // store command line values to vm.
    po::store(po::parse_config_file(infile, desc, true), vm);    // store config file values to vm.
    po::notify(vm);
    //// set operations' actions.
    //if (vm.count("input")) {
    //	// input options with file.
    //	infile.open(infile_name.c_str(), ifstream::in);
    //    po::store(po::parse_config_file(infile, desc, true), vm);    // store config file values to vm.	
    //}
    
    if (vm.count("help")) {
	// help.
	cout << "********************\n"
	     << "**      mapk      **\n"
	     << "********************\n\n"
	     << "mapk is used for calculation of mapk pathway.\n"; 
	cout << "\nYou can use a config file named as mapk.cfg to set options or use the command line. In the config file, use the long option name and = to set values, and multivalues for one option should be in different lines each with the option name and =.\n\n";
	cout << "Reaction type:\n"
	     << "1: one stage one phosphorylation. \n\tReactants: map3k, map3kp, map3k_ras, map3kp_m3kp, ras, m3kp.\n"
	     << "2: one stage two phosphorylations. \n\tReactants: map3k, map3kp, map3k2p, map3k_ras, map3kp_m3kp, map3k2p_m3kp, ras, m3kp.\n";
	cout << "\n********************\n";
	cout << desc << endl;
	return 1;
    }
    
    if (!vm["time"].empty() && (vm["time"].as< vector<double> >().size() == 2)) {
        // set time.
	start_time = vm["time"].as< vector<double> >()[0];
	end_time = vm["time"].as< vector<double> >()[1];
    }
    
    if (vm.count("parameter-f") && vm.count("parameter-b") &&
	vm.count("parameter-kf") && vm.count("parameter-kb")) {
	// set parameters.
	num_chem = vm["parameter-f"].as< vector<double> >().size();
	f = new double [num_chem];
	b = new double [num_chem];
	kf = new double [num_chem];
	kb = new double [num_chem];
	for (int i = 0; i != num_chem; ++i) {
	    f[i] = vm["parameter-f"].as< vector<double> >()[i];
	    b[i] = vm["parameter-b"].as< vector<double> >()[i];
	    kf[i] = vm["parameter-kf"].as< vector<double> >()[i];
	    kb[i] = vm["parameter-kb"].as< vector<double> >()[i];
	}
    }

    if (vm.count("reactant")) {
	// set reactant concentration.
	num_reactant = vm["reactant"].as< vector<double> >().size();
	reactant = new double [num_reactant];
	for (int i = 0; i != num_reactant; ++i) {
	    reactant[i] = vm["reactant"].as< vector<double> >()[i];
	}
    }
    ////////////////////
   
    reactant_name = reactantName(reaction_type); 
    Parameter para(num_chem, f, b, kf, kb);
    Parameter *params = &para;
    outfile.open(outfile_name.c_str(), ofstream::out | ofstream::app);
    // initiat ReactantConcentration.
    ReactantConcentration reaction(reaction_type, reactant_name, num_reactant, product_ord);
    // run ordinary differential equation.
    reaction.odeRun(start_time, end_time, params, reactant, pace_len);
    // output results.
    reaction.outputList(outfile);
    
    cout << "Finished!\n"
	 << "Result in: "
	 << outfile_name
	 << endl;
    return 0;
    
}
////////////////////
