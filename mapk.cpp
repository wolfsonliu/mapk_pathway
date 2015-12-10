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
#include <sstream>
#include "mapk_ode.hpp"
#include "parameter.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    ////////////////////
    // define and declaration
    ////////////////////

    int      reaction_type;
    double   start_time;
    // start of time.
    double   end_time;
    // end of time.
    size_t   num_chem;
    // number of chemical equation, dimension of f,b,kf,kb.
    size_t   num_reactant;
    // number of reactions, dimension of Jacobian.
    long     num_sample;
    // number of latin hypercube sampling.
    double*  reactant;
    // initiation of reactant.
    double   input_roof;
    double   input_floor;
    double   pace_len;
    double*  f;
    // forward reaction coefficient for reaction 1.
    double   f_roof;
    double   f_floor;
    double*  b;
    // backward reaction coefficient for reaction 1.
    double   b_roof;
    double   b_floor;
    double*  kf;
    // forward reaction coefficient for reaction 2.
    double   kf_roof;
    double   kf_floor;
    double*  kb;
    // backward reaction coefficient for reaction 2, usually very small.
    double   kb_roof;
    double   kb_floor;
    int      input_ord;
    // the input reactant (ras) number in reactant, start from 0.
    int      output_ord;
    // the output reactant (mapk2p) number in reactant, start from 0.
    string   infile_name;
    ifstream infile("mapk.cfg");
    string   outfile_name;
    // output file name.
    ofstream outfile;
    // output file
    string   detail_outfile_name;
    ofstream detail_outfile;
    string   parameter_outfile_name;
    ofstream parameter_outfile;
    string   interpolation_outfile_name;
    ofstream interpolation_outfile;
    int detail;
    ////////////////////
    
    ////////////////////
    // boost::program_options
    ////////////////////

    //// used to parse arguments inputed by command line or config file.
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");    // set option description.
    desc.add_options()
        ("help,h", "Produce help message.")
        ("reaction-type,r",
         po::value<int>(&reaction_type),
         "Choose reaction type:\n"
         "\t1: one step reaction with one phosphorylation;\n"
         "\t2: one step reaction with two phosphorylation.")
        ("time,t",
         po::value< vector<double> >()->multitoken(),
         "Input start time and end time with space as seperator.")
        ("pace-length,p",
         po::value<double>(&pace_len),
         "Set the pace lenght of reactions.")
        ("sample-number",
         po::value<long>(&num_sample)->default_value(0),
         "Set the number of latin hypercube sampling.")
        ("parameter-f,A",
         po::value< vector<double> >()->multitoken(),
         "Forward reaction coefficient for the first <b==f>,"
         " must be inputed and seperate by space.")
        ("parameter-f-roof",
         po::value<double>(&f_roof)->default_value(0.0),
         "Forward reaction coefficient roof for the first <b==f>, used for sampling.")
        ("parameter-f-floor",
         po::value<double>(&f_floor)->default_value(0.0),
         "Forward reaction coefficient floor for the first <b==f>, used for sampling.")
        ("parameter-b,B",
         po::value< vector<double> >()->multitoken(),
         "Backward reaction coefficient for the first <b==f>,"
         " must be inputed and seperate by space.")
        ("parameter-b-roof",
         po::value<double>(&b_roof)->default_value(0.0),
         "Backward reaction coefficient roof for the first <b==f>, used for sampling.")
        ("parameter-b-floor",
         po::value<double>(&b_floor)->default_value(0.0),
         "Backward reaction coefficient floor for the first <b==f>, used for sampling.")
        ("parameter-kf,C",
         po::value< vector<double> >()->multitoken(),
         "Forward reaction coefficient for the second <kb-=kf>,"
         " must be inputed and seperate by space.")
        ("parameter-kf-roof",
         po::value<double>(&kf_roof)->default_value(0.0),
         "Forward reaction coefficient roof for the second <kb-=kf>, used for sampling.")
        ("parameter-kf-floor",
         po::value<double>(&kf_floor)->default_value(0.0),
         "Forward reaction coefficient floor for the second <kb-=kf>, used for sampling.")
        ("parameter-kb,D",
         po::value< vector<double> >()->multitoken(),
         "Backward reaction coefficient for the second <kb-=kf>,"
         " must be inputed and seperate by space.")
        ("parameter-kb-roof",
         po::value<double>(&kb_roof)->default_value(0.0),
         "Backward reaction coefficient roof for the second <kb-=kf>, used for sampling.")
        ("parameter-kb-floor",
         po::value<double>(&kb_floor)->default_value(0.0),
         "Backward reaction coefficient floor for the second <kb-=kf>, used for sampling.")
        ("reactant,R",
         po::value< vector<double> >()->multitoken(),
         "Initiation of reactants, must be inputed and seperate by space.")
        ("input-roof",
         po::value<double>(&input_roof)->default_value(0.0),
         "Input roof.")
        ("input-floor",
         po::value<double>(&input_floor)->default_value(0.0),
         "Input floor.")
        ("input-order,I",
         po::value<int>(&input_ord),
         "The order of final out put in reactants, from 0.")
        ("output-order,O",
         po::value<int>(&output_ord),
         "The order of final out put in reactants, from 0.")
        ("detail",
         po::value<int>(&detail)->default_value(1),
         "The order of final out put in reactants, from 0.");
    //("input,i", po::value<string>(&infile_name), "The input file name.")
    //("output,f", po::value<string>(&outfile_name), "The output file name.");
    
    po::variables_map vm;    // make variables map.
    po::store(po::parse_command_line(argc, argv, desc), vm);    // store command line values to vm.
    po::store(po::parse_config_file(infile, desc, true), vm);    // store config file values to vm.
    po::notify(vm);
    //// set operations' actions.
    //if (vm.count("input")) {
    //    // input options with file.
    //    infile.open(infile_name.c_str(), ifstream::in);
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
             << "1: one stage one phosphorylation. \n"
             << "        Reactants: map3k, map3kp,\n"
             << "                   map3k_ras, map3kp_m3kp,\n"
             << "                   ras, m3kp.\n"
             << "2: one stage two phosphorylations. \n"
             << "        Reactants: map3k, map3kp, map3k2p,\n"
             << "                   map3k_ras, map3kp_m3kp, map3k2p_m3kp,\n"
             << "                   ras, m3kp.\n"
             << "3: two stages one phosphorylation. \n"
             << "        Reactants: map3k, map3kp, map2k, map2kp,\n"
             << "                   map3k_ras, map3kp_m3kp, map2k_map3kp, map2kp_m2kp,\n"
             << "                   ras, m3kp, m2kp.\n"
             << "4: two stages two phosphorylations. \n"
             << "        Reactants: map3k, map3kp, map2k, map2kp, map2k2p,\n"
             << "                   map3k_ras, map3kp_m3kp, map2k_map3kp, map2kp_m2kp, map2kp_map3kp, map2k2p_m2kp,\n"
             << "                   ras, m3kp, m2kp\n"
             << "5: three stages one phosphorylation. \n"
             << "        Reactants: map3k, map3kp, map2k, map2kp, mapk, mapkp,\n"
             << "                   map3k_ras, map3kp_m3kp, map2k_map3kp, map2kp_m2kp, mapk_map2kp, mapkp_mkp,\n"
             << "                   ras, m3kp, m2kp, mkp\n"
             << "6: three stages two phosphorylations. \n"
             << "        Reactants: map3k, map3kp, map2k, map2kp, map2k2p, mapk, mapkp, mapk2p,\n"
             << "                   map3k_ras, map3kp_m3kp, map2k_map3kp, map2kp_m2kp, map2kp_map3kp, \n"
             << "                   map2k2p_m2kp, mapk_map2k2p, mapkp_mkp, mapkp_map2k2p, mapk2p_mkp, \n"
             << "                   ras, m3kp, m2kp, mkp\n";
        
        cout << "\n********************\n";
        cout << desc << endl;
        return 1;
    }

    if (!vm["time"].empty() &&
        (vm["time"].as< vector<double> >().size() == 2)) {
        // set time.
        start_time = vm["time"].as< vector<double> >()[0];
        end_time   = vm["time"].as< vector<double> >()[1];
    }
    
    if (vm.count("parameter-f") &&
        vm.count("parameter-b") &&
        vm.count("parameter-kf") &&
        vm.count("parameter-kb")) {
        // set parameters.
        num_chem = vm["parameter-f"].as< vector<double> >().size();
        f  = new double [num_chem];
        b  = new double [num_chem];
        kf = new double [num_chem];
        kb = new double [num_chem];
        for (int i = 0; i != num_chem; ++i) {
            f[i]  = vm["parameter-f"].as< vector<double> >()[i];
            b[i]  = vm["parameter-b"].as< vector<double> >()[i];
            kf[i] = vm["parameter-kf"].as< vector<double> >()[i];
            kb[i] = vm["parameter-kb"].as< vector<double> >()[i];
        }
        ////////////////////

        ////////////////////
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

    ////////////////////
    // open file
    ////////////////////
    ostringstream outstr;     // output.
    ostringstream ioutstr;    // interpolation output.
    ostringstream doutstr;    // detail output.
    ostringstream poutstr;    // parameter output.

    // set the name of parameter_outfile_name.
    poutstr << "mapk_type"
            << reaction_type
            << "_parameter.txt";
    parameter_outfile_name = poutstr.str();
    // open parameter_outfile.
    parameter_outfile.open(parameter_outfile_name.c_str(),
                           ofstream::out | ofstream::app);
    outstr <<"mapk_type"
           << reaction_type
           << "_.csv";
    outfile_name = outstr.str();
    outfile.open(outfile_name.c_str(),
                 ofstream::out | ofstream::app);
    doutstr <<"mapk_type"
            << reaction_type
            << "_detail.csv";
    detail_outfile_name = doutstr.str();
    detail_outfile.open(detail_outfile_name.c_str(),
                        ofstream::out | ofstream::app);
    ioutstr <<"mapk_type"
            << reaction_type
            << "_interpolation.csv";
    interpolation_outfile_name = doutstr.str();
    interpolation_outfile.open(interpolation_outfile_name.c_str(),
                               ofstream::out | ofstream::app);    ////////////////////
    double *reactant_ode = new double[num_reactant];
    for (int i = 0; i != num_reactant; ++i) {
        reactant_ode[i] = reactant[i];
    }    
    ////////////////////


    
    if (num_sample == 0) {
        Parameter para(num_chem, f, b, kf, kb);
        Parameter *param = &para;
        // initiat ReactantConcentration.
        ReactantConcentration reaction(reaction_type);
        // run ordinary differential equation.
        reaction.odeRun(start_time,
                        end_time,
                        param,
                        reactant_ode,
                        pace_len);
        // output results.
        //if (detail == 1) {
            reaction.outputList(detail_outfile);
            //}

        para.outputParameter(parameter_outfile);


    } else {
        LatinVector   latin_f(f_floor, f_roof, num_sample);
        LatinVector   latin_b(b_floor, b_roof, num_sample);
        LatinVector   latin_kf(kf_floor, kf_roof, num_sample);
        LatinVector   latin_kb(kb_floor, kb_roof, num_sample);
        Parameter     para_l(num_chem);
        DataInterpolation input_output;
        DataInterpolation input_time;
        DataInterpolation input_dissipation;

        // write the head of parameter_outfile.
        parameter_outfile << "Num,f,b,kf,kb\n";
        outfile << "input,output,time,dissipation\n";
        for (long i = 0; i != num_sample; ++i) {
            
            
            // Reset parameter to the latin hypercube sampling
            // parameters.
            para_l.setParameter(latin_f[i],
                                latin_b[i],
                                latin_kf[i],
                                latin_kb[i]);

            // Write the information of latin hypercube sampling
            // parameters.
            parameter_outfile << "Latin Hypercube Sampling: "
                              << i << "\n";
            para_l.outputParameter(parameter_outfile);
            parameter_outfile << "\n";

            Parameter *          param_l      = & para_l;
            StableList           stable_list(reaction_type);
            std::vector<Stable>::iterator nsit  = stable_list.situation.end();
            // record the notSufficient pointer.
            double               input_value1 = input_floor;
            double               input_value2 = input_roof;
            long                 whilecount   = 0;
            
            
            do {
                for (int i = 0; i != num_reactant; ++i) {
                    // reset the value of reactant.
                    reactant_ode[i] = reactant[i];
                } // for.
                // use a different input.
                reactant_ode[input_ord] = input_value1; 
                ReactantConcentration reaction_l(reaction_type);
                // run ode of reaction.
                reaction_l.odeRun(start_time,
                                  end_time,
                                  param_l,
                                  reactant_ode,
                                  pace_len);
                if (detail == 1) {
                    // output detail data.
                    detail_outfile << "Latin Hypercube Sampling: "
                                   << i
                                   << " input: "
                                   << input_value1
                                   << "\n";
                    reaction_l.outputList(detail_outfile);
                    detail_outfile << "\n";
                } // if.

                
                stable_list.addInter(nsit,
                                     input_value1,
                                     reaction_l.getFinalOutput(),
                                     reaction_l.getFinalTime(),
                                     reaction_l.getFinalDissipation());
                nsit = stable_list.unSufficient();
                if (nsit != stable_list.situation.end()) {
                    // if get unsufficient pointer, do more calculation.
                    input_value1 = (nsit->input + input_value2) / 2;
                    input_value2 = (nsit + 1)->input;
                } else if (whilecount < 2) {
                    // if this is the first two loop.
                    input_value1 = input_roof;
                } else if (whilecount > 2 &&
                           nsit == stable_list.situation.end()) {
                    break;
                } // if.
                ++whilecount;
            } while(whilecount < 1000);
            
            stable_list.sortList(1);    // sort stable_list.
            // interpolate data.
            input_output.interpolate(stable_list.getArray(1),
                                     stable_list.getArray(2),
                                     stable_list.size(),
                                     2000);
            input_time.interpolate(stable_list.getArray(1),
                                   stable_list.getArray(3),
                                   stable_list.size(),
                                   2000);
            input_dissipation.interpolate(stable_list.getArray(1),
                                          stable_list.getArray(4),
                                          stable_list.size(),
                                          2000);
            interpolation_outfile << "input,output,time,dissipation\n";
            for (size_t i = 0; i != input_output.size(); ++i) {
                interpolation_outfile << input_output.x[i]      << ",";
                interpolation_outfile << input_output.y[i]      << ",";
                interpolation_outfile << input_time.y[i]        << ",";
                interpolation_outfile << input_dissipation.y[i] << "\n";
            } // for.
            // output stable_list.
            outfile << "Latin Hypercube Sampling: "
                    << i << "\n";
            stable_list.outputList(outfile);
            outfile << "\n";
        } // for loop for latin hypercube sampling parameter.
    }
    parameter_outfile.close();
    parameter_outfile.clear();
    outfile.close();
    outfile.clear();
    detail_outfile.close();
    detail_outfile.clear();
    interpolation_outfile.close();
    interpolation_outfile.clear();
    cout << "Finished!\n"
         << endl;
    return 0;
    
}
////////////////////
