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
#include <memory>
#include "mapk_ode.hpp"
#include "parameter.hpp"

using namespace std;


struct DeleteDoubleArray { 
    void operator()(double* p) const {
        delete[] p;
    }
};


//// Function: isStable
//////// Judge whether reach stable.
//bool isStable(std::vector<double>& thevector);


int main(int argc, char* argv[])
{
    ////////////////////
    // define and declaration
    ////////////////////

    int      reaction_type;
    double   start_time;
    // start of time.
    double   end_time;
    // end of time.
    int      num_chem;
    // number of chemical equation, dimension of f,b,kf,kb.
    size_t   num_reactant;
    // number of reactions, dimension of Jacobian.
    long     num_sample;
    // number of latin hypercube sampling.
    unique_ptr<double[]>  reactant;
    // initiation of reactant.
    double   input_roof;
    double   input_floor;
    double   pace_len;
    unique_ptr<double[]>  f;
    // forward reaction coefficient for reaction 1.
    double   f_roof;
    double   f_floor;
    unique_ptr<double[]>  b;
    // backward reaction coefficient for reaction 1.
    double   b_roof;
    double   b_floor;
    unique_ptr<double[]>  kf;
    // forward reaction coefficient for reaction 2.
    double   kf_roof;
    double   kf_floor;
    unique_ptr<double[]>  kb;
    // backward reaction coefficient for reaction 2, usually very small.
    double   kb_roof;
    double   kb_floor;
    int      input_ord;
    // the input reactant (ras) number in reactant, start from 0.
    int      output_ord;
    // the output reactant (mapk2p) number in reactant, start from 0.
    string   infile_name;
    ifstream infile("mapk.cfg");
    string   outfile_tag;    // used before outfile name.
    string   outfile_name;    // output file name.
    ofstream outfile;
    //ostringstream outstr;    
    string   detail_outfile_name;
    ofstream detail_outfile;
    //ostringstream detail_outstr;
    string   parameter_outfile_name;
    ofstream parameter_outfile;
    //ostringstream parameter_outstr;
    string   interpolation_outfile_name;
    ofstream interpolation_outfile;
    //ostringstream interpolation_outstr;
    string   calculate_outfile_name;
    ofstream calculate_outfile;
    string   infofile_name;
    ofstream infofile;
    //ostringstream infofile_str;
    int      detail;
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
        ("function-number",
         po::value<int>(&num_chem),
         "Number of chemical functions.\n")
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
         "The order of final out put in reactants, from 0.")
    //("input,i", po::value<string>(&infile_name), "The input file name.")
        ("outfile-tag",
         po::value<string>(&outfile_tag),
         "The tag added before output file name.");
    
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
        f.reset(new double [num_chem]);
        b.reset(new double [num_chem]);
        kf.reset(new double [num_chem]);
        kb.reset(new double [num_chem]);
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
        reactant.reset(new double [num_reactant]);
        for (int i = 0; i != num_reactant; ++i) {
            reactant[i] = vm["reactant"].as< vector<double> >()[i];
        }
    }
    ////////////////////


    ////////////////////
    double* reactant_ode = new double[num_reactant];
   
    ////////////////////
    ostringstream outnamestr;     // output.
    ostringstream ioutnamestr;    // interpolation output.
    ostringstream doutnamestr;    // detail output.
    ostringstream poutnamestr;    // parameter output.
    ostringstream coutnamestr;    // calculate output.
    ostringstream infonamestr;    // infofile name.
    
    poutnamestr << outfile_tag
                << "mapk_t"
                << reaction_type
                << "_parameter.csv";
    parameter_outfile_name = poutnamestr.str();
    parameter_outfile.open(parameter_outfile_name.c_str(),
                           ofstream::out | ofstream::app);
    bool parameter_header = false;
    
    outnamestr << outfile_tag
               <<"mapk_t"
               << reaction_type
               << ".csv";
    outfile_name = outnamestr.str();
    outfile.open(outfile_name.c_str(),
                 ofstream::out | ofstream::app);
    bool outfile_header = false;
    


    doutnamestr << outfile_tag
                <<"mapk_t"
                << reaction_type
                << "_detail.csv";
    detail_outfile_name = doutnamestr.str();
    detail_outfile.open(detail_outfile_name.c_str(),
                        ofstream::out | ofstream::app);
    bool detail_header = false;
    
    ioutnamestr << outfile_tag
                <<"mapk_t"
                << reaction_type
                << "_interpolation.csv";
    interpolation_outfile_name = ioutnamestr.str();
    interpolation_outfile.open(interpolation_outfile_name.c_str(),
                               ofstream::out | ofstream::app);
    bool interpolation_header = false;
    
    coutnamestr << outfile_tag
                <<"mapk_t"
                << reaction_type
                << "_calculate.csv";
    calculate_outfile_name = coutnamestr.str();
    calculate_outfile.open(calculate_outfile_name.c_str(),
                           ofstream::out | ofstream::app);
    bool calculate_header = false;

    infonamestr << outfile_tag
                << "mapk_t"
                << reaction_type
                << "_info.txt";
    infofile_name = infonamestr.str();
    infofile.open(infofile_name.c_str(), ofstream::out | ofstream::app);
    ////////////////////

    
    if (num_sample == 0) {
        Parameter  para(num_chem, f, b, kf, kb);
        Parameter* param = & para;

        for (int i = 0; i != num_reactant; ++i) {
            reactant_ode[i] = reactant[i];
        } 
        // initiat ReactantConcentration.
        ReactantConcentration reaction(reaction_type);
        // run ordinary differential equation.
        reaction.odeRun(start_time,
                        end_time,
                        param,
                        reactant_ode,
                        pace_len);
        // output results.
        detail_outfile << reaction;
        parameter_outfile << *param;

    } else {
        LatinVector   latin_f(log10(f_floor),
                              log10(f_roof),
                              num_sample);
        LatinVector   latin_b(log10(b_floor),
                              log10(b_roof),
                              num_sample);
        LatinVector   latin_kf(log10(kf_floor),
                               log10(kf_roof),
                               num_sample);
        LatinVector   latin_kb(log10(kb_floor),
                               log10(kb_roof),
                               num_sample);
        Parameter         para_l(num_chem);
        Parameter*        param_l = & para_l;
        DataInterpolation input_output;
        DataInterpolation input_time;
        DataInterpolation input_dissipation;

        for (long si = 0; si != num_sample; ++si) {
            
            
            // Reset parameter to the latin hypercube sampling
            // parameters.
            param_l->setParameter(pow(10, latin_f[si]),
                                  pow(10, latin_b[si]),
                                  pow(10, latin_kf[si]),
                                  pow(10, latin_kb[si]));

            // para_l.printParameter();
            // Write the information of latin hypercube sampling
            // parameters.
 
            StableList  stable_list(reaction_type);
            std::vector<Stable>::iterator nsit  = stable_list.situation.end();
            // record the notSufficient pointer.
            double      input_value     = input_floor;
            double      input_valueroof = input_roof;
            double      input_pace      =
                (input_roof - input_floor) / 200;
            long        whilecount      = 0;
            
            
            do {
                for (int ri = 0; ri != num_reactant; ++ri) {
                    // reset the value of reactant.
                    reactant_ode[ri] = reactant[ri];
                } // for.
                // use a different input.
                reactant_ode[input_ord] = input_value;


                
                ReactantConcentration reaction_l(reaction_type);
                //cout << reaction_l.final_out << "\n";
                // Output parameter and reactant.
                if (!parameter_header) {
                    parameter_outfile << "Sample,"
                                      << "Num,"
                                      << "k,b,kf,kb,atp";
                    for (size_t ni = 0;
                         ni != reaction_l.name.size();
                         ++ni) {
                        parameter_outfile << ","
                                          << reaction_l.name[ni];
                    } // for
                    parameter_outfile << "\n";
                    parameter_header = true;
                } // if no header for parameter_outfile.
                parameter_outfile << si << "," << whilecount << ","
                                  << pow(10, latin_f[si]) << ","
                                  << pow(10, latin_b[si]) << ","
                                  << pow(10, latin_kf[si]) << ","
                                  << pow(10, latin_kb[si]) << ","
                                  << param_l->atp;
                for (size_t rcti = 0;
                     rcti != num_reactant;
                     ++rcti) {
                    parameter_outfile << ","
                                      << reactant_ode[rcti];
                } // for
                parameter_outfile << "\n";
                
                // run ode of reaction.

                reaction_l.odeRun(start_time,
                                  end_time,
                                  param_l,
                                  reactant_ode,
                                  pace_len);
                
                if (!reaction_l.reachStable()) {
                    // if not stable output situation to infofile.
                    infofile << "ODE Not reach Stable: ";
                    infofile << "Sample("
                             << si
                             << ") Num("
                             << whilecount
                             << ") ";
                    infofile << "Parameter("
                             << pow(10,latin_f[si])
                             << ","
                             << pow(10, latin_b[si])
                             << ","
                             << pow(10, latin_kf[si])
                             << ","
                             << pow(10, latin_kb[si])
                             << ")\n";
                    infofile << "Reactant: \n";
                    for (size_t rli = 0;
                         rli != reaction_l.time.size();
                         ++rli) {
                        infofile << reaction_l.time[rli]
                                 << ","
                                 << reaction_l.list[rli]
                                 << ","
                                 << reaction_l.dissipation[rli]
                                 << "\n";
                    } // for
                    infofile << "\n";
                    continue;
                } // if ode result stable.
                
                if (detail == 1) {
                    // if output detail.
                    if (!detail_header) {
                        // Add header.
                        detail_outfile << "Sample,Num,Time,Dissipation";
                        for (size_t ni = 0;
                             ni != reaction_l.name.size();
                             ++ni) {
                            detail_outfile << ","
                                           << reaction_l.name[ni];
                        }//for
                        detail_outfile << "\n";
                        detail_header = true;
                    }//if
                    // output detail data.
                    for (size_t rli = 0;
                         rli != reaction_l.time.size();
                         ++rli) {
                        detail_outfile << si
                                       << ","
                                       << whilecount
                                       << ","
                                       << reaction_l.time[rli]
                                       << ","
                                       << reaction_l.dissipation[rli]
                                       << ","
                                       << reaction_l.list[rli]
                                       << "\n";
                    } // for
                } // if output detail.
                //stable_list.addBack(input_value,
                //                    reaction_l.getFinalOutput(),
                //                    reaction_l.getFinalTime(),
                //                    reaction_l.getFinalDissipation());
                //input_value += input_pace;

                stable_list.addInter(nsit,
                                     input_value,
                                     reaction_l.getFinalOutput(),
                                     reaction_l.getFinalTime(),
                                     reaction_l.getFinalDissipation());
                cout << input_value << "\n";
                nsit = stable_list.unSufficient();
                if (whilecount < 1) {
                    input_value = input_roof;
                } else {
                    if (nsit != stable_list.situation.end()) {
                        // if get unsufficient pointer,
                        // do more calculation.
                        input_value = ((nsit - 1)->input +
                                        nsit->input) / 2;
                    } else if (nsit == stable_list.situation.end()) {
                         //// Stable.
                        if (!stable_list.isStable(1, 2)) {
                            input_valueroof +=
                                (input_roof - input_floor);
                            input_value = input_valueroof;
                            //++whilecount;
                            continue;
                        } // if.
                        break;
                    }
                }// if whilecount.
      
            } while(++whilecount < 1000);

            // sort stable_list.
            stable_list.sortList(1);   

            // Judge fluctuation.
            if (!stable_list.isMonotone()) {
                                    // if not stable output situation to infofile.
                infofile << "Have fluctuation: ";
                infofile << "Sample("
                         << si
                         << ") ";
                infofile << "Parameter("
                         << pow(10,latin_f[si])
                         << ","
                         << pow(10, latin_b[si])
                         << ","
                         << pow(10, latin_kf[si])
                         << ","
                         << pow(10, latin_kb[si])
                         << ")\n";
                infofile << "Input-Output Values: \n";
                for (size_t rli = 0;
                     rli != stable_list.size();
                     ++rli) {
                    infofile << stable_list[rli]
                             << "\n";
                } // for
                infofile << "\n";
                continue;
            }
            
            // Output stable_list.
            if (!outfile_header) {
                outfile << "Sample,Input,Output,Time,Dissipation\n";
                outfile_header = true;
            } // if
            cout << "********************\n";
            for (size_t sli = 0;
                 sli != stable_list.size();
                 ++sli) {
                outfile << si << ",";
                outfile << stable_list[sli];
                cout << stable_list[sli] << "\n";
                outfile << "\n";
            }

            if (!interpolation_header) {
                interpolation_outfile << "Sample,Type,Input,Value\n";
                interpolation_header = true;
            }//if            
            // interpolate data.
            input_output.interpolateL(stable_list.getArray(1),
                                      stable_list.getArray(2),
                                      stable_list.size(),
                                      1000);
            for (int ioi = 0; ioi != input_output.size(); ++ioi) {
                interpolation_outfile << si << ",io,"
                                      << input_output.x[ioi]
                                      << ","
                                      << input_output.y[ioi]
                                      << "\n";
            } //for input_output.
            input_time.interpolateL(stable_list.getArray(1),
                                    stable_list.getArray(3),
                                    stable_list.size(),
                                    1000);
            for (int iti = 0; iti != input_time.size(); ++iti) {
                interpolation_outfile << si << ",it,"
                                      << input_time.x[iti]
                                      << ","
                                      << input_time.y[iti]
                                      << "\n";
            } // for input_time.
            input_dissipation.interpolateL(stable_list.getArray(1),
                                           stable_list.getArray(4),
                                           stable_list.size(),
                                           1000);
            for (int idi = 0; idi != input_dissipation.size(); ++idi) {
                interpolation_outfile << si << ",id,"
                                      << input_dissipation.x[idi]
                                      << ","
                                      << input_dissipation.y[idi]
                                      << "\n";
            } // for input_dissipation.
            //// calculate_outfile.
            if (!calculate_header) {
                calculate_outfile << "Sample,Percent,"
                                  << "Input,Output,"
                                  << "Time,Dissipation\n";
                calculate_header = true;
            }
            int num_9_io = input_output.nearPercentNum(0.9, 2);
            double input_9_io = input_output.x[num_9_io];
            int num_1_io = input_output.nearPercentNum(0.1, 2);
            double input_1_io = input_output.x[num_1_io];
            int num_9_it = input_time.nearValueNum(input_9_io, 1);
            int num_1_it = input_time.nearValueNum(input_1_io, 1);
            int num_9_id = input_dissipation.nearValueNum(input_9_io, 1);
            int num_1_id = input_dissipation.nearValueNum(input_1_io, 1);
            calculate_outfile << si << "," << 0.1 << ","
                              << input_output.x[num_1_io] << ","
                              << input_output.y[num_1_io] << ","
                              << input_time.y[num_1_it] << ","
                              << input_dissipation.y[num_1_id] << "\n";
            calculate_outfile << si << "," << 0.9 << ","
                              << input_output.x[num_9_io] << ","
                              << input_output.y[num_9_io] << ","
                              << input_time.y[num_9_it] << ","
                              << input_dissipation.y[num_9_id] << "\n";

        } // for loop for latin hypercube sampling parameter.
    } // if for sample or not.
    
    ////////////////////
    // End program
    ////////////////////

    delete [] reactant_ode;
    parameter_outfile.close();
    parameter_outfile.clear();
    outfile.close();
    outfile.clear();
    detail_outfile.close();
    detail_outfile.clear();
    interpolation_outfile.close();
    interpolation_outfile.clear();
    calculate_outfile.close();
    calculate_outfile.clear();
    infofile.close();
    infofile.clear();
    cout << "\n"
         << "********************\n"
         << "**    Finished    **\n"
         << "********************\n"
         << endl;
    return 0;
    
}
////////////////////


//bool isStable(std::vector<double>& thevector)
////// Function: isStable
////// Whether reach stable,
////// return true when stable, return false when not stable.
//{
//    std::vector<double>::size_type vct_len = thevector.size();
//    if (vct_len < 6) {
//        return false;
//    }
//    int    big_num = 5;
//    // using how much number to judge.
//    double threshold = 0.0001;    // Here need attention.
//    double tmp = 0.0;
//    double last = thevector[vct_len - 1];
//
//    for (int i = 0; i != big_num; ++i) {
//        tmp = tmp +
//            std::pow(thevector[vct_len - i - 1] -
//                     thevector[vct_len - i - 2],2);
//    }
//    if (last == 0) {
//        return true;
//        // check 0!
//    } else if (last != 0) {
//        tmp = std::sqrt(tmp / big_num);
//        return (tmp / last) < threshold ? true : false;
//    } else {
//        return false;
//    }
//
//}


