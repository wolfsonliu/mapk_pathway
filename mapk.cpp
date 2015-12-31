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
#include "mapk_ode2.hpp"
#include "parameter.hpp"

using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::vector;
using std::unique_ptr;


struct DeleteDoubleArray { 
    void operator()(double* p) const {
        delete[] p;
    }
};



//// Function: isFileExist
//////// Check file exist.
bool isFileExist(const char* file_name);




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
    vector<string> reactant_names;
    // initiation of reactant.
    double   input_roof;
    double   input_floor;
    double   pace_len;
    unique_ptr<double[]>  f;
    double   sf;
    // forward reaction coefficient for reaction 1.
    double   f_roof;
    double   f_floor;
    unique_ptr<double[]>  b;
    double   sb;
    // backward reaction coefficient for reaction 1.
    double   b_roof;
    double   b_floor;
    unique_ptr<double[]>  kf;
    double   skf;
    // forward reaction coefficient for reaction 2.
    double   kf_roof;
    double   kf_floor;
    unique_ptr<double[]>  kb;
    double   skb;
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
    ostringstream outstr;    
    string   detail_outfile_name;
    ofstream detail_outfile;
    ostringstream detail_outstr;
    string   parameter_outfile_name;
    ofstream parameter_outfile;
    ostringstream parameter_outstr;
    string   interpolation_outfile_name;
    ofstream interpolation_outfile;
    ostringstream interpolation_outstr;
    string   calculate_outfile_name;
    ofstream calculate_outfile;
    ostringstream calculate_outstr;
    string   infofile_name;
    ofstream infofile;
    ostringstream infostr;
    int      detail;
    bool     nosample = false;
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
    
    if (reactant_names.size() == 0) {
        reactantName(reactant_names, reaction_type);
    }   
    ////////////////////
    ostringstream outnamestr;     // output.
    ostringstream ioutnamestr;    // interpolation output.
    ostringstream doutnamestr;    // detail output.
    ostringstream poutnamestr;    // parameter output.
    ostringstream coutnamestr;    // calculate output.
    ostringstream infonamestr;    // infofile name.
    

    // outfile.
    outnamestr << outfile_tag
               <<"mapk_t"
               << reaction_type
               << ".csv";
    outfile_name = outnamestr.str();
    if (isFileExist(outfile_name.c_str())) {
        outfile.open(outfile_name.c_str(),
                     ofstream::out | ofstream::app);
    } else {
        outfile.open(outfile_name.c_str(),
                 ofstream::out | ofstream::app);
        outfile << "Sample,Input,Output,Time,Dissipation\n";
    } // if
    

    // detail outfile.
    doutnamestr << outfile_tag
                <<"mapk_t"
                << reaction_type
                << "_detail.csv";
    detail_outfile_name = doutnamestr.str();
    if (isFileExist(detail_outfile_name.c_str())) {
        detail_outfile.open(detail_outfile_name.c_str(),
                            ofstream::out | ofstream::app);
    } else {
        detail_outfile.open(detail_outfile_name.c_str(),
                            ofstream::out | ofstream::app);
        detail_outfile << "Sample,Num,Time,Dissipation";
        for (size_t ni = 0;
             ni != reactant_names.size(); ++ni) {
            detail_outfile << ","
                           << reactant_names[ni];
        } // for
        detail_outfile << "\n";
    } // if
    

    // interpolation outfile.
    ioutnamestr << outfile_tag
                <<"mapk_t"
                << reaction_type
                << "_interpolation.csv";
    interpolation_outfile_name = ioutnamestr.str();
    if (isFileExist(interpolation_outfile_name.c_str())) {
        interpolation_outfile.open(interpolation_outfile_name.c_str(),
                                   ofstream::out | ofstream::app);
    } else {
        interpolation_outfile.open(interpolation_outfile_name.c_str(),
                                   ofstream::out | ofstream::app);
        interpolation_outfile << "Sample,Type,Input,Value\n";
    } // if

    
    // calculate outfile.
    coutnamestr << outfile_tag
                <<"mapk_t"
                << reaction_type
                << "_calculate.csv";
    calculate_outfile_name = coutnamestr.str();
    if (isFileExist(calculate_outfile_name.c_str())) {
        calculate_outfile.open(calculate_outfile_name.c_str(),
                               ofstream::out | ofstream::app);
    } else {
        calculate_outfile.open(calculate_outfile_name.c_str(),
                               ofstream::out | ofstream::app);
        calculate_outfile << "Sample,Percent,"
                          << "Input,Output,"
                          << "Time,Dissipation\n";
    } // if.

    
    // parameter outfile.
    poutnamestr << outfile_tag
                << "mapk_t"
                << reaction_type
                << "_parameter.csv";
    parameter_outfile_name = poutnamestr.str();
    if (isFileExist(parameter_outfile_name.c_str())) {
        parameter_outfile.open(parameter_outfile_name.c_str(),
                               ofstream::out | ofstream::app);
    } else {
        parameter_outfile.open(parameter_outfile_name.c_str(),
                               ofstream::out | ofstream::app);
        parameter_outfile << "Sample," << "Num," << "k,b,kf,kb,atp";
        for (size_t ni = 0;
             ni != reactant_names.size();
             ++ni) {
            parameter_outfile << ","
                              << reactant_names[ni];
        } // for
        parameter_outfile << "\n";
    } // if

        
    // information outfile
    infonamestr << outfile_tag
                << "mapk_t"
                << reaction_type
                << "_info.txt";
    infofile_name = infonamestr.str();
    infofile.open(infofile_name.c_str(), ofstream::out | ofstream::app);

    ////////////////////


    
    if (num_sample == 0) {
        nosample = true;
        num_sample = 1;
    }
    Parameter         para_l(num_chem);
    Parameter*        param_l = & para_l;

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

    DataInterpolation input_output;
    DataInterpolation input_time;
    DataInterpolation input_dissipation;

    for (long si = 0; si != num_sample; ++si) {
        
        
        // Reset parameter to the latin hypercube sampling
        // parameters.
        if (!nosample) {
            sf  = pow(10, latin_f[si]);
            sb  = pow(10, latin_b[si]);
            skf = pow(10, latin_kf[si]);
            skb = pow(10, latin_kb[si]);
            param_l->setParameter(sf, sb, skf, skb);
        } else {
            sf  = f[0];
            sb  = b[0];
            skf = kf[0];
            skb = kb[0];
            param_l->setParameter(f,b,kf,kb);
        }
        
        // para_l.printParameter();
        // Write the information of latin hypercube sampling
        // parameters.

        StableList  stable_list(reaction_type);
        vector<Stable>::iterator nsit =
            stable_list.situation.end();
        // record the notSufficient pointer.
        double      input_value     = input_floor;
        double      input_valueroof = input_roof;
        double      input_pace      =
            (input_roof - input_floor) / 200;
        int         whilecount      = -1;
        
        
        do {
            ++whilecount;
            for (int ri = 0; ri != num_reactant; ++ri) {
                // reset the value of reactant.
                reactant_ode[ri] = reactant[ri];
            } // for.
            // use a different input.
            reactant_ode[input_ord] = input_value;


            
            ReactantConcentration reaction_l(reaction_type);
            //cout << reaction_l.final_out << "\n";
            // Output parameter and reactant.

            parameter_outstr << si  << "," << whilecount << ","
                             << sf  << "," << sb         << ","
                             << skf << "," << skb        << ","
                             << param_l->atp;
                 
            for (size_t rcti = 0;
                 rcti != num_reactant;
                 ++rcti) {
                parameter_outstr << ","
                                 << reactant_ode[rcti];
            } // for
            parameter_outstr << "\n";

            //cout << si << " latin: "
            //     << sf << ", "
            //     << sb << ", "
            //     << skf << ", "
            //     << skb << "\n";
            
            // run ode of reaction.

            reaction_l.odeRun(start_time,
                              end_time,
                              param_l,
                              reactant_ode,
                              pace_len);
            
            if (!reaction_l.reachStable()) {
                // if not stable output situation to infofile.
                infostr << "ODE Not reach Stable: \n"
                        << "Sample(" << si << "," << whilecount << ") \n"
                        << "Parameter(" << sf << "," << sb << ","
                        << skf << "," << skb << ")\n"
                        << "Reactant: \n";
                for (size_t rli = 0;
                     rli != reaction_l.time.size();
                     ++rli) {
                    infostr << reaction_l.time[rli] << ","
                            << reaction_l.list[rli] << ","
                            << reaction_l.dissipation[rli] << "\n";
                } // for
                infostr << "\n";
                continue;
            } // if ode result stable.
            // Judge oscillation.
            if (reaction_l.isOscillation()) {
                // if not stable output situation to infofile.
                infostr << "May have oscillation: \n"
                        << "Sample(" << si << "," << whilecount << ")\n"
                        << "Parameter(" << sf << "," << sb << "," << skf
                        << "," << skb << ")\n";
                infostr << "Reactant: \n";
                for (size_t rli = 0;
                     rli != reaction_l.time.size();
                     ++rli) {
                    infostr << reaction_l.time[rli] << ","
                            << reaction_l.list[rli] << ","
                            << reaction_l.dissipation[rli] << "\n";
                } // for
                infostr << "\n";
                continue;
            }                
            if (detail == 1) {
                // if output detail.

                // output detail data.
                for (size_t rli = 0;
                     rli != reaction_l.time.size();
                     ++rli) {
                    detail_outstr << si << "," << whilecount << ","
                                  << reaction_l.time[rli] << ","
                                  << reaction_l.dissipation[rli] << ","
                                  << reaction_l.list[rli] << "\n";
                } // for
            } // if output detail.
            
            stable_list.addInter(nsit,
                                 input_value,
                                 reaction_l.getFinalOutput(),
                                 reaction_l.getFinalTime(),
                                 reaction_l.getFinalDissipation());
            //cout << input_value << "\n";
            nsit = stable_list.unSufficient();
            if (whilecount < 1) {
                input_value = input_roof;
            } else {
                if (nsit != stable_list.situation.end()) {
                    // if get unsufficient pointer,
                    // do more calculation.
                    input_value = std::pow(10,(std::log10((nsit - 1)->input) +
                                               std::log10(nsit->input)) / 2);
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
  
        } while(whilecount < 1000);
        
        detail_outfile << detail_outstr.str(); // output ostringstream.
        detail_outstr.str(string()); // clear the ostringstream.
        
        // sort stable_list.
        stable_list.sortList(1);   
        
        // Output stable_list.

        //cout << "********************\n";
        for (size_t sli = 0;
             sli != stable_list.size();
             ++sli) {
            outstr << si << "," << stable_list[sli] << "\n";
            //cout << stable_list[sli] << "\n";
        }

       
        // interpolate data.
        input_output.interpolateL(stable_list.getArray(1),
                                  stable_list.getArray(2),
                                  stable_list.size(),
                                  1000);
        for (int ioi = 0; ioi != input_output.size(); ++ioi) {
            interpolation_outstr << si << ",io,"
                                 << input_output.x[ioi] << ","
                                 << input_output.y[ioi] << "\n";
        } //for input_output.
        input_time.interpolateL(stable_list.getArray(1),
                                stable_list.getArray(3),
                                stable_list.size(),
                                1000);
        for (int iti = 0; iti != input_time.size(); ++iti) {
            interpolation_outstr << si << ",it,"
                                 << input_time.x[iti] << ","
                                 << input_time.y[iti] << "\n";
        } // for input_time.
        input_dissipation.interpolateL(stable_list.getArray(1),
                                       stable_list.getArray(4),
                                       stable_list.size(),
                                       1000);
        for (int idi = 0; idi != input_dissipation.size(); ++idi) {
            interpolation_outstr << si << ",id,"
                                 << input_dissipation.x[idi] << ","
                                 << input_dissipation.y[idi] << "\n";
        } // for input_dissipation.
        //// calculate_outfile.

        int num_9_io = input_output.nearPercentNum(0.9, 2);
        double input_9_io = input_output.x[num_9_io];
        int num_1_io = input_output.nearPercentNum(0.1, 2);
        double input_1_io = input_output.x[num_1_io];
        int num_9_it = input_time.nearValueNum(input_9_io, 1);
        int num_1_it = input_time.nearValueNum(input_1_io, 1);
        int num_9_id = input_dissipation.nearValueNum(input_9_io, 1);
        int num_1_id = input_dissipation.nearValueNum(input_1_io, 1);
        calculate_outstr << si << "," << 0.1 << ","
                         << input_output.x[num_1_io] << ","
                         << input_output.y[num_1_io] << ","
                         << input_time.y[num_1_it] << ","
                         << input_dissipation.y[num_1_id] << "\n"
                         << si << "," << 0.9 << ","
                         << input_output.x[num_9_io] << ","
                         << input_output.y[num_9_io] << ","
                         << input_time.y[num_9_it] << ","
                         << input_dissipation.y[num_9_id] << "\n";

    } // for loop for latin hypercube sampling parameter.
    
    ////////////////////
    // End program
    ////////////////////

    
    outfile               << outstr.str();

    calculate_outfile     << calculate_outstr.str();
    parameter_outfile     << parameter_outstr.str();
    interpolation_outfile << interpolation_outstr.str();
    infofile              << infostr.str();
    
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

bool isFileExist(const char* file_name)
{
    ifstream in_file(file_name);
    return in_file.good();
}



