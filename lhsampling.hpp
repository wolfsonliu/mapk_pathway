////////////////////////////////////////
// This file is the head file for Latin Hypercube Sampling.
// Author: Wolfson
// Date: Nov.29, 2015
// Edit: Jan.5, 2016
// Algorithm: shuffle the data for each dimension,
//            and then choose the data from dimensions by their number.
////////////////////////////////////////
#ifndef LATIN_HYPERCUBE_SAMPLING_PARAMETER
#define LATIN_HYPERCUBE_SAMPLING_PARAMETER

#include <algorithm>
#include <random>    // std::default_random_engine
#include <chrono>    // std::chrono::system_clock
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <array>
#include <string>
#include <vector>
#include <memory>


class LatinVector;

////////////////////
// Functions
////////////////////
std::ostream& operator<<(std::ostream& out, LatinVector& latin);

////////////////////
// Classes
////////////////////

class LatinVector
// {{{ Class: LatinVector
//// get the list of random shuffled values generated by random number.
{
    double upper_limit;
    double lower_limit;
    size_t  sample_num;
    std::vector<double> latin_vector;
    
    //// generateRandomList:
    //// Generate latin hypercube sampling random numbers.
    void generateRandomList();

public:
    // Constructor
    LatinVector():
        upper_limit(0.0),
        lower_limit(0.0),
        sample_num(0.0) {}//
    LatinVector(const double& lower,
                const double& upper,
                const size_t& sample);
    // Destructor
    virtual ~LatinVector() {}
    // Friend
    friend std::ostream& operator<<(std::ostream& out,
                                    LatinVector& latin);
    // Operator
    //// operator[]:
    //// Used for subscript operator overloading.
    double& operator[](const size_t& i);
    
    // Member-Function
    //// setLatin:
    //// Set and generate random values.
    void setLatin(const double& lower,
                  const double& upper,
                  const size_t& sample);
    //// getLatin:
    //// Return vector of values.
    std::vector<double> getLatin();
    //// shuffleLatin:
    //// Shuffle again.
    void shuffleLatin();
    //// size:
    //// Return the size of vector.
    size_t size();
    //// printLatin:
    //// Print vector values to screen.
    //void printLatin();
    //// outputLatin:
    //// Output latin vector to file.
    //void outputLatin(std::ofstream& outfile);
};
// }}}


#endif // LATIN_HYPERCUBE_SAMPLING_PARAMETER

/////////
// EOF //
/////////