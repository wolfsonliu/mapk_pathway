//-*-coding:utf-8-*-
////////////////////////////////////////
// This file is the function file for Latin Hypercube Sampling.
// Author: Wolfson
// Date: Nov.29, 2015
// Modified: Jan.5, 2016
// Algorithm: shuffle the data for each dimension,
//            and then choose the data from dimensions by their number.
////////////////////////////////////////

#include "lhsampling.hpp"

////////////////////
// Functions
////////////////////

std::ostream& operator<<(std::ostream& out, LatinVector& latin)
// {{{ operator<<
//// Output LatinVector.
{
    for (size_t i = 0; i != latin.sample_num - 1; ++i) {
        out << latin.latin_vector[i] << ",";
    }
    out << latin.latin_vector.back();
    return out;
}
// }}}


////////////////////
// Class: LatinVector
////////////////////
// {{{ Class: LatinVector

LatinVector::LatinVector(const double& lower,
                         const double& upper,
                         const size_t& sample)
// {{{ Constructor: LatinVector
{
    if (lower > upper) {
        lower_limit = upper;
        upper_limit = lower;
    } else {
        lower_limit = lower;
        upper_limit = upper;
    }

    sample_num = sample;
    generateRandomList();
}
// }}}


void LatinVector::generateRandomList()
// {{{ Member-Function: generateRandomList
//// generate latin hypercube sampling random numbers.
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // set the seed for random shuffle.
    std::mt19937 rd_engine(seed);
    double interval = (upper_limit - lower_limit) / sample_num;

    double period_lowerlimit = 0;

    for (size_t i = 0; i != sample_num; ++i) {
        period_lowerlimit = lower_limit +
            static_cast<double>(i) * interval;
        latin_vector.push_back(period_lowerlimit +
                               interval *
                               (
                                static_cast<double>(rand()) /
                                static_cast<double>(RAND_MAX)
                                )
                               );
    }
    
    std::shuffle(latin_vector.begin(),
                 latin_vector.end(),
                 rd_engine);
    // shuffle the list randomly to make the one dimensional latin hypercube sampling.

}
// }}}


double& LatinVector::operator[](const size_t& i)
// {{{ Member-Function: operator[]
//// Used for subscript operator overloading.
{
    if (i >= sample_num) {
        std::cout << "LatinVector: subscript operator out of bounds."
                  << "\n";
        exit(1);
    }
    return latin_vector[i];
}
// }}}


void LatinVector::setLatin(const double& lower,
                           const double& upper,
                           const size_t& sample)
// {{{ Member-Function: setLatin
//// Set and generate random values.
{
       if (lower > upper) {
           lower_limit = upper;
           upper_limit = lower;
    } else {
           lower_limit = lower;
           upper_limit = upper;
    }
    sample_num = sample;
    generateRandomList();
}
// }}}


std::vector<double> LatinVector::getLatin()
// {{{ Member-Function: getLatin
//// Return vector of values.
{
    std::vector<double> tmp = latin_vector;
    return tmp;
}
// }}}


void LatinVector::shuffleLatin()
// {{{ Member-Function: shuffleLatin
//// shuffle again.
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // set the seed for random shuffle.
    std::mt19937 rd_engine(seed);
    // set the seed for random shuffle.
    std::shuffle(latin_vector.begin(),
                 latin_vector.end(),
                 rd_engine);
}
// }}}


size_t LatinVector::size()
// {{{ Member-Function: size
//// Return the size of vector.
{
    return latin_vector.size();
}
// }}}
////////////////////
// }}}

/////////
// EOF //
/////////
