////////////////////////////////////////////////////////////////////////////////
// This file is the function file for Latin Hypercube Sampling.
// Author: Wolfson
// Algorithm: shuffle the data for each dimension,
//            and then choose the data from dimensions by their number.
////////////////////////////////////////////////////////////////////////////////

#include "parameter.hpp"

////////////////////
// Class: LatinVector
////////////////////


LatinVector::LatinVector(double lower, double upper, long sample)
//// Constructor: LatinVector
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


void LatinVector::generateRandomList()
//// Member-Function: generateRandomList
//// generate latin hypercube sampling random numbers.
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // set the seed for random shuffle.
    std::mt19937 rd_engine(seed);
    double interval = (upper_limit - lower_limit) / sample_num;


    for (int i = 0; i != sample_num; ++i) {
        double period_lowerlimit = lower_limit + i * interval;
        double period_upperlimit = period_lowerlimit + interval;
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


double& LatinVector::operator[](int i)
//// Member-Function: operator[]
//// Used for subscript operator overloading.
{
    if (i >= sample_num) {
        std::cout << "LatinVector: subscript operator out of bounds."
                  << "\n";
        exit(1);
    }
    return latin_vector[i];
}


void LatinVector::setLatin(double lower, double upper, long sample)
//// Member-Function: setLatin
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


std::vector<double> LatinVector::getLatin()
//// Member-Function: getLatin
//// Return vector of values.
{
    std::vector<double> tmp = latin_vector;
    return tmp;
}


void LatinVector::shuffleLatin()
//// Member-Function: shuffleLatin
//// shuffle again.
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // set the seed for random shuffle.
    std::mt19937 rd_engine(seed);
    // set the seed for random shuffle.
    std::shuffle(latin_vector.begin(),
                 latin_vector.end(),
                 rd_engine);
}


void LatinVector::printLatin()
//// Member-Function: printLatin
//// print vector values to screen.
{
    std::cout << "The Latin Vector Values: "
              << "roof( " << upper_limit
              << " ) floor( " << lower_limit << " )\n";
    std::cout << "Num\tValue\n";
    for (long i = 0; i != sample_num; ++i) {
        std::cout << i << "\t" << latin_vector[i] << "\n";
    }
    std::cout << std::endl;
}


void LatinVector::outputLatin(std::ofstream & outfile)
//// Member-Function: outputLatin
//// output latin vector to file.
{
    outfile << "Num,Value\n";
    for (long i = 0; i != sample_num; ++i) {
        outfile << i << "," << latin_vector[i] << "\n";
    }
}
////////////////////
