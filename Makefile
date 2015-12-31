all:
	g++  *.cpp -Wall -std=c++11 -o mapk -lm -lgsl -lgslcblas -lboost_program_options
