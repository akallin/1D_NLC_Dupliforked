//c++ class for creating a general Hamiltonian in c++ vectors
//Roger Melko, November 2007
#ifndef GenHam_H
#define GenHam_H

#include <iostream>
#include <vector>
using namespace std;

#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)

typedef long double h_float;  //precision for Hamiltonian storage

class GENHAM{

public:
  int Vdim; //dimenson of reduced Hilbert space
  
  vector<vector<long> > PosHam;
  vector<vector<h_float> > ValHam;
  //vector<double> DiagHam;
  
  vector<long> Basis;
  vector<long> BasPos;
  
  Array<double,2> Ham;  //full hamiltonian
  
  GENHAM(const int N_ ,const h_float J_, const h_float h_, vector < pair<int,int> > BBond_, bool Low_); 
  void printg();
  //double at(const int , const int );
  Array<double,1> apply(const Array<double,1>&);

  
  void SparseHamJQ();
private:
  int Nsite; //number sites
  bool LowField; //High or Low Field expansion
  
  vector< pair < int,int> > Bond;

    h_float JJ; //heisenberg exchange value
    h_float hh; //next-nearest neighbor exchange value
  
  double HdiagPart(const long, int);
    double HOFFdBondX(const int, const long);
    double HOFFdBondY(const int, const long);

};

#endif
