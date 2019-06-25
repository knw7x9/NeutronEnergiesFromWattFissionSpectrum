// Main.cpp
// CS 4499
// Written By: Katherine Wilsdon
// 19 April 2019
// Dr. Kerby
// Description - Randomly sampling the amount of energy of 100,000 new neutrons
// created from fission will have using the Watt distribution

#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <random>
#include <array>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/tools/minima.hpp>
namespace tools = boost::math::tools;
using namespace std;
using boost::math::quadrature::trapezoidal;
using boost::math::tools::brent_find_minima;


// Random number generator at seed 10
default_random_engine generator(2000);

// Member function declaration
long unsigned int GeneratesRandomNumber(int);
pair<double,double> RandomPoint();
double Normalize(long);
double CreateXorYCoordinate(double);
bool AcceptPoint(pair<double,double>);
double GetYCoordinate(double);
void CreateSample();

// field variables
vector<pair<double,double>> randomNumberPair;
array<long,100> binArray;
long minimum = 0;
long maximum = pow(2, 24);
pair<double, double> maxima;
pair<double, double> minima;



// Generates a random number with a given size
long unsigned int GeneratesRandomNumber(int size){
  uniform_int_distribution<int> distribution(0, size - 1);
  long unsigned int randomNumber = distribution(generator);
  return randomNumber;
}

// Creates a random point in the range x=[0,8] and y=[0,maxima]
pair<double,double> RandomPoint(){
  pair<double,double> point;
  point.first = CreateXorYCoordinate(8.0);
  point.second = CreateXorYCoordinate(maxima.second);
  return point;
}

// Normalize a number between 0 and 1 using the minimum and maximum in the random numbers vector
double Normalize(long num) {
	return ((double)num - minimum) / (maximum - minimum + 1.0);
}

// Create an x or y coordinate between 0 and the upperBound
double CreateXorYCoordinate(double upperBound){
  return Normalize(GeneratesRandomNumber(maximum))*upperBound;
}

// Accept the point if the y coordinate is less than the y coordinate on the curve, else reject the point
bool AcceptPoint(pair<double,double> point){
  if(point.second < GetYCoordinate(point.first)){
    return true;
  } else {
    return false;
  }
}

// Get the Y coordinate of fx
double GetYCoordinate(double xCoordinate){
  return (0.4865*sinh(sqrt(2*xCoordinate))*exp(-xCoordinate));
}

// Create a random sample of 100,000 energies from the Watt Fission Spectrum
void CreateSample(){
  for (int i = 0; i < 100000; i++){
    pair<double,double> point;
    // create a point, determine of the point is accepted or not, if not accepted repeat
    do {
      point = RandomPoint();
    } while (AcceptPoint(point) == false);
    // Add the accepted point to the vector
    randomNumberPair.push_back(point);
  }
}

// Instantiates a bin array by multiplying the normalized vector x value (0 to 8) by 12.5 (for splitting into 100 energy bins), casting the value to an int for the index
void Bins() {
	for (int i = 0; i < randomNumberPair.size(); i++) {
		binArray[static_cast<int> (randomNumberPair.at(i).first * 12.5)]++;
	}
}


int main() {
// fx is the probability distribution of the energy and gx is 1/fx
auto fx = [](auto E){return (0.4865*sinh(sqrt(2*E))*exp(-E));};
auto gx = [](auto E){return 1/(0.4865*sinh(sqrt(2*E))*exp(-E));}; // 1/fx

// integration using trapezoidal
double I = trapezoidal(fx,0.0,8.0);
cout << "The integration using trapezoidal from [0,8] is : " << fixed << setprecision(2) << I << endl<<endl;

// find maximum of fx (P(E))
maxima = brent_find_minima(gx, 0.0, 8.0, 6);
minima = brent_find_minima(fx, 0.0, 8.0, 6);
cout << "The maxima is at f(" << fixed << setprecision(6) << maxima.first << ") = " << maxima.second << endl;
cout << "The minima is at f(" << minima.first << ") = " << minima.second << endl << endl;


// Create a random sample of 100,000 energies from the Watt Fission Spectrum
CreateSample();
cout << "The size of the sample is: " << randomNumberPair.size() << endl << endl;

// Print out the number of neutron engeries in each bin
Bins();
for (int i = 0; i < binArray.size(); i++){
  cout << "Bin " << i + 1 << ": " << binArray[i] << endl;
}
/*for (int i = 0; i < binArray.size(); i++){
  cout << binArray[i] << endl;
}*/

  return 0;
}
