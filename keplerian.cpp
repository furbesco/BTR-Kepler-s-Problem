#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <array>
#include <vector> 

// def the constant for ease
const double Pi = 3.14159265358979323846;


// https://en.wikipedia.org/wiki/Mean_anomaly
// === def the Kepler's equation ===

double KeplerEq(double M, double e){
    double E = M; //initial setting, just to start

    for (int i = 0; i < 50; i++){
        double f = E - e*std::sin(E) - M;
        double fp = 1.0 - e*std::cos(E);
        double dE  = - f/fp;

        E += dE;

        if (std::abs(dE) < 1e-10)
            break; //convergence at this point and no need to continue
    }
    return E;
}



// === main piece of the code ===
int main() {

    // Parameters
    double G = 1.0;
    double m1 = 10.0;
    double m2 = 20.0;
    double a = 300.0; // semi-major axis
    double e = 0.5;   // eccentricity 
    double t0 = 0.0;
    double M0 = 0.0; // at periapsis 
    double dt = 1.0;
    double T  = 500000.0; // tells the loop when to stop

    double Mtot = m1 + m2;
    double mu = G*Mtot;
    double n = std::sqrt(mu / (a*a*a));

    // Output file to store the computation
    std::ofstream file("kepler_output.csv");

    //failsafe for errors
    if (!file.is_open()) {
        std::cout << "Error opening file.\n";
        return 1;
    }

    file << "t, x, y, r, phi\n";
    
}

