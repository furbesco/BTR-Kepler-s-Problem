#include <iostream>
#include <fstream>
#include <cmath>


// https://en.wikipedia.org/wiki/Mean_anomaly
// https://en.wikipedia.org/wiki/True_anomaly

// === Define the constant for ease ===
const double Pi = 3.14159265358979323846;
const double c = 3.0e8; //speed of light, m/s
const double G_si = 6.7e-11;
const double M_sun = 2.0e30;

// === Define the Kepler's equation ===

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



// === Main piece of the code ===
int main() {

    // Parameters
    double m1 = 1.0;
    double m2 = 1.0e-7;
    double a = 60e6; // semi-major axis
    double e = 0.2;   // eccentricity 
    double t0 = 0.0;
    double M0 = 0.0; // at periapsis 
    double dt_si = 100;
    double L = 1.0e3; //in meters, for natural units, freedom of choice to do so
    double T_si  = 4.32e6; // tells the loop when to stop

    double T = T_si  * (1 / L) * c;
    double dt = dt_si * (1 / L) * c;
    double Mtot = (m1 * M_sun * G_si / (c*c * L))  + (m2 * M_sun * G_si / (c*c * L));
    double mu = Mtot;
    double n = std::sqrt(mu / (a*a*a));
    std::cout << "n";

    // Output file
    std::ofstream file("kepler_output.csv");

    // Failsafe 
    if (!file.is_open()) {
        std::cout << "Error opening file.\n";
        return 1;
    }

    file << "t,E,M,x,y,r,phi\n";

    // Number of time steps for the loop
    int N = (int)(T / dt);

    for (int i = 0; i<= N; i++) {
        double t = t0 + i*dt;
        // Mean anomaly: used where time evolution is better as it evolves linearly
        double M = M0 + n*(t - t0);
        // Kepler's equation 
        double E = KeplerEq(M, e);
        // Radius
        double r = a*(1.0 - e*std::cos(E));
        // True anomaly: used for orbit simulation as it gives the angle
        double phi = 2.0 * std::atan(
            std::sqrt((1.0 + e)/(1.0 - e)) * std::tan(E / 2.0));
        // Positions
        double x = r * std::cos(phi);
        double y = r * std::sin(phi);

        file << t << "," << E << "," << M << "," << x << "," << y << "," << r << "," << phi << "\n";
    }

    file.close();

    std::cout << "Simulation completed.";

    return 0;
}

