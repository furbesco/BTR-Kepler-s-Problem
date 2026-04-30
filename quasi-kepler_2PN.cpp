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


// === Structure the parameters so that they are stored together for ease of code later on, from chatgpt ===
struct elems2PN {
    double ar;
    double er;
    double et;
    double ephi;
    double n;
    double Phi;
    double t0;
    double phi0;
    double f4t;
    double g4t;
    double f4phi;
    double g4phi;
};


// === Define the Kepler's and 1PN equation ===

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

// changed to atan2 to get full rotations, otherwise stuck bw pi/2, same as for 1PN
// https://stackoverflow.com/questions/283406/what-is-the-difference-between-atan-and-atan2-in-c
double TA1PN(double u, double ephi) {
    return 2.0 * std::atan2(
        std::sqrt(1.0 + ephi) * std::sin(u / 2.0), 
        std::sqrt(1.0 - ephi) * std::cos(u / 2.0)
    );
}

// === 2PN approximation, formula taken from Memmesheimer ===
double twoPNeq(double l, const elems2PN& elems) {
    double u = l;  // initial guess

    for (int i = 0; i < 60; i++) {
        double v = TA1PN(u, elems.ephi);
        double F = u - elems.et * std::sin(u) + (elems.g4t / std::pow(c, 4)) * (v - u) + (elems.f4t / std::pow(c, 4)) * std::sin(v) - l;

        // numerical derivative, suggested by chatGPT as much easier (I got confused how to do this part analytically)
        double h = 1e-10;
        double up = u + h;
        double vp = TA1PN(up, elems.ephi);
        double Fp = up - elems.et * std::sin(up)+ (elems.g4t / std::pow(c, 4)) * (vp - up) + (elems.f4t / std::pow(c, 4)) * std::sin(vp)- l;
        double dF = (Fp - F) / h;
        double du = -F / dF;
        u += du;

        if (std::abs(du) < 1e-12) break;
    }
    return u;
}


// === Main piece of the code ===
int main() {

    // Parameters
    double m1 = 1.0;
    double m2 = 1.0e-7;
    double a = 60e6; // semi-major axis
    double e = 0.4;   // eccentricity 
    double t0 = 0.0;
    double M0 = 0.0; // at periapsis 
    double dt_si = 500;
    double L = 1.0e3; //in meters, for natural units, freedom of choice to do so
    double T_si  = 4.32e9; // tells the loop when to stop
    double phi0 = 0.0;
    double k = 0.05; // for orbital procession, to have the "swirly" shape

    //double T = T_si  * (1.0 / L) * c;
    //double dt = dt_si * (1.0 / L) * c;

    double Mtot = (m1 * M_sun * G_si / (c*c * L))  + (m2 * M_sun * G_si / (c*c * L));
    double mu = Mtot;
    double n_newton = std::sqrt(mu / (a*a*a));
    std::cout << "n";

    // Initial guess by setting the 2PN values close to the Keplerian

    elems2PN elems;
    elems.ar = a;
    elems.er = e;
    elems.et = 0.8;
    elems.ephi = 0.6;
    elems.Phi = 2.0 * Pi * (1.0 + k);
    elems.n = n_newton;
    elems.t0 = t0;
    elems.phi0 = phi0;
    // look into why functions below only work normally for 10^33, crashes at other values
    elems.f4t = 2e33;
    elems.g4t = 2e33;
    elems.f4phi = 2e33;
    elems.g4phi = 2e33;

    double P = 2.0 * Pi / elems.n; //one orbit
    int Norbits = 6;
    double T = Norbits * P;
    double dt = P / 1000.0;

    // Output file
    std::ofstream file("2PN_output.csv");

    // Failsafe 
    if (!file.is_open()) {
        std::cout << "Error opening file.\n";
        return 1;
    }

    file << "t,l,u,R,v,phi,x,y\n";

   // Number of time steps for the loop
    int N = (int)(T / dt);

    for (int i = 0; i<= N; i++) {
        double t = elems.t0 + i * dt;
        // Mean anomaly: used where time evolution is better as it evolves linearly
        double l = elems.n * (t-elems.t0);
        // completed periods (to overcome the orbit crash)
        int Nrad = static_cast<int>(std::floor(l/(2.0 * Pi)));
        // reduced l, suggested by chatgpt to overcome orbital folding at discontinuities in u
        double l_red = std::fmod(l, 2.0 * Pi);
        if (l_red < 0.0) l_red += 2.0 * Pi;
        // 2PN equation 
        double u = twoPNeq(l_red, elems);
        // Radius
        double R = elems.ar*(1.0 - elems.er*std::cos(u));
        // True anomaly: used for orbit simulation as it gives the angle
        double v = TA1PN(u, elems.ephi);
        //Ortbital movement shift
        double phi_local = v + (elems.f4phi/ std::pow(c, 4)) * std::sin(2.0 * v) + (elems.g4phi/ std::pow(c, 4)) * std::sin(3.0 * v);
        // Full azimuth angle with the 2PN, similar to the 1PN but with no longer simple elliptical orbit
        double phi = elems.phi0 + Nrad * elems.Phi + (elems.Phi / (2.0 * Pi)) * phi_local;
        // Positions
        double x = R * std::cos(phi);
        double y = R * std::sin(phi);

        file << t << "," << l << "," << u << "," << R << ","
             << v << "," << phi << "," << x << "," << y << "\n";
    }

    file.close();

    std::cout << "Simulation completed.";

    return 0;
}