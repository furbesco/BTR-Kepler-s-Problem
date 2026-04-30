#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <map>
#include <cstdlib>
#include <ctime>

// https://en.wikipedia.org/wiki/Mean_anomaly
// https://en.wikipedia.org/wiki/True_anomaly

// === Define the constant for ease ===
const double Pi = 3.14159265358979323846;
const double c = 3.0e8; //speed of light, m/s
const double G_si = 6.7e-11;
const double M_sun = 2.0e30;

// == Configuration file ===
struct Config {
    double m1, m2;
    double a, e;

    double t0, phi0;
    double dt;
    int Norbits;

    int PN_order;

    double L;
    double k;
};

Config load_config(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Error: could not open config file\n";
        exit(1);
    }

    std::string line;
    std::map<std::string, double> values;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string key;
        double value;

        if (std::getline(ss, key, '=') && ss >> value) {
            values[key] = value;
        }
    }

    Config cfg;

    cfg.m1 = values["m1"];
    cfg.m2 = values["m2"];
    cfg.a = values["a"];
    cfg.e = values["e"];

    cfg.t0 = values["t0"];
    cfg.phi0 = values["phi0"];
    cfg.dt = values["dt"];
    cfg.Norbits = (int)values["Norbits"];

    cfg.PN_order = (int)values["PN_order"];
    cfg.L = values["L"];
    cfg.k = values["k"];

    return cfg;
}

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

// now for the 1PN approximation

double onePNeq(double M, double e){
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
    Config cfg = load_config("params.cfg");

    double m1 = cfg.m1;
    double m2 = cfg.m2;
    double a  = cfg.a;
    double e  = cfg.e;

    if (a <= 0 || e < 0 || e >= 1) {
        std::cout << "Invalid parameters\n";
        return 1;
    }

    // orbital frequency
    double Mtot = (m1 * M_sun * G_si / (c*c * cfg.L))+ (m2 * M_sun * G_si / (c*c * cfg.L));

    double mu = Mtot;
    double n= std::sqrt(mu / (a*a*a));

    double P = 2.0 * Pi / n;
    double T = cfg.Norbits * P;

    double dt = (cfg.dt > 0.0) ? cfg.dt : P / 1000.0;

    std::cout << "n = " << n << "\n";

    // Output file
    std::ofstream file("kepler_output.csv");

    // Failsafe 
    if (!file.is_open()) {
        std::cout << "Error opening file.\n";
        return 1;
    }

    file << "# ===== Simulation Parameters =====\n";
    file << "# m1=" << cfg.m1 << "\n";
    file << "# m2=" << cfg.m2<< "\n";
    file << "# a=" << cfg.a << "\n";
    file << "# e=" << cfg.e << "\n";
    file << "# t0=" << cfg.t0 << "\n";
    file << "# phi0=" << cfg.phi0 << "\n";
    file << "# dt="<< dt << "\n";
    file << "# Norbits=" << cfg.Norbits << "\n";

    file << "# L=" << cfg.L << "\n";

    file << "# ================================================\n\n";

    file << "t,E,M,x,y,r,phi\n";

    // Number of time steps for the loop
    int N = (int)(T / dt);

    for (int i = 0; i<= N; i++) {
        double t = cfg.t0 + i*dt;
        // Mean anomaly: used where time evolution is better as it evolves linearly
        double M = n*(t - cfg.t0);
        int Nrad = (int)(M / (2.0 * Pi));
        double M_red = std::fmod(M, 2.0 * Pi);
        if (M_red < 0.0) M_red += 2.0 * Pi;
        // Kepler's equation 
        double E = KeplerEq(M_red, e);
        // Radius
        double r = a*(1.0 - e*std::cos(E));
        // True anomaly: used for orbit simulation as it gives the angle
        double phi_r = 2.0 * std::atan2(
            std::sqrt(1.0 + e) * std::sin(E / 2.0),
            std::sqrt(1.0 - e) * std::cos(E / 2.0)
        );
        if (phi_r < 0.0) {
            phi_r += 2.0 * Pi;
        }
        double phi = phi_r + 2.0 * Pi * Nrad;
        // Positions
        double x = r * std::cos(phi_r);
        double y = r * std::sin(phi_r);

        file << t << "," << E << "," << M << "," << x << "," << y << "," << r << "," << phi << "\n";
    }

    file.close();

    std::cout << "Simulation completed.";

    return 0;
}