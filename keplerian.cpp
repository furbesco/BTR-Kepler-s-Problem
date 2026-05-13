#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <map>
#include <cstdlib>
#include <ctime>
#include <iomanip>

// https://en.wikipedia.org/wiki/Mean_anomaly
// https://en.wikipedia.org/wiki/True_anomaly

// === Define the constant for ease ===
const double Pi = 3.14159265358979323846;
const double c = 3.0e8; //speed of light, m/s
const double G_si = 6.7e-11;
const double M_sun = 2.0e30;

// === Structure the parameters so that they are stored together for ease of code later on, from chatgpt ===
struct elemsKepler {
    double ar;     // radial semi-major axis
    double e;      // eccentricity
    double n;      // mean motion
    double Phi;    // azimuthal angle over one radial period
    double t0;
    double phi0;
    double d;
};

// == Configuration file ===
struct Config {
    double m1, m2;
    double a, e;
    double t0, phi0;
    double dt;
    int Norbits;
    int PN_order;
    double T_years;
    double T_months;
    double d;
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
    cfg.T_years = values["T_years"];
    cfg.T_months = values["T_months"];
    cfg.PN_order = (int)values["PN_order"];
    cfg.d = values["d"];

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

struct Vec2 {
    double x, y;
};

struct State {
    Vec2 r1, r2;
    Vec2 v1, v2;
    Vec2 a1, a2;
};

struct SymTensor2d {
    double xx, xy, yy;
};

SymTensor2d second_mass_moment(const State& state, double m1, double m2) {
    return {
        m1*state.r1.x*state.r1.x + m2*state.r2.x*state.r2.x,
        m1*state.r1.x*state.r1.y + m2*state.r2.x*state.r2.y,
        m1*state.r1.y*state.r1.y + m2*state.r2.y*state.r2.y
    };
}

SymTensor2d second_mass_moment_dot(const State& state, double m1, double m2) {
    return {
        2.0*(m1*state.r1.x*state.v1.x + m2*state.r2.x*state.v2.x),
        m1*(state.r1.x*state.v1.y + state.r1.y*state.v1.x)
        + m2*(state.r2.x*state.v2.y + state.r2.y*state.v2.x),
        2.0*(m1*state.r1.y*state.v1.y + m2*state.r2.y*state.v2.y)
    };
}

SymTensor2d second_mass_moment_ddot(const State& state, double m1, double m2) {
    return {
        2.0*(m1*(state.v1.x*state.v1.x + state.r1.x*state.a1.x)
        + m2*(state.v2.x*state.v2.x + state.r2.x*state.a2.x)),
        m1*(2.0*state.v1.x*state.v1.y + state.r1.x*state.a1.y + state.r1.y*state.a1.x)
        + m2*(2.0*state.v2.x*state.v2.y + state.r2.x*state.a2.y + state.r2.y*state.a2.x),
        2.0*(m1*(state.v1.y*state.v1.y + state.r1.y*state.a1.y)
        + m2*(state.v2.y*state.v2.y + state.r2.y*state.a2.y))
    };
}

// === Main piece of the code ===
int main() {

    // Parameters
    Config cfg = load_config("params.cfg");

    std::cout << "m1 = " << cfg.m1 << "\n";
    std::cout << "m2 = " << cfg.m2 << "\n";
    std::cout << "a = " << cfg.a << "\n";
    std::cout << "e = " << cfg.e << "\n";
    std::cout << "dt = " << cfg.dt << "\n";

    // For the unit conversion
    double a = cfg.a;
    double e = cfg.e;

    double m1_input = cfg.m1;
    double m2_input = cfg.m2;
    double m1_si = m1_input * M_sun;
    double m2_si = m2_input * M_sun;
    double Mtot_si = m1_si + m2_si;
    double L = G_si * Mtot_si / (c*c);

    // geometric mass fractions
    double m1 = m1_input / (m1_input + m2_input);
    double m2 = m2_input / (m1_input + m2_input);
    double mu = 1.0;

    // Chirp masses 
    double Mc = std::pow(m1_si * m2_si, 3.0/5.0) / std::pow(m1_si + m2_si, 1.0/5.0);
    double distance = cfg.d;

    double n_newton = std::sqrt(mu / (a*a*a));
    double epsilon = mu/a;

    elemsKepler elems;
    elems.ar = a;
    elems.e = cfg.e;
    double k_real = 3.0 * epsilon / (1.0 -e*e);
    double k_ex = k_real; //* 2000000.0;
    //elems.Phi = 2.0 * Pi * (1.0 + cfg.k);
    elems.Phi = 2.0 * Pi * (1.0 + k_ex);
    elems.n = n_newton;
    elems.t0 = cfg.t0;
    elems.phi0 = cfg.phi0;


    // Simulation in time or orbits
    double P = 2.0 * Pi / elems.n;
    double T;
    double dt;

    // Fixed number of orbits

    if (cfg.Norbits > 0) {
        T = cfg.Norbits * P;
        std::cout << "Simulation mode: orbits\n";
        if (cfg.dt > 0.0) {
            dt = cfg.dt;
        } else {
            dt = P / 1000.0;
        }
    }

    // Fixed number of years
    else if (cfg.T_years > 0.0) {
        double T_seconds =
            cfg.T_years * 365.25 * 24.0 * 3600.0;
        T = T_seconds * c / L;
        std::cout << "Simulation mode: years\n";
        if (cfg.dt > 0.0) {
            dt = cfg.dt;
        } else {
            dt = P / 100.0;
        }
    }

    // Simulation in months
    else if (cfg.T_months > 0.0) {
        double T_seconds =
            cfg.T_months * 30.44 * 24.0 * 3600.0;
        T = T_seconds * c / L;
        std::cout << "Simulation mode: months\n";
        if (cfg.dt > 0.0) {
            dt = cfg.dt;
        } else {
            dt = P / 100.0;
        }
    }

    else {
        std::cerr
            << "Error: specify either Norbits, "
            << "T_years, or T_months\n";
        exit(1);
    }

    // here time is dimentionless

    double total_orbits = T / P;
    double P_seconds = P * L / c; // seconds
 
    std::cout << "Total simulated orbits = "<< total_orbits << "\n";
    std::cout << "n = " << n_newton << "\n";
    std::cout << "Time per orbit is = " << P_seconds << "\n";

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
    file << "# ================================================\n\n";
    file << "t,E,M,x,y,r,v,phi,"
        << "Ixx,Ixy,Iyy,"
        << "Id_xx,Id_xy,Id_yy,"
        << "Idd_xx,Idd_xy,Idd_yy,"
        << "vx,vy\n";

    // Number of time steps for the loop
    int N = (int)(T / dt);

    for (int i = 0; i<= N; i++) {
        double t_dimentionless = cfg.t0 + i*dt;
        double t = t_dimentionless * L / c; //conversion to seconds
        // Mean anomaly: used where time evolution is better as it evolves linearly
        double M = elems.n*(t_dimentionless - cfg.t0);
        int Nrad = (int)(M / (2.0 * Pi));
        double M_red = std::fmod(M, 2.0 * Pi);
        if (M_red < 0.0) M_red += 2.0 * Pi;
        // Kepler's equation 
        double E = KeplerEq(M_red, e);
        // Radius
        double r = a*(1.0 - e*std::cos(E));
        // True anomaly: used for orbit simulation as it gives the angle
        double v = 2.0 * std::atan2(
            std::sqrt(1.0 + e) * std::sin(E / 2.0),
            std::sqrt(1.0 - e) * std::cos(E / 2.0)
        );
        if (v < 0.0) {
            v += 2.0 * Pi;
        }
        double phi = v + 2.0 * Pi * Nrad;
        // Positions
        double x = r * std::cos(v);
        double y = r * std::sin(v);

        // next point for velocity
        double t_n = t_dimentionless + dt;
        double M_n = elems.n*(t_n - cfg.t0);
        int Nrad_n = (int)(M_n / (2.0 * Pi));
        double M_red_n = std::fmod(M_n, 2.0 * Pi);
        if (M_red_n < 0.0) M_red_n += 2.0 * Pi;

        double E_n = KeplerEq(M_red_n, e);
        double r_n = a*(1.0 - e*std::cos(E_n));

        double v_n = 2.0 * std::atan2(
            std::sqrt(1.0 + e) * std::sin(E_n / 2.0),
            std::sqrt(1.0 - e) * std::cos(E_n / 2.0)
        );
        if (v_n < 0.0) v_n += 2.0 * Pi;

        double x_n = r_n * std::cos(v_n);
        double y_n = r_n * std::sin(v_n);


        // second next point for acceleration
        double t_nn = t_dimentionless + 2.0*dt;
        double M_nn = elems.n*(t_nn - cfg.t0);
        int Nrad_nn = (int)(M_nn / (2.0 * Pi));
        double M_red_nn = std::fmod(M_nn, 2.0 * Pi);
        if (M_red_nn < 0.0) M_red_nn += 2.0 * Pi;

        double E_nn = KeplerEq(M_red_nn, e);
        double r_nn = a*(1.0 - e*std::cos(E_nn));

        double v_nn = 2.0 * std::atan2(
            std::sqrt(1.0 + e) * std::sin(E_nn / 2.0),
            std::sqrt(1.0 - e) * std::cos(E_nn / 2.0)
        );
        if (v_nn < 0.0) v_nn += 2.0 * Pi;

        double x_nn = r_nn * std::cos(v_nn);
        double y_nn = r_nn * std::sin(v_nn);


        // velocities and accelerations
        double vx = (x_n - x) / dt;
        double vy = (y_n - y) / dt;

        double ax = (x_nn - 2.0*x_n + x) / (dt*dt);
        double ay = (y_nn - 2.0*y_n + y) / (dt*dt);


        // two-body state
        double m1_kg = m1 * M_sun;
        double m2_kg = m2 * M_sun;
        double Mtot_phys = m1_kg + m2_kg;

        State state;

        state.r1 = { (m2_kg/Mtot_phys)*x,  (m2_kg/Mtot_phys)*y };
        state.r2 = {-(m1_kg/Mtot_phys)*x, -(m1_kg/Mtot_phys)*y };

        state.v1 = { (m2_kg/Mtot_phys)*vx,  (m2_kg/Mtot_phys)*vy };
        state.v2 = {-(m1_kg/Mtot_phys)*vx, -(m1_kg/Mtot_phys)*vy };

        state.a1 = { (m2_kg/Mtot_phys)*ax,  (m2_kg/Mtot_phys)*ay };
        state.a2 = {-(m1_kg/Mtot_phys)*ax, -(m1_kg/Mtot_phys)*ay };

        // second mass moment
        SymTensor2d I = second_mass_moment(state, m1_kg, m2_kg);
        SymTensor2d Id = second_mass_moment_dot(state, m1_kg, m2_kg);
        SymTensor2d Idd = second_mass_moment_ddot(state, m1_kg, m2_kg);

        file << t << "," << E << "," << M << ","
            << x << "," << y << "," << r << "," << v << "," << phi << ","
            << I.xx << "," << I.xy << "," << I.yy << ","
            << Id.xx << "," << Id.xy << "," << Id.yy << ","
            << Idd.xx << "," << Idd.xy << "," << Idd.yy << ","
            << vx << "," << vy << "\n";
        }

    file.close();

    std::cout << "Simulation completed.";

    return 0;
}