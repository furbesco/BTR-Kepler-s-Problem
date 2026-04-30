#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <map>
#include <cstdlib>


// https://en.wikipedia.org/wiki/Mean_anomaly
// https://en.wikipedia.org/wiki/True_anomaly

// === Define the constant for ease ===
const double Pi = 3.14159265358979323846;
const double c = 3.0e8; //speed of light, m/s
const double G_si = 6.7e-11;
const double M_sun = 2.0e30;


// === Structure the parameters so that they are stored together for ease of code later on, from chatgpt ===
struct elems1PN {
    double ar;     // radial semi-major axis
    double er;     // radial eccentricity
    double et;     // time eccentricity
    double ephi;   // angular eccentricity
    double n;      // mean motion
    double Phi;    // azimuthal angle over one radial period
    double t0;
    double phi0;
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

// now for the 1PN approximation

double onePNeq(double l, double et) {
    double u = l;  // initial guess

    for (int i = 0; i < 50; i++) {
        double f  = u - et * std::sin(u) - l;
        double df = 1.0 - et * std::cos(u);
        double du = -f / df;

        u += du;

        if (std::abs(du) < 1e-12)
            break;
    }
    return u;
}

// changed to atan2 to get full roratition, otherwise stuck bw pi/2 
// https://stackoverflow.com/questions/283406/what-is-the-difference-between-atan-and-atan2-in-c
double TA1PN(double u, double ephi) {
    return 2.0 * std::atan2(
        std::sqrt(1.0 + ephi) * std::sin(u / 2.0), 
        std::sqrt(1.0 - ephi) * std::cos(u / 2.0)
    );
}

struct Vec2 {
    double x, y;
};

struct State {
    Vec2 r1, r2;
    Vec2 v1, v2;
    Vec2 a1, a2;
};

struct Config {
    double m1, m2;
    double a, e;

    double er, et, ephi;

    double t0, phi0;
    double dt;
    int Norbits;

    int PN_order;

    double L;
    double k;
};

// === Parameter configuration ===
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
    cfg.er = values["er"];
    cfg.et = values["et"];
    cfg.ephi = values["ephi"];
    cfg.t0 = values["t0"];
    cfg.phi0 = values["phi0"];
    cfg.dt = values["dt"];
    cfg.Norbits = (int)values["Norbits"];
    cfg.PN_order = (int)values["PN_order"];
    cfg.L = values["L"];
    cfg.k = values["k"];

    return cfg;
}

// === Second mass moment ===

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
        2.0*(m1*(state.v1.x*state.v1.x + state.r1.x*state.a1.x) +
             m2*(state.v2.x*state.v2.x + state.r2.x*state.a2.x)),

        m1*(2.0*state.v1.x*state.v1.y +
            state.r1.x*state.a1.y +
            state.r1.y*state.a1.x)
      + m2*(2.0*state.v2.x*state.v2.y +
            state.r2.x*state.a2.y +
            state.r2.y*state.a2.x),

        2.0*(m1*(state.v1.y*state.v1.y + state.r1.y*state.a1.y) +
             m2*(state.v2.y*state.v2.y + state.r2.y*state.a2.y))
    };
}



// === Main piece of the code ===
int main() {
    
    Config cfg = load_config("params.cfg");

    std::cout << "m1 = " << cfg.m1 << "\n";
    std::cout << "m2 = " << cfg.m2 << "\n";
    std::cout << "a = " << cfg.a << "\n";
    std::cout << "er = " << cfg.er << "\n";
    std::cout << "et = " << cfg.et << "\n";
    std::cout << "ephi = " << cfg.ephi << "\n";
    std::cout << "L = " << cfg.L << "\n";
    std::cout << "dt = " << cfg.dt << "\n";

    double m1 = cfg.m1;
    double m2 = cfg.m2;
    double a = cfg.a;
    double e = cfg.e;

    double Mtot = (m1 * M_sun * G_si / (c*c * cfg.L)) + (m2 * M_sun * G_si / (c*c * cfg.L));
    double mu = Mtot;
    double n_newton = std::sqrt(mu / (a*a*a));

    elems1PN elems;
    elems.ar = a;
    elems.er = cfg.er;
    elems.et = cfg.et;
    elems.ephi = cfg.ephi;
    elems.Phi = 2.0 * Pi * (1.0 + cfg.k);
    elems.n = n_newton;
    elems.t0 = cfg.t0;
    elems.phi0 = cfg.phi0;

    double P = 2.0 * Pi / elems.n;
    int Norbits = cfg.Norbits;
    double T = Norbits * P;

    double dt;
    if (cfg.dt > 0.0) {
        dt = cfg.dt;
    } else {
        dt = P / 1000.0;
    }

    //double T = T_si  * (1.0 / L) * c;
    //double dt = dt_si * (1.0 / L) * c;
    std::cout << "n = " << n_newton << "\n";
    
    // Output file
    std::ofstream file("1PN_output.csv");

    file << "# ===== Simulation Parameters =====\n";
    file << "# m1=" << cfg.m1 << "\n";
    file << "# m2=" << cfg.m2 << "\n";
    file << "# a=" << cfg.a << "\n";
    file << "# e=" <<cfg.e  << "\n";
    file << "# er=" << cfg.er<< "\n";
    file << "# et=" << cfg.et << "\n";
    file << "# ephi="<< cfg.ephi << "\n";
    file << "# t0=" << cfg.t0 << "\n";
    file << "# phi0=" << cfg.phi0 << "\n";
    file << "# dt=" <<dt << "\n";
    file << "# Norbits=" << cfg.Norbits << "\n";
    file << "# PN_order=" << cfg.PN_order << "\n";
    file << "# L=" << cfg.L << "\n";
    file << "# k=" << cfg.k << "\n";
    file << "# =======================================================\n\n";

    // Failsafe 
    if (!file.is_open()) {
        std::cout << "Error opening file.\n";
        return 1;
    }

    file << "t,l,u,R,v,phi,x,y,"
         << "Ixx,Ixy,Iyy,"
         << "Id_xx,Id_xy,Id_yy,"
         << "Idd_xx,Idd_xy,Idd_yy\n";

   // Number of time steps for the loop
    int N = (int)(T / dt);

    for (int i = 0; i<= N; i++) {
        double t = elems.t0 + i * dt;
        // Mean anomaly: used where time evolution is better as it evolves linearly
        double l = elems.n * (t-elems.t0);
        // completed periods (to oveercome thwe orbit crash)
        int Nrad = static_cast<int>(std::floor(l/(2.0 * Pi)));
        // reduced l, suggested by chatgpt to overcome orbital folding at discontinuities in u
        double l_red = std::fmod(l, 2.0 * Pi);
        if (l_red < 0.0) l_red += 2.0 * Pi;
        // 1PN equation 
        double u = onePNeq(l_red, elems.et);
        // Radius
        double R = elems.ar*(1.0 - elems.er*std::cos(u));
        // True anomaly: used for orbit simulation as it gives the angle
        double v = TA1PN(u, elems.ephi);
        //Ortbital movement shift
        double phi = elems.phi0 + Nrad * elems.Phi + v * (elems.Phi/(2.0*Pi));
        // Positions
        double x = R * std::cos(phi);
        double y = R * std::sin(phi);
        // for the second mass moment, the diff in the velocity
        double t_n = t+dt;
        double l_n = elems.n * (t_n - elems.t0);
        double l_red_n = std::fmod(l_n, 2.0 * Pi);
        if (l_red_n < 0.0) l_red_n += 2.0 * Pi;

        double u_n = onePNeq(l_red_n, elems.et);
        double R_n = elems.ar*(1.0 - elems.er*std::cos(u_n));
        double v_n = TA1PN(u_n, elems.ephi);
        int Nrad_n = static_cast<int>(std::floor(l_n/(2.0 * Pi)));
        double phi_n = elems.phi0 + Nrad_n * elems.Phi + v_n * (elems.Phi/(2.0*Pi));

        double x_n = R_n * std::cos(phi_n);
        double y_n = R_n * std::sin(phi_n);

        double t_nn = t + 2.0*dt;
        double l_nn = elems.n * (t_nn - elems.t0);
        double l_red_nn = std::fmod(l_nn, 2.0 * Pi);
        if (l_red_nn < 0.0) l_red_nn += 2.0 * Pi;

        double u_nn = onePNeq(l_red_nn, elems.et);
        double R_nn = elems.ar*(1.0 - elems.er*std::cos(u_nn));
        double v_nn = TA1PN(u_nn, elems.ephi);
        int Nrad_nn = static_cast<int>(std::floor(l_nn/(2.0 * Pi)));
        double phi_nn = elems.phi0 + Nrad_nn * elems.Phi + v_nn * (elems.Phi/(2.0*Pi));

        double x_nn = R_nn * std::cos(phi_nn);
        double y_nn = R_nn * std::sin(phi_nn);

        // velocities
        double vx = (x_n - x)/dt;
        double vy = (y_n - y)/dt;

        //accelerations
        double ax = (x_nn - 2.0*x_n + x) / (dt*dt);
        double ay = (y_nn - 2.0*y_n + y) / (dt*dt);

        // two body state
        double Mtot_phys = m1 + m2;

        State state;

        state.r1 = {(m2/Mtot_phys)*x, (m2/Mtot_phys)*y };
        state.r2 = {-(m1/Mtot_phys)*x, -(m1/Mtot_phys)*y };

        state.v1 = {(m2/Mtot_phys)*vx, (m2/Mtot_phys)*vy };
        state.v2 = {-(m1/Mtot_phys)*vx, -(m1/Mtot_phys)*vy };

        state.a1 = {(m2/Mtot_phys)*ax, (m2/Mtot_phys)*ay };
        state.a2 = {-(m1/Mtot_phys)*ax, -(m1/Mtot_phys)*ay };

        // smm
        SymTensor2d I = second_mass_moment(state, m1, m2);
        SymTensor2d Id = second_mass_moment_dot(state, m1, m2);
        SymTensor2d Idd = second_mass_moment_ddot(state, m1, m2);


        file << t << "," << l << "," << u << "," << R << ","
             << v << "," << phi << "," << x << "," << y << "," 
             << I.xx << "," << I.xy << "," << I.yy << ","
             << Id.xx << "," << Id.xy << "," << Id.yy << ","
             << Idd.xx << "," << Idd.xy << "," << Idd.yy << "\n";
    }

    file.close();

    std::cout << "Simulation completed.";

    return 0;
}