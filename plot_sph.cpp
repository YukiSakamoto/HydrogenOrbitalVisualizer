#include <cmath>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/format.hpp>
#include <array>
#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>

// reference: 
//  http://www.chem.ous.ac.jp/~waka/compchem/hydrogen-like_atom/hy-5.html

// Constants
const double pi = 3.14159265358979323846;

typedef std::array<double,3> position_type;

position_type spherical2cartesian(const double r, const double theta, const double phi)
{
    position_type ret;
    ret[0] = r * std::sin(phi) * std::cos(theta);
    ret[1] = r * std::sin(phi) * std::sin(theta);
    ret[2] = r * std::cos(phi);
    return ret;
}

double
radial_func_1s(const double r)
{
    // R_1s = 2*(Z/a0)^1.5 e^(-Z*r/a0)
    //  where   Z = 1. and a0= 1.  under atomic unit
    double Z = 1.0;
    double a0 = 1.0;
    double r1s = 2.0 * std::pow(Z/a0, 1.5) * std::exp(-r);
    return r1s;
}

double
radial_func_2s(const double r)
{
    // R_2s = R_2,0 = (1/(2sqrt(2))) * (Z/a0)^1.5 * (2-Z*r/a0) *exp(-Z*r/a0 / 2)
    double Z = 1.0;
    double a0 = 1.0;
    double rho = Z*r/a0;
    double r2s = 1/(2.*std::sqrt(2.)) * std::pow(Z/a0, 1.5) * (2.-rho)*std::exp(-rho/2.);
    return r2s;
}

double radial_func_2p(const double r)
{
    double Z = 1.;
    double a0 = 1.;
    double rho = Z * r/a0;
    double r2p = 1./(2.*std::sqrt(6)) * std::pow(Z/a0, 1.5) * rho * std::exp(-rho/2.);
    return r2p;
}

double radial_func(const int n, const int l, const double r)
{
    //std::cout << boost::format("n: %d l: %d") % n % l << std::endl;
    if (n <= l) throw;
    if (n == 1) {
        if (l == 0) {
            return radial_func_1s(r);
        } else {
            throw;
        }
    } else if (n == 2) {
        if (l == 0) {
            return radial_func_2s(r);
        } else if (l == 1) {
            return radial_func_2p(r);
        } else {
            throw;
        }
    } else {
        throw;
    }
    // never get here
}

double
generate_angle(std::vector<double> &v, const double factor = 1.0)
{
    size_t num_grid = v.size();
    double dPhi = factor * pi / num_grid;
    for(int i = 0; i < v.size(); i++) {
        v[i] = dPhi * i;
    }
    return dPhi;
}

double
generate_radial(std::vector<double> &r, const double factor = 1.0)
{
    size_t num_grid = r.size();
    double dr = factor * 1./num_grid;
    for(int i=0; i < r.size(); i++) {
        r[i] = i*dr;
    }
    return dr;
}


void 
print_spherical_coordinate(std::vector<double> const &theta, std::vector<double> const &phi)
{
    for(int i = 0; i < theta.size(); i++) {
        for(int j = 0; j < phi.size(); j++) {
            position_type pos = spherical2cartesian(1., theta[i], phi[j]);
            std::cout << boost::format("%f\t%f\t%f") % pos[0] % pos[1] % pos[2] << std::endl;
        }
    }
    return;
}

void plot_radial_function(const int n, const int l, std::vector<double> const &r)
{

    FILE *gp = popen("gnuplot", "w");
    fprintf(gp, "set ylabel \"|R(r)\"\n");
    fprintf(gp, "plot \"-\" with linespoints\n");
    for(int i = 0; i < r.size(); i++) {
        fprintf(gp, "%f\t%f\n", r[i], radial_func(n, l, r[i]) );
    }
    fclose(gp);
}

void plot_radial_distribution(const int n, const int l, std::vector<double> const &r) 
{
    
    FILE *gp = popen("gnuplot", "w");
    fprintf(gp, "set ylabel \"|r*Psi(r)|^2\"\n");
    fprintf(gp, "plot \"-\" with linespoints\n");
    std::vector<double> R(r.size());
    for(int i = 0; i < r.size()-1; i++) {
        double p = std::pow(r[i], 2) * std::pow(radial_func(n, l, r[i]), 2);
        fprintf(gp, "%f\t%f\n", r[i], p );
    }
    fclose(gp);
}

void calculate_radial_func(const int n, const int l, 
        std::vector<double> const &r, std::vector<double> &R)
{
    for(int i = 0; i < r.size(); i++) {
        R[i] = radial_func(n,l,r[i]);
    }
}


void 
do_plot(const int n, const int l, const int m,
        const size_t num_grid_theta, 
        const size_t num_grid_phi, 
        const size_t num_grid_radial, 
        const double contour, const int plot_type = 1)
{
    // allocate variables
    std::vector<double> phi(num_grid_phi);
    std::vector<double> theta(num_grid_theta);
    std::vector<double> r(num_grid_radial);

    // initialize grid
    double d_phi   = generate_angle(phi, 2.0);
    double d_theta = generate_angle(theta, 1.0);
    double d_r     = generate_radial(r, 20);

    if (plot_type == 1) {
        plot_radial_function(n,l,r);
        return;
    } else if (plot_type == 2) {
        plot_radial_distribution(n, l, r);
        return;
    } else {
        std::vector<double> R_r(num_grid_radial);
        calculate_radial_func(n, l, r, R_r);

        for(int i = 0; i < theta.size(); i++) {
            for(int j = 0; j < phi.size(); j++) {
                std::complex<double> Ylm = boost::math::spherical_harmonic(l, m, theta[i], phi[j]);

                for(int k = 0; k < r.size(); k++) {
                    std::complex<double> psi = R_r[k] * Ylm;

                    double rho = std::norm(psi) * std::pow(r[k], 2.) * d_r * std::sin(theta[i]) * d_theta * d_phi;

                    if (std::abs(rho-contour) < 0.001) {
                        position_type pos_cart = spherical2cartesian(r[k], theta[i], phi[j]);
                        if (0. < psi.real() ) {
                            std::cout << boost::format("%f\t%f\t%f\t 1") % pos_cart[0] % pos_cart[1] % pos_cart[2] << std::endl;
                        } else {
                            std::cout << boost::format("%f\t%f\t%f\t-1") % pos_cart[0] % pos_cart[1] % pos_cart[2] << std::endl;
                        }
                    }
                }
            }
        }
    }
    return;
}

int main(int argc, const char **argv)
{
    // parameters
    const size_t num_grid = 30;
    const size_t num_grid_r = 30;

    const int n = 2;
    const int l = 1;
    const int m = 0;
    const double contour = 1.e-5;

    int plot_type = 3;
    do_plot(n, l, m, num_grid, num_grid, num_grid_r, contour, plot_type);

    // test
    //print_spherical_coordinate(theta,phi);
    
    // calculation

    return 0;
}
