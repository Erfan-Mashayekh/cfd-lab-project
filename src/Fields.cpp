#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {

    for (int j = 1; j < grid.jmax(); j++) {
        for (int i = 1; i < grid.imax(); i++) {
            _F(i, j) = _U(i, j) +
                       _dt * (_nu * Discretization::diffusion(_U, i, j) - Discretization::convection_u(_U, _V, i, j));
            _G(i, j) = _V(i, j) +
                       _dt * (_nu * Discretization::diffusion(_V, i, j) - Discretization::convection_v(_U, _V, i, j));
        }
    }
}

void Fields::calculate_rs(Grid &grid) {

    for (int i = 1; i < grid.imax(); i++) {
        for (int j = 1; j < grid.jmax(); j++) {
            _RS(i, j) = (1 / _dt) * ((_F(i + 1, j) - _F(i, j)) / grid.dx() + (_G(i + 1, j) - _G(i, j)) / grid.dy());
        }
    }
}

void Fields::calculate_velocities(Grid &grid) {
    /*
    Explicit euler is used here to discretize the momentum equation, resulting to the equation 7 and 8.

    Equation 7 and Equation 8 give the closed formula to determine the new velocities.
    */

    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();

    // Velocity estimation on all fluid cells excluding right wall and top wall. (Eq 7,8 WS1)
    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            // U (Eq 7)
            _U(i, j) = _F(i, j) - (_dt / dx) * (_P(i + 1, j) - _P(i, j));

            // V (Eq 8)
            _V(i, j) = _G(i, j) - (_dt / dy) * (_P(i, j + 1) - _P(i, j));
        }
    }
}

double Fields::calculate_dt(Grid &grid) {
    
    double dx = grid.dx();
    double dy = grid.dy();
    double Umax = 0.0, Vmax = 0.0;
    _nu = 0.001;
    cout<<"grid dx = "<<dx<<endl;
    cout<<"grid dy = "<<dy<<endl;
    //no problem with dx and dy

    // Find Maximum values of U and V inside their fields: Umax, Vmax
       for (int i=1;i<=grid.imax();i++)
   {
        for (int j=1;j<=grid.jmax();j++)
        {
            if (Umax<abs(_U(i,j)))
            {
                Umax = abs(_U(i,j));
            }
            
            if (Vmax < abs(_V(i,j)))
            {
                Vmax = abs(_V(i,j));
            }
        }
    }

    cout<<"Umax = "<<Umax<<endl;
    cout<<"Vmax = "<<Vmax<<endl;


    double _dt1 = (0.5/_nu)*pow( (1/pow(grid.dx(),2))+ (1/pow(grid.dy(),2)) , -1 );
    std::cout<<"dt1 = "<<_dt1<<std::endl;

    double _dt2 = grid.dx()/Umax;
    std::cout<<"dt2 = "<<_dt2<<std::endl;

    double _dt3 = grid.dy()/Vmax;
    std::cout<<"dt3 = "<<_dt3<<std::endl;

    _dt = std::min( _dt1 , _dt2);
    _dt = std::min( _dt , _dt3);
    _dt = _tau*_dt;

    return _dt;
}



    double &Fields::p(int i, int j) { return _P(i, j); }
    double &Fields::u(int i, int j) { return _U(i, j); }
    double &Fields::v(int i, int j) { return _V(i, j); }
    double &Fields::f(int i, int j) { return _F(i, j); }
    double &Fields::g(int i, int j) { return _G(i, j); }
    double &Fields::rs(int i, int j) { return _RS(i, j); }

    Matrix<double> &Fields::p_matrix() { return _P; }

    double Fields::dt() const { return _dt; }
