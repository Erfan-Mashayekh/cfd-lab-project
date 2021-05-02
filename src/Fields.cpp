#include "Fields.hpp"

#include <algorithm>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {}

void Fields::calculate_rs(Grid &grid) {}

void Fields::calculate_velocities(Grid &grid) {
    /*
    Explicit euler is used here to discretize the momentum equation, resulting to the equation 7 and 8.

    Equation 7 and Equation 8 give the closed formula to determine the new velocities.
    */

    // Velocity estimation on all fluid cells excluding right wall and top wall. (Eq 7,8 WS1)
    for (int i=1;i<=imax;i++)
    {
        for (int j=1;j<=jmax;j++)
        {
            // U (Eq 7)
            _U(i,j) = _F(i,j) - (_dt/_dx)*(_P(i+1,j) - _P(i,j)) ;

            // V (Eq 8)
            _V(i,j) = _G(i,j) - (_dt/_dy)*(_P(i,j+1) - _P(i,j));
        }
    }

   

double Fields::calculate_dt(Grid &grid) { return _dt; }

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }