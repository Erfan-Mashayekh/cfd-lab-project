#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

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

    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();

    for (int i = 0; i < imax + 1; i++) {
        for (int j = 0; j < jmax + 1; j++) {
            _RS(i, j) = (1 / _dt) * ((_F(i + 1, j) - _F(i, j)) / dx + (_G(i + 1, j) - _G(i, j)) / dy);
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

    // Find Maximum values of U and V inside their fields: Umax, Vmax
    for (int j = 1; j < grid.domain().jmax; j++) {
        for (int i = 1; i < grid.domain().imax; i++) {
            if (_U(i, j) > Umax) {
                Umax = _U(i, j);
            }
            if (_V(i, j) > Vmax) {
                Vmax = _V(i, j);
            }
        }
    }

    //std::vector<double> dt_container;
    //dt_container.push_back((dx * dx * dy * dy) / (dx * dx + dy * dy) / (2.0 * _nu));
    //dt_container.push_back(dx / Umax);
    //dt_container.push_back(dy / Vmax);

    //_dt = std::min_element(dt_container.begin(), dt_container.end());    
    _dt = std::min(dx / Umax, dy / Vmax, (dx * dx * dy * dy) / (dx * dx + dy * dy) / (2.0 * _nu));
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
