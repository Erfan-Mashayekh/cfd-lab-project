#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI, double TI, double alpha, double beta, double GX, double GY)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _GX(GX), _GY(GY)  {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);
    _T = Matrix<double>(imax + 2, jmax + 2, TI);
    _T_new = Matrix<double>(imax + 2, jmax + 2, TI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);

}

// Calculate Fn and Gn
void Fields::calculate_fluxes(Grid &grid, bool energy_eq){

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if(currentCell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID){
            _F(i, j) = _U(i, j) +
                       _dt * (_nu * Discretization::diffusion(_U, i, j) - Discretization::convection_u(_U, _V, i, j)) 
                       - int(energy_eq) * _beta * _dt * Discretization::interpolate(_T, i, j, 1, 0) * _GX;
        }
        if(currentCell->neighbour(border_position::TOP)->type() == cell_type::FLUID){
            _G(i, j) = _V(i, j) +
                        _dt * (_nu * Discretization::diffusion(_V, i, j) - Discretization::convection_v(_U, _V, i, j))
                        - int(energy_eq) * _beta * _dt * Discretization::interpolate(_T, i, j, 0, 1) * _GY;
        } 
    }      
}

// calculate right hand side rhs of the pressure eq
void Fields::calculate_rs(Grid &grid) {

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        _RS(i, j) = (1 / _dt) * ((_F(i, j) - _F(i - 1, j)) * (1 / grid.dx()) + (_G(i, j) - _G(i, j - 1)) * (1 / grid.dy()));
    }
}

// Calculate the velocities at the next time step
void Fields::calculate_velocities(Grid &grid) {
    /*
    Explicit euler is used here to discretize the momentum equation, resulting to the equation 7 and 8.
    Equation 7 and Equation 8 give the closed formula to determine the new velocities.
    */

    double dx = grid.dx();
    double dy = grid.dy();

    // Velocity estimation on all fluid cells excluding right wall and top wall
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        _U(i, j) = _F(i, j) - (_dt / dx) * (_P(i + 1, j) - _P(i, j));
        _V(i, j) = _G(i, j) - (_dt / dy) * (_P(i, j + 1) - _P(i, j));
    }
}


// Calculate the velocities at the next time step
void Fields::calculate_temperature(Grid &grid) {
    _T_new = _T;
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        _T_new(i, j) = _T(i, j) + 
                       _dt * (_alpha * Discretization::diffusion(_T, i, j) - Discretization::convection_T(_T, _U, _V, i, j));
    }
    _T = _T_new;
}

static bool abs_compare(const double& a, const double& b)
{
    return (std::abs(a) < std::abs(b));
}

// calculate dt for adaptive time stepping
double Fields::calculate_dt(Grid &grid, bool energy_eq) {

    double dx = grid.dx();
    double dy = grid.dy();

    // Find Maximum values of U and V inside their fields: Umax, Vmax
    double Umax = 0, Vmax = 0;
    Umax = std::abs(*std::max_element(_U.container()->begin(), _U.container()->end(), abs_compare));
    Vmax = std::abs(*std::max_element(_V.container()->begin(), _V.container()->end(), abs_compare));

    // Comparing 3 dt for Courant-Friedrichs-Levi (CFL) conditions in order to ensure stability
    // and avoid oscillations
    std::vector<double> dt;
    dt.push_back((dx * dx * dy * dy) / (dx * dx + dy * dy) / (2.0 * _nu));
    dt.push_back(grid.dx()/Umax);
    dt.push_back(grid.dy()/Vmax);
    if(energy_eq){
        dt.push_back((dx * dx * dy * dy) / (dx * dx + dy * dy) / (2.0 * _alpha));
    }

    _dt = _tau * (*std::min_element(dt.begin(), dt.end()));
            
    return _dt;
}

double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::T(int i, int j) { return _T(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }
Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::T_matrix() { return _T; }
Matrix<double> &Fields::f_matrix() { return _F; }
Matrix<double> &Fields::g_matrix() { return _G; }
Matrix<double> &Fields::rs_matrix() { return _RS;}


double Fields::dt() const { return _dt; }
