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

    std::cout << "Value of nu = " << _nu << std::endl;
}

void Fields::calculate_fluxes(Grid &grid) {

    for (int j = 1; j < grid.jmax() + 1; j++) {
        for (int i = 1; i < grid.imax(); i++) {
            _F(i, j) = _U(i, j) +
                       _dt * (_nu * Discretization::diffusion(_U, i, j) - Discretization::convection_u(_U, _V, i, j));
        }
    }

    for (int j = 1; j < grid.jmax(); j++) {
        for (int i = 1; i < grid.imax() + 1; i++) {
            _G(i, j) = _V(i, j) +
                       _dt * (_nu * Discretization::diffusion(_V, i, j) - Discretization::convection_v(_U, _V, i, j));
        }
    }

    /*
    //Ghost cell flux BC


      for (int i=1;i<grid.imax();i++){
        for (int j=1;j<grid.jmax();j++){
            // Eq 17 WS1 ghost cell flux BC
            _G(i,0) = _V(i,0);
            _G(i,grid.jmax()) = _V(i,grid.jmax());
            // Eq 17 WS1 ghost cell flux BC
            _F(0,j) = _U(0,j);
            _F(grid.imax(),j) = _U(grid.imax(),j);
        }  
      }

    */
}

void Fields::calculate_rs(Grid &grid) {

    for (int i = 1; i < grid.imax() + 1; i++) {
        for (int j = 1; j < grid.jmax() + 1; j++) {
            _RS(i,j) = 1/_dt * ((_F(i, j) - _F(i-1 ,j)) * 1/grid.dx()+ (_G(i, j) - _G(i, j-1)) * 1/grid.dy());
        }
    }
}

void Fields::calculate_velocities(Grid &grid) {
    /*
    Explicit euler is used here to discretize the momentum equation, resulting to the equation 7 and 8.
    Equation 7 and Equation 8 give the closed formula to determine the new velocities.
    */

    double dx = grid.dx();
    double dy = grid.dy();

    // Velocity estimation on all fluid cells excluding right wall and top wall
    for (int i = 1; i < grid.imax(); i++) {
        for (int j = 1; j < grid.jmax() + 1; j++) {
            // U (Eq 7)
            _U(i, j) = _F(i, j) - (_dt / dx) * (_P(i + 1, j) - _P(i, j));
        }
    }

    for (int i = 1; i < grid.imax() + 1; i++) {
        for (int j = 1; j < grid.jmax(); j++) {
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
    // cout<<"grid dx = "<<dx<<endl;
    // cout<<"grid dy = "<<dy<<endl;
    //no problem with dx and dy

    // Find Maximum values of U and V inside their fields: Umax, Vmax
    for (int j = 1; j < grid.jmax(); j++) {
        for (int i = 1; i < grid.imax(); i++) {
            if (abs(_U(i, j)) > Umax) {
                Umax = abs(_U(i, j));
            }
            if (abs(_V(i, j)) > Vmax) {
                Vmax = abs(_V(i, j));
            }
        }
    }
/*
<<<<<<< HEAD
    _dt = 10;
    _nu = 0.001;
    std::vector<double> dt_container = {(dx * dx * dy * dy) / (dx * dx + dy * dy) / (2.0 * _nu) , dx / Umax , dy / Vmax};
    for (int k = 0; k <= 2; k++){
        if (dt_container[k] < _dt){
            _dt = dt_container[k];
        }
    }
    //_dt = std::min_element(dt_container.begin(), dt_container.end())[0];    
    //_dt = std::min(dx / Umax, dy / Vmax, (dx * dx * dy * dy) / (dx * dx + dy * dy) / (2.0 * _nu));
======= 
*/
    std::cout<<"Umax = "<<Umax<<std::endl;



    double _dt1 = (0.5/_nu)*pow( (1/pow(grid.dx(),2))+ (1/pow(grid.dy(),2)) , -1 );
    // std::cout<<"dt1 = "<<_dt1<<std::endl;

    double _dt2 = grid.dx()/Umax;
    // std::cout<<"dt2 = "<<_dt2<<std::endl;

    double _dt3 = grid.dy()/Vmax;
    // std::cout<<"dt3 = "<<_dt3<<std::endl;

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
