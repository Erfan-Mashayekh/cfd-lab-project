#include "Discretization.hpp"

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double result = (     interpolate(U, i, j, 1, 0)  * interpolate(U, i, j, 1, 0)    -     interpolate(U, i - 1, j, 0, 0)   * interpolate(U, i - 1, j, 0, 0)) / _dx +
                    ( abs(interpolate(U, i, j, 1, 0)) * (U(i, j) - U(i + 1, j)) / 2.0 - abs(interpolate(U, i - 1, j, 0, 0))  * (U(i - 1, j) - U(i, j)) / 2.0 ) / _dx * _gamma +
                    (     interpolate(V, i, j, 1, 0)  * interpolate(U, i, j, 0, 1)    -     interpolate(V, i, j - 1, 1, -1)  * interpolate(U, i, j - 1, 0, 0)) / _dy +
                    ( abs(interpolate(U, i, j, 1, 0)) * (U(i, j) - U(i, j + 1)) / 2.0 - abs(interpolate(V, i, j - 1, 1, -1)) * (U(i, j - 1) - U(i, j)) / 2.0 ) / _dy * _gamma;
    return result;                
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double result = (     interpolate(V, i, j, 0, 1)  * interpolate(V, i, j, 0, 1)    -     interpolate(U, i, j - 1, 0, 0)   * interpolate(V, i, j - 1, 0, 0)) / _dy +
                    ( abs(interpolate(V, i, j, 0, 1)) * (V(i, j) - V(i, j + 1)) / 2.0 - abs(interpolate(V, i, j - 1, 0, 0))  * (V(i, j - 1) - V(i, j)) / 2.0 ) / _dy * _gamma +
                    (     interpolate(U, i, j, 0, 1)  * interpolate(V, i, j, 1, 0)    -     interpolate(U, i - 1, j, -1, 1)  * interpolate(V, i - 1, j, 0, 0)) / _dx +
                    ( abs(interpolate(U, i, j, 0, 1)) * (V(i, j) - V(i + 1, j)) / 2.0 - abs(interpolate(U, i - 1, j, -1, 1)) * (V(i - 1, j) - V(i, j)) / 2.0 ) / _dx * _gamma;
    return result;     
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
    double result = (A(i + 1, j) - 2.0 * A(i, j) + A(i - 1, j)) / (_dx * _dx) +
                    (A(i, j + 1) - 2.0 * A(i, j) + A(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    double result = (A(i, j) + A(i + i_offset, j + j_offset)) / 2.0;
    return result;
}
