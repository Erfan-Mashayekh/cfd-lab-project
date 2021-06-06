#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Fields.hpp"

class Communication {

    public:
        Communication() = default;
        Communication(int imax, int jmax, int iproc,  int jproc);
        // Communication(int rank, double jproc, double iproc); //iprod and jproc is the number of process per y or x
        void init_parallel(int argn, char** args, int &my_rank, int &comm_size);
        void finalize();
        void communicate(const Domain domain, Matrix<double> &field, int &my_rank, int &comm_size) ;
        static double reduce_sum(const Domain &domain, double &value);
        static double reduce_min(double value);

    private:
        Matrix<Domain> _subdomain;
        int _imax;
        int _jmax;
        int _iproc;
        int _jproc;

};

