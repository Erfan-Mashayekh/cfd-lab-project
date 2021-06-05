#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Fields.hpp"

class Communication {

    public:
        Communication() = default;
        Communication(int iproc,  int jproc);
        // Communication(int rank, double jproc, double iproc); //iprod and jproc is the number of process per y or x
        int init_parallel(int argn, char** args);
        void finalize();
        void communicate(Fields &field);
        static double reduce_min(double value);

    private:
        Matrix<Domain> _subdomain;
        int _iproc;
        int _jproc;

};