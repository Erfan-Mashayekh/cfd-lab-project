<<<<<<< HEAD
// #pragma once

// #include <memory>
// #include <string>
// #include <vector>

// #include "Fields.hpp" 

// class Communication {
//     public:
//         Communication(int rank); //iprod and jproc is the number of process per y or x
//         static void init_parallel(int processers, int _size, int my_rank);
//         static void finalize();
//         static void communicate(Fields &field);
//         static double reduce_min(double value);
//         static double reduce_sum(double value);

// };
=======
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
        int init_parallel(int argn, char** args);
        void finalize();
        void communicate(const Domain &domain, Matrix<double> &field);
        double reduce_sum(const Domain &domain, double &value);
        static double reduce_min(double value);

    private:
        Matrix<Domain> _subdomain;
        int _imax;
        int _jmax;
        int _iproc;
        int _jproc;

};
>>>>>>> af65db8120edabe2cd94bdd403787538b3a76bd1
