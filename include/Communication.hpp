#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Field.hpp" 

class Case {
    private:
        int my_rank;

    public:
        Communication(int rank, double jproc, double iproc); //iprod and jproc is the number of process per y or x
        static void init_parallel(int processers);
        static void finalize();
        static void communicate(Field &field);
        static double reduce_min(double value);

