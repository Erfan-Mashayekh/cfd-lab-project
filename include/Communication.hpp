#ifndef COMMUNICATION_HPP
#define COMMUNICATION_HPP

#include <memory>
#include <string>
#include <vector>
#include <mpi.h>

#include "Fields.hpp"

class Communication {

public:
    static void init_parallel(int argn, char** args, int &my_rank, int &comm_size);
    static void barrier();
    static void finalize();
    static void broadcast(void *buffer, int count, MPI_Datatype datatype, int root);
    static void communicate(Matrix<double> &field, const Domain &domain, const int &my_rank) ;
    static void reduce_sum(double &input, double &output);
    static void reduce_min(double &input, double &output);

};

#endif // COMMUNICATION_HPP