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
    static void communicate(Matrix<double> &field, const int &imax, const int &jmax, 
                            const int &iproc, const int &jproc, const int &my_rank) ;
    static double reduce_sum(double &value);
    static double reduce_min(double &value);

};

#endif // COMMUNICATION_HPP