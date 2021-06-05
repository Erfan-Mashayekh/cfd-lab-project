#include "Communication.hpp"
#include <iostream>
#include <mpi.h>

Communication::Communication {}

// Initialize communication
static void Communication::init_parallel(int processers, int _size, int my_rank) {
    MPI_Init(_ra) MPI_Comm_size(MPI_COMM_WORLD, _size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
}

// Finalize communication
static void Communication::finalize() { MPI_Finalize(); }

/*Single call to this should send/receive the respective held
    across all parallel boundaries, without expecting any MPI-specic options in its argu-
    ments.
    */
static void Communication::communicate(Field field) {
    MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
                 recvcount, MPI_Datatype recvtype, source, recvtag, MPI_Comm comm, MPI_Status *status)
}

// Find and return the minimum value over all ranks
static double Communication::reduce_min(double value) {
    double reduction_result = 0;
    MPI_Reduce(&my_rank, &reduction_result, 1, MPI_INT, MPI_MIN, root_rank, MPI_COMM_WORLD);
    return reduction_result;
}
// Compute total sum over all ranks
static double Communication::reduce_sum(double value) {
    double sum = 0;
    MPI_Reduce(&my_rank, &sum, 1, MPI_INT, MPI_MIN, root_rank, MPI_COMM_WORLD);
    return sum;
}
