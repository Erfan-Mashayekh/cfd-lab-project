#include "Communication.hpp"
#include <iostream>
#include <mpi.h>

Communication::Communication(int iproc,  int jproc)
              :_iproc(iproc), _jproc(jproc){}

// Initialize communication
int Communication::init_parallel(int argn, char** args){

    int size, rank;

    // start MPI
    MPI_Init(&argn, &args);

    // Fetch number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Fetch process id (rank) in this communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    return rank;
}

// Finalize communication
void Communication::finalize() { MPI_Finalize(); }

// Send/Receive the data using this communicator
void Communication::communicate(Fields &field) {

    int MPI_Sendrecv(const void *sendbuf,     // The address of the send buffer.
                    int sendcount,            // Number of elements in the send buffer.
                    MPI_Datatype sendtype,    // Data type of the send buffer.
                    int dest,                 // Destination process number (rank number).
                    int sendtag,              // Sent message tag.
                    void *recvbuf,            // The memory address of the receiving buffer.
                    int recvcount,            // Number of elements in the receive buffer.
                    MPI_Datatype recvtype,    // Source process rank number.
                    int source,               // Source process rank number.
                    int recvtag,              // Receive message tag.
                    MPI_Comm comm,            // MPI communicator.
                    MPI_Status *status );     // Status of object. 
}

// // Find and return the minimum value over all ranks
// static double Communication::reduce_min(double value) {
//     double reduction_result = 0;
//     MPI_Reduce(&my_rank, &reduction_result, 1, MPI_INT, MPI_MIN, root_rank, MPI_COMM_WORLD);
//     return reduction_result;
// }

// // Compute total sum over all ranks
// static double Communication::reduce_sum(double value) {
//     double sum = 0;
//     MPI_Reduce(&my_rank, &sum, 1, MPI_INT, MPI_MIN, root_rank, MPI_COMM_WORLD);
//     return sum;
// }
