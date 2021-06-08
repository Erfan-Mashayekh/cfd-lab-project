#include "Communication.hpp"
#include <iostream>
#include <mpi.h>


// Initialize communication
void Communication::init_parallel(int argn, char** args, int &my_rank, int &comm_size){

    // start MPI
    MPI_Init(&argn, &args);

    // Fetch number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Fetch process id (rank) in this communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

}

// Barrier method to synchronise
void Communication::barrier() {

    MPI_Barrier(MPI_COMM_WORLD);
}

// Finalize communication
void Communication::finalize() { 
    // Finalize MPI: Wait until all ranks are here to safely exit
    Communication::barrier();

    MPI_Finalize(); 
}

// Broadcast method
void Communication::broadcast(void *buffer, int count, MPI_Datatype datatype, int root){

    MPI_Bcast(buffer, count, datatype, root, MPI_COMM_WORLD);

}

// Send/Receive the data using this communicator
void Communication::communicate(Matrix<double> &field, const Domain &domain, const int &my_rank) {

    int left = (my_rank % domain.iproc == 0) ? MPI_PROC_NULL : my_rank - 1;
    int right = ((my_rank + 1) % domain.iproc == 0) ? MPI_PROC_NULL : my_rank + 1;

    std::vector<double> recv_buf(domain.jmax, 0.0);
    std::vector<double> send_buf = field.get_col(domain.imax-2);

    // Send/Receive Left-right
    MPI_Sendrecv ( send_buf.data(), domain.jmax, MPI_DOUBLE, right, 0,
                   recv_buf.data(), domain.jmax, MPI_DOUBLE,  left, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    field.set_col(recv_buf, 0);
}


// Find and return the minimum value over all ranks
void Communication::reduce_min(double &input, double &output) {

    MPI_Allreduce(&input, &output, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}


// Compute total sum over all ranks
void Communication::reduce_sum(double &input, double &output) {

    MPI_Allreduce(&input, &output, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


