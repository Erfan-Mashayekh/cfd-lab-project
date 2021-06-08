#include "Communication.hpp"
#include <iostream>
#include <mpi.h>
#include <cmath>

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
void Communication::communicate(Matrix<double> &field, const Domain &domain, const int &my_rank, const int &shift) {

    int left   =       (my_rank % domain.iproc == 0) ? MPI_PROC_NULL : my_rank - 1;
    int right  = ((my_rank + 1) % domain.iproc == 0) ? MPI_PROC_NULL : my_rank + 1;



 // ------------------------Send Right---------------------------------
    std::vector<double> recv_l_buf(domain.jmax, 0.0);
    std::vector<double> send_r_buf = field.get_col(domain.imax - 2);

    MPI_Sendrecv ( send_r_buf.data(), domain.jmax, MPI_DOUBLE, right, 0,
                   recv_l_buf.data(), domain.jmax, MPI_DOUBLE,  left, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if( left != MPI_PROC_NULL ){
        field.set_col(recv_l_buf, domain.imin);
    }

    MPI_Barrier(MPI_COMM_WORLD);



 // ----------------------- Send Left ------------------------------
    std::vector<double> recv_r_buf(domain.jmax, 0.0);
    std::vector<double> send_l_buf = field.get_col(domain.imin + 1);

   
    MPI_Sendrecv ( send_l_buf.data(), domain.jmax, MPI_DOUBLE,  left, 0,
                   recv_r_buf.data(), domain.jmax, MPI_DOUBLE, right, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if( right != MPI_PROC_NULL ){
        field.set_col(recv_r_buf, domain.imax - 1 - shift);
    }

    MPI_Barrier(MPI_COMM_WORLD);

/*******************************************************************
 *******************************************************************/


    int up   =               (std::floor(my_rank / domain.iproc) == 0) ? MPI_PROC_NULL : my_rank - domain.iproc;
    int down = (std::floor(my_rank / domain.iproc) == domain.jproc -1) ? MPI_PROC_NULL : my_rank + domain.iproc;


  // ----------------------- Send down ------------------------------
    std::vector<double> recv_u_buf(domain.imax, 0.0);
    std::vector<double> send_d_buf = field.get_row(domain.jmax - 2);

    MPI_Sendrecv ( send_d_buf.data(), domain.imax, MPI_DOUBLE, down, 0,
                   recv_u_buf.data(), domain.imax, MPI_DOUBLE,   up, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if( up != MPI_PROC_NULL ){
        field.set_row(recv_u_buf, domain.jmin);
    }

    MPI_Barrier(MPI_COMM_WORLD);



  // ----------------------- Send up ------------------------------
    std::vector<double> recv_d_buf(domain.imax, 0.0);
    std::vector<double> send_u_buf = field.get_row(domain.jmin + 1);

    MPI_Sendrecv ( send_u_buf.data(), domain.imax, MPI_DOUBLE,   up, 0,
                   recv_d_buf.data(), domain.imax, MPI_DOUBLE, down, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if( down != MPI_PROC_NULL ){
        field.set_row(recv_d_buf, domain.jmax - 1 - shift);
    }

    MPI_Barrier(MPI_COMM_WORLD);

}


// Find and return the minimum value over all ranks
void Communication::reduce_min(double &input, double &output) {

    MPI_Allreduce(&input, &output, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}


// Compute total sum over all ranks
void Communication::reduce_sum(double &input, double &output) {

    MPI_Allreduce(&input, &output, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


