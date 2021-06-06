#include "Communication.hpp"
#include <iostream>
#include <mpi.h>

Communication::Communication(int imax, int jmax, int iproc,  int jproc)
              :_imax(imax), _jmax(jmax), _iproc(iproc), _jproc(jproc){}

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
void Communication::communicate(const Domain &domain, Matrix<double> &field) {
    
    // Send to right subdomain
    if (domain.col < _iproc){

        // sender
        const void* sendbuf = field.get_col(_imax).data();
        int sendcount = field.get_col(_imax).size();
        MPI_Datatype sendtype = MPI_DOUBLE;
        int sendtag = 1;
        int source = domain.rank(domain.row, domain.col);

        // receiver
        int dest = domain.rank(domain.row, domain.col + 1);
        void *recvbuf = field.get_col(0).size();
        void *recvbuf;
        int recvcount = field.get_col(0).size();
        int recvcount = field.get_col(0).size();
        MPI_Datatype recvtype = MPI_DOUBLE;
        int recvtag = 1;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Status *status = MPI_STATUS_IGNORE;

        MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
    }

    // Send to Left subdomain
    if (domain.col > 0){

        // sender
        const void* sendbuf = field.get_col(1).data();
        int sendcount = field.get_col(1).size();
        MPI_Datatype sendtype = MPI_DOUBLE;
        int sendtag = 2;
        int source = domain.rank(domain.row, domain.col);

        // receiver
        int dest = domain.rank(domain.row, domain.col - 1);
        void *recvbuf = field.get_col(_imax).data();
        int recvcount = field.get_col(_imax).size();
        MPI_Datatype recvtype = MPI_DOUBLE;
        int recvtag = 2;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Status *status = MPI_STATUS_IGNORE;

        MPI_Sendrecv(sendbuf, sendcount, sendtype,   dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
    }
   
    // Send to Top subdomain
    if (domain.row < _jproc){

        // sender
        const void* sendbuf = field.get_row(_jmax).data();
        int sendcount = field.get_row(_jmax).size();
        MPI_Datatype sendtype = MPI_DOUBLE;
        int sendtag = 3;
        int source = domain.rank(domain.row, domain.col);

        // receiver
        int dest = domain.rank(domain.row + 1, domain.col);
        void *recvbuf = field.get_row(0).data();
        int recvcount = field.get_row(0).size();
        MPI_Datatype recvtype = MPI_DOUBLE;
        int recvtag = 3;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Status *status = MPI_STATUS_IGNORE;

        MPI_Sendrecv(sendbuf, sendcount, sendtype,   dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
    }
    
    // Send to Bottom subdomain
    if (domain.row > 0){

        // sender
        const void* sendbuf = field.get_row(1).data();
        int sendcount = field.get_row(1).size();
        MPI_Datatype sendtype = MPI_DOUBLE;
        int sendtag = 4;
        int source = domain.rank(domain.row, domain.col);

        // receiver
        int dest = domain.rank(domain.row - 1, domain.col);
        void *recvbuf = field.get_row(_jmax).data();
        int recvcount = field.get_row(_jmax).size();
        MPI_Datatype recvtype = MPI_DOUBLE;
        int recvtag = 4;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Status *status = MPI_STATUS_IGNORE;

        MPI_Sendrecv(sendbuf, sendcount, sendtype,   dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
    }
}


// // Find and return the minimum value over all ranks
// static double Communication::reduce_min(double value) {
//     double reduction_result = 0;
//     MPI_Reduce(&my_rank, &reduction_result, 1, MPI_INT, MPI_MIN, root_rank, MPI_COMM_WORLD);
//     return reduction_result;
// }

// // Compute total sum over all ranks
// double Communication::reduce_sum(const Domain &domain, double &value) {

//     sum_loc = 0;

// 	sum_loc += value;
    

//     /* accumulate global sum (across processes) */

//     MPI_Reduce(
//         &sum_loc,       /* send buffer */
//         &sum,           /* receive buffer (result on rank=0, see below) */
//         1,              /* send count */
//         MPI_INT,        /* data type */
//         MPI_SUM,        /* reduce operation type */
//         0,              /* result on rank=0 process */
//         MPI_COMM_WORLD  /* communicator */
//         );

//     /* note that only rank=0 has the final (correct) result */

//     if (rank == 0) {
//        printf("rank=%d  result=%d\n",rank,sum);
//     }


//     return sum;

// }
