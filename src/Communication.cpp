#include "Communication.hpp"
#include <iostream>
#include <mpi.h>

Communication::Communication(int imax, int jmax, int iproc,  int jproc)
              :_imax(imax), _jmax(jmax), _iproc(iproc), _jproc(jproc){}

// Initialize communication
static void init_parallel(int argn, char** args, int &my_rank, int &comm_size){

    // start MPI
    MPI_Init(&argn, &args);

    // Fetch number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Fetch process id (rank) in this communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

}

// Finalize communication
static void finalize() { 
    // Finalize MPI: Wait until all ranks are here to safely exit
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize(); 
}

// Send/Receive the data using this communicator
void Communication::communicate(const Domain domain, Matrix<double> &field, int &my_rank, int &comm_size) {

    // Send to right subdomain
    if ( (my_rank + 1) % _iproc != 0 ){

        void* data = field.get_col(_imax).data();
        int count = field.get_col(_imax).size();
        MPI_Datatype datatype = MPI_DOUBLE;
        int tag = 1;
        int destination = my_rank + 1;
        MPI_Comm communicator = MPI_COMM_WORLD;

        MPI_Send(data, count, datatype, destination, tag, communicator);
    }

    // Receive from the left
    if ( my_rank % _iproc != 0 ){
        
        void* data; // field.get_col(_imax).data();
        int count = field.get_col(_imax).size();
        MPI_Datatype datatype = MPI_DOUBLE;
        int source = my_rank - 1;
        int tag = 1;
        MPI_Comm communicator = MPI_COMM_WORLD;
        MPI_Status* status = MPI_STATUS_IGNORE;

        MPI_Recv(data, count, datatype, source,tag, communicator, status);
        
        <<vector pointer>> p = reinterpret_cast< <<vector pointer>> >(voidPointerName);
        field.set_col(vec, 0);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Send to Left subdomain
    if ( my_rank % _iproc != 0 ){

        void* data = field.get_col(_imax).data();
        int count = field.get_col(_imax).size();
        MPI_Datatype datatype = MPI_DOUBLE;
        int tag = 1;
        int destination = my_rank + 1;
        MPI_Comm communicator = MPI_COMM_WORLD;

        MPI_Send(data, count, datatype, destination, tag, communicator);
    }
   
    // Send to Top subdomain

    
    // Send to Bottom subdomain
 
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
