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
    
    int iproc = domain.iproc;
    int jproc = domain.jproc;
    int imax = domain.imax;
    int jmax = domain.jmax;

    void* data;
    MPI_Datatype datatype = MPI_DOUBLE;
    int tag = 1;
    MPI_Comm communicator = MPI_COMM_WORLD;

    // Send to the right subdomain
    if ( (my_rank + 1) % iproc != 0 ){

        data = field.get_col(imax).data();
        int count = field.get_col(imax).size();
        int destination = my_rank + 1;   

        MPI_Send(data, count, datatype, destination, tag, communicator);
    }

    // Receive from the left
    if ( my_rank % iproc != 0 ){
        
        int count = field.get_col(0).size();
        int source = my_rank - 1;
        MPI_Status* status = MPI_STATUS_IGNORE;

        MPI_Recv(data, count, datatype, source,tag, communicator, status);

        // TODO: Check the Casting resluts
        std::vector<double> vec;
        data = &vec;
        auto *voidToVector = reinterpret_cast< std::vector<double>* >(data);

        field.set_col((*voidToVector), 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Send to the Left subdomain
    if ( (my_rank) % iproc != 0 ){

        data = field.get_col(1).data();
        int count = field.get_col(1).size();
        int destination = my_rank - 1;   

        MPI_Send(data, count, datatype, destination, tag, communicator);
    }

    // Receive from the right
    if ( (my_rank + 1) % iproc != 0 ){
        
        int count = field.get_col(imax + 1).size();
        int source = my_rank + 1;
        MPI_Status* status = MPI_STATUS_IGNORE;

        MPI_Recv(data, count, datatype, source,tag, communicator, status);

        // TODO: Check the Casting resluts
        std::vector<double> vec;
        data = &vec;
        auto *voidToVector = reinterpret_cast< std::vector<double>* >(data);

        field.set_col((*voidToVector), imax + 1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

   // Send to the Top subdomain
    if ( my_rank - iproc >= 0 ){

        data = field.get_row(jmax).data();
        int count = field.get_row(jmax).size();
        int destination = my_rank - iproc;   

        MPI_Send(data, count, datatype, destination, tag, communicator);
    }

    // Receive from the Bottom subdomain
    if ( my_rank < iproc * jproc - iproc ){
        
        int count = field.get_row(0).size();
        int source = my_rank + iproc;
        MPI_Status* status = MPI_STATUS_IGNORE;

        MPI_Recv(data, count, datatype, source,tag, communicator, status);

        // TODO: Check the Casting resluts
        std::vector<double> vec;
        data = &vec;
        auto *voidToVector = reinterpret_cast< std::vector<double>* >(data);

        field.set_row((*voidToVector), 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

   // Send to the Bottom subdomain
    if ( my_rank < iproc * jproc - iproc ){

        data = field.get_row(1).data();
        int count = field.get_row(1).size();
        int destination = my_rank + iproc;   

        MPI_Send(data, count, datatype, destination, tag, communicator);
    }

    // Receive from the Top subdomain
    if ( my_rank - iproc >= 0 ){
        
        int count = field.get_row(jmax + 1).size();
        int source = my_rank - iproc;
        MPI_Status* status = MPI_STATUS_IGNORE;

        MPI_Recv(data, count, datatype, source,tag, communicator, status);

        // TODO: Check the Casting resluts
        std::vector<double> vec;
        data = &vec;
        auto *voidToVector = reinterpret_cast< std::vector<double>* >(data);

        field.set_row((*voidToVector), jmax + 1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

}

// // Find and return the minimum value over all ranks
void Communication::reduce_min(double &input, double &output) {

    MPI_Allreduce(&input, &output, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}


// // Compute total sum over all ranks
void Communication::reduce_sum(double &input, double &output) {

    MPI_Allreduce(&input, &output, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


