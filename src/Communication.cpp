#include "Communication.hpp"
#include <array>
#include <iostream>
#include <mpi.h>

using namespace std;

// Initialize communication
void Communication::init_parallel(int argn, char **args, int &my_rank, int &comm_size) {

    // start MPI
    MPI_Init(&argn, &args);

    // Fetch number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Fetch process id (rank) in this communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
}

// Barrier method to synchronise
void Communication::barrier() { MPI_Barrier(MPI_COMM_WORLD); }

// Finalize communication
void Communication::finalize() {
    // Finalize MPI: Wait until all ranks are here to safely exit
    Communication::barrier();

    MPI_Finalize();
}

// Broadcast method
void Communication::broadcast(void *buffer, int count, MPI_Datatype datatype, int root) {

    MPI_Bcast(buffer, count, datatype, root, MPI_COMM_WORLD);
}

// Send/Receive the data using this communicator
void Communication::communicate(Matrix<double> &field, const Domain &domain, const int &my_rank) {

    int iproc = domain.iproc;
    int jproc = domain.jproc;
    int imax = domain.imax;
    int jmax = domain.jmax;

    void *data0;
    void *data1;
    void *data2;
    void *data3;
    void *data4;
    void *data5;
    void *data6;
    void *data7;

    MPI_Datatype datatype = MPI_DOUBLE;
    int tag = 1;
    MPI_Comm communicator = MPI_COMM_WORLD;

    // Send to the right subdomain
    if ((my_rank + 1) % iproc != 0) {

        data0 = field.get_col(imax - 2).data();
       // std::cout << std::endl << "Get col after send to right subdomain is done" << std::endl;

        int count = field.get_col(imax - 2).size();
        int destination = my_rank + 1;

        double *data0new = (double *)data0;
        double bufferarray1[count];

        for (int i = 0; i < count; i++) {
            bufferarray1[i] = (double)*(data0new + i);
            //std::cout << *(data0new + i) << std::endl;
        }

        //     std::cout << std::endl << "Count is " << count << " It should be 50!!! " << std::endl;
        // std::cout << std::endl<< "The size of receive is " << *voidToVector.size() << " It should be 50!!! " <<
        // std::endl;

        MPI_Send(&bufferarray1, count, datatype, destination, tag, communicator);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Receive from the left
    if (my_rank % iproc != 0) {

        int count = field.get_col(0).size();
       // std::cout << std::endl << "Get col after receive from the left subdomain is done" << std::endl;

        int source = my_rank - 1;
        MPI_Status *status = MPI_STATUS_IGNORE;

        double bufferarray2[count];
        MPI_Recv(&bufferarray2, count, datatype, source, tag, communicator, status);

        // TODO: Check the Casting resluts
        // std::vector<double> vec;
        // data1 = &vec;
        // auto *voidToVector = static_cast<std::vector<double> *>(data1);

        // std::cout << std::endl << "You are here" << std::endl;

        // for (int i = 0; i < count; i++) {

        //     std::cout << *(data1new + i) << std::endl;
        // }
        //     std::cout << std::endl << "Count is " << count << " It should be 50!!! " << std::endl;

        // double *data1new = (double *)data1;
        vector<double> values(bufferarray2, bufferarray2 + count);
        // std::cout << std::endl << "Are you here tho?" << std::endl;

        // std::cout << std::endl << "The size of receive is " << (values).size() << std::endl;

        field.set_col(values, 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Send to the Left subdomain
    if ((my_rank) % iproc != 0) {

        data2 = field.get_col(1).data();

        //std::cout << std::endl << "Get col after send to the left subdomain is done" << std::endl;

        int count = field.get_col(1).size();
        int destination = my_rank - 1;

        double *data2new = (double *)data2;
        double bufferarray3[count];

        for (int i = 0; i < count; i++) {
            bufferarray3[i] = (double)*(data2new + i);
          //  std::cout << *(data2new + i) << std::endl;
        }

        MPI_Send(&bufferarray3, count, datatype, destination, tag, communicator);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Receive from the right
    if ((my_rank + 1) % iproc != 0) {

        int count = field.get_col(imax - 1).size();
       // std::cout << std::endl << "Get col after receive from the right subdomain is done" << std::endl;

        int source = my_rank + 1;
        MPI_Status *status = MPI_STATUS_IGNORE;


        double bufferarray4[count];
        MPI_Recv(&bufferarray4, count, datatype, source, tag, communicator, status);
       // std::cout << std::endl << "YOU HAVE RECEIVED!" << std::endl;

        // for(auto item: *voidToVector){
        //      std::cout << item << " ";
        // }


        // for (int i = 0; i < count; i++) {

        //     std::cout << *(data3new + i) << std::endl;
        // }

        vector<double> values(bufferarray4, bufferarray4 + count);

        //std::cout << std::endl << "The size of receive is " << values.size() << std::endl;

        field.set_col(values, imax - 1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // // Send to the Top subdomain
    //  if ( my_rank - iproc >= 0 ){

    //      data4 = field.get_row(jmax - 2).data();
    //      int count = field.get_row(jmax - 2).size();
    //      int destination = my_rank - iproc;

    //      MPI_Send(data4, count, datatype, destination, tag, communicator);
    //  }

    //  // Receive from the Bottom subdomain
    //  if ( my_rank < iproc * jproc - iproc ){

    //      int count = field.get_row(0).size();
    //      int source = my_rank + iproc;
    //      MPI_Status* status = MPI_STATUS_IGNORE;

    //      MPI_Recv(data5, count, datatype, source,tag, communicator, status);

    //      // TODO: Check the Casting resluts
    //      std::vector<double> vec;
    //      data5 = &vec;
    //      auto *voidToVector = reinterpret_cast< std::vector<double>* >(data5);

    //      field.set_row((*voidToVector), 0);
    //  }

    //  MPI_Barrier(MPI_COMM_WORLD);

    // // Send to the Bottom subdomain
    //  if ( my_rank < iproc * jproc - iproc ){

    //      data6 = field.get_row(1).data();
    //      int count = field.get_row(1).size();
    //      int destination = my_rank + iproc;

    //      MPI_Send(data6, count, datatype, destination, tag, communicator);
    //  }

    //  // Receive from the Top subdomain
    //  if ( my_rank - iproc >= 0 ){

    //      int count = field.get_row(jmax - 1).size();
    //      int source = my_rank - iproc;
    //      MPI_Status* status = MPI_STATUS_IGNORE;

    //      MPI_Recv(data7, count, datatype, source,tag, communicator, status);

    //      // TODO: Check the Casting resluts
    //      std::vector<double> vec;
    //      data7 = &vec;
    //      auto *voidToVector = reinterpret_cast< std::vector<double>* >(data7);

    //      field.set_row((*voidToVector), jmax - 1);
    //  }

    //  MPI_Barrier(MPI_COMM_WORLD);
}

// // Find and return the minimum value over all ranks
void Communication::reduce_min(double &input, double &output) {

    MPI_Allreduce(&input, &output, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
}

// // Compute total sum over all ranks
void Communication::reduce_sum(double &input, double &output) {

    MPI_Allreduce(&input, &output, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
