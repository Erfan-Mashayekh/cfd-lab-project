#include <iostream>
#include <string>
#include <mpi.h>

#include "Case.hpp"

int main(int argn, char **args) {

    // Initialize MPI
    MPI_Init(&argn, &args);

    // Get the size of the default communicator
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // Get my process id (rank) in this communicator
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (argn > 1) {
        std::string file_name{args[1]};
        Case problem(file_name, comm_size, my_rank);
    } else {
        if(my_rank == 0){
            std::cout << "Error: No input file is provided to fluidchen." << std::endl;
            std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // problem.simulate();

    MPI_Barrier(MPI_COMM_WORLD);

    // Finalize MPI: Wait until all ranks are here to safely exit
    MPI_Finalize();
}
