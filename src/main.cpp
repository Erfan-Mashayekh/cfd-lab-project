#include <iostream>
#include <string>
#include <mpi.h>

#include "Case.hpp"

int main(int argn, char **args) {

    int comm_size;
    int my_rank;

    Communication::init_parallel(argn, args, my_rank, comm_size);

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

    Communication::finalize();
}
