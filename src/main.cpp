#include <iostream>
#include <string>
#include <mpi.h>
#include <chrono>

#include "Case.hpp"
#include "Communication.hpp"

int main(int argn, char **args) {

    int comm_size;
    int my_rank;

    Communication::init_parallel(argn, args, my_rank, comm_size);

    if (argn <= 1) {
        if(my_rank == 0){
            std::cout << "Error: No input file is provided to fluidchen." << std::endl;
            std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
        }

        Communication::finalize();

        exit(EXIT_FAILURE);
    }

    std::string file_name{args[1]};

    Case problem(file_name, my_rank, comm_size);

    Communication::barrier();

    std::chrono::time_point<std::chrono::system_clock> start, end;

    start = std::chrono::system_clock::now();

    problem.simulate();

    end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_time = end - start;

    std::cout << "The computation took  " << elapsed_time.count() << " seconds\n";

    Communication::finalize();
}
