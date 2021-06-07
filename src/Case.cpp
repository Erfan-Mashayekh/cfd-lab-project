#include "Case.hpp"
#include "Enums.hpp"
#include "Parameters.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <cassert>
#include <regex>
#include <limits>
#include <cstddef>
#include <array>
#include <cmath>

namespace filesystem = std::filesystem;

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>


Case::Case(std::string file_name, int &my_rank, int &comm_size)
     : _my_rank(my_rank), _comm_size(comm_size) { 

    Parameters input; 

    // create a type for struct input 
    const int nitems = 31;
    int blocklengths[nitems] = {1,1,1,1,1, 1,1,1,1,1,
                                1,1,1,1,1, 1,1,1,1,1,
                                1,1,1,1,1, 1,1,1,4,4,
                                4};

    MPI_Datatype types[nitems] = {MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, 
                                  MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, 

                                  MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, 
                                  MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, MPI::DOUBLE, 

                                  MPI::DOUBLE, MPI::INT, MPI::INT, MPI::INT, MPI::INT,
                                  MPI::INT, MPI::INT, MPI::BOOL, MPI::INT, MPI::DOUBLE,

                                  MPI::DOUBLE};

    MPI_Datatype mpi_param_type;
    MPI_Aint offsets[nitems];

    offsets[0] = offsetof(Parameters, xlength);
    offsets[1] = offsetof(Parameters, ylength);
    offsets[2] = offsetof(Parameters, nu);
    offsets[3] = offsetof(Parameters, dt);
    offsets[4] = offsetof(Parameters, omg);
    offsets[5] = offsetof(Parameters, eps);
    offsets[6] = offsetof(Parameters, tau);
    offsets[7] = offsetof(Parameters, gamma);
    offsets[8] = offsetof(Parameters, UI);
    offsets[9] = offsetof(Parameters, VI);
    offsets[10] = offsetof(Parameters, PI);
    offsets[11] = offsetof(Parameters, TI);
    offsets[12] = offsetof(Parameters, GX);
    offsets[13] = offsetof(Parameters, GY);
    offsets[14] = offsetof(Parameters, UIN);
    offsets[15] = offsetof(Parameters, VIN);
    offsets[16] = offsetof(Parameters, TIN);
    offsets[17] = offsetof(Parameters, beta); 
    offsets[18] = offsetof(Parameters, alpha);
    offsets[19] = offsetof(Parameters, t_end);
    offsets[20] = offsetof(Parameters, output_freq);
    offsets[21] = offsetof(Parameters, itermax);
    offsets[22] = offsetof(Parameters, imax);
    offsets[23] = offsetof(Parameters, jmax);
    offsets[24] = offsetof(Parameters, iproc);
    offsets[25] = offsetof(Parameters, jproc);
    offsets[26] = offsetof(Parameters, num_walls);
    offsets[27] = offsetof(Parameters, energy_eq);
    offsets[28] = offsetof(Parameters, wall_idx);
    offsets[29] = offsetof(Parameters, wall_vel);
    offsets[30] = offsetof(Parameters, wall_temp);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_param_type);
    MPI_Type_commit(&mpi_param_type);

    /****************************
     * Rank 0 reads the file
     * *************************/
    if(_my_rank == 0){
        // Read input parameters
        const int MAX_LINE_LENGTH = 1024;
        std::ifstream file(file_name);

        if (file.is_open()) {

            std::string var;
            int index = 0;
            while (!file.eof() && file.good()) {
                file >> var;
                if (var[0] == '#') { /* ignore comment line*/
                    file.ignore(MAX_LINE_LENGTH, '\n');
                } else {
                    if (var == "xlength") file >> input.xlength;
                    if (var == "ylength") file >> input.ylength;
                    if (var == "nu") file >> input.nu;
                    if (var == "dt") file >> input.dt;
                    if (var == "omg") file >> input.omg;
                    if (var == "eps") file >> input.eps;
                    if (var == "tau") file >> input.tau;
                    if (var == "gamma") file >> input.gamma;
                    if (var == "UI") file >> input.UI;
                    if (var == "VI") file >> input.VI;
                    if (var == "PI") file >> input.PI;
                    if (var == "TI") file >> input.TI;
                    if (var == "GX") file >> input.GX;
                    if (var == "GY") file >> input.GY;
                    if (var == "UIN") file >> input.UIN;
                    if (var == "VIN") file >> input.VIN;
                    if (var == "TIN") file >> input.TIN;
                    if (var == "beta") file >> input.beta;
                    if (var == "alpha") file >> input.alpha;
                    if (var == "itermax") file >> input.itermax;
                    if (var == "imax") file >> input.imax;
                    if (var == "jmax") file >> input.jmax;
                    if (var == "iproc") file >> input.iproc;
                    if (var == "jproc") file >> input.jproc;
                    if (var == "num_walls") file >> input.num_walls;
                    // Class members
                    if (var == "t_end") file >> input.t_end;
                    if (var == "dt_value") file >> input.output_freq;
                    if (var == "energy_eq"){
                        std::string state;
                        file >> state;
                        if (state == "on") input.energy_eq = true;
                        else input.energy_eq = false;
                    }

                    /************************************************    
                    * In the following code,
                    * - var reads the 'wall_vel_x' or 'wall_temp_x'
                    * - regex_search for any one/two digit index
                    * in the string and extracts the digit to idx
                    * - then it checks if the read string contains
                    * 'wal_vel' or 'wall_temp'
                    * - Depending on the type of variable (vel/temp)
                    * the respective value is stored in the map with 
                    * value and the idx 
                    *************************************************/
                    std::string str_vel = "wall_vel";
                    std::string str_temp = "wall_temp";
                    std::regex match_idx("[0-9][0-9]*");
                    std::smatch idx;
                    double value;
                    std::regex_search(var, idx, match_idx);

                    if (std::search(var.begin(), var.end(), str_vel.begin(), str_vel.end()) != var.end()){
                        file >> value;
                        if(std::find(input.wall_idx.begin(), input.wall_idx.end(), std::stoi(idx[0])) == input.wall_idx.end()){
                            input.wall_idx[index] = std::stoi(idx[0]);
                            input.wall_vel[index] = value;
                            index++;
                        } else {
                            input.wall_vel[* std::find(input.wall_idx.begin(), input.wall_idx.end(), std::stoi(idx[0]))] = value;
                        }
                    }
                    if (std::search(var.begin(), var.end(), str_temp.begin(), str_temp.end()) != var.end()){
                        file >> value;
                        if(std::find(input.wall_idx.begin(), input.wall_idx.end(), std::stoi(idx[0])) == input.wall_idx.end()){
                            input.wall_idx[index] = std::stoi(idx[0]);
                            input.wall_temp[index] = value;
                            index++;
                        } else {
                            input.wall_temp[* std::find(input.wall_idx.begin(), input.wall_idx.end(), std::stoi(idx[0]))] = value;
                        }
                    }
                }
            }
        }
        file.close();

        if(comm_size != input.iproc * input.jproc){
            std::cout << "ERROR: The simulation must be run on " << input.iproc * input.jproc << " number of processes!. Terminating!" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // Constrcut the geometry file name as the same as file_name
    // But as .pgm format
    _geom_name.assign(file_name, 0, file_name.size() - 4);
    _geom_name.append(".pgm");

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Broadcast input data to all processes in the communicator
    Communication::broadcast(&input, 1, mpi_param_type, 0);

    // Set up wall_vel and wall_temp as std::maps
    std::map<int, double> wall_vel;   
    std::map<int, double> wall_temp; 

    for(size_t i = 0; i < input.wall_idx.size(); i++){
        if(input.wall_idx[i] != -1){
            wall_vel.insert( std::pair<int, double>( input.wall_idx[i], input.wall_vel[i] ) );
            wall_temp.insert( std::pair<int, double>( input.wall_idx[i], input.wall_temp[i] ) );
        }
    }

    // Build up the domain for each process
    Domain subdomain;
    build_domain(subdomain, input.xlength, input.ylength, input.imax, input.jmax, input.iproc, input.jproc);

    _grid = Grid(_geom_name, subdomain, _my_rank);
    _field = Fields(input.nu, input.dt, input.tau, _grid.domain().size_x, _grid.domain().size_y, 
                    input.UI, input.VI, input.PI, input.TI, input.alpha, input.beta, input.GX, input.GY);

    _discretization = Discretization(_grid.domain().dx, _grid.domain().dy, input.gamma);
    _pressure_solver = std::make_unique<SOR>(input.omg);
    _max_iter = input.itermax;
    _tolerance = input.eps;
    _output_freq = input.output_freq;
    _t_end = input.t_end;
    _energy_eq = input.energy_eq;

    // Construct boundaries
    // Inflow
    if (not _grid.inflow_cells().empty()) {
        if(_energy_eq){
            _boundaries.push_back(std::make_unique<InflowBoundary>(_grid.inflow_cells(), input.UIN, input.VIN, input.TIN));
        } else {
            _boundaries.push_back(std::make_unique<InflowBoundary>(_grid.inflow_cells(), input.UIN, input.VIN));
        }
    }
    // Outflow
    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutflowBoundary>(_grid.outflow_cells(), input.PI));
    }

    // Fixed wall
    if (not _grid.fixed_wall_cells().empty()) {
        if(_energy_eq){
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells(), wall_temp));
        } else {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
        }
    }
    // Moving wall
    if (not _grid.moving_wall_cells().empty()) {
        if(_energy_eq){
            _boundaries.push_back(std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), wall_vel, wall_temp));
        } else {
            _boundaries.push_back(std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), wall_vel));
        }
    }
    // Free slip
    if (not _grid.free_slip_cells().empty()) {
        if(_energy_eq){
            _boundaries.push_back(std::make_unique<FreeSlipBoundary>(_grid.free_slip_cells(), wall_temp));
        } else {
            _boundaries.push_back(std::make_unique<FreeSlipBoundary>(_grid.free_slip_cells()));
        } 
    }
}

void Case::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - _case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * - Calculate and apply boundary conditions for all the boundaries in _boundaries container
 *   using apply() member function of Boundary class
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 *   Flux consists of diffusion and convection part, which are located in Discretization class
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver class
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - Calculat the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface. No member functions should be defined in abstract classes. You need to define functions in inherited
 * classes such as MovingWallBoundary class.
 *
 * For information about the classes and functions, you can check the header files.
 */
void Case::simulate() {

    assert(_output_freq > 0); //
    std::cout << "Fluidchen is running and will print vtk output every "
              << _output_freq <<"s until " << _t_end << "s..." << std::endl;

    // initialization
    double t = 0.0;
    double dt_proc = 0;
    double dt = _field.dt();
    int timestep = 0;
    double step = 0;
    // time loop
    while (t < _t_end) {

        // Applying velocity boundary condition for every 4 sides of the wall boundary, inflow, and outflow
        for (auto &boundary : _boundaries) {
            boundary->apply(_field);
        }

        // Calculate Temperature if the energy equation is on
        if(_energy_eq){
            for (auto &boundary : _boundaries) {
                boundary->apply_temperature(_field);
            }
            _field.calculate_temperature(_grid);
            _communication.communicate(_field.T_matrix(), _grid.domain().imax, _grid.domain().jmax, _grid.domain().domain_iproc, _grid.domain().domain_jproc, _my_rank);
        }

        // Calculate Fn and Gn
        _field.calculate_fluxes(_grid, _energy_eq);
        _communication.communicate(_field.f_matrix(), _grid.domain().imax, _grid.domain().jmax, _grid.domain().domain_iproc, _grid.domain().domain_jproc, _my_rank);
        _communication.communicate(_field.g_matrix(), _grid.domain().imax, _grid.domain().jmax, _grid.domain().domain_iproc, _grid.domain().domain_jproc, _my_rank);

        // Calculate Right-hand side of the pressure eq.
        _field.calculate_rs(_grid);

        // SOR Loop
        // Initialization of residual and iteration counter
        int it = 0;
        // Set initial tolerance
        double res_sub =  std::numeric_limits<double>::max();  
        double res =  std::numeric_limits<double>::max();  

        while (res > _tolerance){
            
          
            // TODO: Set pressure Neumann Boundary Conditions
            //_field.set_pressure_bc(_grid);
            
            // Perform SOR Solver and retrieve residual for the loop continuity
            res_sub = _pressure_solver->solve(_field, _grid, _boundaries);
            
            _communication.communicate(_field.p_matrix(), _grid.domain().imax, _grid.domain().jmax, _grid.domain().domain_iproc, _grid.domain().domain_jproc, _my_rank);
            MPI_Barrier(MPI_COMM_WORLD);
            res = _communication.reduce_sum(res_sub);


            // Increment the iteration counter
            it++;

            // Check if SOR didn't converge
            if(it > _max_iter) {
               std::cout << "WARNING! SOR reached maximum number of iterations at t = " 
               << t << ". Residual = " << res << std::endl;
               break;
            }
        }


        // Calculate the velocities at the next time step
        _field.calculate_velocities(_grid);
        _communication.communicate(_field.u_matrix(), _grid.domain().imax, _grid.domain().jmax, _grid.domain().domain_iproc, _grid.domain().domain_jproc, _my_rank);
        _communication.communicate(_field.v_matrix(), _grid.domain().imax, _grid.domain().jmax, _grid.domain().domain_iproc, _grid.domain().domain_jproc, _my_rank);

        // Calculate new time
        t = t + dt;

        // Increment the time step counter
        timestep++;

        // Calculate dt for adaptive time stepping
        dt_proc = _field.calculate_dt(_grid, _energy_eq);
        MPI_Barrier(MPI_COMM_WORLD);
        dt = _communication.reduce_min(dt_proc);

        // Output the vtk every 1s
        if (t >= step + _output_freq) {
            step = step + _output_freq;
            std::cout << "Printing vtk file at t = " << step << "s" << std::endl;
            output_vtk(step);
        }
     
    }

    // Output the final VTK file
    output_vtk(step);

    // End message
    std::cout << "Done!\n";

}

void Case::output_vtk(int timestep) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().imin * dx;
    double y = _grid.domain().jmin * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    // Pressure Array
    vtkDoubleArray *Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

     // Temperature Array
    vtkDoubleArray *Temperature = vtkDoubleArray::New();
    Temperature->SetName("temperature");
    Temperature->SetNumberOfComponents(1);



    // Print pressure and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            double temperature = _field.T(i, j);
            Pressure->InsertNextTuple(&pressure);
            Temperature->InsertNextTuple(&temperature);
        }
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Add Temperature to Structured Grid
    if(_energy_eq){
        structuredGrid->GetCellData()->AddArray(Temperature);
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "" + std::to_string(_my_rank) + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();

}

void Case::build_domain(Domain &domain, double xlength, double ylength, int imax_domain, int jmax_domain, int iproc, int jproc) {

    domain.dx = xlength / (double)imax_domain;
    domain.dy = ylength / (double)jmax_domain;
    domain.domain_size_x = imax_domain;
    domain.domain_size_y = jmax_domain;
    domain.domain_iproc  = iproc;
    domain.domain_jproc  = jproc;
    domain.x_proc = _my_rank % iproc;
    domain.y_proc = std::floor(_my_rank / iproc);

    int max_partition_size_x = std::floor(imax_domain / iproc) + 1;
    int max_partition_size_y = std::floor(jmax_domain / jproc) + 1;

    int residual_partition_size_x = imax_domain - max_partition_size_x * (iproc - 1);
    int residual_partition_size_y = jmax_domain - max_partition_size_y * (jproc - 1);

    if(domain.x_proc == iproc - 1){
        domain.size_x = residual_partition_size_x;
    } else {
        domain.size_x = max_partition_size_x;
    }

    if(domain.y_proc == jproc - 1){
        domain.size_y = residual_partition_size_y;
    } else {
        domain.size_y = max_partition_size_y;
    }

    domain.imin = 0;
    domain.jmin = 0;
    domain.imax = domain.size_x + 2;
    domain.jmax = domain.size_y + 2;

}
