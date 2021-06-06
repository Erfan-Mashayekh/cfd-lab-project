#include "Case.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <cassert>
#include <regex>
#include <limits>

namespace filesystem = std::filesystem;

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

Case::Case(std::string file_name, int argn, char **args) {
    (void)argn;
    (void)args; // Remove if additional arguments are used
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu;      /* viscosity */
    double UI;      /* velocity x-direction */
    double VI;      /* velocity y-direction */
    double PI;      /* pressure */
    double GX;      /* gravitation x-direction */
    double GY;      /* gravitation y-direction */
    double xlength; /* length of the domain x-dir. */
    double ylength; /* length of the domain y-dir. */
    double dt;      /* time step */
    int imax;       /* number of cells x-direction */
    int jmax;       /* number of cells y-direction */
    double gamma;   /* uppwind differencing factor */
    double omg;     /* relaxation factor */
    double tau;     /* safety factor for time step */
    int itermax;    /* max. number of iterations for pressure per time step */
    double eps;     /* accuracy bound for pressure */
    double UIN;     /* inlet velocity x-direction */
    double VIN;     /* inlet velocity y-direction */
    double TI;      /* initial temperature */
    double TIN;     /* inlet temperature */
    double beta;    /* thermal expansion coefficient */
    double alpha;   /* thermal diffusivity */
    int num_walls;  /* number of walls */
    std::map<int, double> wall_vel;   /* Wall velocity against the wall index */
    std::map<int, double> wall_temp;  /* Wall temperature against the wall index */
    int iproc;       /* number of subdomains in x direction */
    int jproc;       /* number of subdomains in y direction */

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "geo_file") file >> _geom_name;
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
                if (var == "UIN") file >> UIN;
                if (var == "VIN") file >> VIN;
                if (var == "energy_eq"){
                    std::string state;
                    file >> state;
                    if (state == "on") _energy_eq = true;
                    else _energy_eq = false;
                }
                if (var == "TI") file >> TI;
                if (var == "TIN") file >> TIN;
                if (var == "beta") file >> beta;
                if (var == "alpha") file >> alpha;
                if (var == "num_walls") file >> num_walls;
                if (var == "iproc") file >> iproc;
                if (var == "jproc") file >> jproc;

                // In the following code,
                // - var reads the 'wall_vel_x' or 'wall_temp_x'
                // - regex_search for any one/two digit index
                // in the string and extracts the digit to idx
                // - then it checks if the read string contains
                // 'wal_vel' or 'wall_temp'
                // - Depending on the type of variable (vel/temp)
                // the respective value is stored in the map with 
                // value and the idx 
                std::string str_vel = "wall_vel";
                std::string str_temp = "wall_temp";
                std::regex match_idx("[0-9][0-9]*");
                std::smatch idx;
                double value;
                std::regex_search(var, idx, match_idx);

                if (std::search(var.begin(), var.end(), str_vel.begin(), str_vel.end()) != var.end()){
                    file >> value;
                    wall_vel.insert( std::pair<int, double>( std::stoi(idx[0]), value ) );
                }
                if (std::search(var.begin(), var.end(), str_temp.begin(), str_temp.end()) != var.end()){
                    file >> value;
                    wall_temp.insert( std::pair<int, double>( std::stoi(idx[0]), value ) );
                }
            }
        }
    }
    file.close();
    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Build up the domain
    Domain domain;
    domain.dx = xlength / (double)imax;
    domain.dy = ylength / (double)jmax;
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;

    _communication = Communication(imax, jmax, iproc, jproc);
    int rank = _communication.init_parallel(argn, args);

    build_domain(domain, imax, jmax, iproc, jproc, rank);

    _grid = Grid(_geom_name, domain);
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI, TI, alpha, beta, GX, GY);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;

    // Construct boundaries
    // Inflow
    if (not _grid.inflow_cells().empty()) {
        if(_energy_eq){
            _boundaries.push_back(std::make_unique<InflowBoundary>(_grid.inflow_cells(), UIN, VIN, TIN));
        } else {
            _boundaries.push_back(std::make_unique<InflowBoundary>(_grid.inflow_cells(), UIN, VIN));
        }
    }
    // Outflow
    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutflowBoundary>(_grid.outflow_cells(), PI));
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
            _communication.communicate(_grid.domain(), _field.T_matrix());
        }

        // Calculate Fn and Gn
        _field.calculate_fluxes(_grid, _energy_eq);
        _communication.communicate(_grid.domain(), _field.f_matrix());
        _communication.communicate(_grid.domain(), _field.g_matrix());

        // Calculate Right-hand side of the pressure eq.
        _field.calculate_rs(_grid);
        _communication.communicate(_grid.domain(), _field.rs_matrix());

        // SOR Loop
        // Initialization of residual and iteration counter
        int it = 0;
        // Set initial tolerance
        double res =  std::numeric_limits<double>::max();        
        
        while (res > _tolerance){
            
          
            // TODO: Set pressure Neumann Boundary Conditions
            //_field.set_pressure_bc(_grid);
            
            // Perform SOR Solver and retrieve residual for the loop continuity
            res = _pressure_solver->solve(_field, _grid, _boundaries);
            _communication.communicate(_grid.domain(), _field.p_matrix());
            

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

        // Calculate new time
        t = t + dt;

        // Increment the time step counter
        timestep++;

        // Calculate dt for adaptive time stepping
        dt = _field.calculate_dt(_grid, _energy_eq);


        // Output the vtk every 1s
        if (t >= step + _output_freq) {
            step = step + _output_freq;
            std::cout << "Printing vtk file at t = " << step << "s" << std::endl;
            output_vtk(step, 0);
        }
     
    }

    // Output the final VTK file
    output_vtk(step, 0);

    // End message
    std::cout << "Done!\n";

}

void Case::output_vtk(int timestep, int my_rank) {
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
        _dict_name + '/' + _case_name + "" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();

}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain, int iproc, int jproc, int rank) {
    
    domain.rank = Matrix<int>(jproc, iproc, 0);

    for (int col = 0; col < iproc; col++){
        for (int row = 0; row < jproc; row++){    
            domain.imin = 0;
            domain.jmin = 0;
            domain.imax = int(imax_domain/iproc) + 2;
            domain.jmax = int(jmax_domain/jproc) + 2;
            domain.size_x = int(imax_domain/iproc);
            domain.size_y = int(jmax_domain/jproc);

            domain.row = row;
            domain.col = col;
            domain.rank(row, col) = rank;
        }
    }
}