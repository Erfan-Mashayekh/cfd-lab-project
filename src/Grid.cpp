#include "Grid.hpp"
#include "Communication.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>
#include <cmath>
#include <mpi.h>

Grid::Grid(std::string geom_name, Domain &domain, const int& my_rank) {

    _domain = domain;

    _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);

    if (geom_name.compare("NONE") == 0) {
        std::cout << "Error: Please provide a geometry data file as a .pgm file in the .dat file!. Exiting!\n";
        exit(EXIT_FAILURE);
    } 

    Matrix<int> geometry_data = Matrix<int>(_domain.domain_size_x + 2, _domain.domain_size_y + 2, 0);
    
    if(my_rank == 0){
        parse_geometry_file(geom_name, geometry_data);
    }

    Communication::broadcast((void*)geometry_data.data(), ((_domain.domain_size_x + 2) * (_domain.domain_size_y + 2)) , MPI::INT, 0);
    
    assign_cell_types(geometry_data, my_rank);
}

void Grid::assign_cell_types(Matrix<int> &geometry_data, const int& my_rank) {
    /*
        Geometry ID for each cell:
            FLUID,         ->  0       
            INFLOW,        ->  1       
            OUTFLOW,       ->  2       
            FIXED_WALL,    ->  3-7     
            MOVING_WALL,   ->  8       
            FREE_SLIP_WALL ->  9                               
    */
    int i = 0;
    int j = 0;
    int i_geom;
    int j_geom;

    for (int j_domain = _domain.jmin; j_domain < _domain.jmax; ++j_domain) {
        
        i = 0;

        for (int i_domain = _domain.imin; i_domain < _domain.imax; ++i_domain) {
            // Shift each subdomain to its correct place in a parallel computation
            i_geom = i_domain + _domain.max_size_x * (my_rank % _domain.iproc);
            j_geom = j_domain + _domain.max_size_y * std::floor(my_rank / _domain.iproc);

            if (geometry_data(i_geom, j_geom) == 0) {
                // Fluid
                _cells(i, j) = Cell(i, j, cell_type::FLUID, geometry_data(i_geom, j_geom));
            } else if (geometry_data(i_geom, j_geom) == 1) {
                // Inlet
                _cells(i, j) = Cell(i, j, cell_type::INFLOW, geometry_data(i_geom, j_geom));
            } else if (geometry_data(i_geom, j_geom) == 2) {
                // Outlet
                _cells(i, j) = Cell(i, j, cell_type::OUTFLOW, geometry_data(i_geom, j_geom));
            } else if (geometry_data(i_geom, j_geom) <= 7) { // Numbers 3-7
                // Fixed wall
                _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data(i_geom, j_geom));
            } else if (geometry_data(i_geom, j_geom) == 8) {
                // Moving wall
                _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL, geometry_data(i_geom, j_geom));
            } else if (geometry_data(i_geom, j_geom) == 9) {
                // Free slip
                _cells(i, j) = Cell(i, j, cell_type::FREE_SLIP_WALL, geometry_data(i_geom, j_geom));
            } 

            ++i;
        }

        ++j;
    }

    // Subdomain Ghost cell in parallel computation 
    for (int i = 0; i < _domain.size_x + 2; ++i) {
        if(_cells(i, 0).type() == cell_type::FLUID){
            _cells(i, 0) = Cell(i, 0, cell_type::GHOST);
        }
        if(_cells(i, _domain.size_y + 1).type() == cell_type::FLUID){
            _cells(i, _domain.size_y + 1) = Cell(i, _domain.size_y + 1, cell_type::GHOST);
        }
    }

    // Subdomain Ghost cell in parallel computation 
    for (int j = 0; j < _domain.size_y + 2; ++j) {
        if(_cells(0, j).type() == cell_type::FLUID){
            _cells(0, j) = Cell(0, j, cell_type::GHOST);
        }
        if(_cells(_domain.size_x + 1, j).type() == cell_type::FLUID){
            _cells(_domain.size_x + 1, j) = Cell(_domain.size_x + 1, j, cell_type::GHOST);
        }
    }  


    for (int i = 0; i < _domain.size_x + 2; ++i){
        for (int j = 0; j < _domain.size_y + 2; ++j) {

            if(_cells(i, j).type() == cell_type::FLUID){
                _fluid_cells.push_back(&_cells(i, j));
            } else if(_cells(i, j).type() == cell_type::INFLOW){
                _inflow_cells.push_back(&_cells(i, j));
            } else if(_cells(i, j).type() == cell_type::OUTFLOW){
                _outflow_cells.push_back(&_cells(i, j));
            } else if(_cells(i, j).type() == cell_type::FIXED_WALL){
                _fixed_wall_cells.push_back(&_cells(i, j));
            } else if(_cells(i, j).type() == cell_type::MOVING_WALL){   
                _moving_wall_cells.push_back(&_cells(i, j));
            } else if(_cells(i, j).type() == cell_type::FREE_SLIP_WALL){    
                _free_slip_cells.push_back(&_cells(i, j));
            } else if(_cells(i, j).type() == cell_type::GHOST){    
                _ghost_cells.push_back(&_cells(i, j));
            }
        }
    }

    // Corner cell neighbour assigment
    // Bottom-Left Corner
    i = 0;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID || 
        _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::TOP);
    }
    if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID || 
        _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::RIGHT);
    } 

    // Top-Left Corner
    i = 0;
    j = _domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID || 
        _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::BOTTOM);
    }
    if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID || 
        _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::RIGHT);
    }  

    // Top-Right Corner
    i = _domain.size_x + 1;
    j = _domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID || 
        _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::BOTTOM);
    }
    if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID || 
        _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::LEFT);
    }  

    // Bottom-Right Corner
    i = _domain.size_x + 1;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID || 
        _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::TOP);
    }
    if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID || 
        _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
        _cells(i, j).add_border(border_position::LEFT);
    } 

    // Bottom cells
    j = 0;
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::TOP);
        }
    }

    // Top Cells
    j = _domain.size_y + 1;

    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
    }

    // Left Cells
    i = 0;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::TOP);
        }       
    }
    // Right Cells
    i = _domain.size_x + 1;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID || 
            _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::GHOST) {
            _cells(i, j).add_border(border_position::TOP);
        }        
    }

    // Inner cells
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        for (int j = 1; j < _domain.size_y + 1; ++j) {
            _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
            _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
            _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
            _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);

            if (_cells(i, j).type() != cell_type::FLUID) {
                if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID || 
                    _cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::GHOST) {
                    _cells(i, j).add_border(border_position::LEFT);
                }
                if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID || 
                    _cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::GHOST) {
                    _cells(i, j).add_border(border_position::RIGHT);
                }
                if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID || 
                    _cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::GHOST) {
                    _cells(i, j).add_border(border_position::BOTTOM);
                }
                if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID || 
                    _cells(i, j).neighbour(border_position::TOP)->type() == cell_type::GHOST) {
                    _cells(i, j).add_border(border_position::TOP);
                }
            }
        }
    }

    // if(my_rank == 1){
    //     for (int col =  _domain.jmax - 1;  _domain.jmax > -1; --col) {
    //         for (size_t row = 0; row < _domain.imax; ++row) {
    //             std::cout << (int) _cells(row, col).type() << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }


    /*******************************************
     * Terminate if Forbidden cells are present 
     ******************************************/

    std::vector<int> forbidden_cells;

    for(size_t i; i < _fixed_wall_cells.size(); i++){

        std::vector<border_position> border_positions = _fixed_wall_cells[i]->borders();

        if(border_positions.size() > 2){
            forbidden_cells.push_back(i);
        }
    }


    if(forbidden_cells.size() > 0){
        
        std::cout << "Error! Forbidden cell found at [" 
        << _fixed_wall_cells.at(i)->i() << ", " << _fixed_wall_cells.at(i)->j() 
        << "]. Exiting ... Try converting to fluid cell..." << std::endl;

        exit(EXIT_FAILURE);        
    }
}



void Grid::parse_geometry_file(std::string filedoc, Matrix<int> &geometry_data) {

    size_t numcols, numrows, depth;

    std::ifstream infile(filedoc);
    std::stringstream ss;
    std::string inputLine = "";

    // First line : version
    getline(infile, inputLine);
    if (inputLine.compare("P2") != 0) {
        std::cerr << "First line of the PGM file should be P2" << std::endl;
    }

    // Second line : comment
    getline(infile, inputLine);

    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numrows >> numcols;
    // Fourth line : depth
    ss >> depth;

    assert(geometry_data.size() == numrows * numcols);

    // Following lines : data
    for (int col = numcols - 1; col > -1; --col) {
        for (size_t row = 0; row < numrows; ++row) {
            ss >> geometry_data(row, col);
        }
    } 

    infile.close();

}

int Grid::imax() const { return _domain.size_x; }
int Grid::jmax() const { return _domain.size_y; }

int Grid::imaxb() const { return _domain.size_x + 2; }
int Grid::jmaxb() const { return _domain.size_y + 2; }

Cell Grid::cell(int i, int j) const { return _cells(i, j); }

double Grid::dx() const { return _domain.dx; }

double Grid::dy() const { return _domain.dy; }

const Domain &Grid::domain() const { return _domain; }

const std::vector<Cell *> &Grid::fluid_cells() const { return _fluid_cells; }

const std::vector<Cell *> &Grid::fixed_wall_cells() const { return _fixed_wall_cells; }

const std::vector<Cell *> &Grid::moving_wall_cells() const { return _moving_wall_cells; }

const std::vector<Cell *> &Grid::free_slip_cells() const{ return _free_slip_cells; }

const std::vector<Cell *> &Grid::inflow_cells() const { return _inflow_cells; }

const std::vector<Cell *> &Grid::outflow_cells() const { return _outflow_cells; }

const std::vector<Cell *> &Grid::ghost_cells() const { return _ghost_cells; }
