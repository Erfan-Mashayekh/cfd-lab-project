#include "Boundary.hpp"
#include <cmath>
#include <iostream>


/***************************************
*
*  Constructors for different boundaries
*
***************************************/


InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double inlet_velocity_x, double inlet_velocity_y)
                : _cells(cells), _inlet_velocity_x(inlet_velocity_x), _inlet_velocity_y(inlet_velocity_y) {}

InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double inlet_velocity_x, double inlet_velocity_y, double inlet_temperature)
                : _cells(cells), _inlet_velocity_x(inlet_velocity_x), _inlet_velocity_y(inlet_velocity_y), _inlet_temperature(inlet_temperature) {}

OutflowBoundary::OutflowBoundary(std::vector<Cell *> cells, double initial_pressure) : _cells(cells), _initial_pressure(initial_pressure) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
                : _cells(cells), _wall_temperature(wall_temperature) {}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity) 
                : _cells(cells), _wall_velocity(wall_velocity) {}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity, std::map<int, double> wall_temperature)
                : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

FreeSlipBoundary::FreeSlipBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FreeSlipBoundary::FreeSlipBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
                : _cells(cells), _wall_temperature(wall_temperature) {}


/**************************************
*
*  Overrides for 'apply()' virtual func.
*
**************************************/


void InflowBoundary::apply(Fields &field){

    for (auto const& cell: _cells) {

        std::vector<border_position> border_pos = cell->borders();
        double sum_pressure = 0;
        int shift_left = 0;
        int shift_down = 0;
        // Account for the staggered grid for u & v.
        // Calculate shift based on the position of the boundary cell
        for (auto const& b_pos: border_pos){
            shift_left = (b_pos == border_position::LEFT) ? 1 : shift_left;
            shift_down = (b_pos == border_position::BOTTOM) ? 1 : shift_down;
            sum_pressure += field.p(cell->neighbour(b_pos)->i(), cell->neighbour(b_pos)->j());
        }
        // Set u, v, F & G (Dirichlet)
        field.u(cell->i() - shift_left, cell->j()) = _inlet_velocity_x;
        field.v(cell->i(), cell->j() - shift_down) = _inlet_velocity_y;
        field.f(cell->i() - shift_left, cell->j()) = field.u(cell->i() - shift_left, cell->j());
        field.g(cell->i(), cell->j() - shift_down) = field.v(cell->i(), cell->j() - shift_down); // Shouldn't we calculate F, G?
        // Set pressure (Neumann)
        field.p(cell->i(), cell->j()) = sum_pressure / border_pos.size();
    }
}

void OutflowBoundary::apply(Fields &field){

    for (auto const& cell : _cells) {
        
        std::vector<border_position> border_pos = cell->borders();
        double sum_pressure = 0;
        int shift_left = 0;
        int shift_down = 0;
        for(auto const& b_pos: border_pos){
            shift_left = (b_pos == border_position::LEFT) ? 1 : shift_left;
            shift_down = (b_pos == border_position::BOTTOM) ? 1 : shift_down;
            // Set u and v (Neumann)
            field.u(cell->i() - shift_left, cell->j()) = (b_pos == border_position::LEFT || b_pos == border_position::RIGHT) ? 
                                                         field.u(cell->neighbour(b_pos)->i() - shift_left, cell->neighbour(b_pos)->j()) :
                                                         field.u(cell->i() - shift_left, cell->j());    
            field.v(cell->i(), cell->j() - shift_down) = (b_pos == border_position::TOP || b_pos == border_position::BOTTOM) ? 
                                                         field.v(cell->neighbour(b_pos)->i(), cell->neighbour(b_pos)->j() - shift_down) :
                                                         field.v(cell->i(), cell->j() - shift_down);
            sum_pressure += field.p(cell->neighbour(b_pos)->i(), cell->neighbour(b_pos)->j());
        }
        // Set F and G (Neumann)
        field.f(cell->i() - shift_left, cell->j()) = field.u(cell->i() - shift_left, cell->j());
        field.g(cell->i(), cell->j() - shift_down) = field.v(cell->i(), cell->j() - shift_down); // Shouldn't we calculate F, G?
        // Set pressure (Neumann)
        field.p(cell->i(), cell->j()) = sum_pressure / border_pos.size();
    }
}

void FixedWallBoundary::apply(Fields &field){

    for(auto const& cell: _cells){
        
        std::vector<border_position> border_pos = cell->borders();

        // Set boundary conditions to zero if there are no bordering fluid cells (at inner obstacle cells)
        if(border_pos.size() == 0){
            field.u(cell->i(), cell->j()) = 0.0;
            field.v(cell->i(), cell->j()) = 0.0;
            field.f(cell->i(), cell->j()) = 0.0;
            field.g(cell->i(), cell->j()) = 0.0;
            continue; 
        }
        // Set boudndary conditions when there is only one boundary cell
        else {
            // Set u, v, f, g
            if(cell->is_border(border_position::RIGHT)){          
                field.v(cell->i(), cell->j() - 1) = - field.v(cell->i() + 1, cell->j() - 1);
                field.v(cell->i(), cell->j()) = - field.v(cell->i() + 1, cell->j());
            }
            if(cell->is_border(border_position::LEFT)){   
                field.v(cell->i(), cell->j()) = - field.v(cell->i() - 1, cell->j());
                field.v(cell->i(), cell->j() - 1) = - field.v(cell->i() - 1, cell->j() - 1);
            }
            if(cell->is_border(border_position::TOP)){
                field.u(cell->i() - 1, cell->j()) = - field.u(cell->i() - 1, cell->j() + 1);
                field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->j() + 1);
            }
            if(cell->is_border(border_position::BOTTOM)){   
                field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->j() - 1);
                field.u(cell->i() - 1, cell->j()) = - field.u(cell->i() - 1, cell->j() - 1);
            }
            if(cell->is_border(border_position::RIGHT)){          
                field.u(cell->i(), cell->j()) = 0.0;
                field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
            }
            if(cell->is_border(border_position::LEFT)){   
                field.u(cell->i()-1, cell->j()) = 0.0;
                field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());
            }
            if(cell->is_border(border_position::TOP)){
                field.v(cell->i(), cell->j()) = 0.0;
                field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
            }
            if(cell->is_border(border_position::BOTTOM)){   
                field.v(cell->i(), cell->j()-1) = 0.0;    
                field.g(cell->i(), cell->j()-1) = field.v(cell->i(), cell->j()-1);
            }
            // Set pressure
            if(border_pos.size() == 1){
                    field.p(cell->i(), cell->j()) = field.p(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j());
            } else if (border_pos.size() == 2) {
                field.p(cell->i(), cell->j()) = 0.5 * (field.p(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j())
                                                    +  field.p(cell->neighbour(border_pos.at(1))->i(), cell->neighbour(border_pos.at(1))->j()));
            }
        }
    }
}


void MovingWallBoundary::apply(Fields &field) {

    /* The top wall is a moving wall thus, the velocity on the top wall would be the velocity of the wall itself.*/

    for(auto cell: _cells){

        std::vector<border_position> border_pos = cell->borders();

        // Skip the two cells at the corners
        if(border_pos.size() == 0){
            continue;
        }

        if(border_pos.size() > 1){
            std::cout << "Warning! Moving wall can have only single fluid neighbour cells!" << std::endl;
        }

        field.v(cell->i(), cell->j()) = 0;
        field.u(cell->i(), cell->j()) = 2.0 * _wall_velocity[cell->wall_id()] - field.u(cell->i(), cell->neighbour(border_pos.at(0))->j());
        field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
        field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
        field.p(cell->i(), cell->j()) = field.p(cell->i(), cell->neighbour(border_position::BOTTOM)->j());
    }
}

void FreeSlipBoundary::apply(Fields &field) {(void)field;}


/***************************************************
*
*  Overrides for 'apply_temperature()' virtual func.
*
***************************************************/


void InflowBoundary::apply_temperature(Fields &field){

    for (auto const& cell: _cells) {
        // Set T (Dirichlet)
        field.T(cell->i(), cell->j()) =_inlet_temperature;
    }
}

void OutflowBoundary::apply_temperature(Fields &field){

    for (auto const& cell : _cells) {
        // Assume single fluid cell neighbour
        std::vector<border_position> border_pos = cell->borders();
        // Set T
        field.T(cell->i(), cell->j()) = field.T(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j());
    }
}

void FixedWallBoundary::apply_temperature(Fields &field) {
    for (auto const &cell : _cells) {

        std::vector<border_position> border_pos = cell->borders();

        // Adiabatic wall (Neumann boundary)
        if (_wall_temperature[cell->wall_id()] == -1) {
            // Neumann boundary condition
            if (border_pos.size() == 1) {
                field.T(cell->i(), cell->j()) = field.T(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j());
            } else if (border_pos.size() == 2) {
                field.T(cell->i(), cell->j()) = 0.5 * (field.T(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j())
                                                     + field.T(cell->neighbour(border_pos.at(1))->i(), cell->neighbour(border_pos.at(1))->j()));
            }
        } else {
            // Dirichlet Boundary Condition
            if (border_pos.size() == 1) {
                field.T(cell->i(), cell->j()) = 2 * _wall_temperature[cell->wall_id()] - field.T(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j());
            } else if (border_pos.size() == 2) {
                field.T(cell->i(), cell->j()) = 2 * _wall_temperature[cell->wall_id()] 
                                                - 0.5 * (field.T(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j()) 
                                                       + field.T(cell->neighbour(border_pos.at(1))->i(), cell->neighbour(border_pos.at(1))->j()));
            }
        }
    }
}

void MovingWallBoundary::apply_temperature(Fields &field) {(void)field;}

void FreeSlipBoundary::apply_temperature(Fields &field) {(void)field;}

