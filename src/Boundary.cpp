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

OutflowBoundary::OutflowBoundary(std::vector<Cell *> cells) : _cells(cells) {}

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
        // Set inlet velocities (Neumann boundary)        
        // Set inlet velocities
        field.u(cell->i(), cell->j()) = _inlet_velocity_x;
        
        field.v(cell->i(), cell->j()) = _inlet_velocity_y;

        // Set F and G
        field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());

        field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());

        // TODO: Set inlet temperature

    }
}

void OutflowBoundary::apply(Fields &field){

    for (auto const& cell : _cells) {

        // Assume single fluid cell neighbour
        std::vector<border_position> border_pos = cell->borders();

        if(border_pos.size() != 1){
            std::cout << "Warning! Outflow can have only single fluid neighbour cells!" << std::endl;
        }
        field.u(cell->i()-1, cell->j()) = field.u(cell->neighbour(border_pos.at(0))->i()-1, cell->neighbour(border_pos.at(0))->j());
        
        field.v(cell->i(), cell->j()) = field.v(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j());

        // Set F and G
        field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());

        field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());

        field.p(cell->i(), cell->j()) = field.p(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j());
    }
}

void FixedWallBoundary::apply(Fields &field){

    for(auto const& cell: _cells){
        
        std::vector<border_position> border_positions = cell->borders();

        // Set velocity to zero at inner obstacle cells
        if(border_positions.size() == 0){
            // u 
            field.u(cell->i(), cell->j()) = 0.0;
            // v 
            field.v(cell->i(), cell->j()) = 0.0;
            // F 
            field.f(cell->i(), cell->j()) = 0.0;
            // G
            field.g(cell->i(), cell->j()) = 0.0;
            // p
            field.p(cell->i(), cell->j()) = 0.0;

            continue;
        }

        for(auto const& border_pos: border_positions){
            if(cell->is_border(border_position::RIGHT)){              
                // u 
                field.u(cell->i(), cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()) = - field.v(cell->neighbour(border_pos)->i(), cell->j());
                // F 
                field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
                // p
                field.p(cell->i(), cell->j()) = field.p(cell->neighbour(border_pos)->i(), cell->j());
            }
            if(cell->is_border(border_position::LEFT)){             
                // u 
                field.u(cell->i()-1, cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()) = - field.v(cell->neighbour(border_pos)->i(), cell->j());
                // F 
                field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());
                // p
                field.p(cell->i(), cell->j()) = field.p(cell->neighbour(border_pos)->i(), cell->j());
            }
            if(cell->is_border(border_position::TOP)){              
                // u 
                field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->neighbour(border_pos)->j());
                // v 
                field.v(cell->i(), cell->j()) = 0.0;
                // G
                field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
                // p
                field.p(cell->i(), cell->j()) = field.p(cell->i(), cell->neighbour(border_pos)->j());
            }
            if(cell->is_border(border_position::BOTTOM)){
                // u 
                field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->neighbour(border_pos)->j());
                // v 
                field.v(cell->i(), cell->j()-1) = 0.0;    
                // G
                field.g(cell->i(), cell->j()-1) = field.v(cell->i(), cell->j()-1);
                // p
                field.p(cell->i(), cell->j()) = field.p(cell->i(), cell->neighbour(border_pos)->j());
            }
        }
    }

    for(auto const& cell: _cells){
        
        std::vector<border_position> border_positions = cell->borders();
        
        for(auto const& border_pos: border_positions){
            if(cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT))
                field.p(cell->i(), cell->j()) = 0.0;
            if(cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT))
                field.p(cell->i(), cell->j()) = 0.0;
            if(cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT))
                field.p(cell->i(), cell->j()) = 0.0;
            if(cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT))
                field.p(cell->i(), cell->j()) = 0.0;                                                
        }
    }
    
    for(auto const& cell: _cells){
        
        std::vector<border_position> border_positions = cell->borders();
        
        for(auto const& border_pos: border_positions){
            if(cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)){ 
                // u 
                field.u(cell->i(), cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()) = 0.0;
                // F 
                field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
                // G
                field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
                // p
                field.p(cell->i(), cell->j()) += 0.5 * field.p(cell->neighbour(border_pos)->i(), cell->neighbour(border_pos)->j());
            }
            if(cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)){                              
                // u 
                field.u(cell->i()-1, cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()) = 0.0;
                // F 
                field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());
                // G
                field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
                // p
                field.p(cell->i(), cell->j()) += 0.5 * field.p(cell->neighbour(border_pos)->i(), cell->neighbour(border_pos)->j());
            }
            if(cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)){               
                // u 
                field.u(cell->i(), cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()-1) = 0.0;
                // F 
                field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
                // G
                field.g(cell->i(), cell->j()-1) = field.v(cell->i(), cell->j()-1);
                // p
                field.p(cell->i(), cell->j()) += 0.5 * field.p(cell->neighbour(border_pos)->i(), cell->neighbour(border_pos)->j());
                // std::cout << field.p(cell->i(), cell->j()) << "  " << std::endl;
            }
            if(cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)){               
                // u 
                field.u(cell->i()-1, cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()-1) = 0.0;
                // F 
                field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());
                // G
                field.g(cell->i(), cell->j()-1) = field.v(cell->i(), cell->j()-1);
                // p
                field.p(cell->i(), cell->j()) += 0.5 * field.p(cell->neighbour(border_pos)->i(), cell->neighbour(border_pos)->j());
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

        // v on the top
        field.v(cell->i(), cell->j()) = 0;
        // u on the top 
        field.u(cell->i(), cell->j()) = 2.0 * _wall_velocity[cell->wall_id()] - field.u(cell->i(), cell->neighbour(border_pos.at(0))->j());
        // F on the left 
        field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
    }
}


void FreeSlipBoundary::apply(Fields &field) {
    // Implement
}
