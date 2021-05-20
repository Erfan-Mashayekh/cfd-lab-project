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


void FixedWallBoundary::apply(Fields &field){
    /** 
     * All boundaries except the top are fixed wall or no slip. 
     * 
     * U on the left and right are set to 0 and v at the bottom is set to 0.
     * 
     * No v-values lie on the vertical boundaries and no u-values lie on 
     * the horizontal boundaries, the boundary value zero is achieved by 
     * averaging the values on both sides of the boundary. Thus Vout = -vin.
     * 
     */
    
    // Applying u,v and F boundary conditions on every border
    /*
    Here we have something strange happening at the code. When we try to ensure the border positions at the right, it
    shows the left and so otherwise (Also the TOP and Bottom). 
    Hence, here we using the border position that is correct according to the i() andj() that is reversed, in order to get the correct border.

    */

    // Loop through all cells an applying boundary conditions if the cells in located at the border.
    for(auto cell: _cells){
        // Fix code duplication
        if(cell->is_border(border_position::RIGHT)){
            // u on the left
            field.u(cell->i(), cell->j()) = 0;
            // v on left 
            field.v(cell->i(), cell->j()) = - field.v(cell->neighbour(border_position::RIGHT)->i(), cell->j());
            // F on the left 
            field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
        }
        else if(cell->is_border(border_position::LEFT)){
            // u on the right
            field.u(cell->i(), cell->j()) = 0;
            // v  on right 
            field.v(cell->i(), cell->j()) = - field.v(cell->neighbour(border_position::LEFT)->i(), cell->j());
            // F on the right
            field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
        }
        else if(cell->is_border(border_position::TOP)){
            // u 
            field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->neighbour(border_position::TOP)->j());
            // v 
            field.v(cell->i(), cell->j()) = 0;    
            // g
            field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j()); 
        }
    }
}


void MovingWallBoundary::apply(Fields &field) {

    /* The top wall is a moving wall thus, the velocity on the top wall would be the velocity of the wall itself.*/

    // Same aboveLoop through all cells an applying boundary conditions if the cells in located at the border.
    // However this here, we have a moving wall so the wall velocity is accounted for.

    for(auto cell: _cells){
        // std::cout << "Moving wall velocity  = " << _wall_velocity[LidDrivenCavity::moving_wall_id] << std::endl;
        if(cell->is_border(border_position::BOTTOM) && cell->type() == cell_type::MOVING_WALL){
            // v on the top
            field.v(cell->i(), cell->j()) = 0;
            // std::cout << "cell->i() = " << cell->j() << "    cell->neighbour()->i() = " << cell->neighbour(border_position::BOTTOM)->j() << std::endl;
            // u on the top 
            field.u(cell->i(), cell->j()) = 2.0 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(cell->i(), cell->neighbour(border_position::BOTTOM)->j());
            // F on the left 
            field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
        }
    }
}
