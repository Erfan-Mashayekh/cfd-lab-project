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
        // u
        field.u(cell->i(), cell->j()) = _inlet_velocity_x;
        // v
        field.v(cell->i(), cell->j()) = _inlet_velocity_y;
        // F
        field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
        // G
        field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
        // p
        field.p(cell->i(), cell->j()) = field.p(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j());
    }
}

void OutflowBoundary::apply(Fields &field){

    for (auto const& cell : _cells) {

        // Assume single fluid cell neighbour
        
        std::vector<border_position> border_pos = cell->borders();

        if(border_pos.size() != 1){
            std::cout << "Warning! Outflow can have only single fluid neighbour cells!" << std::endl;
        }
        // u
        field.u(cell->i()-1, cell->j()) = field.u(cell->neighbour(border_pos.at(0))->i()-1, cell->neighbour(border_pos.at(0))->j());
        // v
        field.v(cell->i(), cell->j()) = field.v(cell->neighbour(border_pos.at(0))->i(), cell->neighbour(border_pos.at(0))->j());
        // F
        field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());
        // G
        field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
        // p
        field.p(cell->i(), cell->j()) = _initial_pressure; 
    }
}

void FixedWallBoundary::apply(Fields &field){

    for(auto const& cell: _cells){
        
        std::vector<border_position> border_positions = cell->borders();

        // Set boundary conditions to zero if there are no bordering fluid cells (at inner obstacle cells)
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
            // T
            field.T(cell->i(), cell->j()) = 0.0;

            continue;
        }

        // Set boudndary conditions when there is only one boundary cell
        else if(border_positions.size() == 1){

            switch(border_positions.at(0)){

                case border_position::RIGHT:          
                    // u 
                    field.u(cell->i(), cell->j()) = 0.0;
                    // v 
                    field.v(cell->i(), cell->j()) = - field.v(cell->neighbour(border_position::RIGHT)->i(), cell->j());
                    // F 
                    field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
                    // p
                    field.p(cell->i(), cell->j()) = field.p(cell->neighbour(border_position::RIGHT)->i(), cell->j());

                    break;

                case border_position::LEFT:         
                    // u - Take the staggered grid position into consideration
                    field.u(cell->i()-1, cell->j()) = 0.0;
                    // v 
                    field.v(cell->i(), cell->j()) = - field.v(cell->neighbour(border_position::LEFT)->i(), cell->j());
                    // F - Take the staggered grid position into consideration 
                    field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());
                    // p
                    field.p(cell->i(), cell->j()) = field.p(cell->neighbour(border_position::LEFT)->i(), cell->j());

                    break;

                case border_position::TOP:       
                    // u 
                    field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->neighbour(border_position::TOP)->j());
                    // v 
                    field.v(cell->i(), cell->j()) = 0.0;
                    // G
                    field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
                    // p
                    field.p(cell->i(), cell->j()) = field.p(cell->i(), cell->neighbour(border_position::TOP)->j());

                    break;

                case border_position::BOTTOM:
                    // u 
                    field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->neighbour(border_position::BOTTOM)->j());
                    // v - Take the staggered grid position into consideration 
                    field.v(cell->i(), cell->j()-1) = 0.0;    
                    // G - Take the staggered grid position into consideration
                    field.g(cell->i(), cell->j()-1) = field.v(cell->i(), cell->j()-1);
                    // p
                    field.p(cell->i(), cell->j()) = field.p(cell->i(), cell->neighbour(border_position::BOTTOM)->j());

                    break;
            }
        }

        // Set boudndary conditions when there are two fluid boundary cells
        else if(border_positions.size() == 2){

            if(cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)){ 
                // u 
                field.u(cell->i(), cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()) = 0.0;
                // u
                field.u(cell->i()-1, cell->j()) = - field.u(cell->i()-1, cell->j()+1);
                // v 
                field.v(cell->i(), cell->j()-1) = - field.v(cell->i()+1, cell->j()-1);
                // F 
                field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
                // G
                field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
                // p
                field.p(cell->i(), cell->j()) = 0.5 * ( field.p(cell->neighbour(border_position::TOP)->i(), cell->neighbour(border_position::TOP)->j())
                                                      + field.p(cell->neighbour(border_position::RIGHT)->i(), cell->neighbour(border_position::RIGHT)->j()) );
            }
            else if(cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)){       
                // u
                field.u(cell->i()-1, cell->j()) = 0.0;
                // v
                field.v(cell->i(), cell->j()) = 0.0;
                // u
                field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->j()+1);
                // v 
                field.v(cell->i(), cell->j()-1) = - field.v(cell->i()-1, cell->j()-1);
                // F 
                field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());
                // G
                field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
                // p
                field.p(cell->i(), cell->j()) = 0.5 * ( field.p(cell->neighbour(border_position::TOP)->i(), cell->neighbour(border_position::TOP)->j())
                                                      + field.p(cell->neighbour(border_position::LEFT)->i(), cell->neighbour(border_position::LEFT)->j()) );
            }
            else if(cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)){    
                // u 
                field.u(cell->i(), cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()-1) = 0.0;
                // u
                field.u(cell->i()-1, cell->j()) = - field.u(cell->i()-1, cell->j()-1);
                // v 
                field.v(cell->i(), cell->j()) = - field.v(cell->i()+1, cell->j());
                // F 
                field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
                // G
                field.g(cell->i(), cell->j()-1) = field.v(cell->i(), cell->j()-1);
                // p
                field.p(cell->i(), cell->j()) = 0.5 * ( field.p(cell->neighbour(border_position::BOTTOM)->i(), cell->neighbour(border_position::BOTTOM)->j())
                                                      + field.p(cell->neighbour(border_position::RIGHT)->i(), cell->neighbour(border_position::RIGHT)->j()) );
            }
            else if(cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)){     
                // u 
                field.u(cell->i()-1, cell->j()) = 0.0;
                // v 
                field.v(cell->i(), cell->j()-1) = 0.0;
                // u
                field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->j()-1);
                // v 
                field.v(cell->i(), cell->j()) = - field.v(cell->i()-1, cell->j());
                // F 
                field.f(cell->i()-1, cell->j()) = field.u(cell->i()-1, cell->j());
                // G
                field.g(cell->i(), cell->j()-1) = field.v(cell->i(), cell->j()-1);
                // p
                field.p(cell->i(), cell->j()) = 0.5 * ( field.p(cell->neighbour(border_position::BOTTOM)->i(), cell->neighbour(border_position::BOTTOM)->j())
                                                      + field.p(cell->neighbour(border_position::LEFT)->i(), cell->neighbour(border_position::LEFT)->j()) );
            }

            else{
                // <TOP, BOTTOM> or <LEFT, RIGHT> configurations are not allowed!
                std::cout << "Warning! Forbidden cells found!" << std::endl;
            }         
        }

        else{
            // More than 2 border cells are not allowed
            std::cout << "Warning! Forbidden cells found!" << std::endl; 
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
        field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
        // G on the left 
        field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
        // p
        field.p(cell->i(), cell->j()) = field.p(cell->i(), cell->neighbour(border_position::BOTTOM)->j());
    }
}

void FreeSlipBoundary::apply(Fields &field) {

    for (auto const &cell : _cells) {

        int i = cell->i();
        int j = cell->j();

        std::vector<border_position> border_positions = cell->borders();

        // Set velocity to zero at inner obstacle cells
        // Need help here guys
        if (border_positions.size() == 0) {

            // u
            field.u(cell->i(), cell->j()) = 0.0;
            // v
            field.v(cell->i(), cell->j()) = 0.0;
            // F
            field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
            // G
            field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
            // p
            field.p(cell->i(), cell->j()) = field.p(cell->i(), cell->i() + 1);

            continue;
        }

        // Set boudndary conditions when there is only one boundary cell
        else if (border_positions.size() == 1) {

            switch (border_positions.at(0)) {

            case border_position::RIGHT:
                field.v(i, j) = field.v(i + 1, j);
                field.u(i, j) = 0;
                field.f(i, j) = field.u(i, j);
                field.p(i, j) = field.p(i + 1, j);
                break;

            case border_position::LEFT:
                field.v(i, j) = field.v(i + 1, j);
                field.u(i - 1, j) = 0;
                field.f(i - 1, j) = field.u(i - 1, j);
                field.p(i, j) = field.p(i + 1, j);
                break;

            case border_position::TOP:
                field.u(i, j) = field.u(i, j + 1);
                field.v(i, j) = 0;
                field.g(i, j) = field.v(i, j);
                field.p(i, j) = field.p(i, j + 1);
                break;

            case border_position::BOTTOM:
                field.u(i, j) = field.u(i, j - 1);
                field.v(i, j - 1) = 0;
                field.g(i, j - 1) = field.v(i, j - 1);
                field.p(i, j) = field.p(i, j - 1);
                break;
            }
        }

        // Set boudndary conditions when there are two fluid boundary cells
        else if (border_positions.size() == 2) {

            if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
                field.u(i, j) = 0;
                field.v(i, j) = 0;
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                field.p(i,j) = (field.p(i,j+1)*field.p(i+1,j))*0.5;

            } else if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
                field.u(i, j) = 0;
                field.v(i, j) = 0;
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                field.p(i,j) = (field.p(i,j+1)*field.p(i-1,j))*0.5;

            } else if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
                field.u(i, j) = 0;
                field.v(i, j) = 0;
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                field.p(i,j) = (field.p(i,j-1)*field.p(i+1,j))*0.5;

            } else if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
                field.u(i, j) = 0;
                field.v(i, j) = 0;
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                field.p(i,j) = (field.p(i,j-1)*field.p(i-1,j))*0.5;

            }

            else {
                // <TOP, BOTTOM> or <LEFT, RIGHT> configurations are not allowed!
                std::cout << "Warning! Forbidden cells found!" << std::endl;
            }
        }

        else {
            // More than 2 border cells are not allowed
            std::cout << "Warning! Forbidden cells found!" << std::endl;
        }
    }
}

void FixedWallBoundary::apply_temperature(Fields &field) {
    for (auto const &cell : _cells) {
        int i = cell->i();
        int j = cell->j();

        std::vector<border_position> border_positions = cell->borders();

        if (_wall_temperature[cell->wall_id()] == -1) {
            // Neumann boundary condition
            if (border_positions.size() == 1) {
                switch (border_positions.at(0)) {
                    case border_position::RIGHT:
                        field.T(i, j) = field.T(i + 1, j);
                        break;

                    case border_position::LEFT:
                        field.T(i, j) = field.T(i - 1, j);
                        break;

                    case border_position::TOP:
                        field.T(i, j) = field.T(i, j + 1);
                        break;

                    case border_position::BOTTOM:
                        field.T(i, j) = field.T(i, j - 1);
                        break;
                }
            } else if (border_positions.size() == 2) {

                if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
                    // NorthEast
                    field.T(i, j) = 0.5 * (field.T(i, j + 1) + field.T(i + 1, j));

                } else if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
                    // Northwest
                    field.T(i, j) = 0.5 * (field.T(i, j + 1) + field.T(i - 1, j));

                } else if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
                    // Southeast
                    field.T(i, j) = 0.5 * (field.T(i, j - 1) + field.T(i + 1, j));

                } else if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
                    // Southwest
                    field.T(i, j) = 0.5 * (field.T(i, j - 1) + field.T(i - 1, j));
                }
            }
        }
        else {
            // Dirichlet Boundary Condition
            if (border_positions.size() == 1) {
                switch (border_positions.at(0)) {

                    case border_position::RIGHT:
                        field.T(i, j) = 2 * _wall_temperature[cell->wall_id()] - field.T(i + 1, j);
                        break;

                    case border_position::LEFT:
                        field.T(i, j) = 2 * _wall_temperature[cell->wall_id()] - field.T(i - 1, j);
                        break;

                    case border_position::TOP:
                        field.T(i, j) = 2 * _wall_temperature[cell->wall_id()] - field.T(i, j + 1);
                        break;

                    case border_position::BOTTOM:
                        field.T(i, j) = 2 * _wall_temperature[cell->wall_id()] - field.T(i, j - 1);
                        break;
                }
            } else if (border_positions.size() == 2) {

                if (cell->is_border(border_position::TOP) && cell->is_border(border_position::RIGHT)) {
                    // NorthEast
                    field.T(i, j) = 2 * _wall_temperature[cell->wall_id()] - 0.5 * (field.T(i + 1, j) + field.T(i, j + 1));

                } else if (cell->is_border(border_position::TOP) && cell->is_border(border_position::LEFT)) {
                    // Northwest
                    field.T(i, j) = 2 * _wall_temperature[cell->wall_id()] - 0.5 * (field.T(i - 1, j) + field.T(i, j + 1));

                } else if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::RIGHT)) {
                    // Southeast
                    field.T(i, j) = 2 * _wall_temperature[cell->wall_id()] - 0.5 * (field.T(i + 1, j) + field.T(i, j - 1));
                } else if (cell->is_border(border_position::BOTTOM) && cell->is_border(border_position::LEFT)) {
                    // Southwest
                    field.T(i, j) = 2 * _wall_temperature[cell->wall_id()] - 0.5 * (field.T(i - 1, j) + field.T(i, j - 1));
                }
            }
        }
    }
}

void InflowBoundary::apply_temperature(Fields &field){

    for (auto const& cell: _cells) {
        // Set T    
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

void FreeSlipBoundary::apply_temperature(Fields &field) {}

void MovingWallBoundary::apply_temperature(Fields &field) {}