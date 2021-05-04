#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

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
       

    for(auto cell: _cells){

        if(cell->is_border(border_position::RIGHT)){
            // u on the left
            field.u(cell->i(), cell->j()) = 0;
            // v on left 
            field.v(cell->i(), cell->j()) = - field.v(cell->i(), cell->j());
            // F on the left 
            field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
            
            //std::cout << " (i , j) " << cell->i() << " " << cell->j() << std::endl;

        }
        if(cell->is_border(border_position::LEFT)){
            // u on the right
            field.u(cell->i(), cell->j()) = 0;
            // v  on right 
            field.v(cell->i(), cell->j()) = - field.v(cell->i(), cell->j());
            // F on the right
            field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());

            //std::cout << " (i , j) " << cell->i() << " " << cell->j() << std::endl;
        }
        else if(cell->is_border(border_position::BOTTOM)){
            // u 
            field.u(cell->i(),0) = - field.u(cell->i(),1);
            // v 
            field.v(cell->i(),0) = 0;    
            // g
            field.g(cell->i(),0) = field.v(cell->i(),0); 
        }
    }


    // // Assigning a maximum grid max
    // int imax = _cells[_cells.size()-1]->i()-1;
    // int jmax = _cells[_cells.size()-1]->j()-1;

    // //Left and right
    // for (int j=1;j<=jmax;j++){
    //     // u on the left
    //     field.u(0,j) = 0;
    //     // u on the right
    //     field.u(imax,j) = 0;

    //     // v on left 
    //     field.v(0,j) = - field.v(1,j);
    //     // v  on right 
    //     field.v(imax+1,j) = - field.v(imax,j);

    //     // F on the left 
    //     field.f(0,j) = field.u(0,j);
    //     // F on the right
    //     field.f(imax,j) = field.u(imax,j);
    // }

    // // Bottom
    // for (int i=1;i<=imax;i++){
    //     // u 
    //     field.u(i,0) = - field.u(i,1);
    //     // v 
    //     field.v(i,0) = 0;    
    //     // g
    //     field.g(i,0) = field.v(i,0); 
    // }
}


MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {

    // Assigning a maximum grid max
    int imax = _cells[_cells.size()-1]->i()-1;
    int jmax = _cells[_cells.size()-1]->j()-1;

    /* The top wall is a moving wall thus, the velocity on the top wall would be the velocity of the wall itself.
    
    */
   for (int i=1;i<=imax;i++){
        // u 
        field.u(i,jmax + 1)  = 2.0 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i,jmax); 
        // v  
        field.v(i,jmax) = 0;
        // g
        field.g(i,jmax) = field.v(i,jmax);  
    }

}