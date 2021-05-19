#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

// void FixedWallBoundary::apply(Fields &field){
//     /**
//      * All boundaries except the top are fixed wall or no slip.
//      *
//      * U on the left and right are set to 0 and v at the bottom is set to 0.
//      *
//      * No v-values lie on the vertical boundaries and no u-values lie on
//      * the horizontal boundaries, the boundary value zero is achieved by
//      * averaging the values on both sides of the boundary. Thus Vout = -vin.
//      *
//      */

//     // Applying u,v and F boundary conditions on every border
//     /*
//     Here we have something strange happening at the code. When we try to ensure the border positions at the right, it
//     shows the left and so otherwise (Also the TOP and Bottom).
//     Hence, here we using the border position that is correct according to the i() andj() that is reversed, in order
//     to get the correct border.

//     */

//     // Loop through all cells an applying boundary conditions if the cells in located at the border.
//     for(auto cell: _cells){
//         // Fix code duplication
//         if(cell->is_border(border_position::RIGHT)){
//             // u on the left
//             field.u(cell->i(), cell->j()) = 0;
//             // v on left
//             field.v(cell->i(), cell->j()) = - field.v(cell->neighbour(border_position::RIGHT)->i(), cell->j());
//             // F on the left
//             field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
//         }
//         else if(cell->is_border(border_position::LEFT)){
//             // u on the right
//             field.u(cell->i(), cell->j()) = 0;
//             // v  on right
//             field.v(cell->i(), cell->j()) = - field.v(cell->neighbour(border_position::LEFT)->i(), cell->j());
//             // F on the right
//             field.f(cell->i(), cell->j()) = field.u(cell->i(), cell->j());
//         }
//         else if(cell->is_border(border_position::TOP)){
//             // u
//             field.u(cell->i(), cell->j()) = - field.u(cell->i(), cell->neighbour(border_position::TOP)->j());
//             // v
//             field.v(cell->i(), cell->j()) = 0;
//             // g
//             field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
//         }
//     }
// }

// list of fixed wall cell cells are passed with the function grid.listoffixedwallcells (needs to be synchronize with
// Erfan's)

void FixedWallBoundary::apply(Fields &field) {

    for (auto thiscell : grid.fixed_wall_cells()) {
        no_slip(thiscell, field, 0);
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

// Same reason with fixed wall above, we  need to make this applicable to any wall.

// void MovingWallBoundary::apply(Fields &field) {

//     /* The top wall is a moving wall thus, the velocity on the top wall would be the velocity of the wall itself.*/

//     // Same aboveLoop through all cells an applying boundary conditions if the cells in located at the border.
//     // However this here, we have a moving wall so the wall velocity is accounted for.

//     for(auto cell: _cells){
//         // std::cout << "Moving wall velocity  = " << _wall_velocity[LidDrivenCavity::moving_wall_id] << std::endl;
//         if(cell->is_border(border_position::BOTTOM) && cell->type() == cell_type::MOVING_WALL){
//             // v on the top
//             field.v(cell->i(), cell->j()) = 0;
//             // std::cout << "cell->i() = " << cell->j() << "    cell->neighbour()->i() = " <<
//             cell->neighbour(border_position::BOTTOM)->j() << std::endl;
//             // u on the top
//             field.u(cell->i(), cell->j()) = 2.0 * _wall_velocity[LidDrivenCavity::moving_wall_id] -
//             field.u(cell->i(), cell->neighbour(border_position::BOTTOM)->j());
//             // F on the left
//             field.g(cell->i(), cell->j()) = field.v(cell->i(), cell->j());
//         }
//     }
// }

// list of moving wall cell cells are passed with the function grid.listofmovingwallcells (needs to be synchronize with
// Erfan's)
void MovingWallBoundary::apply(Fields &field) {

    for (auto thiscell : grid.moving_wall_cells()) {
        no_slip(thiscell, field, _wall_velocity);
    }
}

InflowBoundary::InflowBoundary(std::vector<Cell *> cells, std::map<int, double> UIN, std::map<int, double> TIN)
    : _cells(cells), _UIN(UIN), _TIN(TIN) {}

InflowBoundary::InflowBoundary(std::vector<Cell *> cells, std::map<int, double> UIN) : _cells(cells), _UIN(UIN) {}

// Inflow and outflow cells are passed with the function grid.listinflowcells / outflowcells (needs to be synchronize
// with Erfan's)
void InflowBoundary::apply(Fields &field) {
    // Apply Boundary Conditions
    for (auto thiscell : grid.inflow_cells()) {
        int i = thiscell->i();
        int j = thiscell->j();
        field.u(i, j) = _UIN;
        field.v(i, j) = 0;
        field.t(i, j) = _TIN;
    }

    OutflowBoundary::OutflowBoundary(std::vector<Cell *> cells) : _cells(cells) {}

    void OutflowBoundary::apply(Fields & field) {
        // Apply Bondary Conditions
        for (auto thiscell : grid.outflow_cells()) {
            int i = thiscell->i();
            int j = thiscell->j();
            field.u(i - 1, j) = field.u(i, j);
            field.v(i - 1, j) = field.v(i, j);
        }
    }

    NoSlipBoundary::NoSlipBoundary(std::vector<Cell *> cells) : _cells(cells) {}

    NoSlipBoundary::NoSlipBoundary(std::vector<Cell *> cells, double wall_temperature)
        : _cells(cells), _wall_temperature(wall_temperature) {}

    void NoSlipBoundary::apply(Fields & field, Grid & grid) {

        for (auto &thiscell : grid.no_slip_cells()) {
            no_slip(thiscell, field, 0.0);
        }
    }

    FreeSlipBoundary::FreeSlipBoundary(std::vector<Cell *> cells) : _cells(cells) {}

    FreeSlipBoundary::FreeSlipBoundary(std::vector<Cell *> cells, double wall_temperature)
        : _cells(cells), _wall_temperature(wall_temperature) {}

    // grid.listoffreeslipcells needs to be synchronize with Erfan's Grid.cpp it call for the cell that is a free slip
    // Also needs a function to check whether this current cell has a neighbour on the north, west, east or south (that
    // give a bolean output)
    void FreeSlipBoundary::apply(Fields & field, Grid & grid) {

        for (auto &thiscell : grid.free_slipcells()) {
            int i = thiscell->i();
            int j = thiscell->j();
            if (thiscell->CheckNeighbour().North && !thiscell->CheckNeighbour().West &&
                !thiscell->CheckNeighbour().East) {
                field.u(i - 1, j) = field.u(i - 1, j + 1);
                field.u(i, j) = field.u(i, j + 1);
                field.v(i - 1, j) = 0;
                field.g(i, j) = field.v(i, j);
                // field.p(i,j) = field.p(i,j+1);
            }

            if (thiscell->CheckNeighbour().East && !thiscell->CheckNeighbour().South &&
                !thiscell->CheckNeighbour().North) {
                field.v(i, j - 1) = field.v(i + 1, j - 1);
                field.v(i, j) = field.v(i + 1, j);
                field.u(i + 1, j) = 0;
                field.f(i + 1, j) = field.u(i + 1, j);
                // field.p(i,j) = field.p(i+1,j);
            }

            if (thiscell->CheckNeighbour().South && !thiscell->CheckNeighbour().West &&
                !thiscell->CheckNeighbour().East) {
                field.u(i - 1, j) = field.u(i - 1, j - 1);
                field.u(i, j) = field.u(i, j - 1);
                field.v(i - 1, j) = 0;
                field.g(i, j) = field.v(i, j);
                // field.p(i,j) = field.p(i,j-1);
            }

            if (thiscell->CheckNeighbour().West && !thiscell->CheckNeighbour().South &&
                !thiscell->CheckNeighbour().North) {
                field.v(i, j - 1) = field.v(i - 1, j - 1);
                field.v(i, j) = field.v(i - 1, j);
                field.u(i - 1, j) = 0;
                field.f(i - 1, j) = field.u(i - 1, j);
                // field.p(i,j) = field.p(i-1,j);
            }

            if (thiscell->CheckNeighbour().East && thiscell->CheckNeighbour().North) {
                field.u(i, j) = 0;
                field.u(i - 1, j) = field.u(i - 1, j + 1);
                field.v(i, j) = 0;
                field.v(i, j - 1) = field.v(i + 1, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                // field.p(i,j) = (field.p(i,j+1)*field.p(i+1,j))*0.5;
            }

            if (thiscell->CheckNeighbour().East && thiscell->CheckNeighbour().South) {
                field.u(i, j) = 0;
                field.u(i - 1, j) = field.u(i - 1, j - 1);
                field.v(i, j) = 0;
                field.v(i, j - 1) = field.v(i + 1, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                // field.p(i,j) = (field.p(i,j-1)*field.p(i+1,j))*0.5;
            }

            if (thiscell->CheckNeighbour().West && thiscell->CheckNeighbour().South) {
                field.u(i, j) = 0;
                field.u(i - 1, j) = field.u(i - 1, j - 1);
                field.v(i, j) = 0;
                field.v(i, j - 1) = field.v(i - 1, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                // field.p(i,j) = (field.p(i,j-1)*field.p(i-1,j))*0.5;
            }

            if (thiscell->CheckNeighbour().West && thiscell->CheckNeighbour().North) {
                field.u(i, j) = 0;
                field.u(i - 1, j) = field.u(i - 1, j + 1);
                field.v(i, j) = 0;
                field.v(i, j - 1) = field.v(i - 1, j - 1);
                field.f(i, j) = field.u(i, j);
                field.g(i, j) = field.v(i, j);
                // field.p(i,j) = (field.p(i,j+1)*field.p(i-1,j))*0.5;
            }
        }
    }

    // No slip is applied differently depends on the position of the cell.
    void Boundary::no_slip(Cell * thiscell, Fields & field, double wall_velocity) {
        int i = thiscell->i();
        int j = thiscell->j();

        if (thiscell->CheckNeighbour().North && !thiscell->CheckNeighbour().West && !thiscell->CheckNeighbour().East)

        {
            field.u(i, j) = 2 * wall_velocity - field.u(i, j + 1);
            field.v(i, j) = 0;
        }

        if (thiscell->CheckNeighbour().East && !thiscell->CheckNeighbour().South && !thiscell->CheckNeighbour().North) {
            field.v(i, j) = 2 * wall_velocity - field.v(i + 1, j);
            field.u(i, j) = 0;
        }

        if (thiscell->CheckNeighbour().South && !thiscell->CheckNeighbour().West && !thiscell->CheckNeighbour().East) {
            field.u(i, j) = 2 * wall_velocity - field.u(i, j - 1);
            field.v(i, j - 1) = 0;
        }

        if (thiscell->CheckNeighbour().West && !thiscell->CheckNeighbour().South && !thiscell->CheckNeighbour().North) {
            field.v(i, j) = 2 * wall_velocity - field.v(i - 1, j);
            field.u(i - 1, j) = 0;
        }

        if (thiscell->CheckNeighbour().East && thiscell->CheckNeighbour().North) {
            field.u(i, j) = 0;
            field.v(i, j) = 0;
            field.f(i, j) = field.u(i, j);
            field.g(i, j) = field.v(i, j);
        }

        if (thiscell->CheckNeighbour().East && thiscell->CheckNeighbour().South) {
            field.u(i, j) = 0;
            field.v(i, j - 1) = wall_velocity;
            field.f(i, j) = field.u(i, j);
        }

        if (thiscell->CheckNeighbour().West && thiscell->CheckNeighbour().South) {
            field.u(i - 1, j) = 0;
            field.v(i, j - 1) = 0;
        }

        if (thiscell->CheckNeighbour().West && thiscell->CheckNeighbour().North) {
            field.u(i - 1, j) = 0;
            field.v(i, j) = 0;
            field.g(i, j) = field.v(i, j);
        }
    }
