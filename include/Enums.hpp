#pragma once

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention, which is:
// 0: fluid, 3: fixed wall, 4: moving wall
namespace LidDrivenCavity {
const int moving_wall_id = 8;
const int fixed_wall_id = 4;
const double wall_velocity = 1.0;
} // namespace LidDrivenCavity

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

namespace border {
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;
} // namespace border

enum class cell_type {  /** PGM INDEX **/
    FLUID,              /*     0       */
    INFLOW,             /*     1       */
    OUTFLOW,            /*     2       */
    FIXED_WALL,         /*     3-7     */
    MOVING_WALL,        /*     8       */
    FREE_SLIP_WALL,     /*     9       */
    DEFAULT             /*             */
};                      /***************/
