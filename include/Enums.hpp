#pragma once

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
