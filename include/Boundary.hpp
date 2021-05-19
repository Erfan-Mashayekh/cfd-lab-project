#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"

/**
 * @brief Abstact of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
  public:
    /**
     * @brief Main method to patch the boundary conditons to given field and
     * grid
     *
     * @param[in] Field to be applied
     */
    virtual void apply(Fields &field) = 0;
    virtual ~Boundary() = default;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class FixedWallBoundary : public Boundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature);
    virtual ~FixedWallBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_temperature;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                       std::map<int, double> wall_temperature);
    virtual ~MovingWallBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_velocity;
    std::map<int, double> _wall_temperature;
};

class InflowBoundary : public Boundary {
  public:
    InflowBoundary(std::vector<Cell *> cells, double UIN);
    InflowBoundary(std::vector<Cell *> cells, std::map<int, double> UIN, std::map<int, double> TIN);
    virtual ~InflowBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _UIN;
    std::map<int, double> _TIN;
};

class OutflowBoundary : public Boundary {
  public:
    OutflowBoundary(std::vector<Cell *> cells);
    virtual ~OutflowBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
};

class NoSlipBoundary : public Boundary {
  public:
    NoSlipBoundary(std::vector<Cell *> cells);
    NoSlipBoundary(std::vector<Cell *> cells, double wall_temperature);
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    double _wall_temperature
};

class FreeSlipBoundary : public Boundary {
  public:
    FreeSlipBoundary(std::vector<Cell *> cells);
    FreeSlipBoundary(std::vector<Cell *> cells, double wall_temperature);
    virtual ~FreeSlipBoundary() = default;
    virtual void apply(Fields &field);

  private:
    std::vector<Cell *> _cells;
    double _wall_temperature;
};
