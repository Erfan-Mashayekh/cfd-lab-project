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
    virtual void apply_temperature(Fields &field) = 0;
    virtual ~Boundary() = default;
    
};

/**
 * @brief Boundary condition for the fluid inlet 
 * Dirichlet for velocities specified by inlet_velocity_x and inlet_velocity_y
 */
class InflowBoundary : public Boundary {
  public:
    InflowBoundary(std::vector<Cell *> cells, double inlet_velocity_x, double inlet_velocity_y);
    InflowBoundary(std::vector<Cell *> cells, double inlet_velocity_x, double inlet_velocity_y,
                                                                       double inlet_temperature);
    virtual ~InflowBoundary() = default;
    virtual void apply(Fields &field) override;
    virtual void apply_temperature(Fields &field) override;

  private:
    std::vector<Cell *> _cells;
    double _inlet_velocity_x;
    double _inlet_velocity_y;
    double _inlet_temperature;
};

/**
 * @brief Boundary condition for the fluid outlet 
 * Neumann boundary condition for velocities 
 */
class OutflowBoundary : public Boundary {
  public:
    OutflowBoundary(std::vector<Cell *> cells, double initial_pressure);

    virtual ~OutflowBoundary() = default;
    virtual void apply(Fields &field) override;
    virtual void apply_temperature(Fields &field) override;

  private:
    std::vector<Cell *> _cells;
    double _initial_pressure;
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
    virtual void apply(Fields &field) override;
    virtual void apply_temperature(Fields &field) override;

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
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                                  std::map<int, double> wall_temperature);
    virtual ~MovingWallBoundary() = default;
    virtual void apply(Fields &field) override;
    virtual void apply_temperature(Fields &field) override;

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_velocity;
    std::map<int, double> _wall_temperature;
};

/**
 * @brief Free slip boundary 
 * Neumann boundary conditions for velocity
 */
class FreeSlipBoundary : public Boundary {
  public:
    FreeSlipBoundary(std::vector<Cell *> cells);
    FreeSlipBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature);
    virtual ~FreeSlipBoundary() = default;
    virtual void apply(Fields &field) override;
    virtual void apply_temperature(Fields &field) override;

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_temperature;
};


