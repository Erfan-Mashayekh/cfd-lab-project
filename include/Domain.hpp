#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include "Enums.hpp"

/**
 * @brief Data structure that holds geometrical information
 * necessary for decomposition.
 *
 */
struct Domain {
    /// Minimum x index
    int imin{-1};
    /// Maximum x index
    int imax{-1};

    /// Minimum y index
    int jmin{-1};
    /// Maximum y index
    int jmax{-1};

    /// Cell length
    double dx{-1.0};
    /// Cell height
    double dy{-1.0};

    /// Number of cells in x direction
    int size_x{-1};
    /// Number of cells in y direction
    int size_y{-1};

    /// row number of the subdomain
    int x_proc{-1};
    /// column number of the subdomain
    int y_proc{-1};

    /// Number of cells in x direction, not-decomposed
    int domain_size_x{-1};
    /// Number of cells in y direction, not-decomposed
    int domain_size_y{-1};

    /// Number of processes in x direction
    int domain_iproc{-1};
    /// Number of processes in y direction
    int domain_jproc{-1};

};

#endif // DOMAIN_HPP