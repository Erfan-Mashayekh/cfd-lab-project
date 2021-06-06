// #pragma once

// #include <memory>
// #include <string>
// #include <vector>

// #include "Fields.hpp" 

// class Communication {
//     public:
//         Communication(int rank); //iprod and jproc is the number of process per y or x
//         static void init_parallel(int processers, int _size, int my_rank);
//         static void finalize();
//         static void communicate(Fields &field);
//         static double reduce_min(double value);
//         static double reduce_sum(double value);

// };