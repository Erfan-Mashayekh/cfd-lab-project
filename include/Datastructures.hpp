#ifndef DATASTRUCTURES_HPP
#define DATASTRUCTURES_HPP

#include <vector>
#include <stdexcept>

/**
 * @brief General 2D data structure around std::vector, in column
 * major format.
 *
 */
template <typename T> class Matrix {

  public:
    Matrix<T>() = default;

    /**
     * @brief Constructor with initial value
     *
     * @param[in] number of elements in x direction
     * @param[in] number of elements in y direction
     * @param[in] initial value for the elements
     *
     */
    Matrix<T>(int i_max, int j_max, double init_val) : _imax(i_max), _jmax(j_max) {
        _container.resize(i_max * j_max);
        std::fill(_container.begin(), _container.end(), init_val);
    }

    /**
     * @brief Constructor without an initial value.
     *
     * @param[in] number of elements in x direction
     * @param[in] number of elements in y direction
     *
     */
    Matrix<T>(int i_max, int j_max) : _imax(i_max), _jmax(j_max) { _container.resize(i_max * j_max); }

    /**
     * @brief Element access and modify using index
     *
     * @param[in] x index
     * @param[in] y index
     * @param[out] reference to the value
     */
    T &operator()(int i, int j) {
        try{
            return _container.at(_imax * j + i);
        } catch (...) {
            if(_imax * j + i < 0){
                throw std::out_of_range("Lower bound violated " + std::to_string(i) + " " + std::to_string(j));
            } else {
                throw std::out_of_range("Upper bound violated " + std::to_string(i) + " " + std::to_string(j));
            }
        }
    }

    /**
     * @brief Element access using index
     *
     * @param[in] x index
     * @param[in] y index
     * @param[out] value of the element
     */
    T operator()(int i, int j) const {
        try{
            return _container.at(_imax * j + i);
        } catch (...) {
            if(_imax * j + i < 0){
                throw std::out_of_range("Lower bound violated " + std::to_string(i) + " " + std::to_string(j));
            } else {
                throw std::out_of_range("Upper bound violated " + std::to_string(i) + " " + std::to_string(j));
            }
        }
    }

    /**
     * @brief Pointer representation of underlying data
     *
     * @param[out] pointer to the beginning of the vector
     */
    const T *data() const { return _container.data(); }

    /**
     * @brief Pointer representation of underlying data
     *
     * @param[out] pointer to the beginning of the vector
     */
    std::vector<T> *container() { return &_container; }

    /**
     * @brief Access of the size of the structure
     *
     * @param[out] size of the data structure
     */
    size_t size() const { return _container.size(); }

    /// get the given row of the matrix
    std::vector<T> get_row(int row) {
        std::vector<T> row_data(_imax, -1);
        for (int i = 0; i < _imax; ++i) {
            row_data.at(i) = _container.at(i + _imax * row);
        }
        return row_data;
    }

    /// get the given column of the matrix
    std::vector<T> get_col(int col) {
        std::vector<T> col_data(_jmax, -1);
        for (int i = 0; i < _jmax; ++i) {
            col_data.at(i) = _container.at(col + i * _imax);
        }
        return col_data;
    }

    /// set the given column of matrix to given vector
    void set_col(const std::vector<T> &vec, int col) {
        for (int i = 0; i < _jmax; ++i) {
            _container.at(col + i * _imax) = vec.at(i);
        }
    }

    /// set the given row of matrix to given vector
    void set_row(const std::vector<T> &vec, int row) {
        for (int i = 0; i < _imax; ++i) {
            _container.at(i + row * _imax) = vec.at(i);
        }
    }

    /// get the number of elements in x direction
    size_t imax() const { return _imax; }

    /// get the number of elements in y direction
    size_t jmax() const { return _jmax; }

  private:
    /// Number of elements in x direction
    int _imax;
    /// Number of elements in y direction
    int _jmax;

    /// Data container
    std::vector<T> _container;
};

#endif // DATASTRUCTURES_HPP
