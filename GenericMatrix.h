/****************************************************************************
**
** Copyright (C) 2019 Shaoguang. All rights reserved.
** 
** GenericMatrix, a generic template matrix class,
** the matrix elements are managed by a one-dimension pointer.
** 
** Licensed under the Apache License, Version 2.0 (the "License");
** you may not use this file except in compliance with the License.
** You may obtain a copy of the License at
** 
**     https://www.apache.org/licenses/LICENSE-2.0
** 
** Unless required by applicable law or agreed to in writing, software
** distributed under the License is distributed on an "AS IS" BASIS,
** WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
** See the License for the specific language governing permissions and
** limitations under the License.
**
****************************************************************************/

#ifndef __GERNERICMATRIX_H__
#define __GERNERICMATRIX_H__

#include <ostream>
#include <utility>      // std::move,   std::pair
#include <algorithm>    // std::fill_n, std::copy_n

/*!
    \class GenericMatrix
    \brief The GenericMatrix class defines a generic template matrix class,
    the matrix elements are managed by a one-dimension pointer.

    The GenericMatrix template has one parameter:
    \li \b Elem Element type that is visible to users of the class.
*/
template<typename Elem = float>
class GenericMatrix
{
public:
    using size_type = std::size_t;
    using value_type = Elem;
    using iterator = value_type *;
    using const_iterator = const iterator;

    GenericMatrix();
    GenericMatrix(size_type row, size_type col);
    GenericMatrix(size_type row, size_type col, const Elem &initialValue);
    GenericMatrix(size_type row, size_type col, Elem *data);
    GenericMatrix(const GenericMatrix &other);
    GenericMatrix(GenericMatrix &&other) noexcept;
    ~GenericMatrix();
    GenericMatrix &operator=(const GenericMatrix &other);
    GenericMatrix &operator=(GenericMatrix &&other) noexcept;

    iterator begin() noexcept;
    const_iterator begin() const noexcept;
    const_iterator cbegin() const noexcept;
    iterator end() noexcept;
    const_iterator end() const noexcept;
    const_iterator cend() const noexcept;

    size_type rows() const;
    size_type columns() const;
    size_type size() const;

    Elem *data() noexcept;
    const Elem *data() const noexcept;
    const Elem *constData() const noexcept;
    Elem *data(size_type row, size_type col) noexcept;
    const Elem *data(size_type row, size_type col) const noexcept;
    const Elem *constData(size_type row, size_type col) const noexcept;
    Elem *dataAt(size_type row, size_type col) noexcept;
    const Elem *dataAt(size_type row, size_type col) const noexcept;
    const Elem *constDataAt(size_type row, size_type col) const noexcept;

    Elem &operator()(size_type row, size_type col);
    const Elem &operator()(size_type row, size_type col) const;
    Elem &at(size_type row, size_type col);
    const Elem &at(size_type row, size_type col) const;

    bool empty() const;
    bool isValid() const;
    bool isIdentity() const;
    inline bool isHomomorphicTo(const GenericMatrix<Elem> &m);
    static inline bool isHomomorphic(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);

    void resize(size_type row, size_type col);
    void resize(size_type row, size_type col, const Elem &initialValue);
    void reset(size_type row, size_type col, Elem *data);
    void setToIdentity();
    void fill(const Elem &value);
    void swap(GenericMatrix<Elem> &other);
    GenericMatrix<Elem> reverse() const;
    GenericMatrix<Elem> transposed() const;
    void doHadamardProduct(const GenericMatrix<Elem> &m);
    static GenericMatrix<Elem> hadamardProduct(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);
    GenericMatrix<Elem> roi(size_type row1, size_type col1, size_type row2, size_type col2);

    GenericMatrix<Elem> &operator+=(const GenericMatrix<Elem> &m);
    GenericMatrix<Elem> &operator-=(const GenericMatrix<Elem> &m);
    GenericMatrix<Elem> &operator*=(const GenericMatrix<Elem> &m);
    GenericMatrix<Elem> &operator*=(const Elem &factor);
    GenericMatrix<Elem> &operator/=(const Elem &divisor);
    bool operator==(const GenericMatrix<Elem> &m) const;
    bool operator!=(const GenericMatrix<Elem> &m) const;

    template<typename __ElemDTo, typename __Elem>
    friend GenericMatrix<__ElemDTo> matrix_cast(const GenericMatrix<__Elem> &m);
    template<typename __Elem>
    friend GenericMatrix<__Elem> operator+(const GenericMatrix<__Elem> &m1, const GenericMatrix<__Elem> &m2);
    template<typename __Elem>
    friend GenericMatrix<__Elem> operator-(const GenericMatrix<__Elem> &m1, const GenericMatrix<__Elem> &m2);
    template<typename __Elem>
    friend GenericMatrix<__Elem> operator*(const GenericMatrix<__Elem> &m1, const GenericMatrix<__Elem> &m2);
    template<typename __Elem>
    friend GenericMatrix<__Elem> operator*(const GenericMatrix<__Elem> &matrix, const __Elem &factor);
    template<typename __Elem>
    friend GenericMatrix<__Elem> operator*(const __Elem &factor, const GenericMatrix<__Elem> &matrix);
    template<typename __Elem>
    friend GenericMatrix<__Elem> operator/(const GenericMatrix<__Elem> &matrix, const __Elem &divisor);
    template<typename __Elem>
    friend GenericMatrix<__Elem> operator-(const GenericMatrix<__Elem> &matrix);
    template<typename __Elem>
    friend std::ostream &operator<<(std::ostream &os, const GenericMatrix<__Elem> &matrix);

private:
    void alloc();
    void alloc(size_type row, size_type col);
    void free();
    void setSize(size_type row, size_type col);
    size_type computeOffset(size_type row, size_type col) const;
    bool containsIndex(size_type row, size_type col) const;

private:
    size_type m_rows;
    size_type m_cols;
    Elem *m_data;
};


/*!
    \typedef GenericMatrix::size_type

    Type for declaring the matrix's rows and columns.
*/

/*!
    \typedef GenericMatrix::value_type

    The value type of the template parameter \a Elem.
*/

/*!
    \typedef GenericMatrix::iterator

    An STL-style non-const iterator for GenericMatrix. 
*/

/*!
    \typedef GenericMatrix::const_iterator

    An STL-style const iterator for GenericMatrix.
*/

/*!
    Constructs a invalid matrix.

    \sa isValid()
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix()
    : m_rows(0), m_cols(0), m_data(nullptr)
{
}

/*!
    Constructs a \a row x \a col matrix without initializing the contents.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col)
    : m_rows(row), m_cols(col), m_data(nullptr)
{
    alloc();
}

/*!
    Constructs a \a row x \a col matrix and set the element pointer to \a data.

    \note \a data pointer will be free when destroy this matrix.

    \sa reset()
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col, Elem *data)
    : m_rows(row), m_cols(col), m_data(data)
{
}

/*!
    Constructs a \a row x \a col matrix and initialize all values with \a initialValue.

    \sa fill()
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col, const Elem &initialValue)
    : m_rows(row), m_cols(col), m_data(nullptr)
{
    alloc();
    fill(initialValue);
}

/*!
    Destroys the matrix.
*/
template<typename Elem>
GenericMatrix<Elem>::~GenericMatrix()
{
    free();
}
    
/*!
    \internal

    Alloc the matrix memory from \a m_rows and \a m_cols.
*/
template<typename Elem>
void GenericMatrix<Elem>::alloc()
{
    if (m_rows > 0 && m_cols > 0) {
        m_data = new Elem[m_rows * m_cols];
    }
}

/*!
    \internal

    Alloc the matrix memory from \a row and \a col.
*/
template<typename Elem>
void GenericMatrix<Elem>::alloc(size_type row, size_type col)
{
    if (row > 0 && col > 0) {
        m_data = new Elem[row * col];
    }
}

/*!
    \internal

    Frees the matrix memory.
*/
template<typename Elem>
void GenericMatrix<Elem>::free()
{
    if (m_data)  {
        delete[] m_data; 
        m_data = nullptr;
    }
}

/*!
    \internal

    Sets this matrix's rows and columns with parameter \a row and \a col.
*/
template<typename Elem>
void GenericMatrix<Elem>::setSize(size_type row, size_type col)
{
    m_rows = row;
    m_cols = col;
}

/*!
    \internal

    Compute internal data offset by \a row index and \a col index.
*/
template<typename Elem>
typename GenericMatrix<Elem>::size_type GenericMatrix<Elem>::computeOffset(size_type row, size_type col) const
{
    return row * m_cols + col;
}

/*!
    \internal

    Returns true if \a row index is less than this matrix rows 
    and \a col index is less than this matrix cols, otherwise return false.
*/
template<typename Elem>
bool GenericMatrix<Elem>::containsIndex(size_type row, size_type col) const
{
    return (row < m_rows && col < m_cols);
}

/*!
    Constructs a copy of \a other.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(const GenericMatrix &other)
    : m_rows(other.m_rows), m_cols(other.m_cols), m_data(nullptr)
{
    alloc();
    std::copy_n(other.data(), other.size(), m_data);
}

/*!
    Move-constructs a GenericMatrix instance, making it point at the same object that \a other was pointing to.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(GenericMatrix &&other) noexcept
    : m_rows(other.m_rows), m_cols(other.m_cols), m_data(other.m_data)
{
    other.m_data = nullptr;
    other.setSize(0, 0);
}

/*!
    Assigns \a other to this matrix and returns a reference to this matrix.
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator=(const GenericMatrix &other)
{
    if (this == &other) {
        return *this;
    }

    resize(other.rows(), other.columns());
    std::fill_n(other.m_data, other.size(), m_data);

    return *this;
}

/*!
    Move-assigns \a other to this GenericMatrix instance.
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator=(GenericMatrix &&other) noexcept
{
    if (this == &other) {
        return *this;
    }
    free();

    setSize(other.m_rows, other.m_cols);
    other.setSize(0, 0);
    std::swap(m_data. other.m_data);

    return *this;
}

template<typename Elem>
typename GenericMatrix<Elem>::iterator GenericMatrix<Elem>::begin() noexcept
{
    return m_data;
}

template<typename Elem>
typename GenericMatrix<Elem>::const_iterator GenericMatrix<Elem>::begin() const noexcept
{
    return m_data;
}

template<typename Elem>
typename GenericMatrix<Elem>::const_iterator GenericMatrix<Elem>::cbegin() const noexcept
{
    return m_data;
}

template<typename Elem>
typename GenericMatrix<Elem>::iterator GenericMatrix<Elem>::end() noexcept
{
    return m_data + size();
}

template<typename Elem>
typename GenericMatrix<Elem>::const_iterator GenericMatrix<Elem>::end() const noexcept
{
    return m_data + size();
}

template<typename Elem>
typename GenericMatrix<Elem>::const_iterator GenericMatrix<Elem>::cend() const noexcept
{
    return m_data + size();
}

/*!
    Returns the number of matrix rows.

    \sa columns(), size()
*/
template<typename Elem>
typename GenericMatrix<Elem>::size_type GenericMatrix<Elem>::rows() const
{
    return m_rows;
}

/*!
    Returns the number of matrix columns.

    \sa rows(), size()
*/
template<typename Elem>
typename GenericMatrix<Elem>::size_type GenericMatrix<Elem>::columns() const
{
    return m_cols;
}

/*!
    Returns the number of matrix elements.

    \sa rows(), columns()
*/
template<typename Elem>
typename GenericMatrix<Elem>::size_type GenericMatrix<Elem>::size() const
{
    return m_rows * m_cols;
}

/*!
    Returns a pointer to the raw data of this matrix.

    \sa constData()
*/
template<typename Elem>
Elem *GenericMatrix<Elem>::data() noexcept
{
    return m_data;
}

/*!
    Returns a constant pointer to the raw data of this matrix.

    \sa constData()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::data() const noexcept
{
    return m_data;
}

/*!
    Returns a constant pointer to the raw data of this matrix.

    \sa data()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::constData() const noexcept
{
    return m_data;
}

/*!
    Returns a offset \a row and \a col pointer to the raw data of this matrix.

    \note No bounds checking is performed.

    \sa dataAt()
*/
template<typename Elem>
Elem *GenericMatrix<Elem>::data(size_type row, size_type col) noexcept
{
    return m_data + computeOffset(row, col);
}

/*!
    Returns a constant offset \a row and \a col pointer to the raw data of this matrix.

    \note No bounds checking is performed.

    \sa dataAt()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::data(size_type row, size_type col) const noexcept
{
    return m_data + computeOffset(row, col);
}

/*!
    Returns a constant offset \a row and \a col pointer to the raw data of this matrix.

    \note No bounds checking is performed.

    \sa constDataAt()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::constData(size_type row, size_type col) const noexcept
{
    return m_data + computeOffset(row, col);
}

/*!
    Returns a offset \a row and \a col pointer to the raw data of this matrix.

    If \a row or \a col out of range, returns end();

    \sa constDataAt()
*/
template<typename Elem>
Elem *GenericMatrix<Elem>::dataAt(size_type row, size_type col) noexcept
{
    if (row > m_rows || col > m_cols)
        return end();
    return m_data + computeOffset(row, col);
}

/*!
    Returns a constant offset \a row and \a col pointer to the raw data of this matrix.

    If \a row or \a col out of range, returns end();

    \sa constDataAt()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::dataAt(size_type row, size_type col) const noexcept
{
    if (row > m_rows || col > m_cols)
        return end();
    return m_data + computeOffset(row, col);
}

/*!
    Returns a constant offset \a row and \a col pointer to the raw data of this matrix.
    If \a row or \a col out of range, returns end();

    \sa dataAt()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::constDataAt(size_type row, size_type col) const noexcept
{
    if (row > m_rows || col > m_cols)
        return end();
    return m_data + computeOffset(row, col);
}

/*!
    Returns a constant reference to the element at position (\a row, \a column) in this matrix.

    \note No bounds checking is performed.

    \sa at()
*/
template<typename Elem>
const Elem &GenericMatrix<Elem>::operator()(size_type row, size_type column) const
{
    return m_data[computeOffset(row, column)];
}

/*!
    Returns a reference to the element at position (\a row, \a column)
    in this matrix so that the element can be assigned to.

    \note No bounds checking is performed.

    \sa at()
*/
template<typename Elem>
Elem &GenericMatrix<Elem>::operator()(size_type row, size_type column)
{
    return m_data[computeOffset(row, column)];
}

/*!
    Returns a constant reference to the element at position (\a row, \a column)
    in this matrix, with bounds checking.

    \exception std::out_of_range if position (\a row, \a col) is not within the range of the matrix.

    \sa operator()()
*/
template<typename Elem>
const Elem &GenericMatrix<Elem>::at(size_type row, size_type col) const
{
    if (row >= 0 && row < m_rows && col >= 0 && col < m_cols)
        return m_data[computeOffset(row, col)];
    throw std::out_of_range("out of range.");
}

/*!
    Returns a reference to the element at position (\a row, \a column)
    in this matrix, with bounds checking.

    \exception std::out_of_range if position (\a row, \a col) is not within the range of the matrix.

    \sa operator()()
*/
template<typename Elem>
Elem &GenericMatrix<Elem>::at(size_type row, size_type col)
{
    if (row >= 0 && row < m_rows && col >= 0 && col < m_cols)
        return m_data[computeOffset(row, col)];
    throw std::out_of_range("out of range.");
}

/*!
    Returns \c true if this matrix element's size equal to 0,
    otherwise returns \c false.
*/
template<typename Elem>
bool GenericMatrix<Elem>::empty() const
{
    return 0 == size();
}

/*!
    Returns \c true if this matrix internal data is not null pointer 
    and matrix rows/cols greater than 0, otherwise returns \c false.
*/
template<typename Elem>
bool GenericMatrix<Elem>::isValid() const
{
    return m_data && size() > 0;
}

/*!
    Returns \c true if this matrix is the identity; false otherwise.

    \sa setToIdentity()
*/
template<typename Elem>
bool GenericMatrix<Elem>::isIdentity() const
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            if (i == j) {
                if (m_data[computeOffset(i, j)] != 1)
                    return false;
            } else {
                if (m_data[computeOffset(i, j)] != 0)
                    return false;
            }
        }
    }
    return true;
}

/*!
    Reconstructs a \a row x \a col identity matrix without initializing the contents.
*/
template<typename Elem>
void GenericMatrix<Elem>::resize(size_type row, size_type col)
{
    if (this->size() != row * col) {
        free();
        alloc(row, col);
    } 
    setSize(row, col);
}

/*!
    Reconstructs a \a row x \a col identity matrix and initialize all values with \a initialValue.
*/
template<typename Elem>
void GenericMatrix<Elem>::resize(size_type row, size_type col, const Elem &initialValue)
{
    resize(row, col);
    fill(initialValue);
}

/*!
    Reset this matrix to \a row x \a col matrix and set the elments pointer to \a data.

    \note \a data pointer will be free when destroy this matrix.
*/
template<typename Elem>
void GenericMatrix<Elem>::reset(size_type row, size_type col, Elem *data)
{
    free();
    m_rows = row;
    m_cols = col;
    m_data = data;
}

/*!
    Sets this matrix to the identity.

    \sa isIdentity()
*/
template<typename Elem>
void GenericMatrix<Elem>::setToIdentity()
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            if (j == i) {
                m_data[computeOffset(i, j)] = 1;
            } else {
                m_data[computeOffset(i, j)] = 0;
            }
        }
    }
}

/*!
    Returns a new matrix reversed by this matrix.
*/
template<typename Elem>
GenericMatrix<Elem> GenericMatrix<Elem>::reverse() const
{
    GenericMatrix<Elem> result(m_rows, m_cols);
    auto reverse_it = m_data + size() - 1;
    auto reverse_it_end = m_data - 1;
    auto its = result.begin();
    while (reverse_it != reverse_it_end) {
        *its++ = *reverse_it--;
    }
    return result;
}

/*!
    Returns this matrix, transposed about its diagonal.
*/
template<typename Elem>
GenericMatrix<Elem> GenericMatrix<Elem>::transposed() const
{
    GenericMatrix<Elem> result(m_cols, m_rows);
    for (size_type i = 0; i < m_cols; ++i) {
        for (size_type j = 0; j < m_rows; ++j) {
            result(i, j) = m_data[computeOffset(j, i)];
        }
    }  
    return result;
}

/*!
    Fills all elements of this matrix with \a value.
*/
template<typename Elem>
void GenericMatrix<Elem>::fill(const Elem &value)
{
    std::fill_n(m_data, m_rows * m_cols, value);   
}

/*!
    Do the hadamard product of this matrix and the matrix \a m.

    \exception std::invalid_argument if \a m1 and \a m2 are not homomorphic.
    \sa hadamardProduct()
*/
template<typename Elem>
void GenericMatrix<Elem>::doHadamardProduct(const GenericMatrix<Elem> &m)
{
    if (!GenericMatrix<Elem>::isHomomorphic(*this, m))
        throw std::invalid_argument("invalid argument.");

    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; i < m_cols; ++j) {
            m_data[computeOffset(i, j)] *= m(i, j);
        }
    }
}

/*!
    Returns the hadamard product of the matrix \a m1 and the matrix \a m2.

    \exception std::invalid_argument if \a m1 and \a m2 are not homomorphic.

    \sa doHadamardProduct()
*/
template<typename Elem>
GenericMatrix<Elem> GenericMatrix<Elem>::hadamardProduct(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
{
    if (!GenericMatrix<Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("invalid argument.");

    GenericMatrix<Elem> result(m1.rows(), m1.columns());
    for (size_type i = 0; i < m1.rows(); ++i) {
        for (size_type j = 0; i < m2.columns(); ++j) {
            result(i, j) = m1(i, j) * m2(i, j);
        }
    }
    return result;
}

/*!
    Extracts region of interest of this matrix.

    \note No bounds checking is performed.
*/
template<typename Elem>
GenericMatrix<Elem> GenericMatrix<Elem>::roi(size_type row1, size_type col1, size_type row2, size_type col2)
{
    if (!containsIndex(row1, col1) || !containsIndex(row2, col2)) {
        return GenericMatrix<Elem>(); 
    }
    std::pair<size_type, size_type> rowMinMax = std::minmax(row1, row2);
    std::pair<size_type, size_type> colMinMax = std::minmax(col1, col2);
    GenericMatrix<Elem> result(rowMinMax.second - rowMinMax.first + 1, colMinMax.second - colMinMax.first + 1);
    for (size_type i = rowMinMax.first; i <= rowMinMax.second; ++i) {
        std::copy_n(constData(i, colMinMax.first), result.columns(), result.data((i - rowMinMax.first), 0));
    }
    return result;
}

/*!
    Returns \c true if this matrix and matrix \a m are homomorphic; false otherwise.

    \sa isHomomorphic()
*/
template<typename Elem>
inline bool GenericMatrix<Elem>::isHomomorphicTo(const GenericMatrix<Elem> &m)
{
    return((m_rows == m.m_rows) && (m_cols == m.m_cols));
}

/*!
    Returns \c true if matrix \a m1 and matrix \a m2 are homomorphic; false otherwise.

    \sa isHomomorphicTo()
*/
template<typename Elem>
inline bool GenericMatrix<Elem>::isHomomorphic(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
{
    return ((m1.m_rows == m2.m_rows) && (m1.m_cols == m2.m_cols));
}

/*!
    Swaps matrix \a other with this matrix.
*/
template<typename Elem>
void GenericMatrix<Elem>::swap(GenericMatrix<Elem> &other)
{
    std::swap(*this, other);
}

/*!
    Adds the contents of \a m to this matrix.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.

    \sa operator-=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator+=(const GenericMatrix<Elem> &m)
{
    if (!isHomomorphicTo(m))
        throw std::invalid_argument("invalid argument.");

    for (size_type i = 0; i < m.size(); ++i) {
        m_data[i] += m.m_data[i];
    }
    return *this;
}

/*!
    Subtracts the contents of \a m from this matrix.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.

    \sa operator+=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator-=(const GenericMatrix<Elem> &m)
{
    if (!isHomomorphicTo(m))
        throw std::invalid_argument("invalid argument.");

    for (size_type i = 0; i < m.size(); ++i) {
        m_data[i] -= m.m_data[i];
    }
    return *this;
}

/*!
    Multiplies this matrix and matrix \a m, in format: this = this * m.

    \exception std::invalid_argument if this matrix's columns is not equal to matrix \a m's rows.
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator*=(const GenericMatrix<Elem> &m)
{
    if (m_cols != m.m_rows)
        throw std::invalid_argument("invalid argument.");

    GenericMatrix<Elem> result(m_rows, m.m_cols, 0);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            for (size_type k = 0; k < m_cols; ++k) {
                result(i, j) += m_data[computeOffset(i, k)] * m(k, j);
            }
        }
    }
    return (*this = std::move(result));
}

/*!
    Multiplies all elements of this matrix by \a factor.

    \sa operator/=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator*=(const Elem &factor)
{
    for (size_type i = 0; i < size(); ++i) {
        m_data[i] *= factor;
    }
    return *this;
}

/*!
    Divides all elements of this matrix by \a divisor.

    \sa operator*=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator/=(const Elem &divisor)
{
    for (size_type i = 0; i < size(); ++i) {
        m_data[i] /= divisor;
    }
    return *this;
}

/*!
    Returns \c true if this matrix is identical to \a m; false otherwise.

    \sa operator!=()
*/
template<typename Elem>
bool GenericMatrix<Elem>::operator==(const GenericMatrix<Elem> &m) const
{
    for (size_type i = 0; i < size(); ++i) {
        if (m_data[i] != m.m_data[i])
            return false;
    }
    return true;
}

/*!
    Returns \c true if this matrix is not identical to \a m; false otherwise.

    \sa operator==()
*/
template<typename Elem>
bool GenericMatrix<Elem>::operator!=(const GenericMatrix<Elem> &m) const
{
    for (size_type i = 0; i < size(); ++i) {
        if (m_data[i] != m.m_data[i])
            return true;
    }
    return false;
}

/*****************************************************************************
  friend functions
 *****************************************************************************/

/*!
    \relates GenericMatrix

    Converts a matrix to a different type matrix.

    \note The \b __ElemDst must can do \c static_cast<__ElemDst>(__Elem).
*/
template<typename __ElemDTo, typename __Elem>
GenericMatrix<__ElemDTo> matrix_cast(const GenericMatrix<__Elem> &m)
{
    GenericMatrix<__ElemDTo> result(m.rows(), m.columns());
    auto itm = m.begin();
    auto its = result.begin();
    while (itm != m.end()) {
        *its++ = static_cast<__ElemDTo>(*itm++);
    }
    return result;
}

/*!
    \relates GenericMatrix

    Returns the sum of \a m1 and \a m2.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.
*/

template<typename __Elem>
GenericMatrix<__Elem> operator+(const GenericMatrix<__Elem> &m1, const GenericMatrix<__Elem> &m2)
{
    if (!GenericMatrix<__Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("invalid argument.");

    using size_type = typename GenericMatrix<__Elem>::size_type;
    GenericMatrix<__Elem> result(m1.m_rows, m1.m_cols);
    for (size_type i = 0; i < result.size(); ++i) {
        result.m_data[i] = m1.m_data[i] + m2.m_data[i];
    }
    return result;
}

/*!
    \relates GenericMatrix

    Returns the difference of \a m1 and \a m2.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.
*/
template<typename __Elem>
GenericMatrix<__Elem> operator-(const GenericMatrix<__Elem> &m1, const GenericMatrix<__Elem> &m2)
{
    if (!GenericMatrix<__Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("invalid argument.");

    using size_type = typename GenericMatrix<__Elem>::size_type;
    GenericMatrix<__Elem> result(m1.m_rows, m1.m_cols);
    for (size_type i = 0; i < result.size(); ++i) {
        result.m_data[i] = m1.m_data[i] - m2.m_data[i];
    }
    return result;
}

/*!
    \relates GenericMatrix

    Returns the product of the \c M1xNN matrix \a m1 and the \c NNxM2 matrix \a m2
    to produce a \c M1xM2 matrix result.

    \exception std::invalid_argument if \a m1's columns is not equal to matrix \a m2's rows.
*/
template<typename __Elem>
GenericMatrix<__Elem> operator*(const GenericMatrix<__Elem> &m1, const GenericMatrix<__Elem> &m2)
{
    if (m1.m_cols != m2.m_rows)
        throw std::invalid_argument("invalid argument.");

    using size_type = typename GenericMatrix<__Elem>::size_type;
    GenericMatrix<__Elem> result(m1.m_rows, m2.m_cols, 0);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            for (size_type k = 0; k < m1.m_cols; ++k) {
                result(i, j) += m1(i, k) * m2(k, j);
            }
        }
    }
    return result;
}

/*!
    \relates GenericMatrix

    Returns the result of multiplying all elements of \a matrix by \a factor.
*/
template<typename __Elem>
GenericMatrix<__Elem> operator*(const GenericMatrix<__Elem> &matrix, const __Elem &factor)
{
    using size_type = typename GenericMatrix<__Elem>::size_type;
    GenericMatrix<__Elem> result(matrix.m_rows, matrix.m_cols);
    for (size_type i = 0; i < matrix.size(); ++i) {
        result.m_data[i] = matrix.m_data[i] * factor;
    }
    return result;
}

/*!
    \relates GenericMatrix

    Returns the result of multiplying all elements of \a matrix by \a factor.
*/
template<typename __Elem>
GenericMatrix<__Elem> operator*(const __Elem &factor, const GenericMatrix<__Elem> &matrix)
{
    using size_type = typename GenericMatrix<__Elem>::size_type;
    GenericMatrix<__Elem> result(matrix.m_rows, matrix.m_cols);
    for (size_type i = 0; i < matrix.size(); ++i) {
        result.m_data[i] = matrix.m_data[i] * factor;
    }
    return result;
}

/*!
    \relates GenericMatrix

    Returns the result of dividing all elements of \a matrix by \a divisor.
*/
template<typename __Elem>
GenericMatrix<__Elem> operator/(const GenericMatrix<__Elem> &matrix, const __Elem &divisor)
{
    using size_type = typename GenericMatrix<__Elem>::size_type;
    GenericMatrix<__Elem> result(matrix.m_rows, matrix.m_cols);
    for (size_type i = 0; i < matrix.size(); ++i) {
        result.m_data[i] = matrix.m_data[i] / divisor;
    }
    return result;
}

/*!
    \relates GenericMatrix

    Returns the negation of \a matrix.
*/
template<typename __Elem>
GenericMatrix<__Elem> operator-(const GenericMatrix<__Elem> &matrix)
{
    using size_type = typename GenericMatrix<__Elem>::size_type;
    GenericMatrix<__Elem> result(matrix.m_rows, matrix.m_cols);
    for (size_type i = 0; i < matrix.size(); ++i) {
        result.m_data[i] = -matrix.m_data[i];
    }
    return result;
}

/*!
    \relates GenericMatrix

    Writes the given \a matrix to the given \a stream and returns a
    reference to the stream.
*/
template<typename __Elem>
std::ostream &operator<<(std::ostream &stream, const GenericMatrix<__Elem> &matrix)
{
    using size_type = typename GenericMatrix<__Elem>::size_type;
    stream << "GenericMatrix<" << matrix.m_rows << ", " << matrix.m_cols << ", " << typeid(__Elem).name() << ">(" << std::endl;
    for (size_type i = 0; i < matrix.m_rows; ++i) {
        for (size_type j = 0; j < matrix.m_cols; ++j) {
            stream << stream.width(10) << matrix(i, j);
        }
        stream << std::endl;
    }
    stream << ')';
    return stream;
}

#endif // __GERNERICMATRIX_H__
