/****************************************************************************
**
** Copyright (C) 2019 Shaoguang. All rights reserved.
** 
** GenericMatrix2, a generic template matrix class, 
** the matrix elements are managed by a two-dimension pointer.
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

#ifndef __GERNERICMATRIX2_H__
#define __GERNERICMATRIX2_H__

#include <ostream>
#include <utility>

/*!
    \class GenericMatrix2
    \brief The GenericMatrix2 class defines a generic template matrix class,
    the matrix elements are managed by a two-dimension pointer.

    The GenericMatrix2 template has one parameter:
    \li \b Elem Element type that is visible to users of the class.
*/
template<typename Elem = float>
class GenericMatrix2
{
public:
    using size_type = std::size_t;
    using value_type = Elem;

    GenericMatrix2();
    GenericMatrix2(size_type row, size_type col);
    GenericMatrix2(size_type row, size_type col, const Elem &initialValue);
    GenericMatrix2(const GenericMatrix2 &other);
    GenericMatrix2(GenericMatrix2 &&other);
    ~GenericMatrix2();
    GenericMatrix2 &operator=(const GenericMatrix2 &other);
    GenericMatrix2 &operator=(GenericMatrix2 &&other);

    size_type rows() const noexcept;
    size_type columns() const noexcept;

    Elem **data();
    const Elem * const *data() const;
    const Elem * const *constData() const;

    Elem *operator[](size_type row);
    const Elem *operator[](size_type row) const;
    Elem *rowData(size_type row);
    const Elem *rowData(size_type row) const;
    const Elem *constRowData(size_type row) const;
    Elem *at(size_type row);
    const Elem *at(size_type row) const;

    Elem &operator()(size_type row, size_type col);
    const Elem &operator()(size_type row, size_type col) const;
    Elem &at(size_type row, size_type col);
    const Elem &at(size_type row, size_type col) const;

    bool isValid() const noexcept;
    bool isIdentity() const;
    inline bool isHomomorphicTo(const GenericMatrix2<Elem> &m);
    static inline bool isHomomorphic(const GenericMatrix2<Elem> &m1, const GenericMatrix2<Elem> &m2);

    void resize(size_type row, size_type col);
    void resize(size_type row, size_type col, const Elem &initialValue);
    void setToIdentity();
    void fill(const Elem &value);
    void swap(GenericMatrix2<Elem> &other);
    GenericMatrix2<Elem> transposed() const;
    void doHadamardProduct(const GenericMatrix2<Elem> &m);
    static GenericMatrix2<Elem> hadamardProduct(const GenericMatrix2<Elem> &m1, const GenericMatrix2<Elem> &m2);

    template<typename ElemDst>
    GenericMatrix2<ElemDst> cast() const;

    GenericMatrix2<Elem> &operator+=(const GenericMatrix2<Elem> &m);
    GenericMatrix2<Elem> &operator-=(const GenericMatrix2<Elem> &m);
    GenericMatrix2<Elem> &operator*=(const GenericMatrix2<Elem> &m);
    GenericMatrix2<Elem> &operator*=(const Elem &factor);
    GenericMatrix2<Elem> &operator/=(const Elem &divisor);
    bool operator==(const GenericMatrix2<Elem> &m) const;
    bool operator!=(const GenericMatrix2<Elem> &m) const;

    template<typename __Elem>
    friend GenericMatrix2<__Elem> operator+(const GenericMatrix2<__Elem> &m1, const GenericMatrix2<__Elem> &m2);
    template<typename __Elem>
    friend GenericMatrix2<__Elem> operator-(const GenericMatrix2<__Elem> &m1, const GenericMatrix2<__Elem> &m2);
    template<typename __Elem>
    friend GenericMatrix2<__Elem> operator*(const GenericMatrix2<__Elem> &m1, const GenericMatrix2<__Elem> &m2);
    template<typename __Elem>
    friend GenericMatrix2<__Elem> operator*(const GenericMatrix2<__Elem> &matrix, const __Elem &factor);
    template<typename __Elem>
    friend GenericMatrix2<__Elem> operator*(const __Elem &factor, const GenericMatrix2<__Elem> &matrix);
    template<typename __Elem>
    friend GenericMatrix2<__Elem> operator/(const GenericMatrix2<__Elem> &matrix, const __Elem &divisor);
    template<typename __Elem>
    friend GenericMatrix2<__Elem> operator-(const GenericMatrix2<__Elem> &matrix);
    template<typename __Elem>
    friend std::ostream &operator<<(std::ostream &os, const GenericMatrix2<__Elem> &matrix);

private:
    void alloc();
    void free();

private:
    size_type m_rows;
    size_type m_cols;
    Elem **m_data;
};


/*!
    \typedef GenericMatrix2::size_type

    Type for declaring the matrix's rows and columns.
*/

/*!
    \typedef GenericMatrix2::value_type

    The value type of the template parameter \a Elem.
*/

/*!
    Constructs a invalid matrix.

    \sa isValid()
*/
template<typename Elem>
GenericMatrix2<Elem>::GenericMatrix2()
    : m_rows(0), m_cols(0), m_data(nullptr)
{
}

/*!
    Constructs a \a row x \a col matrix without initializing the contents.
*/
template<typename Elem>
GenericMatrix2<Elem>::GenericMatrix2(size_type row, size_type col)
    : m_rows(row), m_cols(col), m_data(nullptr)
{
    alloc();
}

/*!
    Constructs a \a row x \a col matrix and initialize all values with \a initialValue.

    \sa fill()
*/
template<typename Elem>
GenericMatrix2<Elem>::GenericMatrix2(size_type row, size_type col, const Elem &initialValue)
    : m_rows(row), m_cols(col), m_data(nullptr)
{
    alloc();
    fill(initialValue);
}

/*!
    Destroys the matrix.
*/
template<typename Elem>
GenericMatrix2<Elem>::~GenericMatrix2()
{
    free();
}
    
/*!
    \internal

    Alloc the matrix memory from \a m_rows and \a m_cols.
*/
template<typename Elem>
void GenericMatrix2<Elem>::alloc()
{
    if (m_rows > 0 && m_cols > 0) {
        m_data = new Elem*[m_rows];
        for (size_type i = 0; i < m_rows; ++i)
            m_data[i] = new Elem[m_cols];
    }
}

/*!
    \internal

    Frees the matrix memory.
*/
template<typename Elem>
void GenericMatrix2<Elem>::free()
{
    if (m_data)  {
        for (size_type i = 0; i < m_rows; ++i)  {
            if (m_data[i])  {
                delete[] m_data[i];
                m_data[i] = nullptr;
            }
        }
        delete[] m_data; 
        m_data = nullptr;
    }
}

/*!
    Constructs a copy of \a other.
*/
template<typename Elem>
GenericMatrix2<Elem>::GenericMatrix2(const GenericMatrix2 &other)
    : m_rows(other.m_rows), m_cols(other.m_cols), m_data(nullptr)
{
    alloc();
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_data[i][j] = other.m_data[i][j];
        }
    }
}

/*!
    Move-constructs a GenericMatrix2 instance, making it point at the same object that \a other was pointing to.
*/
template<typename Elem>
GenericMatrix2<Elem>::GenericMatrix2(GenericMatrix2 &&other)
    : m_rows(other.m_rows), m_cols(other.m_cols), m_data(other.m_data)
{
    other.m_data = nullptr;
    other.m_rows = 0;
    other.m_cols = 0;
}

/*!
    Assigns \a other to this matrix and returns a reference to this matrix.
*/
template<typename Elem>
GenericMatrix2<Elem> &GenericMatrix2<Elem>::operator=(const GenericMatrix2 &other)
{
    if (this == &other) {
        return *this;
    }

    if (m_rows != other.m_rows || m_cols != other.m_cols) {
        free();
        m_rows = other.m_rows;
        m_cols = other.m_cols;
        alloc();
    }

    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_data[i][j] = other.m_data[i][j];
        }
    }
    return *this;
}

/*!
    Move-assigns \a other to this GenericMatrix2 instance.
*/
template<typename Elem>
GenericMatrix2<Elem> &GenericMatrix2<Elem>::operator=(GenericMatrix2 &&other)
{
    if (this == &other) {
        return *this;
    }
    free();

    m_rows = other.m_rows;  
    m_cols = other.m_cols;
    other.m_rows = 0; 
    other.m_cols = 0;
    m_data = other.m_data; 
    other.m_data = nullptr;

    return *this;
}

/*!
    Returns the number of matrix rows.

    \sa columns()
*/
template<typename Elem>
typename GenericMatrix2<Elem>::size_type GenericMatrix2<Elem>::rows() const noexcept
{
    return m_rows;
}

/*!
    Returns the number of matrix columns.

    \sa rows()
*/
template<typename Elem>
typename GenericMatrix2<Elem>::size_type GenericMatrix2<Elem>::columns() const noexcept
{
    return m_cols;
}

/*!
    Returns a pointer to the raw data of this matrix.

    \sa constData()
*/
template<typename Elem>
Elem **GenericMatrix2<Elem>::data()
{
    return m_data;
}

/*!
    Returns a constant pointer to the raw data of this matrix.

    \sa constData()
*/
template<typename Elem>
const Elem * const *GenericMatrix2<Elem>::data() const
{
    return m_data;
}

/*!
    Returns a constant pointer to the raw data of this matrix.

    \sa data()
*/
template<typename Elem>
const Elem * const *GenericMatrix2<Elem>::constData() const
{
    return m_data;
}

/*!
    Returns a pointer to the \a row data of this matrix.

    \note No bounds checking is performed.

    \sa rowData()
*/
template<typename Elem>
Elem *GenericMatrix2<Elem>::operator[](size_type row)
{
    return m_data[row];
}

/*!
    Returns a pointer to the \a row data of this matrix.

    \note No bounds checking is performed.

    \sa constRowData()
*/
template<typename Elem>
const Elem *GenericMatrix2<Elem>::operator[](size_type row) const
{
    return m_data[row];
}

/*!
    Returns a pointer to the \a row data of this matrix.

    \note No bounds checking is performed.

    \sa constRowData()
*/
template<typename Elem>
Elem *GenericMatrix2<Elem>::rowData(size_type row)
{
    return m_data[row];
}

/*!
    Returns a constant pointer to the \a row data of this matrix.

    \note No bounds checking is performed.

    \sa constRowData()
*/
template<typename Elem>
const Elem *GenericMatrix2<Elem>::rowData(size_type row) const
{
    return m_data[row];
}

/*!
    Returns a constant pointer to the \a row data of this matrix.

    \note No bounds checking is performed.

    \sa rowData()
*/
template<typename Elem>
const Elem *GenericMatrix2<Elem>::constRowData(size_type row) const
{
    return m_data[row];
}

/*!
    Returns a pointer to the \a row data of this matrix, with bounds checking.

    \exception std::out_of_range if position \a row is not within the range of the matrix.

    \sa rowData(), constRowDataAt()
*/
template<typename Elem>
Elem *GenericMatrix2<Elem>::at(size_type row)
{
    if (row < m_rows)
        return m_data[row];
    throw std::out_of_range("out of range.");
}

/*!
    Returns a constant pointer to the \a row data of this matrix, with bounds checking.

    \exception std::out_of_range if position \a row is not within the range of the matrix.

    \sa constRowData(), rowDataAt()
*/
template<typename Elem>
const Elem *GenericMatrix2<Elem>::at(size_type row) const
{
    if (row < m_rows)
        return m_data[row];
    throw std::out_of_range("out of range.");
}

/*!
    Returns a constant reference to the element at position (\a row, \a column) in this matrix.

    \note No bounds checking is performed.

    \sa at()
*/
template<typename Elem>
const Elem &GenericMatrix2<Elem>::operator()(size_type row, size_type column) const
{
    return m_data[row][column];
}

/*!
    Returns a reference to the element at position (\a row, \a column)
    in this matrix so that the element can be assigned to.

    \note No bounds checking is performed.

    \sa at()
*/
template<typename Elem>
Elem &GenericMatrix2<Elem>::operator()(size_type row, size_type column)
{
    return m_data[row][column];
}

/*!
    Returns a constant reference to the element at position (\a row, \a column)
    in this matrix, with bounds checking.

    \exception std::out_of_range if position (\a row, \a col) is not within the range of the matrix.

    \sa operator()()
*/
template<typename Elem>
const Elem &GenericMatrix2<Elem>::at(size_type row, size_type col) const
{
    if (row >= 0 && row < m_rows && col >= 0 && col < m_cols)
        return m_data[row][col];
    throw std::out_of_range("out of range.");
}

/*!
    Returns a reference to the element at position (\a row, \a column)
    in this matrix, with bounds checking.

    \exception std::out_of_range if position (\a row, \a col) is not within the range of the matrix.

    \sa operator()()
*/
template<typename Elem>
Elem &GenericMatrix2<Elem>::at(size_type row, size_type col)
{
    if (row >= 0 && row < m_rows && col >= 0 && col < m_cols)
        return m_data[row][col];
    throw std::out_of_range("out of range.");
}

/*!
    Returns \c true if this matrix internal data is not null pointer and matrix rows/cols greater than 0, otherwise returns \c false.
*/
template<typename Elem>
bool GenericMatrix2<Elem>::isValid() const noexcept
{
    return m_data && m_rows > 0 && m_cols > 0;
}

/*!
    Returns \c true if this matrix is the identity; false otherwise.

    \sa setToIdentity()
*/
template<typename Elem>
bool GenericMatrix2<Elem>::isIdentity() const
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            if (i == j) {
                if (m_data[i][j] != 1)
                    return false;
            } else {
                if (m_data[i][j] != 0)
                    return false;
            }
        }
    }
    return true;
}

/*!
    Reconstructs a \a row x \a col identity matrix without initializing the contents.

    \sa isIdentity()
*/
template<typename Elem>
void GenericMatrix2<Elem>::resize(size_type row, size_type col)
{
    if (row != m_rows && col != m_cols) {
        m_rows = row;
        m_cols = col;
        alloc();
    }
}

/*!
    Reconstructs a \a row x \a col identity matrix and initialize all values with \a initialValue.
*/
template<typename Elem>
void GenericMatrix2<Elem>::resize(size_type row, size_type col, const Elem &initialValue)
{
    resize(row, col);
    fill(initialValue);
}

/*!
    Sets this matrix to the identity.

    \sa isIdentity()
*/
template<typename Elem>
void GenericMatrix2<Elem>::setToIdentity()
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            if (j == i) {
                m_data[i][j] = 1;
            } else {
                m_data[i][j] = 0;
            }
        }
    }
}

/*!
    Returns this matrix, transposed about its diagonal.
*/
template<typename Elem>
GenericMatrix2<Elem> GenericMatrix2<Elem>::transposed() const
{
    GenericMatrix2<Elem> result(m_cols, m_rows);
    for (size_type i = 0; i < m_cols; ++i) {
        for (size_type j = 0; j < m_rows; ++j) {
            result.m_data[i][j] = m_data[j][i];
        }
    }  
    return result;
}

/*!
    Fills all elements of this matrix with \a value.
*/
template<typename Elem>
void GenericMatrix2<Elem>::fill(const Elem &value)
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_data[i][j] = value;
        }
    }        
}

/*!
    Do the hadamard product of this matrix and the matrix \a m.

    \exception std::invalid_argument if \a m1 and \a m2 are not homomorphic.
    \sa hadamardProduct()
*/
template<typename Elem>
void GenericMatrix2<Elem>::doHadamardProduct(const GenericMatrix2<Elem> &m)
{
    if (!GenericMatrix2<Elem>::isHomomorphic(*this, m))
        throw std::invalid_argument("invalid argument.");

    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; i < m_cols; ++j) {
            m_data[i][j] *= m.m_data[i][j];
        }
    }
}

/*!
    Returns the hadamard product of the matrix \a m1 and the matrix \a m2.

    \exception std::invalid_argument if \a m1 and \a m2 are not homomorphic.

    \sa doHadamardProduct()
*/
template<typename Elem>
GenericMatrix2<Elem> GenericMatrix2<Elem>::hadamardProduct(const GenericMatrix2<Elem> &m1, const GenericMatrix2<Elem> &m2)
{
    if (!GenericMatrix2<Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("invalid argument.");

    GenericMatrix2<Elem> result(m1.m_rows, m1.m_cols);
    for (size_type i = 0; i < m1.m_rows; ++i) {
        for (size_type j = 0; i < m2.m_cols; ++j) {
            result.m_data[i][j] = m1.m_data[i][j] * m2.m_data[i][j];
        }
    }
    return result;
}

/*!
    Returns a new element type matrix.

    \note The \b ElemDst must can do \c static_cast<ElemDst>(Elem).
*/
template<typename Elem>
template<typename ElemDst>
GenericMatrix2<ElemDst> GenericMatrix2<Elem>::cast() const
{
    GenericMatrix2<ElemDst> result(m_rows, m_cols);
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            result(i, j) = static_cast<ElemDst>(m_data[i][j]);
        }
    }      
    return result;
}

/*!
    Returns \c true if this matrix and matrix \a m are homomorphic; false otherwise.

    \sa isHomomorphic()
*/
template<typename Elem>
inline bool GenericMatrix2<Elem>::isHomomorphicTo(const GenericMatrix2<Elem> &m)
{
    return((m_rows == m.m_rows) && (m_cols == m.m_cols));
}

/*!
    Returns \c true if matrix \a m1 and matrix \a m2 are homomorphic; false otherwise.

    \sa isHomomorphicTo()
*/
template<typename Elem>
inline bool GenericMatrix2<Elem>::isHomomorphic(const GenericMatrix2<Elem> &m1, const GenericMatrix2<Elem> &m2)
{
    return ((m1.m_rows == m2.m_rows) && (m1.m_cols == m2.m_cols));
}

/*!
    Swaps matrix \a other with this matrix.
*/
template<typename Elem>
void GenericMatrix2<Elem>::swap(GenericMatrix2<Elem> &other)
{
    std::swap(m_data, other.m_data);
    std::swap(m_rows, other.m_rows);
    std::swap(m_cols, other.m_cols);
}

/*!
    Adds the contents of \a m to this matrix.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.

    \sa operator-=()
*/
template<typename Elem>
GenericMatrix2<Elem> &GenericMatrix2<Elem>::operator+=(const GenericMatrix2<Elem> &m)
{
    if (!isHomomorphicTo(m))
        throw std::invalid_argument("invalid argument.");

    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_data[i][j] += m.m_data[i][j];
        }
    }
    return *this;
}

/*!
    Subtracts the contents of \a m from this matrix.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.

    \sa operator+=()
*/
template<typename Elem>
GenericMatrix2<Elem> &GenericMatrix2<Elem>::operator-=(const GenericMatrix2<Elem> &m)
{
    if (!isHomomorphicTo(m))
        throw std::invalid_argument("invalid argument.");

    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_data[i][j] -= m.m_data[i][j];
        }
    }
    return *this;
}

/*!
    Multiplies this matrix and matrix \a m, in format: this = this * m.

    \exception std::invalid_argument if this matrix's columns is not equal to matrix \a m's rows.
*/
template<typename Elem>
GenericMatrix2<Elem> &GenericMatrix2<Elem>::operator*=(const GenericMatrix2<Elem> &m)
{
    if (m_cols != m.m_rows)
        throw std::invalid_argument("invalid argument.");

    GenericMatrix2<Elem> result(m_rows, m.m_cols, 0);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            for (size_type k = 0; k < m_cols; ++k) {
                result.m_data[i][j] += (m_data[i][k] * m.m_data[k][j]);
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
GenericMatrix2<Elem> &GenericMatrix2<Elem>::operator*=(const Elem &factor)
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_data[i][j] *= factor;
        }
    }
    return *this;
}

/*!
    Divides all elements of this matrix by \a divisor.

    \sa operator*=()
*/
template<typename Elem>
GenericMatrix2<Elem> &GenericMatrix2<Elem>::operator/=(const Elem &divisor)
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_data[i][j] /= divisor;
        }
    }
    return *this;
}

/*!
    Returns \c true if this matrix is identical to \a m; false otherwise.

    \sa operator!=()
*/
template<typename Elem>
bool GenericMatrix2<Elem>::operator==(const GenericMatrix2<Elem> &m) const
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            if (m_data[i][j] != m.m_data[i][j])
                return false;
        }
    }
    return true;
}

/*!
    Returns \c true if this matrix is not identical to \a m; false otherwise.

    \sa operator==()
*/
template<typename Elem>
bool GenericMatrix2<Elem>::operator!=(const GenericMatrix2<Elem> &m) const
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            if (m_data[i][j] != m.m_data[i][j])
                return true;
        }
    }
    return false;
}

/*****************************************************************************
  friend functions
 *****************************************************************************/

/*!
    \relates GenericMatrix2

    Returns the sum of \a m1 and \a m2.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.
*/

template<typename __Elem>
GenericMatrix2<__Elem> operator+(const GenericMatrix2<__Elem> &m1, const GenericMatrix2<__Elem> &m2)
{
    if (!GenericMatrix2<__Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("invalid argument.");

    using size_type = typename GenericMatrix2<__Elem>::size_type;
    GenericMatrix2<__Elem> result(m1.m_rows, m1.m_cols);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            result.m_data[i][j] = m1.m_data[i][j] + m2.m_data[i][j];
        }
    }
    return result;
}

/*!
    \relates GenericMatrix2

    Returns the difference of \a m1 and \a m2.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.
*/
template<typename __Elem>
GenericMatrix2<__Elem> operator-(const GenericMatrix2<__Elem> &m1, const GenericMatrix2<__Elem> &m2)
{
    if (!GenericMatrix2<__Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("invalid argument.");

    using size_type = typename GenericMatrix2<__Elem>::size_type;
    GenericMatrix2<__Elem> result(m1.m_rows, m1.m_cols);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            result.m_data[i][j] = m1.m_data[i][j] - m2.m_data[i][j];
        }
    }
    return result;
}

/*!
    \relates GenericMatrix2

    Returns the product of the \c M1xNN matrix \a m1 and the \c NNxM2 matrix \a m2
    to produce a \c M1xM2 matrix result.

    \exception std::invalid_argument if \a m1's columns is not equal to matrix \a m2's rows.
*/
template<typename __Elem>
GenericMatrix2<__Elem> operator*(const GenericMatrix2<__Elem> &m1, const GenericMatrix2<__Elem> &m2)
{
    if (m1.m_cols != m2.m_rows)
        throw std::invalid_argument("invalid argument.");

    using size_type = typename GenericMatrix2<__Elem>::size_type;
    GenericMatrix2<__Elem> result(m1.m_rows, m2.m_cols, 0);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            for (size_type k = 0; k < m1.m_cols; ++k) {
                result.m_data[i][j] += m1.m_data[i][k] * m2.m_data[k][j];
            }
        }
    }
    return result;
}

/*!
    \relates GenericMatrix2

    Returns the result of multiplying all elements of \a matrix by \a factor.
*/
template<typename __Elem>
GenericMatrix2<__Elem> operator*(const GenericMatrix2<__Elem> &matrix, const __Elem &factor)
{
    using size_type = typename GenericMatrix2<__Elem>::size_type;
    GenericMatrix2<__Elem> result(matrix.m_rows, matrix.m_cols);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            result.m_data[i][j] = matrix.m_data[i][j] * factor;
        }
    }
    return result;
}

/*!
    \relates GenericMatrix2

    Returns the result of multiplying all elements of \a matrix by \a factor.
*/
template<typename __Elem>
GenericMatrix2<__Elem> operator*(const __Elem &factor, const GenericMatrix2<__Elem> &matrix)
{
    using size_type = typename GenericMatrix2<__Elem>::size_type;
    GenericMatrix2<__Elem> result(matrix.m_rows, matrix.m_cols);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            result.m_data[i][j] = matrix.m_data[i][j] * factor;
        }
    }
    return result;
}

/*!
    \relates GenericMatrix2

    Returns the result of dividing all elements of \a matrix by \a divisor.
*/
template<typename __Elem>
GenericMatrix2<__Elem> operator/(const GenericMatrix2<__Elem> &matrix, const __Elem &divisor)
{
    using size_type = typename GenericMatrix2<__Elem>::size_type;
    GenericMatrix2<__Elem> result(matrix.m_rows, matrix.m_cols);
    for (size_type i = 0; i < result.m_rows; ++i) {
        for (size_type j = 0; j < result.m_cols; ++j) {
            result.m_data[i][j] = matrix.m_data[i][j] / divisor;
        }
    }
    return result;
}

/*!
    \relates GenericMatrix2

    Returns the negation of \a matrix.
*/
template<typename __Elem>
GenericMatrix2<__Elem> operator-(const GenericMatrix2<__Elem> &matrix)
{
    using size_type = typename GenericMatrix2<__Elem>::size_type;
    GenericMatrix2<__Elem> result(matrix.m_rows, matrix.m_cols);
    for (size_type i = 0; i < matrix.m_rows; ++i) {
        for (size_type j = 0; j < matrix.m_cols; ++j) {
            result.m_data[i][j] = -matrix.m_data[i][j];
        }
    }
    return result;
}

/*!
    \relates GenericMatrix2

    Writes the given \a matrix to the given \a stream and returns a
    reference to the stream.
*/
template<typename __Elem>
std::ostream &operator<<(std::ostream &stream, const GenericMatrix2<__Elem> &matrix)
{
    using size_type = typename GenericMatrix2<__Elem>::size_type;
    stream << "GenericMatrix2<" << matrix.m_rows << ", " << matrix.m_cols << ", " << typeid(__Elem).name() << ">(" << std::endl;
    for (size_type i = 0; i < matrix.m_rows; ++i) {
        for (size_type j = 0; j < matrix.m_cols; ++j) {
            stream << stream.width(10) << matrix.m_data[i][j];
        }
        stream << std::endl;
    }
    stream << ')';
    return stream;
}

#endif // __GERNERICMATRIX2_H__
