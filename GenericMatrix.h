/*!
 * \file	GenericMatrix.h
 * \brief	The header file of the GenericMatrix class.
 * \date	2019-09-19
 * \author	shaoguang
 */

#ifndef __GENERICMATRIX_H__
#define __GENERICMATRIX_H__

#include <ostream>
#include <utility>

/*!
	\enum Initialization
	\brief Indicates uninitialized.
*/
enum class Initialization {
	Uninitialized				///< Indicates uninitialized.
};
static constexpr Initialization Uninitialized = Initialization::Uninitialized;

/*!
	\class GenericMatrix
	\brief The GenericMatrix class defines a generic template matrix class.

	The GenericMatrix template has one parameter:

	\li \b Elem Element type that is visible to users of the class.
*/
template<typename Elem = float>
class GenericMatrix
{
public:
    using size_type = std::size_t;
    using value_type = typename Elem;

    GenericMatrix();
    GenericMatrix(size_type row, size_type col);
    GenericMatrix(size_type col, size_type row, Initialization);
    GenericMatrix(const GenericMatrix &other);
    GenericMatrix(GenericMatrix &&other);
    ~GenericMatrix();
    GenericMatrix &operator=(const GenericMatrix &other);
    GenericMatrix &operator=(GenericMatrix &&other);

    size_type rows() const noexcept;
    size_type columns() const noexcept;

    Elem **data();
    const Elem **data() const;
    const Elem **constData() const;

    Elem *rowData(size_type row);
    const Elem *rowData(size_type row) const;
    const Elem *constRowData(size_type row) const;
    Elem *rowDataAt(size_type row);
    const Elem *rowDataAt(size_type row) const;
    const Elem *constRowDataAt(size_type row) const;

    bool isValid() const noexcept;
    bool isIdentity() const;
    inline bool isHomomorphicTo(const GenericMatrix<Elem> &m);
    static inline bool isHomomorphic(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);

    void resize(size_type row, size_type col);
    void resize(size_type row, size_type col, Initialization);
    void setToIdentity();
    void fill(const Elem &value);
    void swap(GenericMatrix<Elem> &other);
    GenericMatrix<Elem> transposed() const;
    void doHadamardProduct(const GenericMatrix<Elem> &m);
    static GenericMatrix<Elem> hadamardProduct(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);

    template<typename ElemDst>
    GenericMatrix<ElemDst> cast() const;

    Elem &operator()(size_type row, size_type col);
    const Elem &operator()(size_type row, size_type col) const;

    Elem &at(size_type row, size_type col);
    const Elem &at(size_type row, size_type col) const;

    GenericMatrix<Elem> &operator+=(const GenericMatrix<Elem> &m);
    GenericMatrix<Elem> &operator-=(const GenericMatrix<Elem> &m);
    GenericMatrix<Elem> &operator*=(const GenericMatrix<Elem> &m);
    GenericMatrix<Elem> &operator*=(const Elem &factor);
    GenericMatrix<Elem> &operator/=(const Elem &divisor);
    bool operator==(const GenericMatrix<Elem> &m) const;
    bool operator!=(const GenericMatrix<Elem> &m) const;

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
    void free();

private:
    size_type m_rows;
    size_type m_cols;
    Elem **m_p;
};

/*****************************************************************************
    definitions
 *****************************************************************************/

/*!
    \typedef GenericMatrix::size_type

    Type for declaring the matrix's rows and columns.
*/

/*!
    \typedef GenericMatrix::value_type

    The value type of the template parameter \a Elem.
*/

/*!
    \fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix()

    Constructs a invalid matrix.

    \sa isValid()
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix()
    : m_rows(0), m_cols(0), m_p(nullptr)
{
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col)

    Constructs a \a row x \a col identity matrix.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col)
    : m_rows(row), m_cols(col), m_p(nullptr)
{
    alloc();
    setToIdentity();
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col, Initialization)

    Constructs a \a row x \a col matrix without initializing the contents.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col, Initialization)
    : m_rows(row), m_cols(col), m_p(nullptr)
{
    alloc();
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem>::~GenericMatrix()

    Destroys the matrix.
*/
template<typename Elem>
GenericMatrix<Elem>::~GenericMatrix()
{
    free();
}
    
/*!
    \fn template<typename Elem> void GenericMatrix<Elem>::alloc()
    \internal

    Alloc the matrix memory from \a _rows and \a _cols.
*/
template<typename Elem>
void GenericMatrix<Elem>::alloc()
{
    if (m_rows > 0)  {
        m_p = new Elem*[m_rows];
        if (m_cols > 0)  {
            for (size_type i = 0; i < m_rows; ++i)
                m_p[i] = new Elem[m_cols];
        } else {
            for (size_type i = 0; i < m_rows; ++i)
                m_p[i] = nullptr;
        }
    }
}

/*!
    \fn template<typename Elem> void GenericMatrix<Elem>::free()
    \internal

    Frees the matrix memory.
*/
template<typename Elem>
void GenericMatrix<Elem>::free()
{
    if (m_p)  {
        for (size_type i = 0; i < m_rows; ++i)  {
            if (m_p[i]) {
				delete[] m_p[i];
				m_p[i] = nullptr;
			}
        }
        delete[] m_p; 
        m_p = nullptr;
    }
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix(const GenericMatrix &other)

    Constructs a copy of \a other.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(const GenericMatrix &other)
    : m_rows(other._rows), m_cols(other._cols), m_p(nullptr)
{
    alloc();
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_p[i][j] = other._p[i][j];
        }
    }
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix(GenericMatrix &&other)

    Move-constructs a GenericMatrix instance, making it point at the same object that \a other was pointing to.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(GenericMatrix &&other)
    : m_rows(other._rows), m_cols(other._cols), m_p(other._p)
{
    other._p = nullptr;
    other._rows = 0;
    other._cols = 0;
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator=(const GenericMatrix &other)

    Assigns \a other to this matrix and returns a reference to this matrix.
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator=(const GenericMatrix &other)
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
            m_p[i][j] = other._p[i][j];
        }
    }
    return *this;
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator=(GenericMatrix &&other)

    Move-assigns \a other to this GenericMatrix instance.
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator=(GenericMatrix &&other)
{
    if (this == &other) {
        return *this;
    }
    free();

    m_rows = other._rows;  
    m_cols = other._cols;
    other._rows = 0; 
    other._cols = 0;
    m_p = other._p; 
    other._p = nullptr;

    return *this;
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem>::size_type GenericMatrix<Elem>::rows() const

    Returns the number of matrix rows.

    \sa columns()
*/
template<typename Elem>
typename GenericMatrix<Elem>::size_type GenericMatrix<Elem>::rows() const noexcept
{
    return m_rows;
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem>::size_type GenericMatrix<Elem>::columns() const

    Returns the number of matrix columns.

    \sa rows()
*/
template<typename Elem>
typename GenericMatrix<Elem>::size_type GenericMatrix<Elem>::columns() const noexcept
{
    return m_cols;
}

/*!
    \fn template<typename Elem> Elem **GenericMatrix<Elem>::data()

    Returns a pointer to the raw data of this matrix.

    \sa constData()
*/
template<typename Elem>
Elem **GenericMatrix<Elem>::data()
{
    return m_p;
}

/*!
    \fn template<typename Elem> const Elem **GenericMatrix<Elem>::data() const

    Returns a constant pointer to the raw data of this matrix.

    \sa constData()
*/
template<typename Elem>
const Elem **GenericMatrix<Elem>::data() const
{
    return m_p;
}

/*!
    \fn template<typename Elem> const Elem **GenericMatrix<Elem>::constData() const

    Returns a constant pointer to the raw data of this matrix.

    \sa data()
*/
template<typename Elem>
const Elem **GenericMatrix<Elem>::constData() const
{
    return m_p;
}

/*!
    \fn template<typename Elem> Elem *GenericMatrix<Elem>::rowData(size_type row)

    Returns a pointer to the \a row data of this matrix.

    \note No bounds checking is performed.

    \sa constRowData()
*/
template<typename Elem>
Elem *GenericMatrix<Elem>::rowData(size_type row)
{
    return m_p[row];
}

/*!
    \fn template<typename Elem> Elem *GenericMatrix<Elem>::rowData(size_type row)

    Returns a constant pointer to the \a row data of this matrix.

    \note No bounds checking is performed.

    \sa constRowData()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::rowData(size_type row) const
{
    return m_p[row];
}

/*!
    \fn template<typename Elem> const Elem *GenericMatrix<Elem>::constRowData(size_type row) const

    Returns a constant pointer to the \a row data of this matrix.

    \note No bounds checking is performed.

    \sa rowData()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::constRowData(size_type row) const
{
    return m_p[row];
}

/*!
    \fn template<typename Elem> Elem *GenericMatrix<Elem>::rowDataAt(size_type row)

    Returns a pointer to the \a row data of this matrix, with bounds checking.

    \exception std::out_of_range if position \a row is not within the range of the matrix.

    \sa rowData(), constRowDataAt()
*/
template<typename Elem>
Elem *GenericMatrix<Elem>::rowDataAt(size_type row)
{
    if (row < m_rows)
        return m_p[row];
    throw std::out_of_range("\fn[GenericMatrix::rowDataAt()] out of range.");
}

/*!
    \fn template<typename Elem> const Elem *GenericMatrix<Elem>::rowDataAt(size_type row) const

    Returns a constant pointer to the \a row data of this matrix, with bounds checking.

    \exception std::out_of_range if position \a row is not within the range of the matrix.

    \sa rowData(), constRowDataAt()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::rowDataAt(size_type row) const
{
    if (row < m_rows)
        return m_p[row];
    throw std::out_of_range("\fn[GenericMatrix::rowDataAt() const] out of range.");
}

/*!
    \fn const Elem *GenericMatrix<Elem>::constRowDataAt(size_type row) const

    Returns a constant pointer to the \a row data of this matrix, with bounds checking.

    \exception std::out_of_range if position \a row is not within the range of the matrix.

    \sa constRowData(), rowDataAt()
*/
template<typename Elem>
const Elem *GenericMatrix<Elem>::constRowDataAt(size_type row) const
{
    if (row < m_rows)
        return m_p[row];
    throw std::out_of_range("\fn[GenericMatrix::constRowDataAt()] out of range.");
}

/*!
    \fn template<typename Elem> bool GenericMatrix<Elem>::isValid() const

    Returns \c true if this matrix internal data is null pointer, otherwise returns \c false.
*/
template<typename Elem>
bool GenericMatrix<Elem>::isValid() const noexcept
{
    return m_p;
}

/*!
    \fn template<typename Elem> bool GenericMatrix<Elem>::isIdentity() const

    Returns \c true if this matrix is the identity; false otherwise.

    \sa setToIdentity()
*/
template<typename Elem>
bool GenericMatrix<Elem>::isIdentity() const
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            if (i == j) {
                if (m_p[i][j] != 1)
                    return false;
            } else {
                if (m_p[i][j] != 0)
                    return false;
            }
        }
    }
    return true;
}

/*!
    \fn template<typename Elem> void GenericMatrix<Elem>::resize(size_type row, size_type col)

    Reconstructs a \a row x \a col identity matrix.

    \sa isIdentity()
*/
template<typename Elem>
void GenericMatrix<Elem>::resize(size_type row, size_type col)
{
    resize(row, col, Uninitialized);
    setToIdentity();
}

/*!
    \fn template<typename Elem> void GenericMatrix<Elem>::resize(size_type row, size_type col, Initialization)

    Reconstructs a \a row x \a col identity matrix  without initializing the contents.
*/
template<typename Elem>
void GenericMatrix<Elem>::resize(size_type row, size_type col, Initialization)
{
    free();
    m_rows = row;
    m_cols = col;
    alloc();
}

/*!
    \fn template<typename Elem> void GenericMatrix<Elem>::setToIdentity()

    Sets this matrix to the identity.

    \sa isIdentity()
*/
template<typename Elem>
void GenericMatrix<Elem>::setToIdentity()
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            if (j == i)
                m_p[i][j] = 1;
            else
                m_p[i][j] = 0;
        }
    }
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> GenericMatrix<Elem>::transposed() const

    Returns this matrix, transposed about its diagonal.
*/
template<typename Elem>
GenericMatrix<Elem> GenericMatrix<Elem>::transposed() const
{
    GenericMatrix<Elem> result(m_cols, m_rows, Uninitialized);
    for (int i = 0; i < m_cols; ++i)
        for (int j = 0; j < m_rows; ++j)
            result._p[i][j] = m_p[j][i];
    return result;
}

/*!
    \fn template<typename Elem> void GenericMatrix<Elem>::fill(const Elem &value)

    Fills all elements of this matrix with \a value.
*/
template<typename Elem>
void GenericMatrix<Elem>::fill(const Elem &value)
{
    for (size_type i = 0; i < m_rows; ++i)
        for (size_type j = 0; j < m_cols; ++j)
            m_p[i][j] = value;
}

/*!
    \fn template<typename Elem> void GenericMatrix<Elem>::doHadamardProduct(const GenericMatrix<Elem> &m)

    Do the hadamard product of this matrix and the matrix \a m.

    \exception std::invalid_argument if \a m1 and \a m2 are not homomorphic.
    \sa hadamardProduct()
*/
template<typename Elem>
void GenericMatrix<Elem>::doHadamardProduct(const GenericMatrix<Elem> &m)
{
    if (!GenericMatrix<Elem>::isHomomorphic(*this, m))
        throw std::invalid_argument("\fn[GenericMatrix::hadamardProduct()] only homomorphic matrix can do hadamard product operation.");

    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; i < m_cols; ++j) {
            m_p[i][j] *= m._m[i][j];
        }
    }
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> GenericMatrix<Elem>::hadamardProduct(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)

    Returns the hadamard product of the matrix \a m1 and the matrix \a m2.

    \exception std::invalid_argument if \a m1 and \a m2 are not homomorphic.
    \sa doHadamardProduct()
*/
template<typename Elem>
GenericMatrix<Elem> GenericMatrix<Elem>::hadamardProduct(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
{
    if (!GenericMatrix<Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("\fn[GenericMatrix::hadamardProduct()] only homomorphic matrix can do hadamard product operation.");

    GenericMatrix<Elem> result(m1._rows, m1._cols, Uninitialized);
    for (size_type i = 0; i < m1._rows; ++i) {
        for (size_type j = 0; i < m2._cols; ++j) {
            result._p[i][j] = m1._p[i][j] * m2._p[i][j];
        }
    }
    return result;
}

/*!
    \fn template<typename Elem> template<typename ElemDst> GenericMatrix<ElemDst> GenericMatrix<Elem>::cast() const

    Returns a new element type matrix.
*/
template<typename Elem>
template<typename ElemDst>
GenericMatrix<ElemDst> GenericMatrix<Elem>::cast() const
{
    GenericMatrix<ElemDst> result(m_rows, m_cols, Uninitialized);
    for (size_type i = 0; i < m_rows; ++i)
        for (size_type j = 0; j < m_cols; ++j)
            result(i, j) = static_cast<ElemDst>(m_p[i][j]);
    return result;
}

/*!
    \fn template<typename Elem> bool GenericMatrix<Elem>::isHomomorphicTo(const GenericMatrix<Elem> &m)

    Returns \c true if this matrix and matrix \a m are homomorphic; false otherwise.

    \sa isHomomorphic()
*/
template<typename Elem>
inline bool GenericMatrix<Elem>::isHomomorphicTo(const GenericMatrix<Elem> &m)
{
    return((m_rows == m._rows) && (m_cols == m._cols));
}

/*!
    \fn template<typename Elem> bool GenericMatrix<Elem>::isHomomorphic(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)

    Returns \c true if matrix \a m1 and matrix \a m2 are homomorphic; false otherwise.

    \sa isHomomorphicTo()
*/
template<typename Elem>
inline bool GenericMatrix<Elem>::isHomomorphic(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
{
    return ((m1._rows == m2._rows) && (m1._cols == m2._cols));
}

/*!
    \fn template<typename Elem> void GenericMatrix<Elem>::swap(GenericMatrix<Elem> &m)

    Swaps matrix \a other with this matrix.
*/
template<typename Elem>
void GenericMatrix<Elem>::swap(GenericMatrix<Elem> &other)
{
    std::swap(m_p, other._p);
    std::swap(m_rows, other._rows);
    std::swap(m_cols, other._cols);
}

/*****************************************************************************
  operators overloading
 *****************************************************************************/

/*!
    \fn template<typename Elem> const Elem &GenericMatrix<Elem>::operator()(size_type row, size_type column) const

    Returns a constant reference to the element at position (\a row, \a column) in this matrix.

    \note No bounds checking is performed.

    \sa at()
*/
template<typename Elem>
const Elem &GenericMatrix<Elem>::operator()(size_type row, size_type column) const
{
    return m_p[row][column];
}

/*!
    \fn template<typename Elem> Elem &GenericMatrix<Elem>::operator()(size_type row, size_type column)

    Returns a reference to the element at position (\a row, \a column) 
    in this matrix so that the element can be assigned to.

    \note No bounds checking is performed.

    \sa at()
*/
template<typename Elem>
Elem &GenericMatrix<Elem>::operator()(size_type row, size_type column)
{
    return m_p[row][column];
}

/*!
    \fn template<typename Elem> const Elem &GenericMatrix<Elem>::at(size_type row, size_type col) const

    Returns a constant reference to the element at position (\a row, \a column) 
    in this matrix, with bounds checking.

    \exception std::out_of_range if position (\a row, \a col) is not within the range of the matrix.

    \sa operator()()
*/
template<typename Elem>
const Elem &GenericMatrix<Elem>::at(size_type row, size_type col) const
{
    if (row >= 0 && row < m_rows && col >= 0 && col < m_cols)
        return m_p[row][col];
    throw std::out_of_range("\fn[GenericMatrix::at()] out of range.");
}

/*!
    \fn template<typename Elem> Elem &GenericMatrix<Elem>::at(size_type row, size_type col)

    Returns a reference to the element at position (\a row, \a column) 
    in this matrix, with bounds checking.

    \exception std::out_of_range if position (\a row, \a col) is not within the range of the matrix.

    \sa operator()()
*/
template<typename Elem>
Elem &GenericMatrix<Elem>::at(size_type row, size_type col)
{
    if (row >= 0 && row < m_rows && col >= 0 && col < m_cols)
        return m_p[row][col];
    throw std::out_of_range("\fn[GenericMatrix::at()] out of range.");
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator+=(const GenericMatrix<Elem> &m)

    Adds the contents of \a m to this matrix.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.

    \sa operator-=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator+=(const GenericMatrix<Elem> &m)
{
    if (!isHomomorphicTo(m))
        throw std::invalid_argument("\fn[GenericMatrix::operator+=()] only homomorphic matrix can do add operation.");

    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_p[i][j] += m._p[i][j];
        }
    }
    return *this;
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator-=(const GenericMatrix<Elem> &m)

    Subtracts the contents of \a m from this matrix.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.

    \sa operator+=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator-=(const GenericMatrix<Elem> &m)
{
    if (!isHomomorphicTo(m))
        throw std::invalid_argument("\fn[GenericMatrix::operator-=()] only homomorphic matrix can do subtract operation.");

    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_p[i][j] -= m._p[i][j];
        }
    }
    return *this;
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator*=(const GenericMatrix<Elem> &m)

    Multiplies this matrix and matrix \a m, in format: this = this * m.

    \exception std::invalid_argument if this matrix's columns is not equal to matrix \a m's rows.
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator*=(const GenericMatrix<Elem> &m)
{
    if (m_cols != m._rows)
        throw std::invalid_argument("\fn[GenericMatrix::operator*=()] can't do multiplication operation.");

    GenericMatrix<Elem> result(m_rows, m._cols, Uninitialized);
    result.fill(0);
    for (size_type i = 0; i < result._rows; ++i) {
        for (size_type j = 0; j < result._cols; ++j) {
            for (size_type k = 0; k < m_cols; ++k) {
                result._p[i][j] += (m_p[i][k] * m._p[k][j]);
            }
        }
    }
    return (*this = std::move(result));
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator*=(const Elem &factor)

    Multiplies all elements of this matrix by \a factor.

    \sa operator/=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator*=(const Elem &factor)
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_p[i][j] *= factor;
        }
    }
    return *this;
}

/*!
    \fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator/=(const Elem &divisor)

    Divides all elements of this matrix by \a divisor.

    \sa operator*=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator/=(const Elem &divisor)
{
    for (size_type i = 0; i < m_rows; ++i) {
        for (size_type j = 0; j < m_cols; ++j) {
            m_p[i][j] /= divisor;
        }
    }
    return *this;
}

/*!
    \fn template<typename Elem> bool GenericMatrix<Elem>::operator==(const GenericMatrix<Elem> &m) const

    Returns \c true if this matrix is identical to \a m; false otherwise.

    \sa operator!=()
*/
template<typename Elem>
bool GenericMatrix<Elem>::operator==(const GenericMatrix<Elem> &m) const
{
    for (size_type i = 0; i < m_rows; ++i)
        for (int j = 0; j < m_cols; ++j) {
            if (m_p[i][j] != m._p[i][j])
                return false;
        }
    return true;
}

/*!
    \fn template<typename Elem> bool GenericMatrix<Elem>::operator!=(const GenericMatrix<Elem> &m) const

    Returns \c true if this matrix is not identical to \a m; false otherwise.

    \sa operator==()
*/
template<typename Elem>
bool GenericMatrix<Elem>::operator!=(const GenericMatrix<Elem> &m) const
{
    for (size_type i = 0; i < m_rows; ++i)
        for (int j = 0; j < m_cols; ++j) {
            if (m_p[i][j] != m._p[i][j])
                return true;
        }
    return false;
}

/*****************************************************************************
  friend functions
 *****************************************************************************/

/*!
    \fn GenericMatrix<Elem> operator+(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
    \relates GenericMatrix

    Returns the sum of \a m1 and \a m2.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.
*/

template<typename Elem>
GenericMatrix<Elem> operator+(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
{
    if (!GenericMatrix<Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("\fn[GenericMatrix friend operator+()] only homomorphic matrix can do add operation.");

    using size_type = typename GenericMatrix<Elem>::size_type;
    GenericMatrix<Elem> result(m1.m_rows, m1.m_cols, Uninitialized);
    for (size_type i = 0; i < result._rows; ++i) {
        for (size_type j = 0; j < result._cols; ++j) {
            result._p[i][j] = m1.m_p[i][j] + m2.m_p[i][j];
        }
    }
    return result;
}

/*!
    \fn GenericMatrix<Elem> operator-(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
    \relates GenericMatrix

    Returns the difference of \a m1 and \a m2.

    \exception std::invalid_argument if this matrix and matrix \a m are not homomorphic.
*/
template<typename Elem>
GenericMatrix<Elem> operator-(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
{
    if (!GenericMatrix<Elem>::isHomomorphic(m1, m2))
        throw std::invalid_argument("\fn[GenericMatrix friend operator+()] only homomorphic matrix can do subtract operation.");

    using size_type = typename GenericMatrix<Elem>::size_type;
    GenericMatrix<Elem> result(m1.m_rows, m1.m_cols, Uninitialized);
    for (size_type i = 0; i < result._rows; ++i) {
        for (size_type j = 0; j < result._cols; ++j) {
            result._p[i][j] = m1.m_p[i][j] - m2.m_p[i][j];
        }
    }
    return result;
}

/*!
    \fn GenericMatrix<Elem> operator*(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
    \relates GenericMatrix

    Returns the product of the \c M1xNN matrix \a m1 and the \c NNxM2 matrix \a m2
    to produce a \c M1xM2 matrix result.

    \exception std::invalid_argument if \a m1's columns is not equal to matrix \a m2's rows.
*/
template<typename Elem>
GenericMatrix<Elem> operator*(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2)
{
    if (m1.m_cols != m2.m_rows)
        throw std::invalid_argument("input matrix can not do multiplication operation.");

    using size_type = typename GenericMatrix<Elem>::size_type;
    GenericMatrix<Elem> result(m1.m_rows, m2.m_cols, Uninitialized);
    result.fill(0);
    for (size_type i = 0; i < result._rows; ++i) {
        for (size_type j = 0; j < result._cols; ++j) {
            for (size_type k = 0; k < m1.m_cols; ++k) {
                result._p[i][j] += m1.m_p[i][k] * m2.m_p[k][j];
            }
        }
    }
    return result;
}

/*!
    \fn GenericMatrix<Elem> operator*(const GenericMatrix<Elem> &matrix, const Elem &factor)
    \relates GenericMatrix

    Returns the result of multiplying all elements of \a matrix by \a factor.
*/
template<typename Elem>
GenericMatrix<Elem> operator*(const GenericMatrix<Elem> &matrix, const Elem &factor)
{
    using size_type = typename GenericMatrix<Elem>::size_type;
    GenericMatrix<Elem> result(matrix.m_rows, matrix.m_cols, Uninitialized);
    for (size_type i = 0; i < result._rows; ++i) {
        for (size_type j = 0; j < result._cols; ++j) {
            result._p[i][j] = matrix.m_p[i][j] * factor;
        }
    }
    return result;
}

/*!
    \fn GenericMatrix<Elem> operator*(const Elem &factor, const GenericMatrix<Elem> &matrix)
    \relates GenericMatrix

    Returns the result of multiplying all elements of \a matrix by \a factor.
*/
template<typename Elem>
GenericMatrix<Elem> operator*(const Elem &factor, const GenericMatrix<Elem> &matrix)
{
    using size_type = typename GenericMatrix<Elem>::size_type;
    GenericMatrix<Elem> result(matrix.m_rows, matrix.m_cols, Uninitialized);
    for (size_type i = 0; i < result._rows; ++i) {
        for (size_type j = 0; j < result._cols; ++j) {
            result._p[i][j] = matrix.m_p[i][j] * factor;
        }
    }
    return result;
}

/*!
    \fn GenericMatrix<Elem> operator/(const GenericMatrix<Elem> &matrix, const Elem &divisor)
    \relates GenericMatrix

    Returns the result of dividing all elements of \a matrix by \a divisor.
*/
template<typename Elem>
GenericMatrix<Elem> operator/(const GenericMatrix<Elem> &matrix, const Elem &divisor)
{
    using size_type = typename GenericMatrix<Elem>::size_type;
    GenericMatrix<Elem> result(matrix.m_rows, matrix.m_cols, Uninitialized);
    for (size_type i = 0; i < result._rows; ++i) {
        for (size_type j = 0; j < result._cols; ++j) {
            result._p[i][j] = matrix.m_p[i][j] / divisor;
        }
    }
    return result;
}

/*!
    \fn GenericMatrix<Elem> operator-(const GenericMatrix<Elem> &matrix)
    \relates GenericMatrix

    Returns the negation of \a matrix.
*/
template<typename Elem>
GenericMatrix<Elem> operator-(const GenericMatrix<Elem> &matrix)
{
    using size_type = typename GenericMatrix<Elem>::size_type;
    GenericMatrix<Elem> result(matrix.m_rows, matrix.m_cols, Uninitialized);
    for (size_type i = 0; i < matrix.m_rows; ++i)
        for (size_type j = 0; j < matrix.m_cols; ++j)
            result._p[i][j] = -matrix.m_p[i][j];
    return result;
}

/*!
    \fn std::ostream &operator<<(std::ostream &os, const GenericMatrix<Elem> &matrix)
    \relates GenericMatrix

    Writes the given \a matrix to the given \a stream and returns a
    reference to the stream.
*/
template<typename Elem>
std::ostream &operator<<(std::ostream &stream, const GenericMatrix<Elem> &matrix)
{
    using size_type = typename GenericMatrix<Elem>::size_type;
    stream << "GenericMatrix<" << matrix.m_rows << ", " << matrix.m_cols << ", " << typeid(Elem).name() << ">(" << std::endl;
    for (size_type i = 0; i < matrix.m_rows; ++i) {
        for (size_type j = 0; j < matrix.m_cols; ++j)
            stream << stream.width(10) << matrix.m_p[i][j];
        stream << std::endl;
    }
    stream << ')';
    return stream;
}

#endif // !__GENERICMATRIX_H__
