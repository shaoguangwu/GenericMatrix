/*!
 * \file	GenericMatrix.h
 * \brief	The header file of the GenericMatrix class.
 * \date	2019-08-28
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
	\ingroup DataStructure
	\inmodule Utils

	The GenericMatrix template has one parameter:

	\li \b Elem Element type that is visible to users of the class.
*/
template<typename Elem = float>
class GenericMatrix
{
public:
	using size_type		= typename unsigned int;
	using value_type	= typename Elem;

	GenericMatrix();
	explicit GenericMatrix(Initialization);
	GenericMatrix(size_type row, size_type col);
	explicit GenericMatrix(size_type col, size_type row, Initialization);
	GenericMatrix(const GenericMatrix &other);
	GenericMatrix(GenericMatrix &&other);
	~GenericMatrix();
	GenericMatrix &operator=(const GenericMatrix &other);
	GenericMatrix &operator=(GenericMatrix &&other);

	size_type rows() const noexcept;
	size_type columns() const noexcept;

	bool isValid() const noexcept;
	bool isIdentity() const;
	inline bool isHomomorphicTo(const GenericMatrix<Elem> &m);
	static inline bool isHomomorphic(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);

	void setToIdentity();
	void fill(const Elem &value);
	void swap(GenericMatrix<Elem> &other);
	GenericMatrix<Elem> transposed() const;
	void doHadamardProduct(const GenericMatrix<Elem> &m);
	static GenericMatrix<Elem> hadamardProduct(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);

	Elem &operator()(size_type row, size_type col);
	const Elem &operator()(size_type row, size_type col) const;

	Elem &at(size_type row, size_type col);
	const Elem &at(size_type row, size_type col) const;

	GenericMatrix<Elem> &operator+=(const GenericMatrix<Elem> &m);
	GenericMatrix<Elem> &operator-=(const GenericMatrix<Elem> &m);
	GenericMatrix<Elem> &operator*=(const GenericMatrix<Elem> &m);
	GenericMatrix<Elem> &operator*=(Elem factor);
	GenericMatrix<Elem> &operator/=(Elem divisor);
	bool operator==(const GenericMatrix<Elem> &m) const;
	bool operator!=(const GenericMatrix<Elem> &m) const;

	template<typename Elem>
	friend GenericMatrix<Elem> operator+(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);
	template<typename Elem>
	friend GenericMatrix<Elem> operator-(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);
	template<typename Elem>
	friend GenericMatrix<Elem> operator*(const GenericMatrix<Elem> &m1, const GenericMatrix<Elem> &m2);
	template<typename Elem>
	friend GenericMatrix<Elem> operator*(const GenericMatrix<Elem> &matrix, Elem factor);
	template<typename Elem>
	friend GenericMatrix<Elem> operator*(Elem factor, const GenericMatrix<Elem> &matrix);
	template<typename Elem>
	friend GenericMatrix<Elem> operator/(const GenericMatrix<Elem> &matrix, Elem divisor);
	template<typename Elem>
	friend GenericMatrix<Elem> operator-(const GenericMatrix<Elem> &matrix);
	template<typename Elem>
	friend std::ostream &operator<<(std::ostream &os, const GenericMatrix<Elem> &matrix);

private:
	void alloc();

private:
	size_type _rows;
	size_type _cols;
	Elem **_p;
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

	Constructs a 1x1 null matrix.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix()
	: _rows(1), _cols(1)
{
	alloc();
	_p[0][0] = 0;
}

/*!
	\fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix(Initialization)

	Constructs a 1x1 matrix without initializing the contents.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(Initialization)
	: _rows(1), _cols(1)
{
	alloc();
}

/*!
	\fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col)

	Constructs a \a row x \a col identity matrix.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(size_type row, size_type col)
	: _rows(row), _cols(col)
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
	: _rows(row), _cols(col)
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
	if (_p) {
		for (size_type i = 0; i < _rows; ++i) {
			if (_p[i]) {
				delete[] _p[i];
			}
		}
		delete[] _p; _p = nullptr;
	}
}
	
/*!
	\fn template<typename Elem> GenericMatrix<Elem>::~GenericMatrix()
	\internal

	Alloc the matrix memory from \a _rows and \a _cols.
*/
template<typename Elem>
void GenericMatrix<Elem>::alloc()
{
	_p = new Elem*[_rows];
	for (size_type i = 0; i < _rows; ++i)
		_p[i] = new Elem[_cols];
}

/*!
	\fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix(const GenericMatrix &other)

	Constructs a copy of \a other.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(const GenericMatrix &other)
	: _rows(other._rows), _cols(other._cols)
{
	alloc();
	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; j < _cols; ++j) {
			_p[i][j] = other._p[i][j];
		}
	}
}

/*!
	\fn template<typename Elem> GenericMatrix<Elem>::GenericMatrix(GenericMatrix &&other)

	Move-constructs a GenericMatrix instance, making it point at the same object that \a other was pointing to.
*/
template<typename Elem>
GenericMatrix<Elem>::GenericMatrix(GenericMatrix &&other)
	: _rows(other._rows), _cols(other._cols), _p(other._p)
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

	if (_rows != other._rows || _cols != other._cols) {
		for (size_type i = 0; i < _rows; ++i) {
			delete[] _p[i];
		}
		delete[] _p;

		_rows = other._rows; _cols = other._cols;
		alloc();
	}

	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; j < _cols; ++j) {
			_p[i][j] = other._p[i][j];
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
	if (_p) {
		for (size_type i = 0; i < _rows; ++i) {
			if (_p[i])
				delete[] _p[i];
		}
		delete[] _p; _p = nullptr;
	}

	_rows = other._rows;  _cols = other._cols;
	other._rows = 0; other._cols = 0;
	_p = other._p; other._p = nullptr;

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
	return _rows;
}

/*!
	\fn template<typename Elem> GenericMatrix<Elem>::size_type GenericMatrix<Elem>::columns() const

	Returns the number of matrix columns.

	\sa rows()
*/
template<typename Elem>
typename GenericMatrix<Elem>::size_type GenericMatrix<Elem>::columns() const noexcept
{
	return _cols;
}

/*!
	\fn template<typename Elem> bool GenericMatrix<Elem>::isValid() const

	Returns \c true if this matrix internal data is null pointer, otherwise returns \c false.
*/
template<typename Elem>
bool GenericMatrix<Elem>::isValid() const noexcept
{
	return _p;
}

/*!
	\fn template<typename Elem> bool GenericMatrix<Elem>::isIdentity() const

	Returns \c true if this matrix is the identity; false otherwise.

	\sa setToIdentity()
*/
template<typename Elem>
bool GenericMatrix<Elem>::isIdentity() const
{
	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; j < _cols; ++j) {
			if (i == j) {
				if (_p[i][j] != 1)
					return false;
			} else {
				if (_p[i][j] != 0)
					return false;
			}
		}
	}
	return true;
}

/*!
	\fn template<typename Elem> void GenericMatrix<Elem>::setToIdentity()

	Sets this matrix to the identity.

	\sa isIdentity()
*/
template<typename Elem>
void GenericMatrix<Elem>::setToIdentity()
{
	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; j < _cols; ++j) {
			if (j == i)
				_p[i][j] = 1;
			else
				_p[i][j] = 0;
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
	GenericMatrix<Elem> result(_cols, _rows, Uninitialized);
	for (int i = 0; i < _cols; ++i)
		for (int j = 0; j < _rows; ++j)
			result._p[i][j] = _p[j][i];
	return result;
}

/*!
	\fn template<typename Elem> void GenericMatrix<Elem>::fill(const Elem &value)

	Fills all elements of this matrix with \a value.
*/
template<typename Elem>
void GenericMatrix<Elem>::fill(const Elem &value)
{
	for (size_type i = 0; i < _rows; ++i)
		for (size_type j = 0; j < _cols; ++j)
			_p[i][j] = value;
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

	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; i < _cols; ++j) {
			_p[i][j] *= m._m[i][j];
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
	\fn template<typename Elem> bool GenericMatrix<Elem>::isHomomorphicTo(const GenericMatrix<Elem> &m)

	Returns \c true if this matrix and matrix \a m are homomorphic; false otherwise.

	\sa isHomomorphic()
*/
template<typename Elem>
inline bool GenericMatrix<Elem>::isHomomorphicTo(const GenericMatrix<Elem> &m)
{
	return((_rows == m._rows) && (_cols == m._cols));
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
	std::swap(_p, other._p);
	std::swap(_rows, other._rows);
	std::swap(_cols, other._cols);
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
	return _p[row][column];
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
	return _p[row][column];
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
	if (row >= 0 && row < _rows && col >= 0 && col < _cols)
		return _p[row][col];
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
	if (row >= 0 && row < _rows && col >= 0 && col < _cols)
		return _p[row][col];
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

	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; j < _cols; ++j) {
			_p[i][j] += m._p[i][j];
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

	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; j < _cols; ++j) {
			_p[i][j] -= m._p[i][j];
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
	if (_cols != m._rows)
		throw std::invalid_argument("\fn[GenericMatrix::operator*=()] can't do multiplication operation.");

	GenericMatrix<Elem> result(_rows, m._cols, Uninitialized);
	result.fill(0);
	for (size_type i = 0; i < result._rows; ++i) {
		for (size_type j = 0; j < result._cols; ++j) {
			for (size_type k = 0; k < _cols; ++k) {
				result._p[i][j] += (_p[i][k] * m._p[k][j]);
			}
		}
	}
	return (*this = std::move(result));
}

/*!
	\fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator*=(Elem factor)

	Multiplies all elements of this matrix by \a factor.

	\sa operator/=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator*=(Elem factor)
{
	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; j < _cols; ++j) {
			_p[i][j] *= factor;
		}
	}
	return *this;
}

/*!
	\fn template<typename Elem> GenericMatrix<Elem> &GenericMatrix<Elem>::operator/=(Elem divisor)

	Divides all elements of this matrix by \a divisor.

	\sa operator*=()
*/
template<typename Elem>
GenericMatrix<Elem> &GenericMatrix<Elem>::operator/=(Elem divisor)
{
	for (size_type i = 0; i < _rows; ++i) {
		for (size_type j = 0; j < _cols; ++j) {
			_p[i][j] /= divisor;
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
	for (size_type i = 0; i < _rows; ++i)
		for (int j = 0; j < _cols; ++j) {
			if (_p[i][j] != m._p[i][j])
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
	for (size_type i = 0; i < _rows; ++i)
		for (int j = 0; j < _cols; ++j) {
			if (_p[i][j] != m._p[i][j])
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
	GenericMatrix<Elem> result(m1._rows, m1._cols, Uninitialized);
	for (size_type i = 0; i < result._rows; ++i) {
		for (size_type j = 0; j < result._cols; ++j) {
			result._p[i][j] = m1._p[i][j] + m2._p[i][j];
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
	GenericMatrix<Elem> result(m1._rows, m1._cols, Uninitialized);
	for (size_type i = 0; i < result._rows; ++i) {
		for (size_type j = 0; j < result._cols; ++j) {
			result._p[i][j] = m1._p[i][j] - m2._p[i][j];
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
	if (m1._cols != m2._rows)
		throw std::invalid_argument("input matrix can not do multiplication operation.");

	using size_type = typename GenericMatrix<Elem>::size_type;
	GenericMatrix<Elem> result(m1._rows, m2._cols, Uninitialized);
	result.fill(0);
	for (size_type i = 0; i < result._rows; ++i) {
		for (size_type j = 0; j < result._cols; ++j) {
			for (size_type k = 0; k < m1._cols; ++k) {
				result._p[i][j] += m1._p[i][k] * m2._p[k][j];
			}
		}
	}
	return result;
}

/*!
	\fn GenericMatrix<Elem> operator*(const GenericMatrix<Elem> &matrix, Elem factor)
	\relates GenericMatrix

	Returns the result of multiplying all elements of \a matrix by \a factor.
*/
template<typename Elem>
GenericMatrix<Elem> operator*(const GenericMatrix<Elem> &matrix, Elem factor)
{
	using size_type = typename GenericMatrix<Elem>::size_type;
	GenericMatrix<Elem> result(matrix._rows, matrix._cols, Uninitialized);
	for (size_type i = 0; i < result._rows; ++i) {
		for (size_type j = 0; j < result._cols; ++j) {
			result._p[i][j] = matrix._p[i][j] * factor;
		}
	}
	return result;
}

/*!
	\fn GenericMatrix<Elem> operator*(Elem factor, const GenericMatrix<Elem> &matrix)
	\relates GenericMatrix

	Returns the result of multiplying all elements of \a matrix by \a factor.
*/
template<typename Elem>
GenericMatrix<Elem> operator*(Elem factor, const GenericMatrix<Elem> &matrix)
{
	using size_type = typename GenericMatrix<Elem>::size_type;
	GenericMatrix<Elem> result(matrix._rows, matrix._cols, Uninitialized);
	for (size_type i = 0; i < result._rows; ++i) {
		for (size_type j = 0; j < result._cols; ++j) {
			result._p[i][j] = matrix._p[i][j] * factor;
		}
	}
	return result;
}

/*!
	\fn GenericMatrix<Elem> operator/(const GenericMatrix<Elem> &matrix, Elem divisor)
	\relates GenericMatrix

	Returns the result of dividing all elements of \a matrix by \a divisor.
*/
template<typename Elem>
GenericMatrix<Elem> operator/(const GenericMatrix<Elem> &matrix, Elem divisor)
{
	using size_type = typename GenericMatrix<Elem>::size_type;
	GenericMatrix<Elem> result(matrix._rows, matrix._cols, Uninitialized);
	for (size_type i = 0; i < result._rows; ++i) {
		for (size_type j = 0; j < result._cols; ++j) {
			result._p[i][j] = matrix._p[i][j] / divisor;
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
	GenericMatrix<Elem> result(matrix._rows, matrix._cols, Uninitialized);
	for (size_type i = 0; i < matrix._rows; ++i)
		for (size_type j = 0; j < matrix._cols; ++j)
			result._p[i][j] = -matrix._p[i][j];
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
	stream << "GenericMatrix<" << matrix._rows << ", " << matrix._cols << ", " << typeid(Elem).name() << ">(" << std::endl;
	for (size_type i = 0; i < matrix._rows; ++i) {
		for (size_type j = 0; j < matrix._cols; ++j)
			stream << stream.width(10) << matrix._p[i][j];
		stream << std::endl;
	}
	stream << ')';
	return stream;
}

#endif // !__GENERICMATRIX_H__
