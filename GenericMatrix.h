/*!
 * \file	GenericMatrix.h
 * \brief	The header file of the GenericMatrix class.
 * \date	2019-08-22
 * \author	shaoguang
 */

#ifndef GENERIC_MATRIX_H
#define GENERIC_MATRIX_H
 
/*!
	\enum Initialization
	\brief Indicates uninitialized.
	\note enum class
*/
enum class Initialization {
	Uninitialized				///< Indicates uninitialized.
};
static constexpr Initialization Uninitialized = Initialization::Uninitialized;

/*!
	\class GenericMatrix
*/
template <int COL, int ROW, typename T>
class GenericMatrix
{
public:
	using size_t	= typename int;
	using value_t	= typename T;

	GenericMatrix() { setToIdentity(); }
	explicit GenericMatrix(Initialization) {}
	explicit GenericMatrix(const T *values);

	const T& operator()(size_t row, size_t column) const;
	T& operator()(size_t row, size_t column);

	bool isIdentity() const;
	void setToIdentity();

	void fill(T value);

	GenericMatrix<ROW, COL, T> transposed() const;

	GenericMatrix<COL, ROW, T>& operator+=(const GenericMatrix<COL, ROW, T>& other);
	GenericMatrix<COL, ROW, T>& operator-=(const GenericMatrix<COL, ROW, T>& other);
	GenericMatrix<COL, ROW, T>& operator*=(T factor);
	GenericMatrix<COL, ROW, T>& operator/=(T divisor);
	bool operator==(const GenericMatrix<COL, ROW, T>& other) const;
	bool operator!=(const GenericMatrix<COL, ROW, T>& other) const;

	void copyDataTo(T *values) const;

	T *data() { return *m; }
	const T *data() const { return *m; }
	const T *constData() const { return *m; }

	template<int NN, int MM, typename TT>
	friend GenericMatrix<NN, MM, TT> operator+(const GenericMatrix<NN, MM, TT>& m1, const GenericMatrix<NN, MM, TT>& m2);
	template<int NN, int MM, typename TT>
	friend GenericMatrix<NN, MM, TT> operator-(const GenericMatrix<NN, MM, TT>& m1, const GenericMatrix<NN, MM, TT>& m2);
	template<int NN, int M1, int M2, typename TT>
	friend GenericMatrix<M1, M2, TT> operator*(const GenericMatrix<NN, M2, TT>& m1, const GenericMatrix<M1, NN, TT>& m2);
	template<int NN, int MM, typename TT>
	friend GenericMatrix<NN, MM, TT> operator-(const GenericMatrix<NN, MM, TT>& matrix);
	template<int NN, int MM, typename TT>
	friend GenericMatrix<NN, MM, TT> operator*(TT factor, const GenericMatrix<NN, MM, TT>& matrix);
	template<int NN, int MM, typename TT>
	friend GenericMatrix<NN, MM, TT> operator*(const GenericMatrix<NN, MM, TT>& matrix, TT factor);
	template<int NN, int MM, typename TT>
	friend GenericMatrix<NN, MM, TT> operator/(const GenericMatrix<NN, MM, TT>& matrix, TT divisor);

private:
	T m[COL][ROW];    // Column-major order to match OpenGL.
};

// Define aliases for the useful variants of GenericMatrix.
typedef GenericMatrix<2, 2, float> Matrix2x2;
typedef GenericMatrix<2, 3, float> Matrix2x3;
typedef GenericMatrix<2, 4, float> Matrix2x4;
typedef GenericMatrix<3, 2, float> Matrix3x2;
typedef GenericMatrix<3, 3, float> Matrix3x3;
typedef GenericMatrix<3, 4, float> Matrix3x4;
typedef GenericMatrix<4, 2, float> Matrix4x2;
typedef GenericMatrix<4, 3, float> Matrix4x3;

//===================================================================================

template<int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T>::GenericMatrix(const T *values)
{
	for (int col = 0; col < COL; ++col)
		for (int row = 0; row < ROW; ++row)
			m[col][row] = values[row * COL + col];
}

template <int COL, int ROW, typename T>
const T& GenericMatrix<COL, ROW, T>::operator()(size_t row, size_t column) const
{
	if (!(row >= 0 && row < ROW && column >= 0 && column < COL))
		throw std::out_of_range("\fn[GenericMatrix::operator()] out of range.");
	return m[column][row];
}

template <int COL, int ROW, typename T>
T& GenericMatrix<COL, ROW, T>::operator()(size_t row, size_t column)
{
	if (!(row >= 0 && row < ROW && column >= 0 && column < COL))
		throw std::out_of_range("\fn[GenericMatrix::operator()] out of range.");
	return m[column][row];
}

template <int COL, int ROW, typename T>
bool GenericMatrix<COL, ROW, T>::isIdentity() const
{
	for (int col = 0; col < COL; ++col) {
		for (int row = 0; row < ROW; ++row) {
			if (row == col) {
				if (m[col][row] != 1.0f)
					return false;
			}
			else {
				if (m[col][row] != 0.0f)
					return false;
			}
		}
	}
	return true;
}

template <int COL, int ROW, typename T>
void GenericMatrix<COL, ROW, T>::setToIdentity()
{
	for (int col = 0; col < COL; ++col) {
		for (int row = 0; row < ROW; ++row) {
			if (row == col)
				m[col][row] = 1.0f;
			else
				m[col][row] = 0.0f;
		}
	}
}

template <int COL, int ROW, typename T>
void GenericMatrix<COL, ROW, T>::fill(T value)
{
	for (int col = 0; col < COL; ++col)
		for (int row = 0; row < ROW; ++row)
			m[col][row] = value;
}

template <int COL, int ROW, typename T>
GenericMatrix<ROW, COL, T> GenericMatrix<COL, ROW, T>::transposed() const
{
	GenericMatrix<ROW, COL, T> result(dip::Uninitialized);
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			result.m[row][col] = m[col][row];
	return result;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T>& GenericMatrix<COL, ROW, T>::operator+=(const GenericMatrix<COL, ROW, T>& other)
{
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			m[col][row] += other.m[col][row];
	return *this;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T>& GenericMatrix<COL, ROW, T>::operator-=(const GenericMatrix<COL, ROW, T>& other)
{
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			m[col][row] -= other.m[col][row];
	return *this;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T>& GenericMatrix<COL, ROW, T>::operator*=(T factor)
{
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			m[col][row] *= factor;
	return *this;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T>& GenericMatrix<COL, ROW, T>::operator/=(T divisor)
{
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			m[col][row] /= divisor;
	return *this;
}

template <int COL, int ROW, typename T>
bool GenericMatrix<COL, ROW, T>::operator==(const GenericMatrix<COL, ROW, T>& other) const
{
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col) {
			if (m[col][row] != other.m[col][row])
				return false;
		}
	return true;
}

template <int COL, int ROW, typename T>
bool GenericMatrix<COL, ROW, T>::operator!=(const GenericMatrix<COL, ROW, T>& other) const
{
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col) {
			if (m[col][row] != other.m[col][row])
				return true;
		}
	return false;
}

template <int COL, int ROW, typename T>
void GenericMatrix<COL, ROW, T>::copyDataTo(T *values) const
{
	for (int col = 0; col < COL; ++col)
		for (int row = 0; row < ROW; ++row)
			values[row * COL + col] = T(m[col][row]);
}

//\\====================================================================================//\\

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T> operator+(const GenericMatrix<COL, ROW, T>& m1, const GenericMatrix<COL, ROW, T>& m2)
{
	GenericMatrix<COL, ROW, T> result(dip::Uninitialized);
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			result.m[col][row] = m1.m[col][row] + m2.m[col][row];
	return result;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T> operator-(const GenericMatrix<COL, ROW, T>& m1, const GenericMatrix<COL, ROW, T>& m2)
{
	GenericMatrix<COL, ROW, T> result(dip::Uninitialized);
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			result.m[col][row] = m1.m[col][row] - m2.m[col][row];
	return result;
}

template <int COL, int M1, int M2, typename T>
GenericMatrix<M1, M2, T> operator*(const GenericMatrix<COL, M2, T>& m1, const GenericMatrix<M1, COL, T>& m2)
{
	GenericMatrix<M1, M2, T> result(dip::Uninitialized);
	for (int row = 0; row < M2; ++row) {
		for (int col = 0; col < M1; ++col) {
			T sum(0.0f);
			for (int j = 0; j < COL; ++j)
				sum += m1.m[j][row] * m2.m[col][j];
			result.m[col][row] = sum;
		}
	}
	return result;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T> operator-(const GenericMatrix<COL, ROW, T>& matrix)
{
	GenericMatrix<COL, ROW, T> result(dip::Uninitialized);
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			result.m[col][row] = -matrix.m[col][row];
	return result;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T> operator*(T factor, const GenericMatrix<COL, ROW, T>& matrix)
{
	GenericMatrix<COL, ROW, T> result(dip::Uninitialized);
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			result.m[col][row] = matrix.m[col][row] * factor;
	return result;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T> operator*(const GenericMatrix<COL, ROW, T>& matrix, T factor)
{
	GenericMatrix<COL, ROW, T> result(dip::Uninitialized);
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			result.m[col][row] = matrix.m[col][row] * factor;
	return result;
}

template <int COL, int ROW, typename T>
GenericMatrix<COL, ROW, T> operator/(const GenericMatrix<COL, ROW, T>& matrix, T divisor)
{
	GenericMatrix<COL, ROW, T> result(dip::Uninitialized);
	for (int row = 0; row < ROW; ++row)
		for (int col = 0; col < COL; ++col)
			result.m[col][row] = matrix.m[col][row] / divisor;
	return result;
}

template <int COL, int ROW, typename T>
std::ostream &operator<<(std::ostream &os, const GenericMatrix<COL, ROW, T> &m)
{
	os << "GenericMatrix<" << COL << ", " << ROW << ", " << typeid(T).name() << ">(" << std::endl;
	for (int row = 0; row < ROW; ++row) {
		for (int col = 0; col < COL; ++col)
			os << os.width(10) << m(row, col);

		os << std::endl;
	}
	os << ')';
	return os;
}

#endif // !GENERIC_MATRIX_H
