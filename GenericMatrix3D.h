/****************************************************************************
**
** Copyright (C) 2019 Shaoguang. All rights reserved.
** 
** GenericMatrix3D, a generic template matrix class for three-dimension,
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

#ifndef __GERNERICMATRIX3D_H__
#define __GERNERICMATRIX3D_H__

#include <ostream>
#include <utility>      // std::move,   std::pair
#include <algorithm>    // std::fill_n, std::copy_n

#include "GenericMatrix.h"

/*!
    \class GenericMatrix3D
    \brief The GenericMatrix3D class defines a generic template matrix class for three-dimension,
    the matrix elements are managed by a one-dimension pointer.

    The GenericMatrix3D template has one parameter:
    \li \b Elem Element type that is visible to users of the class.
*/
template<typename Elem = float>
class GenericMatrix3D
{
public:
    using size_type = std::size_t;
    using value_type = Elem;
    using iterator = value_type *;
    using const_iterator = const iterator;

    GenericMatrix3D();
    GenericMatrix3D(size_type layer, size_type row, size_type col);
    GenericMatrix3D(size_type layer, size_type row, size_type col, const Elem &initialValue);
    GenericMatrix3D(size_type layer, size_type row, size_type col, Elem *data);
    GenericMatrix3D(const GenericMatrix3D &other);
    GenericMatrix3D(GenericMatrix3D &&other) noexcept;
    ~GenericMatrix3D();
    GenericMatrix3D &operator=(const GenericMatrix3D &other);
    GenericMatrix3D &operator=(GenericMatrix3D &&other) noexcept;
    iterator begin() noexcept;
    const_iterator begin() const noexcept;
    const_iterator cbegin() const noexcept;
    iterator end() noexcept;
    const_iterator end() const noexcept;
    const_iterator cend() const noexcept;
    Elem *data() noexcept;
    const Elem *data() const noexcept;
    const Elem *constData() const noexcept;
    Elem *data(size_type layer, size_type row, size_type col) noexcept;
    const Elem *data(size_type layer, size_type row, size_type col) const noexcept;
    const Elem *constData(size_type layer, size_type row, size_type col) const noexcept;
    size_type size() const;
    size_type oneSize() const;
    bool empty() const;
    bool isValid() const;
    Elem &operator()(size_type layer, size_type row, size_type col);
    const Elem &operator()(size_type layer, size_type row, size_type col) const;
    Elem &at(size_type layer, size_type row, size_type col);
    const Elem &at(size_type layer, size_type row, size_type col) const;
    void fill(const Elem &initialValue);
    void swap(GenericMatrix3D &other);
    void resize(size_type layer, size_type row, size_type col);
    void resize(size_type layer, size_type row, size_type col, const Elem &initialValue);
    void reset(size_type layer, size_type row, size_type col, Elem *data);
    size_type layers() const;
    size_type rows() const;
    size_type cols() const;
    GenericMatrix3D<Elem> roi(size_type layer1, size_type row1, size_type col1, size_type layer2, size_type row2, size_type col2);

    // GenericMatrix<Elem> oneOf(size_type layer);
    // GenericMatrix<Elem> oneOfAt(size_type layer);
    GenericMatrix<Elem> copyOne(size_type layer);
    GenericMatrix<Elem> copyOneAt(size_type layer);

private:
    void alloc();

    void alloc(size_type layers, size_type rows, size_type cols);
    void free();
    void setSize(size_type layers, size_type rows, size_type cols);

    size_type computeOffset(size_type layer, size_type row, size_type col) const;

    bool containsIndex(size_type layer, size_type row, size_type col) const;

    template<typename __Elem>
    friend std::ostream &operator<<(std::ostream &stream, const GenericMatrix3D<__Elem> &matrix);

private:
    size_type m_layers;
    size_type m_rows;
    size_type m_cols;
    value_type *m_data;
};

/*!
    \var GenericMatrix3D::m_layers

    The number of layers of this three-dimension matrix.
*/

/*!
    \var GenericMatrix3D::m_rows

    The number of rows of this three-dimension matrix.
*/

/*!
    \var GenericMatrix3D::m_cols

    The number of cols of this three-dimension matrix.
*/

template<typename Elem>
GenericMatrix3D<Elem>::GenericMatrix3D()
    : m_layers(0), m_rows(0), m_cols(0), m_data(nullptr)
{
}

template<typename Elem>
GenericMatrix3D<Elem>::GenericMatrix3D(size_type layer, size_type row, size_type col)
    : m_layers(layer), m_rows(row), m_cols(col), m_data(nullptr)
{
    alloc();
}

template<typename Elem>
GenericMatrix3D<Elem>::GenericMatrix3D(size_type layer, size_type row, size_type col, const Elem &initialValue)
    : m_layers(layer), m_rows(row), m_cols(col), m_data(nullptr)
{
    alloc();
    fill(initialValue);
}

template<typename Elem>
GenericMatrix3D<Elem>::GenericMatrix3D(size_type layer, size_type row, size_type col, Elem *data)
    : m_layers(layer), m_rows(row), m_cols(col), m_data(data)
{   
}
template<typename Elem>
GenericMatrix3D<Elem>::GenericMatrix3D(const GenericMatrix3D &other)
    : m_layers(other.m_layers), m_rows(other.m_rows), m_cols(other.m_cols), m_data(nullptr)
{
    alloc();
    std::copy_n(other.data(), other.size(), m_data);
}
template<typename Elem>
GenericMatrix3D<Elem>::GenericMatrix3D(GenericMatrix3D &&other) noexcept
    : m_layers(other.m_layers), m_rows(other.m_rows), m_cols(other.m_cols), m_data(other.data())
{
    other.m_data = nullptr;
    other.setSize(0, 0, 0);
}
template<typename Elem>
GenericMatrix3D<Elem>::~GenericMatrix3D()
{
    free();
}
template<typename Elem>
GenericMatrix3D<Elem> &GenericMatrix3D<Elem>::operator=(const GenericMatrix3D &other)
{
    if (this == &other){
        return *this;
    }
    resize(other.m_layers, other.m_rows, other.m_cols);
    std::copy_n(other.data(), other.size(), m_data);
}
template<typename Elem>
GenericMatrix3D<Elem> &GenericMatrix3D<Elem>::operator=(GenericMatrix3D &&other) noexcept
{
    if (this == &other){
        return *this;
    }
    setSize(other.m_layers, other.m_rows, other.m_cols);
    m_data = other.data();
    other.m_data = nullptr;
    other.setSize(0, 0, 0);
}
template<typename Elem>
typename GenericMatrix3D<Elem>::iterator GenericMatrix3D<Elem>::begin() noexcept
{
    return m_data;
}
template<typename Elem>
typename GenericMatrix3D<Elem>::const_iterator GenericMatrix3D<Elem>::begin() const noexcept
{
    return m_data;
}
template<typename Elem>
typename GenericMatrix3D<Elem>::const_iterator GenericMatrix3D<Elem>::cbegin() const noexcept
{
    return m_data;
}
template<typename Elem>
typename GenericMatrix3D<Elem>::iterator GenericMatrix3D<Elem>::end() noexcept
{
    return m_data + size();
}
template<typename Elem>
typename GenericMatrix3D<Elem>::const_iterator GenericMatrix3D<Elem>::end() const noexcept
{
    return m_data + size();
}
template<typename Elem>
typename GenericMatrix3D<Elem>::const_iterator GenericMatrix3D<Elem>::cend() const noexcept
{
    return m_data + size();
}
template<typename Elem>
Elem *GenericMatrix3D<Elem>::data() noexcept
{
    return m_data;
}
template<typename Elem>
const Elem *GenericMatrix3D<Elem>::data() const noexcept
{
    return m_data;
}
template<typename Elem>
const Elem *GenericMatrix3D<Elem>::constData() const noexcept
{
    return m_data;
}
template<typename Elem>
Elem *GenericMatrix3D<Elem>::data(size_type layer, size_type row, size_type col) noexcept
{
    return m_data + computeOffset(layer, row, col);
}
template<typename Elem>
const Elem *GenericMatrix3D<Elem>::data(size_type layer, size_type row, size_type col) const noexcept
{
    return m_data + computeOffset(layer, row, col);
}
template<typename Elem>
const Elem *GenericMatrix3D<Elem>::constData(size_type layer, size_type row, size_type col) const noexcept
{
    return m_data + computeOffset(layer, row, col);
}
template<typename Elem>
typename GenericMatrix3D<Elem>::size_type GenericMatrix3D<Elem>::size() const
{
    return m_layers * m_rows * m_cols;
}
template<typename Elem>
typename GenericMatrix3D<Elem>::size_type GenericMatrix3D<Elem>::oneSize() const
{
    return m_rows * m_cols;
}
template<typename Elem>
bool GenericMatrix3D<Elem>::empty() const
{
    return 0 == size();
}

template<typename Elem>
bool GenericMatrix3D<Elem>::isValid() const
{
    return m_data && size() > 0;
}
template<typename Elem>
Elem &GenericMatrix3D<Elem>::operator()(size_type layer, size_type row, size_type col)
{
    return m_data[computeOffset(layer, row, col)];
}
template<typename Elem>
const Elem &GenericMatrix3D<Elem>::operator()(size_type layer, size_type row, size_type col) const
{
    return m_data[computeOffset(layer, row, col)];
}
template<typename Elem>
Elem &GenericMatrix3D<Elem>::at(size_type layer, size_type row, size_type col)
{
    if (layer < m_layers && row < m_rows && col < m_cols)
        return m_data[computeOffset(layer, row, col)];
    throw std::out_of_range("out of range.");
}
template<typename Elem>
const Elem &GenericMatrix3D<Elem>::at(size_type layer, size_type row, size_type col) const
{
    if (layer < m_layers && row < m_rows && col < m_cols)
        return m_data[computeOffset(layer, row, col)];
    throw std::out_of_range("out of range.");
}
template<typename Elem>
void GenericMatrix3D<Elem>::fill(const Elem &initialValue)
{
    std::fill_n(m_data, size(), initialValue);
}
template<typename Elem>
void GenericMatrix3D<Elem>::swap(GenericMatrix3D &other)
{
    std::swap(*this, other);
}
template<typename Elem>
void GenericMatrix3D<Elem>::resize(size_type layer, size_type row, size_type col)
{
    if (size() != layer * row * col) {
        free();
        alloc(layer, row, col);
    }
    setSize(layer, row, col);
}
template<typename Elem>
void GenericMatrix3D<Elem>::resize(size_type layer, size_type row, size_type col, const Elem &initialValue)
{
    resize(layer, row, col);
    fill(initialValue);
}
template<typename Elem>
void GenericMatrix3D<Elem>::reset(size_type layer, size_type row, size_type col, Elem *data)
{
    free();
    setSize(layer, row, col);
    m_data = data;
}
template<typename Elem>
typename GenericMatrix3D<Elem>::size_type GenericMatrix3D<Elem>::layers() const
{
    return m_layers;
}
template<typename Elem>
typename GenericMatrix3D<Elem>::size_type GenericMatrix3D<Elem>::rows() const
{
    return m_rows;
}
template<typename Elem>
typename GenericMatrix3D<Elem>::size_type GenericMatrix3D<Elem>::cols() const
{
    return m_cols;
}

template<typename Elem>
GenericMatrix3D<Elem> GenericMatrix3D<Elem>::roi(size_type layer1, size_type row1, size_type col1, size_type layer2, size_type row2, size_type col2)
{
    if (!containsIndex(layer1, row1, col1) || !containsIndex(layer2, row2, col2)) {
        return GenericMatrix3D<Elem>(); 
    }
    std::pair<size_type, size_type> layerMinMax = std::minmax(layer1, layer2);
    std::pair<size_type, size_type> rowMinMax = std::minmax(row1, row2);
    std::pair<size_type, size_type> colMinMax = std::minmax(col1, col2);
    GenericMatrix3D<Elem> result(layerMinMax.second - layerMinMax.first + 1, rowMinMax.second - rowMinMax.first + 1, colMinMax.second - colMinMax.first + 1);

    for (size_type i = layerMinMax.first; i <= layerMinMax.second; ++i) {
        for (size_type j = rowMinMax.first; j <= rowMinMax.second; ++j) {
            std::copy_n(constData(i, j, colMinMax.first), result.cols(), result.data(i - layerMinMax.first, j - rowMinMax.first, 0));
        }
    }
    return result;
}

template<typename Elem>
GenericMatrix<Elem> GenericMatrix3D<Elem>::copyOne(size_type layer)
{
    GenericMatrix<Elem> result(m_rows, m_cols);
    std::copy_n(begin() + computeOffset(layer, 0, 0), result.size(), result.data());
    return result;
}

template<typename Elem>
GenericMatrix<Elem> GenericMatrix3D<Elem>::copyOneAt(size_type layer)
{
    if (layer < m_layers) {
        return GenericMatrix<Elem>();
    }
    GenericMatrix<Elem> result(m_rows, m_cols);
    std::copy_n(begin() + computeOffset(layer, 0, 0), result.size(), result.data());
    return result;
}

template<typename Elem>
void GenericMatrix3D<Elem>::alloc()
{
    if (size() > 0)
    {
        m_data = new value_type[size()];
    }
}
template<typename Elem>
void GenericMatrix3D<Elem>::alloc(size_type layers, size_type rows, size_type cols)
{
    size_type size = layers * rows * cols;
    if (size > 0) {
        m_data = new value_type[size];
    }
}
template<typename Elem>
void GenericMatrix3D<Elem>::free()
{
    if (m_data)
    {
        delete[] m_data;
        m_data = nullptr;
    }
}
template<typename Elem>
void GenericMatrix3D<Elem>::setSize(size_type layers, size_type rows, size_type cols)
{
    m_layers = layers;
    m_rows = rows;
    m_cols = cols;
}
template<typename Elem>
typename GenericMatrix3D<Elem>::size_type GenericMatrix3D<Elem>::computeOffset(size_type layer, size_type row, size_type col) const
{
    return ((layer > 0 ? layer * oneSize() : 0) + row * m_cols + col);
}
template<typename Elem>
bool GenericMatrix3D<Elem>::containsIndex(size_type layer, size_type row, size_type col) const
{
    return (layer < m_layers && row < m_rows && col < m_cols);
}

template<typename __Elem>
std::ostream &operator<<(std::ostream &stream, const GenericMatrix3D<__Elem> &matrix)
{
    using size_type = typename GenericMatrix3D<__Elem>::size_type;
    stream << "GenericMatrix3D<" << matrix.m_layers << ", " << matrix.m_rows << ", " << matrix.m_cols << ", " << typeid(__Elem).name() << ">(" << std::endl;
    for (size_type l = 0; l < matrix.m_layers; ++l) {
        for (size_type r = 0; r < matrix.m_rows; ++r) {
            for (size_type c = 0; c < matrix.m_cols; ++c) {
                stream << stream.width(10) << matrix(l, r, c);
            }
            stream << std::endl;
        }
        stream << std::endl;
    }
    stream << ')';
    return stream;
}

#endif // __GERNERICMATRIX3D_H__