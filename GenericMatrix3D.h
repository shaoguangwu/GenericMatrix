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

    GenericMatrix3D()
        : m_layers(0), m_rows(0), m_cols(0), m_data(nullptr)
    {
    }

    GenericMatrix3D(size_type layer, size_type row, size_type col)
        : m_layers(layer), m_rows(row), m_cols(col), m_data(nullptr)
    {
        alloc();
    }

    GenericMatrix3D(size_type layer, size_type row, size_type col, const Elem &initialValue)
        : m_layers(layer), m_rows(row), m_cols(col), m_data(nullptr)
    {
        alloc();
        fill(initialValue);
    }

    GenericMatrix3D(size_type layer, size_type row, size_type col, Elem *data)
        : m_layers(layer), m_rows(row), m_cols(col), m_data(data)
    {   
    }

    GenericMatrix3D(const GenericMatrix3D &other)
        : m_layers(other.m_layers), m_rows(other.m_rows), m_cols(other.m_cols), m_data(nullptr)
    {
        alloc();
        std::copy_n(other.data(), other.size(), m_data);
    }

    GenericMatrix3D(GenericMatrix3D &&other) noexcept
        : m_layers(other.m_layers), m_rows(other.m_rows), m_cols(other.m_cols), m_data(other.data())
    {
        other.data() = nullptr;
        other.m_layers = 0;
        other.rows = 0;
        other.cols = 0;
    }

    ~GenericMatrix3D()
    {
        free();
    }

    GenericMatrix3D &operator=(const GenericMatrix3D &other)
    {
        if (*this == other){
            return *this;
        }
        resize(other.m_layers, other.m_rows, other.m_cols);
        std::copy_n(other.data(), other.size(), m_data);
    }

    GenericMatrix3D &operator=(GenericMatrix3D &&other) noexcept
    {
        if (*this == other){
            return *this;
        }
        setSize(other.m_layers, other.m_rows, other.m_cols);
        m_data = other.data();
        other.m_data = nullptr;
        other.setSize(0, 0, 0);
    }

    iterator begin() noexcept
    {
        return m_data;
    }
    const_iterator begin() const noexcept
    {
        return m_data;
    }
    const_iterator cbegin() const noexcept
    {
        return m_data;
    }
    iterator end() noexcept
    {
        return m_data + size();
    }
    const_iterator end() const noexcept
    {
        return m_data + size();
    }
    const_iterator cend() const noexcept
    {
        return m_data + size();
    }

    Elem *data() noexcept
    {
        return m_data;
    }

    const Elem *data() const noexcept
    {
        return m_data;
    }

    const Elem *constData() const noexcept
    {
        return m_data;
    }

    size_type size() const
    {
        return m_layers * m_rows * m_cols;
    }

    size_type oneSize() const
    {
        return m_rows * m_cols;
    }

    bool empty() const
    {
        return 0 == size();
    }
    bool isValid() const
    {
        return m_data && size() > 0;
    }

    Elem &operator()(size_type layer, size_type row, size_type col)
    {
        return m_data[computeOffset(layer, row, col)];
    }
    const Elem &operator()(size_type layer, size_type row, size_type col) const
    {
        return m_data[computeOffset(layer, row, col)];
    }
    Elem &at(size_type layer, size_type row, size_type col)
    {
        if (layer < m_layers && row < m_rows && col < m_cols)
            return m_data[computeOffset(layer, row, col)];
        throw std::out_of_range("out of range.");
    }
    const Elem &at(size_type layer, size_type row, size_type col) const
    {
        if (layer < m_layers && row < m_rows && col < m_cols)
            return m_data[computeOffset(layer, row, col)];
        throw std::out_of_range("out of range.");
    }

    void fill(const Elem &initialValue)
    {
        std::fill_n(m_data, size(), initialValue);
    }

    void swap(GenericMatrix3D &other)
    {
        std::swap(*this, other);
    }

    void resize(size_type layer, size_type row, size_type col)
    {
        if (size() != layer * row * col) {
            free();
            alloc(layer, row, col);
        }
        setSize(layer, row, col);
    }

    void resize(size_type layer, size_type row, size_type col, const Elem &initialValue)
    {
        resize(layer, row, col);
        fill(initialValue);
    }

    void reset(size_type layer, size_type row, size_type col, Elem *data)
    {
        free();
        setSize(layer, row, col);
        m_data = data;
    }

    size_type layers() const
    {
        return m_layers;
    }
    size_type rows() const
    {
        return m_rows;
    }
    size_type cols() const
    {
        return m_cols;
    }


private:
    void alloc()
    {
        if (size() > 0)
        {
            m_data = new value_type[size()];
        }
    }
    void alloc(size_type layers, size_type rows, size_type cols)
    {
        size_type size = layers * rows * cols;
        if (size > 0) {
            m_data = new value_type[size];
        }
    }
    void free()
    {
        if (m_data)
        {
            delete[] m_data;
            m_data = nullptr;
        }
    }
    void setSize(size_type layers, size_type rows, size_type cols)
    {
        m_layers = layers;
        m_rows = rows;
        m_cols = cols;
    }

    size_type computeOffset(size_type layer, size_type row, size_type col) const
    {
        return ((layer > 0 ? layer * oneSize() : 0) + row * m_cols + col);
    }

    bool containsIndex(size_type layer, size_type row, size_type col) const
    {
        return (layer < m_layers && row < m_rows && col < m_cols);
    }

private:
    size_type m_layers;
    size_type m_rows;
    size_type m_cols;
    value_type *m_data;
};

#endif // __GERNERICMATRIX3D_H__