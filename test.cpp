#include <iostream>
#include "GenericMatrix3D.h"

int main()
{
    using size = std::size_t;
    size layer = 4;
    size row = 4;
    size col = 4;

    int v = 1;
    GenericMatrix3D<int> matrix(layer, row, col);
    for (size l = 0; l < layer; ++l)
    {
        for (size r = 0; r < row; ++r)
        {
            for (size c = 0; c < col; ++c)
            {
                matrix(l, r, c) = v++;
            }
        }
    }

    std::cout << matrix << std::endl;

    auto matrix_roi = matrix.roi(0, 1, 0, 2, 3, 3);
    std::cout << matrix_roi << std::endl;

    matrix.swap(matrix_roi);
    std::cout << "after swap: " << matrix << std::endl;
    std::cout << matrix_roi << std::endl;

    return 0;
}