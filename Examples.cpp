/**
 *  Copyright (C) 2019 Shaoguang. All rights reserved.
 * 
 *  Examples for the GernericMatrix class and the GernericMatrix2 class.
 *  
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *  
 *      https://www.apache.org/licenses/LICENSE-2.0
 *  
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include <iostream>
#include "GenericMatrix.h"
#include "GenericMatrix2.h"

int main(int argc, char *argv[])
{
    GenericMatrix<> matrix(2, 4);
    matrix.fill(1);
    
    GenericMatrix2<> matrix2(4, 4, 0);
    matrix2.setToIdentity();

    std::cout << matrix << std::endl;
    std::cout << matrix2 << std::endl;

    return 0;
}