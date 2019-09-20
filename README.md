# GenericMatrix

C++ generic matrix template class, C++11 Required Support.

C++通用矩阵类设计，需编译器支持C++11.

## Quick Example

```cpp
    GenericMatrix<> matrix(10, 5);
    std::cout << matrix << std::endl;
```

The Result is:

```cpp
GenericMatrix<10, 5, f>(
         01         00         00         00         00
         00         01         00         00         00
         00         00         01         00         00
         00         00         00         01         00
         00         00         00         00         01
         00         00         00         00         00
         00         00         00         00         00
         00         00         00         00         00
         00         00         00         00         00
         00         00         00         00         00
)
```
