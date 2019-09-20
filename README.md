# GenericMatrix

C++ cross platform generic matrix template class, C++11 Required Support.

C++跨平台通用模板矩阵类，需编译器支持C++11.

## Quick Example

linux or windows:

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

## License

[Apache License 2.0](LICENSE)
