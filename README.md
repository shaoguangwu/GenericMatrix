# GenericMatrix

C++ cross platform generic matrix template class, C++11 Required Support.

More interface information see "doc/html/index.html".

C++跨平台通用模板矩阵类，需编译器支持C++11.

接口详细文档见“doc/html/index.html”。

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
