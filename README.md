# GenericMatrix & GenericMatrix2

C++ cross platform generic matrix template class, C++11 Required Support.

Two classes:
|Class||Brief|
|---|---|
|GenericMatrix|A generic template matrix class, the matrix elements are managed by a one-dimension pointer.
|GenericMatrix2|A generic template matrix class, the matrix elements are managed by a two-dimension pointer.

More interface information see "doc/html/index.html".

C++跨平台通用模板矩阵类，需编译器支持C++11.

接口详细文档见“doc/html/index.html”。

## Quick Example

linux or windows:

```cpp
    GenericMatrix<> matrix(2, 4);
    matrix.fill(1);

    GenericMatrix2<> matrix2(4, 4, 0);
    matrix2.setToIdentity();

    std::cout << matrix << std::endl;
    std::cout << matrix2 << std::endl;
```

(Compiled by gcc7.4.0) The Result is:

```cpp
GenericMatrix<2, 4, f>(
         01         01         01         01
         01         01         01         01
)
GenericMatrix2<4, 4, f>(
         01         00         00         00
         00         01         00         00
         00         00         01         00
         00         00         00         01
)
```

## License

[Apache License 2.0](LICENSE)
