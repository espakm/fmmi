## Fast Matrix Multiplication and Inversion

This project contains an implementation of matrix multiplication and inversion
based on the algorithms described in the following article:

  - https://www.lehigh.edu/~gi02/m242/08linstras.pdf


### Underlying data types

Since the algorithm is recursive and it relies on matrix partitioning, it was
important that selecting a matrix partition happens at low cost.

The project contains two alternative matrix types:

  - `smatrix`
  - `dmatrix`

The former implementation provides zero-cost abstraction for selecting a matrix
partition. This means that selecting a sub-block of a matrix has absolutely no
run-time cost. This comes at the price of reduced flexibility. The size of the
original matrix and the partitions to select have to be known at compile time.

The latter type (`dmatrix`) allows specifying the matrix size using any run-time
expression. This allows more flexibility, but it comes at a price that partition
selection has a certain run-time cost. For any sub-matrix, a new object needs to
be created that stores the memory address of the first (top-left) element of the
matrix, along with the height, width and stride of the matrix. The stride is the
distance of the matrix rows as stored in the memory. Using this information lets
us work with blocks of matrices at their original location, eliminating the need
for copying elements to a contiguous memory area.

Not surprisingly, benchmarking showed better performance using `smatrix` instead
of `dmatrix`.


### Differences to the original algorithms

The implementations differ from the algorithms that are described in the referred
article, in a couple of ways.

For matrices whose height or width is not a power of 2, the article suggests to
add a "padding" of zeros on the right and at the bottom, so that the "padded"
matrix is the smallest square matrix whose height and width is the same power of
2. It is shown that this does not alter the results of the multiplication and it
does not increase the computation time in terms of algorithmic complexity.

In a practical implementation, this padding can be achieved either by allocating
a larger space for the matrix and filling up the superseding area with zeros, or
by modifying the indexing operator to return zeros when an "outsider" element is
referred to. Both solutions have significant run-time costs.

To avoid this extra cost, I chose to partition the matrices for multiplication
instead of extending them to the square of power of 2 size, even when the size
of the matrix was not divisible by 2. In this case (when the height or width was
odd), the first column or row vector of the matrix has been split off, the
operation has been applied first on these vectors, then on the other (even sized)
partition of the matrix and then the results have been combined.

This implementation performed much faster for (2^n-1)x(2^n-1) sized matrices
than for (2^n)x(2^n) sized ones what approves that using padding would have been
impractical.

Three implementations are provided for matrix multiplication:

  - `mul()`: classic, iterative matrix multiplication algorithm
  - `mul_rec()`: recursive matrix multiplication as described in the article,
      extending the matrix to power of 2 size and filling up with zeros before
      the recursion starts
  - `mul_rec_part()`: recursive matrix multiplication as described above,
      splitting off left or top vectors from odd-sized matrices.

This method of partitioning has been applied only for the multiplication of odd
sized matrices. For matrix inversion, the matrices have been extended with the
unity matrix from the bottom right and zeros from the bottom and the right, as
suggested by the paper.

The other main difference is related to the matrix inversion algorithm. In the
paper an assumption is made that matrix A (whose inverse is to be determined) is
symmetric and positive definite. The described algorithm uses this assumption.

However, not every invertible matrix has these properties. For instance, the
following matrix is invertible but not symmetric:

```
5  6  6  8
2  2  2  8
6  6  2  8
2  3  6  7
```

So that we can compute the inverse of any invertible matrix, the original
blockwise inversion algorithm has been implemented, which the algorithm of the
article is based on.

In contrast with the algorithm of the article, the implemented algorithm uses
six matrix multiplications and two matrix inversions. (The algorithm of the
article makes four matrix multiplications and two matrix inversions.)

The implemented algorithm uses the recursive matrix multiplication algorithm.

For comparison, a classic, iterative algorithm has also been implemented, using
Gauss-Jordan elimination.

The two implementations provided for matrix inversion are:

  - `inv()`: classic, iterative matrix inversion with Gauss-Jordan elimination.
  - `inv_rec()`: recursive inversion as describe in the article with the
      modification discussed above to make it work for not positive-definite
      and symmetric matrices.


### Testing and benchmarking

Requirements for building and testing the project is a C++17 compatible compiler
and CMake 3.13.0 or newer.

To build the project:

```
mkdir build
cd build
cmake ..
cmake --build .
```

To run the tests and benchmarks (still from the `build` subdirectory):

```
./test/test
```

By default, tests and benchmarks are executed for matrix sizes up to 16, but
this can be changed by re-configuring the project and changing the following
variables:

  - `DMATRIX_REC_TEST_MAX_MATRIX_SIZE`
  - `SMATRIX_REC_TEST_MAX_MATRIX_SIZE`

Valid values are powers of 2 from 2 to 4096. The project can be re-configured
by the `ccmake .` or `cmake-gui .` command in the build directory. (Latter on
Windows.)

Note that increasing the max matrix size for the `smatrix` tests increases the
time and the memory consumption of the compilation significantly. Depending on
compiler and the available memory, the compilation might be aborted due to the
lack of resources.

By default, all the tests and benchmarks are executed. The following tags can be
used to filter out the required tests only:

  - `[dmatrix]`
  - `[smatrix]`
  - `[mul]`
  - `[mul_rec]`
  - `[inv]`
  - `[inv_rec]`
  - `[equals]`
  - `[benchmark]`

For instance, the benchmarks of the recursive inversion for static sized
matrices can be performed by the following command:

```
./test/test "[smatrix][inv_rec][benchmark]"
```

## Evaluation

Performance has been measured in terms of speed for six matrix element types
(`int16_t`, `int32_t`, `int64_t`, `float`, `double` and `long double`) and 16
matrix sizes (1x1, 2x2, 3x3, 4x4, 7x7, 8x8, 15x15, 16x16, 31x31, 32x32, 63x63,
64x64, 127x127, 128x128, 255x255 and 256x256) and for both matrix type
implementations (`dmatrix` and `smatrix`).

First, the iterative and the recursive *multiplication* algorithms have been
executed alternatingly, for each element type and size combinations with the
`dmatrix` implementation. Then the same has been repeated with the interative
and recursive *inversion* algorithms, also with `dmatrix` arguments. After the
process has been completed, it has been repeated with the `smatrix` functions.

The results obtained with AppleClang 12.0.0 and GCC 11.2.0 can be found in the
following file:

  - `benchmark.ods`

The cells of the spreadsheet have been highlighted in green where a recursive
algorithm executed faster.

The first observation is that the smatrix algorithms were almost always much
faster than their dmatrix equivalents. For large matrices and for the iterative
algorithms, the difference was less significant between the smatrix and dmatrix
versions, but for the recursive algorithms, the smatrix variants often performed
over 80% faster than the dmatrix ones.

The most interesting question, however, is whether there are cases when the
recursive algorithms can outperform the iterative ones. It can be observed that
the advantage of the recursive algorithms came out in particular situations and
mainly with smatrix arguments:

  - multiplication of 2x2 matrices (both smatrix and dmatrix)
  - multiplication of `int64_t` smatrix objects from size 15x15 up
  - multiplication of any smatrix objects from size 31x31 up
  - multiplication of any 127x127 and 255x255 dmatrix objects without padding
  - inversion of most dmatrix and smatrix objects up to size 4x4
  - inversion of most smatrix objects up to size 16x16
  - inversion of `int64_t` smatrix objects from size 31x31 up.

The rational explanation of the observations is that the recursive algorithms
performed better when the cost of recursion could be eliminated or outweighed by
the benefit of reduced algorithmic complexity.
