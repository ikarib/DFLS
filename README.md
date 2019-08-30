# C version of NEWUOA_H

This provides a C implementation of DFLS algorithm for minimizing a function of many variables.
The method is *derivatives free* (only the function values are needed).

These algorithms are modifications and based on the Software Newuoa, authored by M. J. D. Powell,
to minimize sum of squares by taking advantage of the problem structure.
   min_{x \in R^n }  F(x) := Sum_{i=1}^{mv}  v_err_i(x)^2
where v_err(x) : R^n \to R^{mv} is a vector function.

This subroutine seeks the least value of sum of the squres of the components of v_err(x)
by combing trust region method and Levenberg-Marquardt method.

References:

1.  M.J.D. Powell, "The NEWUOA software for unconstrained minimization
    without derivatives", in Large-Scale Nonlinear Optimization, editors
    G. Di Pillo and M. Roma, Springer (2006), pages 255-297.
2.  H. Zhang, A. R. Conn, and K. Scheinberg, A Derivative-Free Algorithm for Least-Squares
    Minimization, SIAM Journal on Optimization, 2010, Vol. 20, No. 6, pages 3555-3576.

The present code is based on the original FORTRAN version written by Hongchao Zhang,
https://www.math.lsu.edu/~hozhang/
and has been converted to C by Iskander Karibzhanov.

In addition to being usable from C code, this version of NEWUOA_H has a few improvements over the
FORTRAN version:

* any objective function can be used (the function is an argument of the
  method) so different objective functions can co-exist in the same code
  (*e.g.* hierarchical optimization with NEWUOA_H is possible);

* a return code indicates the success of the method or the reason of the
  failure;

* no warnings about variables being uninitialized.

Of course, this version produces the same output for the tests and it has a
more C-like interface compared to a version produced by `f2c`.

## Usage

The code consists in two files [`newuoa_h.c`](./newuoa_h.c) and [`newuoa_h.h`](./newuoa_h.h)
and is documented in the header [`newuoa_h.h`](./newuoa_h.h).

To build the library, edit the `Makefile` to suit your preferences and
do:
```
make
```
Then copy `libnewuoa_h.a` and `newuoa_h.h` to appropriate directories.

A number of macros can defined for the compiler to adapt the type of variables
used in the code (see [`newuoa_h.h`](./newuoa_h.h) for details).

You may check the code (this requires the original FORTRAN code available
from Hongchao Zhang on request at hozhang@math.lsu.edu):
```
make test
```


## Legal notice

The present code is based on the original FORTRAN version written by Hongchao Zhang
who kindly provides his code on demand (at hozhang@math.lsu.edu) and has
been converted to C by Iskander Karibzhanov.

Copyright (c) 2009, Mike Powell (FORTRAN version).
Copyright (c) 2010, Hongchao Zhang (FORTRAN version).
Copyright (c) 2019, Iskander Karibzhanov (C version).

Read the accompanying [`LICENSE`](../LICENSE) file for details.

