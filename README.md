# DFLS algorithm for minimizing a function of many variables.

The method is "derivatives free" (only the function values are needed).

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
