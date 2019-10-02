/*
 * newuoa_h.h -
 *
 *  Definitions for DFLS algorithm for minimizing a function of many variables.
 *  The method is "derivatives free" (only the function values are needed).
 *  
 *  These algorithms are modifications and based on the Software Newuoa, authored by M. J. D. Powell,
 *  to minimize sum of squares by taking advantage of the problem structure.
 *   min_{x \in R^n }  F(x) := Sum_{i=1}^{mv}  v_err_i(x)^2
 *  where v_err(x) : R^n \to R^{mv} is a vector function.
 *  This subroutine seeks the least value of sum of the squares of the components of v_err(x)
 *  by combing trust region method and Levenberg-Marquardt method 
 *
 *  References:
 *
 *  1.  M.J.D. Powell, "The NEWUOA software for unconstrained minimization
 *      without derivatives", in Large-Scale Nonlinear Optimization, editors
 *      G. Di Pillo and M. Roma, Springer (2006), pages 255-297.
 *  2.  H. Zhang, A. R. Conn, and K. Scheinberg, A Derivative-Free Algorithm for Least-Squares
 *      Minimization, SIAM Journal on Optimization, 2010, Vol. 20, No. 6, pages 3555-3576.
 *      
 * The present code is based on the original FORTRAN version written by Hongchao Zhang,
 * https://www.math.lsu.edu/~hozhang/
 * and has been converted to C by Iskander Karibzhanov.
 *
 */

#ifndef _NEWUOA_H
#define  _NEWUOA_H 1

#ifndef LOGICAL
# define LOGICAL int
#endif

#ifndef INTEGER
# define INTEGER long
#endif

#undef REAL
#ifdef SINGLE_PRECISION
# define REAL float
#else
# define REAL double
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Prototype of the objective function assumed by the NEWUOA_H routine:
           F(X) = Sum_{i=1}^{MV} V_ERR_i(X)^2
   N is the number of variables X(1),X(2),...,X(N),
   MV is the number of computed values of the vector function V_ERR.
   DATA is anything else needed by the objective function.
   Here: N, MV, X \in R^N and DATA are input, V_ERR \in R^{MV} are output. */
typedef void newuoa_dfovec(const INTEGER n, const INTEGER mv, const REAL* x,
	REAL* v_err, const void* data);

/* This subroutine seeks the least value of a function of many variables, by
   a trust region method that forms quadratic models by interpolation.  There
   can be some freedom in the interpolation conditions, which is taken up by
   minimizing the Frobenius norm of the change to the second derivative of
   the quadratic model, beginning with a zero matrix. The arguments of the
   subroutine are as follows.

   N must be set to the number of variables and must be at least two.
   MV must be set to the lengh of the vector function V_ERR(X):  R^N \to R^{MV}.
   The maximum number variables in this codes: NMAX =  100
   The maximum lengh of the vector function V_ERR(X): MMAX = 400
   If N > 100 or M > 400, the parameter NMAX and MMAX need to be creased in
   subroutines NEWUOB_H and TRSAPP_H

   NPT is the number of interpolation conditions. Its value must be in the
   interval [N+2,2N+1].  Recommended: NPT = 2*N+1

   Function DFOVEC(N, MV, X, V_ERR, DATA) must be provided by the user to compute
   the values of the vector function V_ERR(X) : R^N to R^{MV} at the variables X(1),X(2),...,X(N). 
   DATA is anything else needed by the objective function (unused by NEWUOA_H itself).

   Initial values of the variables must be set in X(1),X(2),...,X(N). They
   will be changed to the values that give the least calculated F(X) = Sum_{i=1}^{MV} V_ERR_i(X)^2.

   RHOBEG and RHOEND must be set to the initial and final values of a trust
   region radius, so both must be positive with RHOEND<=RHOBEG. Typically
   RHOBEG should be about one tenth of the greatest expected change to a
   variable, and RHOEND should indicate the accuracy that is required in
   the final values of the variables. Default: RHOBEG = 1.0, RHOEND = 10^{-8}

   The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
   amount of printing. Specifically, there is no output if IPRINT=0 and
   there is output only at the return if IPRINT=1. Otherwise, each new
   value of RHO is printed, with the best vector of variables so far and
   the corresponding value of the objective function. Further, each new
   value of F with its variables are output if IPRINT=3.

   MAXFUN must be set to an upper bound on the number of calls of OBJFUN.
   Default:  MAXFUN= 100(n+1), i.e 100 (simplex) gradients for reasonable accuracy.
             MAXFUN= infinity, to let the algorithm explore the lowest function value
                     as much as it could.

   The array W will be used for working space. Its length must be at least
   (NPT+11)*(NPT+N)+N*(5*N+11)/2+MV*(NPT+N*(N+7)/2+7)

   The returned value should be NEWUOA_SUCCESS, but a different value can be
   returned upon error (see `newuoa_reason` for an explanatory message). */
extern int newuoa_h(const INTEGER n, const INTEGER npt,
		newuoa_dfovec* dfovec, const void* data,
		REAL* x, const REAL rhobeg, const REAL rhoend,
		const INTEGER iprint, const INTEGER maxfun,
		REAL* w, const INTEGER mv);

/* Possible values returned by NEWUOA. */
#define NEWUOA_INITIAL_ITERATE       (2) /* only used internaly */
#define NEWUOA_ITERATE               (1) /* caller is requested to evaluate
                                            the objective function and call
                                            newuoa_iterate */
#define NEWUOA_SUCCESS               (0) /* algorithm converged */
#define NEWUOA_BAD_NPT              (-1) /* NPT is not in the required
                                            interval */
#define NEWUOA_ROUNDING_ERRORS      (-2) /* too much cancellation in a
                                            denominator */
#define NEWUOA_TOO_MANY_EVALUATIONS (-3) /* maximum number of function
                                            evaluations exceeded */
#define NEWUOA_STEP_FAILED          (-4) /* trust region step has failed to
                                            reduce quadratic approximation */
#define NEWUOA_BAD_ADDRESS          (-5) /* illegal NULL address */
#define NEWUOA_CORRUPTED            (-6) /* corrupted or misused workspace */

/* Get a textual explanation of the status returned by NEWUOA. */
extern const char* newuoa_reason(int status);

#ifdef __cplusplus
}
#endif

#endif /* _NEWUOA_H */

/*---------------------------------------------------------------------------*/
