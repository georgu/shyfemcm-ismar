
!--------------------------------------------------------------------------
!
!                   S P A R S K I T   V E R S I O N  2.
!
! Welcome  to SPARSKIT  VERSION  2.  SPARSKIT is  a  package of  FORTRAN
! subroutines  for working  with  sparse matrices.  It includes  general
! sparse  matrix  manipulation  routines  as  well as  a  few  iterative
! solvers, see detailed description of contents below.
!
!    Copyright (C) 2005  the Regents of the University of Minnesota
!
! SPARSKIT is  free software; you  can redistribute it and/or  modify it
! under the terms of the  GNU Lesser General Public License as published
! by the  Free Software Foundation [version  2.1 of the  License, or any
! later version.]
!
! A copy of  the licencing agreement is attached in  the file LGPL.  For
! additional information  contact the Free Software  Foundation Inc., 59
! Temple Place - Suite 330, Boston, MA 02111, USA or visit the web-site
!
!  http://www.gnu.org/copyleft/lesser.html
!
! DISCLAIMER
! ----------
!
! SPARSKIT  is distributed  in  the hope  that  it will  be useful,  but
! WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
! MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
! Lesser General Public License for more details.
!
! For more information contact saad@cs.umn.edu
! or see https://www-users.cs.umn.edu/~saad/software/SPARSKIT/
!
!    This file is part of SHYFEM. (m)
!
!    The original file is called ITSOL/iters.f
!
!--------------------------------------------------------------------------

! revision log :
!
! 29.10.2023    ggu     imported into SHYFEM

!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!         Basic Iterative Solvers with Reverse Communication           c
!----------------------------------------------------------------------c
!     This file currently has several basic iterative linear system    c
!     solvers. They are:                                               c
!     CG       -- Conjugate Gradient Method                            c
!     CGNR     -- Conjugate Gradient Method (Normal Residual equation) c
!     BCG      -- Bi-Conjugate Gradient Method                         c
!     DBCG     -- BCG with partial pivoting                            c
!     BCGSTAB  -- BCG stabilized                                       c
!     TFQMR    -- Transpose-Free Quasi-Minimum Residual method         c
!     FOM      -- Full Orthogonalization Method                        c
!     GMRES    -- Generalized Minimum RESidual method                  c
!     FGMRES   -- Flexible version of Generalized Minimum              c
!                 RESidual method                                      c
!     DQGMRES  -- Direct versions of Quasi Generalize Minimum          c
!                 Residual method                                      c
!----------------------------------------------------------------------c
!     They all have the following calling sequence:
!      subroutine solver(n, rhs, sol, ipar, fpar, w)
!      integer n, ipar(16)
!      real*8 rhs(n), sol(n), fpar(16), w(*)
!     Where
!     (1) 'n' is the size of the linear system,
!     (2) 'rhs' is the right-hand side of the linear system,
!     (3) 'sol' is the solution to the linear system,
!     (4) 'ipar' is an integer parameter array for the reverse
!     communication protocol,
!     (5) 'fpar' is an floating-point parameter array storing
!     information to and from the iterative solvers.
!     (6) 'w' is the work space (size is specified in ipar)
!
!     They are preconditioned iterative solvers with reverse
!     communication. The preconditioners can be applied from either
!     from left or right or both (specified by ipar(2), see below).
!
!     Author: Kesheng John Wu (kewu@mail.cs.umn.edu) 1993
!
!     NOTES:
!
!     (1) Work space required by each of the iterative solver
!     routines is as follows:
!       CG      == 5 * n
!       CGNR    == 5 * n
!       BCG     == 7 * n
!       DBCG    == 11 * n
!       BCGSTAB == 8 * n
!       TFQMR   == 11 * n
!       FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!       GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!       FGMRES  == 2*n*(m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
!                  default m=15)
!       DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
!
!     (2) ALL iterative solvers require a user-supplied DOT-product
!     routine named DISTDOT. The prototype of DISTDOT is
!
!     real*8 function distdot(n,x,ix,y,iy)
!     integer n, ix, iy
!     real*8 x(1+(n-1)*ix), y(1+(n-1)*iy)
!
!     This interface of DISTDOT is exactly the same as that of
!     DDOT (or SDOT if real == real*8) from BLAS-1. It should have
!     same functionality as DDOT on a single processor machine. On a
!     parallel/distributed environment, each processor can perform
!     DDOT on the data it has, then perform a summation on all the
!     partial results.
!
!     (3) To use this set of routines under SPMD/MIMD program paradigm,
!     several things are to be noted: (a) 'n' should be the number of
!     vector elements of 'rhs' that is present on the local processor.
!     (b) if RHS(i) is on processor j, it is expected that SOL(i)
!     will be on the same processor, i.e. the vectors are distributed
!     to each processor in the same way. (c) the preconditioning and
!     stopping criteria specifications have to be the same on all
!     processor involved, ipar and fpar have to be the same on each
!     processor. (d) DISTDOT should be replaced by a distributed
!     dot-product function.
!
!     ..................................................................
!     Reverse Communication Protocols
!
!     When a reverse-communication routine returns, it could be either
!     that the routine has terminated or it simply requires the caller
!     to perform one matrix-vector multiplication. The possible matrices
!     that involve in the matrix-vector multiplications are:
!     A       (the matrix of the linear system),
!     A^T     (A transposed),
!     Ml^{-1} (inverse of the left preconditioner),
!     Ml^{-T} (inverse of the left preconditioner transposed),
!     Mr^{-1} (inverse of the right preconditioner),
!     Mr^{-T} (inverse of the right preconditioner transposed).
!     For all the matrix vector multiplication, v = A u. The input and
!     output vectors are supposed to be part of the work space 'w', and
!     the starting positions of them are stored in ipar(8:9), see below.
!
!     The array 'ipar' is used to store the information about the solver.
!     Here is the list of what each element represents:
!
!     ipar(1) -- status of the call/return.
!     A call to the solver with ipar(1) == 0 will initialize the
!     iterative solver. On return from the iterative solver, ipar(1)
!     carries the status flag which indicates the condition of the
!     return. The status information is divided into two categories,
!     (1) a positive value indicates the solver requires a matrix-vector
!     multiplication,
!     (2) a non-positive value indicates termination of the solver.
!     Here is the current definition:
!       1 == request a matvec with A,
!       2 == request a matvec with A^T,
!       3 == request a left preconditioner solve (Ml^{-1}),
!       4 == request a left preconditioner transposed solve (Ml^{-T}),
!       5 == request a right preconditioner solve (Mr^{-1}),
!       6 == request a right preconditioner transposed solve (Mr^{-T}),
!      10 == request the caller to perform stopping test,
!       0 == normal termination of the solver, satisfied the stopping
!            criteria,
!      -1 == termination because iteration number is greater than the
!            preset limit,
!      -2 == return due to insufficient work space,
!      -3 == return due to anticipated break-down / divide by zero,
!            in the case where Arnoldi procedure is used, additional
!            error code can be found in ipar(12), where ipar(12) is
!            the error code of orthogonalization procedure MGSRO:
!               -1: zero input vector
!               -2: input vector contains abnormal numbers
!               -3: input vector is a linear combination of others
!               -4: trianguler system in GMRES/FOM/etc. has nul rank
!      -4 == the values of fpar(1) and fpar(2) are both <= 0, the valid
!            ranges are 0 <= fpar(1) < 1, 0 <= fpar(2), and they can
!            not be zero at the same time
!      -9 == while trying to detect a break-down, an abnormal number is
!            detected.
!     -10 == return due to some non-numerical reasons, e.g. invalid
!            floating-point numbers etc.
!
!     ipar(2) -- status of the preconditioning:
!       0 == no preconditioning
!       1 == left preconditioning only
!       2 == right preconditioning only
!       3 == both left and right preconditioning
!
!     ipar(3) -- stopping criteria (details of this will be
!     discussed later).
!
!     ipar(4) -- number of elements in the array 'w'. if this is less
!     than the desired size, it will be over-written with the minimum
!     requirement. In which case the status flag ipar(1) = -2.
!
!     ipar(5) -- size of the Krylov subspace (used by GMRES and its
!     variants), e.g. GMRES(ipar(5)), FGMRES(ipar(5)),
!     DQGMRES(ipar(5)).
!
!     ipar(6) -- maximum number of matrix-vector multiplies, if not a
!     positive number the iterative solver will run till convergence
!     test is satisfied.
!
!     ipar(7) -- current number of matrix-vector multiplies. It is
!     incremented after each matrix-vector multiplication. If there
!     is preconditioning, the counter is incremented after the
!     preconditioning associated with each matrix-vector multiplication.
!
!     ipar(8) -- pointer to the input vector to the requested matrix-
!     vector multiplication.
!
!     ipar(9) -- pointer to the output vector of the requested matrix-
!     vector multiplication.
!
!     To perform v = A * u, it is assumed that u is w(ipar(8):ipar(8)+n-1)
!     and v is stored as w(ipar(9):ipar(9)+n-1).
!
!     ipar(10) -- the return address (used to determine where to go to
!     inside the iterative solvers after the caller has performed the
!     requested services).
!
!     ipar(11) -- the result of the external convergence test
!     On final return from the iterative solvers, this value
!     will be reflected by ipar(1) = 0 (details discussed later)
!
!     ipar(12) -- error code of MGSRO, it is
!                  1 if the input vector to MGSRO is linear combination
!                    of others,
!                  0 if MGSRO was successful,
!                 -1 if the input vector to MGSRO is zero,
!                 -2 if the input vector contains invalid number.
!
!     ipar(13) -- number of initializations. During each initilization
!                 residual norm is computed directly from M_l(b - A x).
!
!     ipar(14) to ipar(16) are NOT defined, they are NOT USED by
!     any iterative solver at this time.
!
!     Information about the error and tolerance are stored in the array
!     FPAR. So are some internal variables that need to be saved from
!     one iteration to the next one. Since the internal variables are
!     not the same for each routine, we only define the common ones.
!
!     The first two are input parameters:
!     fpar(1) -- the relative tolerance,
!     fpar(2) -- the absolute tolerance (details discussed later),
!
!     When the iterative solver terminates,
!     fpar(3) -- initial residual/error norm,
!     fpar(4) -- target residual/error norm,
!     fpar(5) -- current residual norm (if available),
!     fpar(6) -- current residual/error norm,
!     fpar(7) -- convergence rate,
!
!     fpar(8:10) are used by some of the iterative solvers to save some
!     internal information.
!
!     fpar(11) -- number of floating-point operations. The iterative
!     solvers will add the number of FLOPS they used to this variable,
!     but they do NOT initialize it, nor add the number of FLOPS due to
!     matrix-vector multiplications (since matvec is outside of the
!     iterative solvers). To insure the correct FLOPS count, the
!     caller should set fpar(11) = 0 before invoking the iterative
!     solvers and account for the number of FLOPS from matrix-vector
!     multiplications and preconditioners.
!
!     fpar(12:16) are not used in current implementation.
!
!     Whether the content of fpar(3), fpar(4) and fpar(6) are residual
!     norms or error norms depends on ipar(3). If the requested
!     convergence test is based on the residual norm, they will be
!     residual norms. If the caller want to test convergence based the
!     error norms (estimated by the norm of the modifications applied
!     to the approximate solution), they will be error norms.
!     Convergence rate is defined by (Fortran 77 statement)
!     fpar(7) = log10(fpar(3) / fpar(6)) / (ipar(7)-ipar(13))
!     If fpar(7) = 0.5, it means that approximately every 2 (= 1/0.5)
!     steps the residual/error norm decrease by a factor of 10.
!
!     ..................................................................
!     Stopping criteria,
!
!     An iterative solver may be terminated due to (1) satisfying
!     convergence test; (2) exceeding iteration limit; (3) insufficient
!     work space; (4) break-down. Checking of the work space is
!     only done in the initialization stage, i.e. when it is called with
!     ipar(1) == 0. A complete convergence test is done after each
!     update of the solutions. Other conditions are monitored
!     continuously.
!
!     With regard to the number of iteration, when ipar(6) is positive,
!     the current iteration number will be checked against it. If
!     current iteration number is greater the ipar(6) than the solver
!     will return with status -1. If ipar(6) is not positive, the
!     iteration will continue until convergence test is satisfied.
!
!     Two things may be used in the convergence tests, one is the
!     residual 2-norm, the other one is 2-norm of the change in the
!     approximate solution. The residual and the change in approximate
!     solution are from the preconditioned system (if preconditioning
!     is applied). The DQGMRES and TFQMR use two estimates for the
!     residual norms. The estimates are not accurate, but they are
!     acceptable in most of the cases. Generally speaking, the error
!     of the TFQMR's estimate is less accurate.
!
!     The convergence test type is indicated by ipar(3). There are four
!     type convergence tests: (1) tests based on the residual norm;
!     (2) tests based on change in approximate solution; (3) caller
!     does not care, the solver choose one from above two on its own;
!     (4) caller will perform the test, the solver should simply continue.
!     Here is the complete definition:
!      -2 == || dx(i) || <= rtol * || rhs || + atol
!      -1 == || dx(i) || <= rtol * || dx(1) || + atol
!       0 == solver will choose test 1 (next)
!       1 == || residual || <= rtol * || initial residual || + atol
!       2 == || residual || <= rtol * || rhs || + atol
!     999 == caller will perform the test
!     where dx(i) denote the change in the solution at the ith update.
!     ||.|| denotes 2-norm. rtol = fpar(1) and atol = fpar(2).
!
!     If the caller is to perform the convergence test, the outcome
!     should be stored in ipar(11).
!     ipar(11) = 0 -- failed the convergence test, iterative solver
!     should continue
!     ipar(11) = 1 -- satisfied convergence test, iterative solver
!     should perform the clean up job and stop.
!
!     Upon return with ipar(1) = 10,
!     ipar(8)  points to the starting position of the change in
!              solution Sx, where the actual solution of the step is
!              x_j = x_0 + M_r^{-1} Sx.
!              Exception: ipar(8) < 0, Sx = 0. It is mostly used by
!              GMRES and variants to indicate (1) Sx was not necessary,
!              (2) intermediate result of Sx is not computed.
!     ipar(9)  points to the starting position of a work vector that
!              can be used by the caller.
!
!     NOTE: the caller should allow the iterative solver to perform
!     clean up job after the external convergence test is satisfied,
!     since some of the iterative solvers do not directly
!     update the 'sol' array. A typical clean-up stage includes
!     performing the final update of the approximate solution and
!     computing the convergence information (e.g. values of fpar(3:7)).
!
!     NOTE: fpar(4) and fpar(6) are not set by the accelerators (the
!     routines implemented here) if ipar(3) = 999.
!
!     ..................................................................
!     Usage:
!
!     To start solving a linear system, the user needs to specify
!     first 6 elements of the ipar, and first 2 elements of fpar.
!     The user may optionally set fpar(11) = 0 if one wants to count
!     the number of floating-point operations. (Note: the iterative
!     solvers will only add the floating-point operations inside
!     themselves, the caller will have to add the FLOPS from the
!     matrix-vector multiplication routines and the preconditioning
!     routines in order to account for all the arithmetic operations.)
!
!     Here is an example:
!     ipar(1) = 0	! always 0 to start an iterative solver
!     ipar(2) = 2	! right preconditioning
!     ipar(3) = 1	! use convergence test scheme 1
!     ipar(4) = 10000	! the 'w' has 10,000 elements
!     ipar(5) = 10	! use *GMRES(10) (e.g. FGMRES(10))
!     ipar(6) = 100	! use at most 100 matvec's
!     fpar(1) = 1.0E-6	! relative tolerance 1.0E-6
!     fpar(2) = 1.0E-10 ! absolute tolerance 1.0E-10
!     fpar(11) = 0.0	! clearing the FLOPS counter
!
!     After the above specifications, one can start to call an iterative
!     solver, say BCG. Here is a piece of pseudo-code showing how it can
!     be done,
!
! 10   call bcg(n,rhs,sol,ipar,fpar,w)
!      if (ipar(1).eq.1) then
!         call amux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!         goto 10
!      else if (ipar(1).eq.2) then
!         call atmux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!         goto 10
!      else if (ipar(1).eq.3) then
!         left preconditioner solver
!         goto 10
!      else if (ipar(1).eq.4) then
!         left preconditioner transposed solve
!         goto 10
!      else if (ipar(1).eq.5) then
!         right preconditioner solve
!         goto 10
!      else if (ipar(1).eq.6) then
!         right preconditioner transposed solve
!         goto 10
!      else if (ipar(1).eq.10) then
!         call my own stopping test routine
!         goto 10
!      else if (ipar(1).gt.0) then
!         ipar(1) is an unspecified code
!      else
!         the iterative solver terminated with code = ipar(1)
!      endif
!
!     This segment of pseudo-code assumes the matrix is in CSR format,
!     AMUX and ATMUX are two routines from the SPARSKIT MATVEC module.
!     They perform matrix-vector multiplications for CSR matrices,
!     where w(ipar(8)) is the first element of the input vectors to the
!     two routines, and w(ipar(9)) is the first element of the output
!     vectors from them. For simplicity, we did not show the name of
!     the routine that performs the preconditioning operations or the
!     convergence tests.
!-----------------------------------------------------------------------
      subroutine cg(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
!-----------------------------------------------------------------------
!     This is a implementation of the Conjugate Gradient (CG) method
!     for solving linear system.
!
!     NOTE: This is not the PCG algorithm. It is a regular CG algorithm.
!     To be consistent with the other solvers, the preconditioners are
!     applied by performing Ml^{-1} A Mr^{-1} P in place of A P in the
!     CG algorithm.  PCG uses the preconditioner differently.
!
!     fpar(7) is used here internally to store <r, r>.
!     w(:,1) -- residual vector
!     w(:,2) -- P, the conjugate direction
!     w(:,3) -- A P, matrix multiply the conjugate direction
!     w(:,4) -- temporary storage for results of preconditioning
!     w(:,5) -- change in the solution (sol) is stored here until
!               termination of this solver
!-----------------------------------------------------------------------
!     external functions used
!
      real*8 distdot
      logical stopbis, brkdn
      external distdot, stopbis, brkdn, bisinit
!
!     local variables
!
      integer i
      real*8 alpha
      logical lp,rp
      save
!
!     check the status of the call
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80), ipar(10)
!
!     initialization
!
      call bisinit(ipar,fpar,5*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
!
!     request for matrix vector multiplication A*x in the initialization
!
      ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = ipar(8) + n
      ipar(10) = 1
      do i = 1, n
         w(i,2) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = 1
      do i = 1, n
         w(i,2) = rhs(i) - w(i,3)
      enddo
      fpar(11) = fpar(11) + n
!
!     if left preconditioned
!
      if (lp) then
         ipar(1) = 3
         ipar(9) = 1
         ipar(10) = 2
         return
      endif
!
 20   if (lp) then
         do i = 1, n
            w(i,2) = w(i,1)
         enddo
      else
         do i = 1, n
            w(i,1) = w(i,2)
         enddo
      endif
!
      fpar(7) = distdot(n,w,1,w,1)
      fpar(11) = fpar(11) + 2 * n
      fpar(3) = sqrt(fpar(7))
      fpar(5) = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
!
!     before iteration can continue, we need to compute A * p, which
!     includes the preconditioning operations
!
 30   if (rp) then
         ipar(1) = 5
         ipar(8) = n + 1
         if (lp) then
            ipar(9) = ipar(8) + n
         else
            ipar(9) = 3*n + 1
         endif
         ipar(10) = 3
         return
      endif
!
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n + 1
      endif
      if (lp) then
         ipar(9) = 3*n+1
      else
         ipar(9) = n+n+1
      endif
      ipar(10) = 4
      return
!
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = n+n+1
         ipar(10) = 5
         return
      endif
!
!     continuing with the iterations
!
 60   ipar(7) = ipar(7) + 1
      alpha = distdot(n,w(1,2),1,w(1,3),1)
      fpar(11) = fpar(11) + 2*n
      if (brkdn(alpha,ipar)) goto 900
      alpha = fpar(7) / alpha
      do i = 1, n
         w(i,5) = w(i,5) + alpha * w(i,2)
         w(i,1) = w(i,1) - alpha * w(i,3)
      enddo
      fpar(11) = fpar(11) + 4*n
!
!     are we ready to terminate ?
!
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 4*n + 1
         ipar(9) = 3*n + 1
         ipar(10) = 6
         return
      endif
 70   if (ipar(3).eq.999) then
         if (ipar(11).eq.1) goto 900
      else if (stopbis(n,ipar,1,fpar,w,w(1,2),alpha)) then
         goto 900
      endif
!
!     continue the iterations
!
      alpha = fpar(5)*fpar(5) / fpar(7)
      fpar(7) = fpar(5)*fpar(5)
      do i = 1, n
         w(i,2) = w(i,1) + alpha * w(i,2)
      enddo
      fpar(11) = fpar(11) + 2*n
      goto 30
!
!     clean up -- necessary to accommodate the right-preconditioning
!
 900  if (rp) then
         if (ipar(1).lt.0) ipar(12) = ipar(1)
         ipar(1) = 5
         ipar(8) = 4*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 7
         return
      endif
 80   if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,4))
      else
         call tidycg(n,ipar,fpar,sol,w(1,5))
      endif
!
      return
      end
!-----end-of-cg
!-----------------------------------------------------------------------
      subroutine cgnr(n,rhs,sol,ipar,fpar,wk)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n),sol(n),fpar(16),wk(n,*)
!-----------------------------------------------------------------------
!     CGNR -- Using CG algorithm solving A x = b by solving
!     Normal Residual equation: A^T A x = A^T b
!     As long as the matrix is not singular, A^T A is symmetric
!     positive definite, therefore CG (CGNR) will converge.
!
!     Usage of the work space:
!     wk(:,1) == residual vector R
!     wk(:,2) == the conjugate direction vector P
!     wk(:,3) == a scratch vector holds A P, or A^T R
!     wk(:,4) == a scratch vector holds intermediate results of the
!                preconditioning
!     wk(:,5) == a place to hold the modification to SOL
!
!     size of the work space WK is required = 5*n
!-----------------------------------------------------------------------
!     external functions used
!
      real*8 distdot
      logical stopbis, brkdn
      external distdot, stopbis, brkdn, bisinit
!
!     local variables
!
      integer i
      real*8 alpha, zz, zzm1
      logical lp, rp
      save
!
!     check the status of the call
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110), ipar(10)
!
!     initialization
!
      call bisinit(ipar,fpar,5*n,1,lp,rp,wk)
      if (ipar(1).lt.0) return
!
!     request for matrix vector multiplication A*x in the initialization
!
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      ipar(10) = 1
      do i = 1, n
         wk(i,1) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      do i = 1, n
         wk(i,1) = rhs(i) - wk(i,2)
      enddo
      fpar(11) = fpar(11) + n
!
!     if left preconditioned, precondition the initial residual
!
      if (lp) then
         ipar(1) = 3
         ipar(10) = 2
         return
      endif
!
 20   if (lp) then
         do i = 1, n
            wk(i,1) = wk(i,2)
         enddo
      endif
!
      zz = distdot(n,wk,1,wk,1)
      fpar(11) = fpar(11) + 2 * n
      fpar(3) = sqrt(zz)
      fpar(5) = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
!
!     normal iteration begins here, first half of the iteration
!     computes the conjugate direction
!
 30   continue
!
!     request the caller to perform a A^T r --> wk(:,3)
!
      if (lp) then
         ipar(1) = 4
         ipar(8) = 1
         if (rp) then
            ipar(9) = n + n + 1
         else
            ipar(9) = 3*n + 1
         endif
         ipar(10) = 3
         return
      endif
!
 40   ipar(1) = 2
      if (lp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 1
      endif
      if (rp) then
         ipar(9) = 3*n + 1
      else
         ipar(9) = n + n + 1
      endif
      ipar(10) = 4
      return
!
 50   if (rp) then
         ipar(1) = 6
         ipar(8) = ipar(9)
         ipar(9) = n + n + 1
         ipar(10) = 5
         return
      endif
!
 60   ipar(7) = ipar(7) + 1
      zzm1 = zz
      zz = distdot(n,wk(1,3),1,wk(1,3),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(zz,ipar)) goto 900
      if (ipar(7).gt.3) then
         alpha = zz / zzm1
         do i = 1, n
            wk(i,2) = wk(i,3) + alpha * wk(i,2)
         enddo
         fpar(11) = fpar(11) + 2 * n
      else
         do i = 1, n
            wk(i,2) = wk(i,3)
         enddo
      endif
!
!     before iteration can continue, we need to compute A * p
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = n + 1
         if (lp) then
            ipar(9) = ipar(8) + n
         else
            ipar(9) = 3*n + 1
         endif
         ipar(10) = 6
         return
      endif
!
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n + 1
      endif
      if (lp) then
        ipar(9) = 3*n+1
      else
         ipar(9) = n+n+1
      endif
      ipar(10) = 7
      return
!
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = n+n+1
         ipar(10) = 8
         return
      endif
!
!     update the solution -- accumulate the changes in w(:,5)
!
 90   ipar(7) = ipar(7) + 1
      alpha = distdot(n,wk(1,3),1,wk(1,3),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha,ipar)) goto 900
      alpha = zz / alpha
      do i = 1, n
         wk(i,5) = wk(i,5) + alpha * wk(i,2)
         wk(i,1) = wk(i,1) - alpha * wk(i,3)
      enddo
      fpar(11) = fpar(11) + 4 * n
!
!     are we ready to terminate ?
!
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 4*n + 1
         ipar(9) = 3*n + 1
         ipar(10) = 9
         return
      endif
 100  if (ipar(3).eq.999) then
         if (ipar(11).eq.1) goto 900
      else if (stopbis(n,ipar,1,fpar,wk,wk(1,2),alpha)) then
         goto 900
      endif
!
!     continue the iterations
!
      goto 30
!
!     clean up -- necessary to accommodate the right-preconditioning
!
 900  if (rp) then
         if (ipar(1).lt.0) ipar(12) = ipar(1)
         ipar(1) = 5
         ipar(8) = 4*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      endif
 110  if (rp) then
         call tidycg(n,ipar,fpar,sol,wk(1,4))
      else
         call tidycg(n,ipar,fpar,sol,wk(1,5))
      endif
      return
      end
!-----end-of-cgnr
!-----------------------------------------------------------------------
      subroutine bcg(n,rhs,sol,ipar,fpar,w)
      implicit none
      integer n, ipar(16)
      real*8  fpar(16), rhs(n), sol(n), w(n,*)
!-----------------------------------------------------------------------
!     BCG: Bi Conjugate Gradient method. Programmed with reverse
!     communication, see the header for detailed specifications
!     of the protocol.
!
!     in this routine, before successful return, the fpar's are
!     fpar(3) == initial residual norm
!     fpar(4) == target residual norm
!     fpar(5) == current residual norm
!     fpar(7) == current rho (rhok = <r, s>)
!     fpar(8) == previous rho (rhokm1)
!
!     w(:,1) -- r, the residual
!     w(:,2) -- s, the dual of the 'r'
!     w(:,3) -- p, the projection direction
!     w(:,4) -- q, the dual of the 'p'
!     w(:,5) -- v, a scratch vector to store A*p, or A*q.
!     w(:,6) -- a scratch vector to store intermediate results
!     w(:,7) -- changes in the solution
!-----------------------------------------------------------------------
!     external routines used
!
      real*8 distdot
      logical stopbis,brkdn
      external distdot, stopbis, brkdn
!
      real*8 one
      parameter(one=1.0D0)
!
!     local variables
!
      integer i
      real*8 alpha
      logical rp, lp
      save
!
!     status of the program
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110), ipar(10)
!
!     initialization, initial residual
!
      call bisinit(ipar,fpar,7*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
!
!     compute initial residual, request a matvecc
!
      ipar(1) = 1
      ipar(8) = 3*n+1
      ipar(9) = ipar(8) + n
      do i = 1, n
         w(i,4) = sol(i)
      enddo
      ipar(10) = 1
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      do i = 1, n
         w(i,1) = rhs(i) - w(i,5)
      enddo
      fpar(11) = fpar(11) + n
      if (lp) then
         ipar(1) = 3
         ipar(8) = 1
         ipar(9) = n+1
         ipar(10) = 2
         return
      endif
!
 20   if (lp) then
         do i = 1, n
            w(i,1) = w(i,2)
            w(i,3) = w(i,2)
            w(i,4) = w(i,2)
         enddo
      else
         do i = 1, n
            w(i,2) = w(i,1)
            w(i,3) = w(i,1)
            w(i,4) = w(i,1)
         enddo
      endif
!
      fpar(7) = distdot(n,w,1,w,1)
      fpar(11) = fpar(11) + 2 * n
      fpar(3) = sqrt(fpar(7))
      fpar(5) = fpar(3)
      fpar(8) = one
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
      if (ipar(3).ge.0.and.fpar(5).le.fpar(4)) then
         fpar(6) = fpar(5)
         goto 900
      endif
!
!     end of initialization, begin iteration, v = A p
!
 30   if (rp) then
         ipar(1) = 5
         ipar(8) = n + n + 1
         if (lp) then
            ipar(9) = 4*n + 1
         else
            ipar(9) = 5*n + 1
         endif
         ipar(10) = 3
         return
      endif
!
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n + n + 1
      endif
      if (lp) then
         ipar(9) = 5*n + 1
      else
         ipar(9) = 4*n + 1
      endif
      ipar(10) = 4
      return
!
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 4*n + 1
         ipar(10) = 5
         return
      endif
!
 60   ipar(7) = ipar(7) + 1
      alpha = distdot(n,w(1,4),1,w(1,5),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha,ipar)) goto 900
      alpha = fpar(7) / alpha
      do i = 1, n
         w(i,7) = w(i,7) + alpha * w(i,3)
         w(i,1) = w(i,1) - alpha * w(i,5)
      enddo
      fpar(11) = fpar(11) + 4 * n
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 6*n + 1
         ipar(9) = 5*n + 1
         ipar(10) = 6
         return
      endif
 70   if (ipar(3).eq.999) then
         if (ipar(11).eq.1) goto 900
      else if (stopbis(n,ipar,1,fpar,w,w(1,3),alpha)) then
         goto 900
      endif
!
!     A^t * x
!
      if (lp) then
         ipar(1) = 4
         ipar(8) = 3*n + 1
         if (rp) then
            ipar(9) = 4*n + 1
         else
            ipar(9) = 5*n + 1
         endif
         ipar(10) = 7
         return
      endif
!
 80   ipar(1) = 2
      if (lp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 3*n + 1
      endif
      if (rp) then
         ipar(9) = 5*n + 1
      else
         ipar(9) = 4*n + 1
      endif
      ipar(10) = 8
      return
!
 90   if (rp) then
         ipar(1) = 6
         ipar(8) = ipar(9)
         ipar(9) = 4*n + 1
         ipar(10) = 9
         return
      endif
!
 100  ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,2) = w(i,2) - alpha * w(i,5)
      enddo
      fpar(8) = fpar(7)
      fpar(7) = distdot(n,w,1,w(1,2),1)
      fpar(11) = fpar(11) + 4 * n
      if (brkdn(fpar(7), ipar)) return
      alpha = fpar(7) / fpar(8)
      do i = 1, n
         w(i,3) = w(i,1) + alpha * w(i,3)
         w(i,4) = w(i,2) + alpha * w(i,4)
      enddo
      fpar(11) = fpar(11) + 4 * n
!
!     end of the iterations
!
      goto 30
!
!     some clean up job to do
!
 900  if (rp) then
         if (ipar(1).lt.0) ipar(12) = ipar(1)
         ipar(1) = 5
         ipar(8) = 6*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      endif
 110  if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,6))
      else
         call tidycg(n,ipar,fpar,sol,w(1,7))
      endif
      return
!-----end-of-bcg
      end
!-----------------------------------------------------------------------
      subroutine bcgstab(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,8)
!-----------------------------------------------------------------------
!     BCGSTAB --- Bi Conjugate Gradient stabilized (BCGSTAB)
!     This is an improved BCG routine. (1) no matrix transpose is
!     involved. (2) the convergence is smoother.
!
!
!     Algorithm:
!     Initialization - r = b - A x, r0 = r, p = r, rho = (r0, r),
!     Iterate -
!     (1) v = A p
!     (2) alpha = rho / (r0, v)
!     (3) s = r - alpha v
!     (4) t = A s
!     (5) omega = (t, s) / (t, t)
!     (6) x = x + alpha * p + omega * s
!     (7) r = s - omega * t
!     convergence test goes here
!     (8) beta = rho, rho = (r0, r), beta = rho * alpha / (beta * omega)
!         p = r + beta * (p - omega * v)
!
!     in this routine, before successful return, the fpar''s are
!     fpar(3) == initial (preconditionied-)residual norm
!     fpar(4) == target (preconditionied-)residual norm
!     fpar(5) == current (preconditionied-)residual norm
!     fpar(6) == current residual norm or error
!     fpar(7) == current rho (rhok = <r, r0>)
!     fpar(8) == alpha
!     fpar(9) == omega
!
!     Usage of the work space W
!     w(:, 1) = r0, the initial residual vector
!     w(:, 2) = r, current residual vector
!     w(:, 3) = s
!     w(:, 4) = t
!     w(:, 5) = v
!     w(:, 6) = p
!     w(:, 7) = tmp, used in preconditioning, etc.
!     w(:, 8) = delta x, the correction to the answer is accumulated
!               here, so that the right-preconditioning may be applied
!               at the end
!-----------------------------------------------------------------------
!     external routines used
!
      real*8 distdot
      logical stopbis, brkdn
      external distdot, stopbis, brkdn
!
      real*8 one
      parameter(one=1.0D0)
!
!     local variables
!
      integer i
      real*8 alpha,beta,rho,omega
      logical lp, rp
      save lp, rp
!
!     where to go
!
      if (ipar(1).gt.0) then
         goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110) ipar(10)
      else if (ipar(1).lt.0) then
         goto 900
      endif
!
!     call the initialization routine
!
      call bisinit(ipar,fpar,8*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
!
!     perform a matvec to compute the initial residual
!
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      do i = 1, n
         w(i,1) = sol(i)
      enddo
      ipar(10) = 1
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      do i = 1, n
         w(i,1) = rhs(i) - w(i,2)
      enddo
      fpar(11) = fpar(11) + n
      if (lp) then
         ipar(1) = 3
         ipar(10) = 2
         return
      endif
!
 20   if (lp) then
         do i = 1, n
            w(i,1) = w(i,2)
            w(i,6) = w(i,2)
         enddo
      else
         do i = 1, n
            w(i,2) = w(i,1)
            w(i,6) = w(i,1)
         enddo
      endif
!
      fpar(7) = distdot(n,w,1,w,1)
      fpar(11) = fpar(11) + 2 * n
      fpar(5) = sqrt(fpar(7))
      fpar(3) = fpar(5)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2 * n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * fpar(3) + fpar(2)
      endif
      if (ipar(3).ge.0) fpar(6) = fpar(5)
      if (ipar(3).ge.0 .and. fpar(5).le.fpar(4) .and. &
     &     ipar(3).ne.999) then
         goto 900
      endif
!
!     beginning of the iterations
!
!     Step (1), v = A p
 30   if (rp) then
         ipar(1) = 5
         ipar(8) = 5*n+1
         if (lp) then
            ipar(9) = 4*n + 1
         else
            ipar(9) = 6*n + 1
         endif
         ipar(10) = 3
         return
      endif
!
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 5*n+1
      endif
      if (lp) then
         ipar(9) = 6*n + 1
      else
         ipar(9) = 4*n + 1
      endif
      ipar(10) = 4
      return
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 4*n + 1
         ipar(10) = 5
         return
      endif
!
 60   ipar(7) = ipar(7) + 1
!
!     step (2)
      alpha = distdot(n,w(1,1),1,w(1,5),1)
      fpar(11) = fpar(11) + 2 * n
      if (brkdn(alpha, ipar)) goto 900
      alpha = fpar(7) / alpha
      fpar(8) = alpha
!
!     step (3)
      do i = 1, n
         w(i,3) = w(i,2) - alpha * w(i,5)
      enddo
      fpar(11) = fpar(11) + 2 * n
!
!     Step (4): the second matvec -- t = A s
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = n+n+1
         if (lp) then
            ipar(9) = ipar(8)+n
         else
            ipar(9) = 6*n + 1
         endif
         ipar(10) = 6
         return
      endif
!
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = n+n+1
      endif
      if (lp) then
         ipar(9) = 6*n + 1
      else
         ipar(9) = 3*n + 1
      endif
      ipar(10) = 7
      return
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 3*n + 1
         ipar(10) = 8
         return
      endif
 90   ipar(7) = ipar(7) + 1
!
!     step (5)
      omega = distdot(n,w(1,4),1,w(1,4),1)
      fpar(11) = fpar(11) + n + n
      if (brkdn(omega,ipar)) goto 900
      omega = distdot(n,w(1,4),1,w(1,3),1) / omega
      fpar(11) = fpar(11) + n + n
      if (brkdn(omega,ipar)) goto 900
      fpar(9) = omega
      alpha = fpar(8)
!
!     step (6) and (7)
      do i = 1, n
         w(i,7) = alpha * w(i,6) + omega * w(i,3)
         w(i,8) = w(i,8) + w(i,7)
         w(i,2) = w(i,3) - omega * w(i,4)
      enddo
      fpar(11) = fpar(11) + 6 * n + 1
!
!     convergence test
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 7*n + 1
         ipar(9) = 6*n + 1
         ipar(10) = 9
         return
      endif
      if (stopbis(n,ipar,2,fpar,w(1,2),w(1,7),one))  goto 900
 100  if (ipar(3).eq.999.and.ipar(11).eq.1) goto 900
!
!     step (8): computing new p and rho
      rho = fpar(7)
      fpar(7) = distdot(n,w(1,2),1,w(1,1),1)
      omega = fpar(9)
      beta = fpar(7) * fpar(8) / (fpar(9) * rho)
      do i = 1, n
         w(i,6) = w(i,2) + beta * (w(i,6) - omega * w(i,5))
      enddo
      fpar(11) = fpar(11) + 6 * n + 3
      if (brkdn(fpar(7),ipar)) goto 900
!
!     end of an iteration
!
      goto 30
!
!     some clean up job to do
!
 900  if (rp) then
         if (ipar(1).lt.0) ipar(12) = ipar(1)
         ipar(1) = 5
         ipar(8) = 7*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      endif
 110  if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,7))
      else
         call tidycg(n,ipar,fpar,sol,w(1,8))
      endif
!
      return
!-----end-of-bcgstab
      end
!-----------------------------------------------------------------------
      subroutine tfqmr(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
!-----------------------------------------------------------------------
!     TFQMR --- transpose-free Quasi-Minimum Residual method
!     This is developed from BCG based on the principle of Quasi-Minimum
!     Residual, and it is transpose-free.
!
!     It uses approximate residual norm.
!
!     Internally, the fpar's are used as following:
!     fpar(3) --- initial residual norm squared
!     fpar(4) --- target residual norm squared
!     fpar(5) --- current residual norm squared
!
!     w(:,1) -- R, residual
!     w(:,2) -- R0, the initial residual
!     w(:,3) -- W
!     w(:,4) -- Y
!     w(:,5) -- Z
!     w(:,6) -- A * Y
!     w(:,7) -- A * Z
!     w(:,8) -- V
!     w(:,9) -- D
!     w(:,10) -- intermediate results of preconditioning
!     w(:,11) -- changes in the solution
!-----------------------------------------------------------------------
!     external functions
!
      real*8 distdot
      logical brkdn
      external brkdn, distdot
!
      real*8 one,zero
      parameter(one=1.0D0,zero=0.0D0)
!
!     local variables
!
      integer i
      logical lp, rp
      real*8 eta,sigma,theta,te,alpha,rho,tao
      save
!
!     status of the call (where to go)
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10,20,40,50,60,70,80,90,100,110), ipar(10)
!
!     initializations
!
      call bisinit(ipar,fpar,11*n,2,lp,rp,w)
      if (ipar(1).lt.0) return
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + 6*n
      do i = 1, n
         w(i,1) = sol(i)
      enddo
      ipar(10) = 1
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      do i = 1, n
         w(i,1) = rhs(i) - w(i,7)
         w(i,9) = zero
      enddo
      fpar(11) = fpar(11) + n
!
      if (lp) then
         ipar(1) = 3
         ipar(9) = n+1
         ipar(10) = 2
         return
      endif
 20   continue
      if (lp) then
         do i = 1, n
            w(i,1) = w(i,2)
            w(i,3) = w(i,2)
         enddo
      else
         do i = 1, n
            w(i,2) = w(i,1)
            w(i,3) = w(i,1)
         enddo
      endif
!
      fpar(5) = sqrt(distdot(n,w,1,w,1))
      fpar(3) = fpar(5)
      tao = fpar(5)
      fpar(11) = fpar(11) + n + n
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + n + n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * tao + fpar(2)
      endif
      te = zero
      rho = zero
!
!     begin iteration
!
 30   sigma = rho
      rho = distdot(n,w(1,2),1,w(1,3),1)
      fpar(11) = fpar(11) + n + n
      if (brkdn(rho,ipar)) goto 900
      if (ipar(7).eq.1) then
         alpha = zero
      else
         alpha = rho / sigma
      endif
      do i = 1, n
         w(i,4) = w(i,3) + alpha * w(i,5)
      enddo
      fpar(11) = fpar(11) + n + n
!
!     A * x -- with preconditioning
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = 3*n + 1
         if (lp) then
            ipar(9) = 5*n + 1
         else
            ipar(9) = 9*n + 1
         endif
         ipar(10) = 3
         return
      endif
!
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 3*n + 1
      endif
      if (lp) then
         ipar(9) = 9*n + 1
      else
         ipar(9) = 5*n + 1
      endif
      ipar(10) = 4
      return
!
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 5*n + 1
         ipar(10) = 5
         return
      endif
 60   ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,8) = w(i,6) + alpha * (w(i,7) + alpha * w(i,8))
      enddo
      sigma = distdot(n,w(1,2),1,w(1,8),1)
      fpar(11) = fpar(11) + 6 * n
      if (brkdn(sigma,ipar)) goto 900
      alpha = rho / sigma
      do i = 1, n
         w(i,5) = w(i,4) - alpha * w(i,8)
      enddo
      fpar(11) = fpar(11) + 2*n
!
!     the second A * x
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = 4*n + 1
         if (lp) then
            ipar(9) = 6*n + 1
         else
            ipar(9) = 9*n + 1
         endif
         ipar(10) = 6
         return
      endif
!
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 4*n + 1
      endif
      if (lp) then
         ipar(9) = 9*n + 1
      else
         ipar(9) = 6*n + 1
      endif
      ipar(10) = 7
      return
!
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 6*n + 1
         ipar(10) = 8
         return
      endif
 90   ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,3) = w(i,3) - alpha * w(i,6)
      enddo
!
!     update I
!
      theta = distdot(n,w(1,3),1,w(1,3),1) / (tao*tao)
      sigma = one / (one + theta)
      tao = tao * sqrt(sigma * theta)
      fpar(11) = fpar(11) + 4*n + 6
      if (brkdn(tao,ipar)) goto 900
      eta = sigma * alpha
      sigma = te / alpha
      te = theta * eta
      do i = 1, n
         w(i,9) = w(i,4) + sigma * w(i,9)
         w(i,11) = w(i,11) + eta * w(i,9)
         w(i,3) = w(i,3) - alpha * w(i,7)
      enddo
      fpar(11) = fpar(11) + 6 * n + 6
      if (ipar(7).eq.1) then
         if (ipar(3).eq.-1) then
            fpar(3) = eta * sqrt(distdot(n,w(1,9),1,w(1,9),1))
            fpar(4) = fpar(1)*fpar(3) + fpar(2)
            fpar(11) = fpar(11) + n + n + 4
         endif
      endif
!
!     update II
!
      theta = distdot(n,w(1,3),1,w(1,3),1) / (tao*tao)
      sigma = one / (one + theta)
      tao = tao * sqrt(sigma * theta)
      fpar(11) = fpar(11) + 8 + 2*n
      if (brkdn(tao,ipar)) goto 900
      eta = sigma * alpha
      sigma = te / alpha
      te = theta * eta
      do i = 1, n
         w(i,9) = w(i,5) + sigma * w(i,9)
         w(i,11) = w(i,11) + eta * w(i,9)
      enddo
      fpar(11) = fpar(11) + 4*n + 3
!
!     this is the correct over-estimate
!      fpar(5) = sqrt(real(ipar(7)+1)) * tao
!     this is an approximation
      fpar(5) = tao
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 10*n + 1
         ipar(9) = 9*n + 1
         ipar(10) = 9
         return
      else if (ipar(3).lt.0) then
         fpar(6) = eta * sqrt(distdot(n,w(1,9),1,w(1,9),1))
         fpar(11) = fpar(11) + n + n + 2
      else
         fpar(6) = fpar(5)
      endif
      if (fpar(6).gt.fpar(4) .and. (ipar(7).lt.ipar(6) &
     &     .or. ipar(6).le.0)) goto 30
 100  if (ipar(3).eq.999.and.ipar(11).eq.0) goto 30
!
!     clean up
!
 900  if (rp) then
         if (ipar(1).lt.0) ipar(12) = ipar(1)
         ipar(1) = 5
         ipar(8) = 10*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      endif
 110  if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,10))
      else
         call tidycg(n,ipar,fpar,sol,w(1,11))
      endif
!
      return
      end
!-----end-of-tfqmr
!-----------------------------------------------------------------------
      subroutine fom(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
!-----------------------------------------------------------------------
!     This a version of The Full Orthogonalization Method (FOM) 
!     implemented with reverse communication. It is a simple restart 
!     version of the FOM algorithm and is implemented with plane 
!     rotations similarly to GMRES.
!
!  parameters:
!  ----------- 
!     ipar(5) == the dimension of the Krylov subspace
!     after every ipar(5) iterations, the FOM will restart with
!     the updated solution and recomputed residual vector.
!
!     the work space in `w' is used as follows:
!     (1) the basis for the Krylov subspace, size n*(m+1);
!     (2) the Hessenberg matrix, only the upper triangular
!     portion of the matrix is stored, size (m+1)*m/2 + 1
!     (3) three vectors, all are of size m, they are
!     the cosine and sine of the Givens rotations, the third one holds
!     the residuals, it is of size m+1.
!
!     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
!     Note: m == ipar(5). The default value for this is 15 if
!     ipar(5) <= 1.
!-----------------------------------------------------------------------
!     external functions used
!
      real*8 distdot
      external distdot
!
      real*8 one, zero
      parameter(one=1.0D0, zero=0.0D0)
!
!     local variables, ptr and p2 are temporary pointers,
!     hes points to the Hessenberg matrix,
!     vc, vs point to the cosines and sines of the Givens rotations
!     vrn points to the vectors of residual norms, more precisely
!     the right hand side of the least square problem solved.
!
      integer i,ii,idx,k,m,ptr,p2,prs,hes,vc,vs,vrn
      real*8 alpha, c, s 
      logical lp, rp
      save
!
!     check the status of the call
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
!
!     initialization
!
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      endif
      idx = n * (m+1)
      hes = idx + n
      vc = hes + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
!
!     request for matrix vector multiplication A*x in the initialization
!
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do i = 1, n
         w(n+i) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      if (lp) then
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         enddo
         ipar(1) = 3
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
      endif
      fpar(11) = fpar(11) + n
!
 20   alpha = sqrt(distdot(n,w,1,w,1))
      fpar(11) = fpar(11) + 2*n + 1
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
         if (abs(ipar(3)).eq.2) then
            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
            fpar(11) = fpar(11) + 2*n
         else
            fpar(4) = fpar(1) * alpha + fpar(2)
         endif
         fpar(3) = alpha
      endif
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 300
      endif
      alpha = one / alpha
      do ii = 1, n
         w(ii) = alpha * w(ii)
      enddo
      fpar(11) = fpar(11) + n
!
!     request for (1) right preconditioning
!     (2) matrix vector multiplication
!     (3) left preconditioning
!
 110  k = k + 1
      if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         if (lp) then
            ipar(9) = k*n + 1
         else
            ipar(9) = idx + 1
         endif
         ipar(10) = 3
         return
      endif
!
 30   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (k-1)*n + 1
      endif
      if (lp) then
         ipar(9) = idx + 1
      else
         ipar(9) = 1 + k*n
      endif
      ipar(10) = 4
      return
!
 40   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = k*n + 1
         ipar(10) = 5
         return
      endif
!
!     Modified Gram-Schmidt orthogonalization procedure
!     temporary pointer 'ptr' is pointing to the current column of the
!     Hessenberg matrix. 'p2' points to the new basis vector
!
 50   ipar(7) = ipar(7) + 1
      ptr = k * (k - 1) / 2 + hes
      p2 = ipar(9)
      call mgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),     &
     &     ipar(12))
      if (ipar(12).lt.0) goto 200
!
!     apply previous Givens rotations to column.
!
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      enddo
!
!     end of one Arnoldi iteration, alpha will store the estimated
!     residual norm at current stage
!
      fpar(11) = fpar(11) + 6*k

      prs = vrn+k
      alpha = fpar(5) 
      if (w(p2) .ne. zero) alpha = abs(w(p2+1)*w(prs)/w(p2)) 
      fpar(5) = alpha
!
      if ((k.ge.m) .or. (ipar(3).ge.0 .and. alpha.le.fpar(4)) &
     &     .or. (ipar(6).gt.0 .and. ipar(7).ge.ipar(6))) &
     &     goto 200
!
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      alpha = - s * w(prs)
      w(prs) = c * w(prs)
      w(prs+1) = alpha
!
      if (w(p2).ne.zero) goto 110
!
!     update the approximate solution, first solve the upper triangular
!     system, temporary pointer ptr points to the Hessenberg matrix,
!     prs points to the right-hand-side (also the solution) of the system.
!
 200  ptr = hes + k * (k + 1) / 2
      prs = vrn + k
      if (w(ptr).eq.zero) then
!
!     if the diagonal elements of the last column is zero, reduce k by 1
!     so that a smaller trianguler system is solved
!
         k = k - 1
         if (k.gt.0) then
            goto 200
         else
            ipar(1) = -3
            ipar(12) = -4
            goto 300
         endif
      endif
      w(prs) = w(prs) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i - 1
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(prs) * w(ptr+ii)
         enddo
         prs = prs - 1
         w(prs) = w(prs) / w(ptr)
      enddo
!
      do ii = 1, n
         w(ii) = w(ii) * w(prs)
      enddo
      do i = 1, k-1
         prs = prs + 1
         ptr = i*n
         do ii = 1, n
            w(ii) = w(ii) + w(prs) * w(ptr+ii)
         enddo
      enddo
      fpar(11) = fpar(11) + 2*(k-1)*n + n + k*(k+1)
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = 1
         ipar(9) = idx + 1
         ipar(10) = 6
         return
      endif
!
 60   if (rp) then
         do i = 1, n
            sol(i) = sol(i) + w(idx+i)
         enddo
      else
         do i = 1, n
            sol(i) = sol(i) + w(i)
         enddo
      endif
      fpar(11) = fpar(11) + n
!
!     process the complete stopping criteria
!
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = idx + 1
         ipar(10) = 7
         return
      else if (ipar(3).lt.0) then
         if (ipar(7).le.m+1) then
            fpar(3) = abs(w(vrn+1))
            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         endif
         alpha = abs(w(vrn+k))
      endif
      fpar(6) = alpha
!
!     do we need to restart ?
!
 70   if (ipar(12).ne.0) then
         ipar(1) = -3
         goto 300
      endif
      if (ipar(7).lt.ipar(6) .or. ipar(6).le.0) then
         if (ipar(3).ne.999) then
            if (fpar(6).gt.fpar(4)) goto 100
         else
            if (ipar(11).eq.0) goto 100
         endif
      endif
!
!     termination, set error code, compute convergence rate
!
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero .and. &
     &     ipar(7).gt.ipar(13)) then
         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      else
         fpar(7) = zero
      endif
      return
      end
!-----end-of-fom-------------------------------------------------------- 
!-----------------------------------------------------------------------
      subroutine gmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
!-----------------------------------------------------------------------
!     This a version of GMRES implemented with reverse communication.
!     It is a simple restart version of the GMRES algorithm.
!
!     ipar(5) == the dimension of the Krylov subspace
!     after every ipar(5) iterations, the GMRES will restart with
!     the updated solution and recomputed residual vector.
!
!     the space of the `w' is used as follows:
!     (1) the basis for the Krylov subspace, size n*(m+1);
!     (2) the Hessenberg matrix, only the upper triangular
!     portion of the matrix is stored, size (m+1)*m/2 + 1
!     (3) three vectors, all are of size m, they are
!     the cosine and sine of the Givens rotations, the third one holds
!     the residuals, it is of size m+1.
!
!     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
!     Note: m == ipar(5). The default value for this is 15 if
!     ipar(5) <= 1.
!-----------------------------------------------------------------------
!     external functions used
!
      real*8 distdot
      external distdot
!
      real*8 one, zero
      parameter(one=1.0D0, zero=0.0D0)
!
!     local variables, ptr and p2 are temporary pointers,
!     hess points to the Hessenberg matrix,
!     vc, vs point to the cosines and sines of the Givens rotations
!     vrn points to the vectors of residual norms, more precisely
!     the right hand side of the least square problem solved.
!
      integer i,ii,idx,k,m,ptr,p2,hess,vc,vs,vrn
      real*8 alpha, c, s
      logical lp, rp
      save
!
!     check the status of the call
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
!
!     initialization
!
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      endif
      idx = n * (m+1)
      hess = idx + n
      vc = hess + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
!
!     request for matrix vector multiplication A*x in the initialization
!
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do i = 1, n
         w(n+i) = sol(i)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      if (lp) then
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         enddo
         ipar(1) = 3
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
      endif
      fpar(11) = fpar(11) + n
!
 20   alpha = sqrt(distdot(n,w,1,w,1))
      fpar(11) = fpar(11) + 2*n
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
         if (abs(ipar(3)).eq.2) then
            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
            fpar(11) = fpar(11) + 2*n
         else
            fpar(4) = fpar(1) * alpha + fpar(2)
         endif
         fpar(3) = alpha
      endif
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 300
      endif
      alpha = one / alpha
      do ii = 1, n
         w(ii) = alpha * w(ii)
      enddo
      fpar(11) = fpar(11) + n
!
!     request for (1) right preconditioning
!     (2) matrix vector multiplication
!     (3) left preconditioning
!
 110  k = k + 1
      if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         if (lp) then
            ipar(9) = k*n + 1
         else
            ipar(9) = idx + 1
         endif
         ipar(10) = 3
         return
      endif
!
 30   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (k-1)*n + 1
      endif
      if (lp) then
         ipar(9) = idx + 1
      else
         ipar(9) = 1 + k*n
      endif
      ipar(10) = 4
      return
!
 40   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = k*n + 1
         ipar(10) = 5
         return
      endif
!
!     Modified Gram-Schmidt orthogonalization procedure
!     temporary pointer 'ptr' is pointing to the current column of the
!     Hessenberg matrix. 'p2' points to the new basis vector
!
 50   ipar(7) = ipar(7) + 1
      ptr = k * (k - 1) / 2 + hess
      p2 = ipar(9)
      call mgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1), ipar(12))
      if (ipar(12).lt.0) goto 200
!
!     apply previous Givens rotations and generate a new one to eliminate
!     the subdiagonal element.
!
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      enddo
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      p2 = vrn + k
      alpha = - s * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = alpha
!
!     end of one Arnoldi iteration, alpha will store the estimated
!     residual norm at current stage
!
      fpar(11) = fpar(11) + 6*k + 2
      alpha = abs(alpha)
      fpar(5) = alpha
      if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.fpar(4)) &
     &     .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
!
!     update the approximate solution, first solve the upper triangular
!     system, temporary pointer ptr points to the Hessenberg matrix,
!     p2 points to the right-hand-side (also the solution) of the system.
!
 200  ptr = hess + k * (k + 1) / 2
      p2 = vrn + k
      if (w(ptr).eq.zero) then
!
!     if the diagonal elements of the last column is zero, reduce k by 1
!     so that a smaller trianguler system is solved [It should only
!     happen when the matrix is singular, and at most once!]
!
         k = k - 1
         if (k.gt.0) then
            goto 200
         else
            ipar(1) = -3
            ipar(12) = -4
            goto 300
         endif
      endif
      w(p2) = w(p2) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i - 1
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
         enddo
         p2 = p2 - 1
         w(p2) = w(p2) / w(ptr)
      enddo
!
      do ii = 1, n
         w(ii) = w(ii) * w(p2)
      enddo
      do i = 1, k-1
         ptr = i*n
         p2 = p2 + 1
         do ii = 1, n
            w(ii) = w(ii) + w(p2) * w(ptr+ii)
         enddo
      enddo
      fpar(11) = fpar(11) + 2*k*n - n + k*(k+1)
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = 1
         ipar(9) = idx + 1
         ipar(10) = 6
         return
      endif
!
 60   if (rp) then
         do i = 1, n
            sol(i) = sol(i) + w(idx+i)
         enddo
      else
         do i = 1, n
            sol(i) = sol(i) + w(i)
         enddo
      endif
      fpar(11) = fpar(11) + n
!
!     process the complete stopping criteria
!
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = idx + 1
         ipar(10) = 7
         return
      else if (ipar(3).lt.0) then
         if (ipar(7).le.m+1) then
            fpar(3) = abs(w(vrn+1))
            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         endif
         fpar(6) = abs(w(vrn+k))
      else
         fpar(6) = fpar(5)
      endif
!
!     do we need to restart ?
!
 70   if (ipar(12).ne.0) then
         ipar(1) = -3
         goto 300
      endif
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0) .and. &
     &     ((ipar(3).eq.999.and.ipar(11).eq.0) .or. &
     &     (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
!
!     termination, set error code, compute convergence rate
!
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero .and. &
     &     ipar(7).gt.ipar(13)) then
         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      else
         fpar(7) = zero
      endif
      return
      end
!-----end-of-gmres
!-----------------------------------------------------------------------
      subroutine dqgmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
!-----------------------------------------------------------------------
!     DQGMRES -- Flexible Direct version of Quasi-General Minimum
!     Residual method. The right preconditioning can be varied from
!     step to step.
!
!     Work space used = n + lb * (2*n+4)
!     where lb = ipar(5) + 1 (default 16 if ipar(5) <= 1)
!-----------------------------------------------------------------------
!     local variables
!
      real*8 one,zero,deps
      parameter(one=1.0D0,zero=0.0D0)
      parameter(deps=1.0D-33)
!
      integer i,ii,j,jp1,j0,k,ptrw,ptrv,iv,iw,ic,is,ihm,ihd,lb,ptr
      real*8 alpha,beta,psi,c,s,distdot
      logical lp,rp,full
      external distdot,bisinit
      save
!
!     where to go
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 40, 50, 60, 70) ipar(10)
!
!     locations of the work arrays. The arrangement is as follows:
!     w(1:n) -- temporary storage for the results of the preconditioning
!     w(iv+1:iw) -- the V's
!     w(iw+1:ic) -- the W's
!     w(ic+1:is) -- the COSINEs of the Givens rotations
!     w(is+1:ihm) -- the SINEs of the Givens rotations
!     w(ihm+1:ihd) -- the last column of the Hessenberg matrix
!     w(ihd+1:i) -- the inverse of the diagonals of the Hessenberg matrix
!
      if (ipar(5).le.1) then
         lb = 16
      else
         lb = ipar(5) + 1
      endif
      iv = n
      iw = iv + lb * n
      ic = iw + lb * n
      is = ic + lb
      ihm = is + lb
      ihd = ihm + lb
      i = ihd + lb
!
!     parameter check, initializations
!
      full = .false.
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
      ipar(1) = 1
      if (lp) then
         do ii = 1, n
            w(iv+ii) = sol(ii)
         enddo
         ipar(8) = iv+1
         ipar(9) = 1
      else
         do ii = 1, n
            w(ii) = sol(ii)
         enddo
         ipar(8) = 1
         ipar(9) = iv+1
      endif
      ipar(10) = 1
      return
!
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      if (lp) then
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
         ipar(1) = 3
         ipar(8) = 1
         ipar(9) = iv+1
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(iv+i) = rhs(i) - w(iv+i)
         enddo
      endif
      fpar(11) = fpar(11) + n
!
 20   alpha = sqrt(distdot(n, w(iv+1), 1, w(iv+1), 1))
      fpar(11) = fpar(11) + (n + n)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2*n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * alpha + fpar(2)
      endif
      fpar(3) = alpha
      fpar(5) = alpha
      psi = alpha
      if (alpha.le.fpar(4)) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 80
      endif
      alpha = one / alpha
      do i = 1, n
         w(iv+i) = w(iv+i) * alpha
      enddo
      fpar(11) = fpar(11) + n
      j = 0
!
!     iterations start here
!
 30   j = j + 1
      if (j.gt.lb) j = j - lb
      jp1 = j + 1
      if (jp1.gt.lb) jp1 = jp1 - lb
      ptrv = iv + (j-1)*n + 1
      ptrw = iv + (jp1-1)*n + 1
      if (.not.full) then
         if (j.gt.jp1) full = .true.
      endif
      if (full) then
         j0 = jp1+1
         if (j0.gt.lb) j0 = j0 - lb
      else
         j0 = 1
      endif
!
!     request the caller to perform matrix-vector multiplication and
!     preconditioning
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = ptrv
         ipar(9) = ptrv + iw - iv
         ipar(10) = 3
         return
      else
         do i = 0, n-1
            w(ptrv+iw-iv+i) = w(ptrv+i)
         enddo
      endif
!
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = ptrv
      endif
      if (lp) then
         ipar(9) = 1
      else
         ipar(9) = ptrw
      endif
      ipar(10) = 4
      return
!
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = ptrw
         ipar(10) = 5
         return
      endif
!
!     compute the last column of the Hessenberg matrix
!     modified Gram-schmidt procedure, orthogonalize against (lb-1)
!     previous vectors
!
 60   continue
      call mgsro(full,n,n,lb,jp1,fpar(11),w(iv+1),w(ihm+1),     &
     &     ipar(12))
      if (ipar(12).lt.0) then
         ipar(1) = -3
         goto 80
      endif
      beta = w(ihm+jp1)
!
!     incomplete factorization (QR factorization through Givens rotations)
!     (1) apply previous rotations [(lb-1) of them]
!     (2) generate a new rotation
!
      if (full) then
         w(ihm+jp1) = w(ihm+j0) * w(is+jp1)
         w(ihm+j0) = w(ihm+j0) * w(ic+jp1)
      endif
      i = j0
      do while (i.ne.j)
         k = i+1
         if (k.gt.lb) k = k - lb
         c = w(ic+i)
         s = w(is+i)
         alpha = w(ihm+i)
         w(ihm+i) = c * alpha + s * w(ihm+k)
         w(ihm+k) = c * w(ihm+k) - s * alpha
         i = k
      enddo
      call givens(w(ihm+j), beta, c, s)
      if (full) then
         fpar(11) = fpar(11) + 6 * lb
      else
         fpar(11) = fpar(11) + 6 * j
      endif
!
!     detect whether diagonal element of this column is zero
!
      if (abs(w(ihm+j)).lt.deps) then
         ipar(1) = -3
         goto 80
      endif
      w(ihd+j) = one / w(ihm+j)
      w(ic+j) = c
      w(is+j) = s
!
!     update the W's (the conjugate directions) -- essentially this is one
!     step of triangular solve.
!
      ptrw = iw+(j-1)*n + 1
      if (full) then
         do i = j+1, lb
            alpha = -w(ihm+i)*w(ihd+i)
            ptr = iw+(i-1)*n+1
            do ii = 0, n-1
               w(ptrw+ii) = w(ptrw+ii) + alpha * w(ptr+ii)
            enddo
         enddo
      endif
      do i = 1, j-1
         alpha = -w(ihm+i)*w(ihd+i)
         ptr = iw+(i-1)*n+1
         do ii = 0, n-1
            w(ptrw+ii) = w(ptrw+ii) + alpha * w(ptr+ii)
         enddo
      enddo
!
!     update the solution to the linear system
!
      alpha = psi * c * w(ihd+j)
      psi = - s * psi
      do i = 1, n
         sol(i) = sol(i) + alpha * w(ptrw-1+i)
      enddo
      if (full) then
         fpar(11) = fpar(11) + lb * (n+n)
      else
         fpar(11) = fpar(11) + j * (n+n)
      endif
!
!     determine whether to continue,
!     compute the desired error/residual norm
!
      ipar(7) = ipar(7) + 1
      fpar(5) = abs(psi)
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = 1
         ipar(10) = 6
         return
      endif
      if (ipar(3).lt.0) then
         alpha = abs(alpha)
         if (ipar(7).eq.2 .and. ipar(3).eq.-1) then
            fpar(3) = alpha*sqrt(distdot(n, w(ptrw), 1, w(ptrw), 1))
            fpar(4) = fpar(1) * fpar(3) + fpar(2)
            fpar(6) = fpar(3)
         else
            fpar(6) = alpha*sqrt(distdot(n, w(ptrw), 1, w(ptrw), 1))
         endif
         fpar(11) = fpar(11) + 2 * n
      else
         fpar(6) = fpar(5)
      endif
      if (ipar(1).ge.0 .and. fpar(6).gt.fpar(4) .and. (ipar(6).le.0 &
     &     .or. ipar(7).lt.ipar(6))) goto 30
 70   if (ipar(3).eq.999 .and. ipar(11).eq.0) goto 30
!
!     clean up the iterative solver
!
 80   fpar(7) = zero
      if (fpar(3).ne.zero .and. fpar(6).ne.zero .and. &
     &     ipar(7).gt.ipar(13)) &
     &     fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).ne.0) then
            ipar(1) = 0
         else if (fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
      return
      end
!-----end-of-dqgmres
!-----------------------------------------------------------------------
      subroutine fgmres(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
!-----------------------------------------------------------------------
!     This a version of FGMRES implemented with reverse communication.
!
!     ipar(5) == the dimension of the Krylov subspace
!
!     the space of the `w' is used as follows:
!     >> V: the bases for the Krylov subspace, size n*(m+1);
!     >> W: the above bases after (left-)multiplying with the
!     right-preconditioner inverse, size m*n;
!     >> a temporary vector of size n;
!     >> the Hessenberg matrix, only the upper triangular portion
!     of the matrix is stored, size (m+1)*m/2 + 1
!     >> three vectors, first two are of size m, they are the cosine
!     and sine of the Givens rotations, the third one holds the
!     residuals, it is of size m+1.
!
!     TOTAL SIZE REQUIRED == n*(2m+1) + (m+1)*m/2 + 3*m + 2
!     Note: m == ipar(5). The default value for this is 15 if
!     ipar(5) <= 1.
!-----------------------------------------------------------------------
!     external functions used
!
      real*8 distdot
      external distdot
!
      real*8 one, zero
      parameter(one=1.0D0, zero=0.0D0)
!
!     local variables, ptr and p2 are temporary pointers,
!     hess points to the Hessenberg matrix,
!     vc, vs point to the cosines and sines of the Givens rotations
!     vrn points to the vectors of residual norms, more precisely
!     the right hand side of the least square problem solved.
!
      integer i,ii,idx,iz,k,m,ptr,p2,hess,vc,vs,vrn
      real*8 alpha, c, s
      logical lp, rp
      save
!
!     check the status of the call
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20, 30, 40, 50, 60) ipar(10)
!
!     initialization
!
      if (ipar(5).le.1) then
         m = 15
      else
         m = ipar(5)
      endif
      idx = n * (m+1)
      iz = idx + n
      hess = iz + n*m
      vc = hess + (m+1) * m / 2 + 1
      vs = vc + m
      vrn = vs + m
      i = vrn + m + 1
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
!
!     request for matrix vector multiplication A*x in the initialization
!
 100  ipar(1) = 1
      ipar(8) = n+1
      ipar(9) = 1
      ipar(10) = 1
      k = 0
      do ii = 1, n
         w(ii+n) = sol(ii)
      enddo
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      fpar(11) = fpar(11) + n
      if (lp) then
         do i = 1, n
            w(n+i) = rhs(i) - w(i)
         enddo
         ipar(1) = 3
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i) = rhs(i) - w(i)
         enddo
      endif
!
 20   alpha = sqrt(distdot(n,w,1,w,1))
      fpar(11) = fpar(11) + n + n
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
         if (abs(ipar(3)).eq.2) then
            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
            fpar(11) = fpar(11) + 2*n
         else
            fpar(4) = fpar(1) * alpha + fpar(2)
         endif
         fpar(3) = alpha
      endif
      fpar(5) = alpha
      w(vrn+1) = alpha
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
         ipar(1) = 0
         fpar(6) = alpha
         goto 300
      endif
      alpha = one / alpha
      do ii = 1, n
         w(ii) = w(ii) * alpha
      enddo
      fpar(11) = fpar(11) + n
!
!     request for (1) right preconditioning
!     (2) matrix vector multiplication
!     (3) left preconditioning
!
 110  k = k + 1
      if (rp) then
         ipar(1) = 5
         ipar(8) = k*n - n + 1
         ipar(9) = iz + ipar(8)
         ipar(10) = 3
         return
      else
         do ii = 0, n-1
            w(iz+k*n-ii) = w(k*n-ii)
         enddo
      endif
!
 30   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (k-1)*n + 1
      endif
      if (lp) then
         ipar(9) = idx + 1
      else
         ipar(9) = 1 + k*n
      endif
      ipar(10) = 4
      return
!
 40   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = k*n + 1
         ipar(10) = 5
         return
      endif
!
!     Modified Gram-Schmidt orthogonalization procedure
!     temporary pointer 'ptr' is pointing to the current column of the
!     Hessenberg matrix. 'p2' points to the new basis vector
!
 50   ptr = k * (k - 1) / 2 + hess
      p2 = ipar(9)
      ipar(7) = ipar(7) + 1
      call mgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),     &
     &     ipar(12))
      if (ipar(12).lt.0) goto 200
!
!     apply previous Givens rotations and generate a new one to eliminate
!     the subdiagonal element.
!
      p2 = ptr + 1
      do i = 1, k-1
         ptr = p2
         p2 = p2 + 1
         alpha = w(ptr)
         c = w(vc+i)
         s = w(vs+i)
         w(ptr) = c * alpha + s * w(p2)
         w(p2) = c * w(p2) - s * alpha
      enddo
      call givens(w(p2), w(p2+1), c, s)
      w(vc+k) = c
      w(vs+k) = s
      p2 = vrn + k
      alpha = - s * w(p2)
      w(p2) = c * w(p2)
      w(p2+1) = alpha
      fpar(11) = fpar(11) + 6 * k
!
!     end of one Arnoldi iteration, alpha will store the estimated
!     residual norm at current stage
!
      alpha = abs(alpha)
      fpar(5) = alpha
      if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.fpar(4)) &
     &      .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
!
!     update the approximate solution, first solve the upper triangular
!     system, temporary pointer ptr points to the Hessenberg matrix,
!     p2 points to the right-hand-side (also the solution) of the system.
!
 200  ptr = hess + k * (k + 1 ) / 2
      p2 = vrn + k
      if (w(ptr).eq.zero) then
!
!     if the diagonal elements of the last column is zero, reduce k by 1
!     so that a smaller trianguler system is solved [It should only
!     happen when the matrix is singular!]
!
         k = k - 1
         if (k.gt.0) then
            goto 200
         else
            ipar(1) = -3
            ipar(12) = -4
            goto 300
         endif
      endif
      w(p2) = w(p2) / w(ptr)
      do i = k-1, 1, -1
         ptr = ptr - i - 1
         do ii = 1, i
            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
         enddo
         p2 = p2 - 1
         w(p2) = w(p2) / w(ptr)
      enddo
!
      do i = 0, k-1
         ptr = iz+i*n
         do ii = 1, n
            sol(ii) = sol(ii) + w(p2)*w(ptr+ii)
         enddo
         p2 = p2 + 1
      enddo
      fpar(11) = fpar(11) + 2*k*n + k*(k+1)
!
!     process the complete stopping criteria
!
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = -1
         ipar(9) = idx + 1
         ipar(10) = 6
         return
      else if (ipar(3).lt.0) then
         if (ipar(7).le.m+1) then
            fpar(3) = abs(w(vrn+1))
            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
         endif
         fpar(6) = abs(w(vrn+k))
      else if (ipar(3).ne.999) then
         fpar(6) = fpar(5)
      endif
!
!     do we need to restart ?
!
 60   if (ipar(12).ne.0) then
         ipar(1) = -3
         goto 300
      endif
      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0).and. &
     &     ((ipar(3).eq.999.and.ipar(11).eq.0) .or. &
     &     (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
!
!     termination, set error code, compute convergence rate
!
      if (ipar(1).gt.0) then
         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
            ipar(1) = 0
         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero .and. &
     &     ipar(7).gt.ipar(13)) then
         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      else
         fpar(7) = zero
      endif
      return
      end
!-----end-of-fgmres
!-----------------------------------------------------------------------
      subroutine dbcg (n,rhs,sol,ipar,fpar,w)
      implicit none
      integer n,ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(n,*)
!-----------------------------------------------------------------------
! Quasi GMRES method for solving a linear
! system of equations a * sol = y.  double precision version.
! this version is without restarting and without preconditioning.
! parameters :
! -----------
! n     = dimension of the problem
!
! y     = w(:,1) a temporary storage used for various operations
! z     = w(:,2) a work vector of length n.
! v     = w(:,3:4) size n x 2
! w     = w(:,5:6) size n x 2
! p     = w(:,7:9) work array of dimension n x 3
! del x = w(:,10)  accumulation of the changes in solution
! tmp   = w(:,11)  a temporary vector used to hold intermediate result of
!                  preconditioning, etc.
!
! sol   = the solution of the problem . at input sol must contain an
!         initial guess to the solution.
!    ***  note:   y is destroyed on return.
!
!-----------------------------------------------------------------------
! subroutines and functions called:
! 1) matrix vector multiplication and preconditioning through reverse
!     communication
!
! 2) implu, uppdir, distdot (blas)
!-----------------------------------------------------------------------
! aug. 1983  version.    author youcef saad. yale university computer
! science dept. some  changes made july 3, 1986.
! references: siam j. sci. stat. comp., vol. 5, pp. 203-228 (1984)
!-----------------------------------------------------------------------
!     local variables
!
      real*8 one,zero
      parameter(one=1.0D0,zero=0.0D0)
!
      real*8 t,sqrt,distdot,ss,res,beta,ss1,delta,x,zeta,umm
      integer k,j,i,i2,ip2,ju,lb,lbm1,np,indp
      logical lp,rp,full, perm(3)
      real*8 ypiv(3),u(3),usav(3)
      external tidycg
      save
!
!     where to go
!
      if (ipar(1).le.0) ipar(10) = 0
      goto (110, 120, 130, 140, 150, 160, 170, 180, 190, 200) ipar(10)
!
!     initialization, parameter checking, clear the work arrays
!
      call bisinit(ipar,fpar,11*n,1,lp,rp,w)
      if (ipar(1).lt.0) return
      perm(1) = .false.
      perm(2) = .false.
      perm(3) = .false.
      usav(1) = zero
      usav(2) = zero
      usav(3) = zero
      ypiv(1) = zero
      ypiv(2) = zero
      ypiv(3) = zero
!-----------------------------------------------------------------------
!     initialize constants for outer loop :
!-----------------------------------------------------------------------
      lb = 3
      lbm1 = 2
!
!     get initial residual vector and norm
!
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + n
      do i = 1, n
         w(i,1) = sol(i)
      enddo
      ipar(10) = 1
      return
 110  ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      if (lp) then
         do i = 1, n
            w(i,1) = rhs(i) - w(i,2)
         enddo
         ipar(1) = 3
         ipar(8) = 1
         ipar(9) = n+n+1
         ipar(10) = 2
         return
      else
         do i = 1, n
            w(i,3) = rhs(i) - w(i,2)
         enddo
      endif
      fpar(11) = fpar(11) + n
!
 120  fpar(3) = sqrt(distdot(n,w(1,3),1,w(1,3),1))
      fpar(11) = fpar(11) + n + n
      fpar(5) = fpar(3)
      fpar(7) = fpar(3)
      zeta = fpar(3)
      if (abs(ipar(3)).eq.2) then
         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
         fpar(11) = fpar(11) + 2*n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * zeta + fpar(2)
      endif
      if (ipar(3).ge.0.and.fpar(5).le.fpar(4)) then
         fpar(6) = fpar(5)
         goto 900
      endif
!
!     normalize first arnoldi vector
!
      t = one/zeta
      do 22 k=1,n
         w(k,3) = w(k,3)*t
         w(k,5) = w(k,3)
 22   continue
      fpar(11) = fpar(11) + n
!
!     initialize constants for main loop
!
      beta = zero
      delta = zero
      i2 = 1
      indp = 0
      i = 0
!
!     main loop: i = index of the loop.
!
!-----------------------------------------------------------------------
 30   i = i + 1
!
      if (rp) then
         ipar(1) = 5
         ipar(8) = (1+i2)*n+1
         if (lp) then
            ipar(9) = 1
         else
            ipar(9) = 10*n + 1
         endif
         ipar(10) = 3
         return
      endif
!
 130  ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (1+i2)*n + 1
      endif
      if (lp) then
         ipar(9) = 10*n + 1
      else
         ipar(9) = 1
      endif
      ipar(10) = 4
      return
!
 140  if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 1
         ipar(10) = 5
         return
      endif
!
!     A^t * x
!
 150  ipar(7) = ipar(7) + 1
      if (lp) then
         ipar(1) = 4
         ipar(8) = (3+i2)*n + 1
         if (rp) then
            ipar(9) = n + 1
         else
            ipar(9) = 10*n + 1
         endif
         ipar(10) = 6
         return
      endif
!
 160  ipar(1) = 2
      if (lp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = (3+i2)*n + 1
      endif
      if (rp) then
         ipar(9) = 10*n + 1
      else
         ipar(9) = n + 1
      endif
      ipar(10) = 7
      return
!
 170  if (rp) then
         ipar(1) = 6
         ipar(8) = ipar(9)
         ipar(9) = n + 1
         ipar(10) = 8
         return
      endif
!-----------------------------------------------------------------------
!     orthogonalize current v against previous v's and
!     determine relevant part of i-th column of u(.,.) the
!     upper triangular matrix --
!-----------------------------------------------------------------------
 180  ipar(7) = ipar(7) + 1
      u(1) = zero
      ju = 1
      k = i2
      if (i .le. lbm1) ju = 0
      if (i .lt. lb) k = 0
 31   if (k .eq. lbm1) k=0
      k=k+1
!
      if (k .ne. i2) then
         ss  = delta
         ss1 = beta
         ju = ju + 1
         u(ju) = ss
      else
         ss = distdot(n,w(1,1),1,w(1,4+k),1)
         fpar(11) = fpar(11) + 2*n
         ss1= ss
         ju = ju + 1
         u(ju) = ss
      endif
!
      do 32  j=1,n
         w(j,1) = w(j,1) - ss*w(j,k+2)
         w(j,2) = w(j,2) - ss1*w(j,k+4)
 32   continue
      fpar(11) = fpar(11) + 4*n
!
      if (k .ne. i2) goto 31
!
!     end of Mod. Gram. Schmidt loop
!
      t = distdot(n,w(1,2),1,w(1,1),1)
!
      beta   = sqrt(abs(t))
      delta  = t/beta
!
      ss = one/beta
      ss1 = one/ delta
!
!     normalize and insert new vectors
!
      ip2 = i2
      if (i2 .eq. lbm1) i2=0
      i2=i2+1
!
      do 315 j=1,n
         w(j,i2+2)=w(j,1)*ss
         w(j,i2+4)=w(j,2)*ss1
 315  continue
      fpar(11) = fpar(11) + 4*n
!-----------------------------------------------------------------------
!     end of orthogonalization.
!     now compute the coefficients u(k) of the last
!     column of the  l . u  factorization of h .
!-----------------------------------------------------------------------
      np = min0(i,lb)
      full = (i .ge. lb)
      call implu(np, umm, beta, ypiv, u, perm, full)
!-----------------------------------------------------------------------
!     update conjugate directions and solution
!-----------------------------------------------------------------------
      do 33 k=1,n
         w(k,1) = w(k,ip2+2)
 33   continue
      call uppdir(n, w(1,7), np, lb, indp, w, u, usav, fpar(11))
!-----------------------------------------------------------------------
      if (i .eq. 1) goto 34
      j = np - 1
      if (full) j = j-1
      if (.not.perm(j)) zeta = -zeta*ypiv(j)
 34   x = zeta/u(np)
      if (perm(np))goto 36
      do 35 k=1,n
         w(k,10) = w(k,10) + x*w(k,1)
 35   continue
      fpar(11) = fpar(11) + 2 * n
!-----------------------------------------------------------------------
 36   if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 9*n + 1
         ipar(9) = 10*n + 1
         ipar(10) = 9
         return
      endif
      res = abs(beta*zeta/umm)
      fpar(5) = res * sqrt(distdot(n, w(1,i2+2), 1, w(1,i2+2), 1))
      fpar(11) = fpar(11) + 2 * n
      if (ipar(3).lt.0) then
         fpar(6) = x * sqrt(distdot(n,w,1,w,1))
         fpar(11) = fpar(11) + 2 * n
         if (ipar(7).le.3) then
            fpar(3) = fpar(6)
            if (ipar(3).eq.-1) then
               fpar(4) = fpar(1) * sqrt(fpar(3)) + fpar(2)
            endif
         endif
      else
         fpar(6) = fpar(5)
      endif
!---- convergence test -----------------------------------------------
 190  if (ipar(3).eq.999.and.ipar(11).eq.0) then
         goto 30
      else if (fpar(6).gt.fpar(4) .and. (ipar(6).gt.ipar(7) .or. &
     &        ipar(6).le.0)) then
         goto 30
      endif
!-----------------------------------------------------------------------
!     here the fact that the last step is different is accounted for.
!-----------------------------------------------------------------------
      if (.not. perm(np)) goto 900
      x = zeta/umm
      do 40 k = 1,n
         w(k,10) = w(k,10) + x*w(k,1)
 40   continue
      fpar(11) = fpar(11) + 2 * n
!
!     right preconditioning and clean-up jobs
!
 900  if (rp) then
         if (ipar(1).lt.0) ipar(12) = ipar(1)
         ipar(1) = 5
         ipar(8) = 9*n + 1
         ipar(9) = ipar(8) + n
         ipar(10) = 10
         return
      endif
 200  if (rp) then
         call tidycg(n,ipar,fpar,sol,w(1,11))
      else
         call tidycg(n,ipar,fpar,sol,w(1,10))
      endif
      return
      end
!-----end-of-dbcg-------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine implu(np,umm,beta,ypiv,u,permut,full)
      real*8 umm,beta,ypiv(*),u(*),x, xpiv
      logical full, perm, permut(*)
      integer np,k,npm1
!-----------------------------------------------------------------------
!     performs implicitly one step of the lu factorization of a
!     banded hessenberg matrix.
!-----------------------------------------------------------------------
      npm1  = 0
      if (np .le. 1) goto 12
      npm1 = np - 1
!
!     -- perform  previous step of the factorization-
!
      do 6 k=1,npm1
         if (.not. permut(k)) goto 5
         x=u(k)
         u(k) = u(k+1)
         u(k+1) = x
 5       u(k+1) = u(k+1) - ypiv(k)*u(k)
 6    continue
!-----------------------------------------------------------------------
!     now determine pivotal information to be used in the next call
!-----------------------------------------------------------------------
 12   umm = u(np)
      perm = (beta .gt. abs(umm))
      if (.not. perm) goto 4
      xpiv = umm / beta
      u(np) = beta
      goto 8
 4    xpiv = beta/umm
 8    permut(np) = perm
      ypiv(np) = xpiv
      if (.not. full) return
!     shift everything up if full...
      do 7 k=1,npm1
         ypiv(k) = ypiv(k+1)
         permut(k) = permut(k+1)
 7    continue
      return
!-----end-of-implu
      end
!-----------------------------------------------------------------------
      subroutine uppdir(n,p,np,lbp,indp,y,u,usav,flops)
      real*8 p(n,lbp), y(*), u(*), usav(*), x, flops
      integer k,np,n,npm1,j,ju,indp,lbp
!-----------------------------------------------------------------------
!     updates the conjugate directions p given the upper part of the
!     banded upper triangular matrix u.  u contains the non zero
!     elements of the column of the triangular matrix..
!-----------------------------------------------------------------------
      real*8 zero
      parameter(zero=0.0D0)
!
      npm1=np-1
      if (np .le. 1) goto 12
      j=indp
      ju = npm1
 10   if (j .le. 0) j=lbp
      x = u(ju) /usav(j)
      if (x .eq. zero) goto 115
      do 11 k=1,n
         y(k) = y(k) - x*p(k,j)
 11   continue
      flops = flops + 2*n
 115  j = j-1
      ju = ju -1
      if (ju .ge. 1) goto 10
 12   indp = indp + 1
      if (indp .gt. lbp) indp = 1
      usav(indp) = u(np)
      do 13 k=1,n
         p(k,indp) = y(k)
 13   continue
      return
!-----------------------------------------------------------------------
!-------end-of-uppdir---------------------------------------------------
      end
      subroutine givens(x,y,c,s)
      real*8 x,y,c,s
!-----------------------------------------------------------------------
!     Given x and y, this subroutine generates a Givens' rotation c, s.
!     And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
!     (See P 202 of "matrix computation" by Golub and van Loan.)
!-----------------------------------------------------------------------
      real*8 t,one,zero
      parameter (zero=0.0D0,one=1.0D0)
!
      if (x.eq.zero .and. y.eq.zero) then
         c = one
         s = zero
      else if (abs(y).gt.abs(x)) then
         t = x / y
         x = sqrt(one+t*t)
         s = sign(one / x, y)
         c = t*s
      else if (abs(y).le.abs(x)) then
         t = y / x
         y = sqrt(one+t*t)
         c = sign(one / y, x)
         s = t*c
      else
!
!     X or Y must be an invalid floating-point number, set both to zero
!
         x = zero
         y = zero
         c = one
         s = zero
      endif
      x = abs(x*y)
!
!     end of givens
!
      return
      end
!-----end-of-givens
!-----------------------------------------------------------------------
      logical function stopbis(n,ipar,mvpi,fpar,r,delx,sx)
      implicit none
      integer n,mvpi,ipar(16)
      real*8 fpar(16), r(n), delx(n), sx, distdot
      external distdot
!-----------------------------------------------------------------------
!     function for determining the stopping criteria. return value of
!     true if the stopbis criteria is satisfied.
!-----------------------------------------------------------------------
      if (ipar(11) .eq. 1) then
         stopbis = .true.
      else
         stopbis = .false.
      endif
      if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
         ipar(1) = -1
         stopbis = .true.
      endif
      if (stopbis) return
!
!     computes errors
!
      fpar(5) = sqrt(distdot(n,r,1,r,1))
      fpar(11) = fpar(11) + 2 * n
      if (ipar(3).lt.0) then
!
!     compute the change in the solution vector
!
         fpar(6) = sx * sqrt(distdot(n,delx,1,delx,1))
         fpar(11) = fpar(11) + 2 * n
         if (ipar(7).lt.mvpi+mvpi+1) then
!
!     if this is the end of the first iteration, set fpar(3:4)
!
            fpar(3) = fpar(6)
            if (ipar(3).eq.-1) then
               fpar(4) = fpar(1) * fpar(3) + fpar(2)
            endif
         endif
      else
         fpar(6) = fpar(5)
      endif
!
!     .. the test is struct this way so that when the value in fpar(6)
!       is not a valid number, STOPBIS is set to .true.
!
      if (fpar(6).gt.fpar(4)) then
         stopbis = .false.
         ipar(11) = 0
      else
         stopbis = .true.
         ipar(11) = 1
      endif
!
      return
      end
!-----end-of-stopbis
!-----------------------------------------------------------------------
      subroutine tidycg(n,ipar,fpar,sol,delx)
      implicit none
      integer i,n,ipar(16)
      real*8 fpar(16),sol(n),delx(n)
!-----------------------------------------------------------------------
!     Some common operations required before terminating the CG routines
!-----------------------------------------------------------------------
      real*8 zero
      parameter(zero=0.0D0)
!
      if (ipar(12).ne.0) then
         ipar(1) = ipar(12) 
      else if (ipar(1).gt.0) then
         if ((ipar(3).eq.999 .and. ipar(11).eq.1) .or. &
     &        fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
      if (fpar(3).gt.zero .and. fpar(6).gt.zero .and. &
     &     ipar(7).gt.ipar(13)) then
         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      else
         fpar(7) = zero
      endif
      do i = 1, n
         sol(i) = sol(i) + delx(i)
      enddo
      return
      end
!-----end-of-tidycg
!-----------------------------------------------------------------------
      logical function brkdn(alpha, ipar)
      implicit none
      integer ipar(16)
      real*8 alpha, beta, zero, one
      parameter (zero=0.0D0, one=1.0D0)
!-----------------------------------------------------------------------
!     test whether alpha is zero or an abnormal number, if yes,
!     this routine will return .true.
!
!     If alpha == 0, ipar(1) = -3,
!     if alpha is an abnormal number, ipar(1) = -9.
!-----------------------------------------------------------------------
      brkdn = .false.
      if (alpha.gt.zero) then
         beta = one / alpha
         if (.not. beta.gt.zero) then
            brkdn = .true.
            ipar(1) = -9
         endif
      else if (alpha.lt.zero) then
         beta = one / alpha
         if (.not. beta.lt.zero) then
            brkdn = .true.
            ipar(1) = -9
         endif
      else if (alpha.eq.zero) then
         brkdn = .true.
         ipar(1) = -3
      else
         brkdn = .true.
         ipar(1) = -9
      endif
      return
      end
!-----end-of-brkdn
!-----------------------------------------------------------------------
      subroutine bisinit(ipar,fpar,wksize,dsc,lp,rp,wk)
      implicit none
      integer i,ipar(16),wksize,dsc
      logical lp,rp
      real*8  fpar(16),wk(*)
!-----------------------------------------------------------------------
!     some common initializations for the iterative solvers
!-----------------------------------------------------------------------
      real*8 zero, one
      parameter(zero=0.0D0, one=1.0D0)
!
!     ipar(1) = -2 inidcate that there are not enough space in the work
!     array
!
      if (ipar(4).lt.wksize) then
         ipar(1) = -2
         ipar(4) = wksize
         return
      endif
!
      if (ipar(2).gt.2) then
         lp = .true.
         rp = .true.
      else if (ipar(2).eq.2) then
         lp = .false.
         rp = .true.
      else if (ipar(2).eq.1) then
         lp = .true.
         rp = .false.
      else
         lp = .false.
         rp = .false.
      endif
      if (ipar(3).eq.0) ipar(3) = dsc
!     .. clear the ipar elements used
      ipar(7) = 0
      ipar(8) = 0
      ipar(9) = 0
      ipar(10) = 0
      ipar(11) = 0
      ipar(12) = 0
      ipar(13) = 0
!
!     fpar(1) must be between (0, 1), fpar(2) must be positive,
!     fpar(1) and fpar(2) can NOT both be zero
!     Normally returns ipar(1) = -4 to indicate any of above error
!
      if (fpar(1).lt.zero .or. fpar(1).ge.one .or. fpar(2).lt.zero .or.     &
     &     (fpar(1).eq.zero .and. fpar(2).eq.zero)) then
         if (ipar(1).eq.0) then
            ipar(1) = -4
            return
         else
            fpar(1) = 1.0D-6
            fpar(2) = 1.0D-16
         endif
      endif
!     .. clear the fpar elements
      do i = 3, 10
         fpar(i) = zero
      enddo
      if (fpar(11).lt.zero) fpar(11) = zero
!     .. clear the used portion of the work array to zero
      do i = 1, wksize
         wk(i) = zero
      enddo
!
      return
!-----end-of-bisinit
      end
!-----------------------------------------------------------------------
      subroutine mgsro(full,lda,n,m,ind,ops,vec,hh,ierr)
      implicit none
      logical full
      integer lda,m,n,ind,ierr
      real*8  ops,hh(m),vec(lda,m)
!-----------------------------------------------------------------------
!     MGSRO  -- Modified Gram-Schmidt procedure with Selective Re-
!               Orthogonalization
!     The ind'th vector of VEC is orthogonalized against the rest of
!     the vectors.
!
!     The test for performing re-orthogonalization is performed for
!     each indivadual vectors. If the cosine between the two vectors
!     is greater than 0.99 (REORTH = 0.99**2), re-orthogonalization is
!     performed. The norm of the 'new' vector is kept in variable NRM0,
!     and updated after operating with each vector.
!
!     full   -- .ture. if it is necessary to orthogonalize the ind'th
!               against all the vectors vec(:,1:ind-1), vec(:,ind+2:m)
!               .false. only orthogonalize againt vec(:,1:ind-1)
!     lda    -- the leading dimension of VEC
!     n      -- length of the vector in VEC
!     m      -- number of vectors can be stored in VEC
!     ind    -- index to the vector to be changed
!     ops    -- operation counts
!     vec    -- vector of LDA X M storing the vectors
!     hh     -- coefficient of the orthogonalization
!     ierr   -- error code
!               0 : successful return
!               -1: zero input vector
!               -2: input vector contains abnormal numbers
!               -3: input vector is a linear combination of others
!
!     External routines used: real*8 distdot
!-----------------------------------------------------------------------
      integer i,k
      real*8  nrm0, nrm1, fct, thr, distdot, zero, one, reorth
      parameter (zero=0.0D0, one=1.0D0, reorth=0.98D0)
      external distdot
!
!     compute the norm of the input vector
!
      nrm0 = distdot(n,vec(1,ind),1,vec(1,ind),1)
      ops = ops + n + n
      thr = nrm0 * reorth
      if (nrm0.le.zero) then
         ierr = - 1
         return
      else if (nrm0.gt.zero .and. one/nrm0.gt.zero) then
         ierr = 0
      else
         ierr = -2
         return
      endif
!
!     Modified Gram-Schmidt loop
!
      if (full) then
         do 40 i = ind+1, m
            fct = distdot(n,vec(1,ind),1,vec(1,i),1)
            hh(i) = fct
            do 20 k = 1, n
               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 20         continue
            ops = ops + 4 * n + 2
            if (fct*fct.gt.thr) then
               fct = distdot(n,vec(1,ind),1,vec(1,i),1)
               hh(i) = hh(i) + fct
               do 30 k = 1, n
                  vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 30            continue
               ops = ops + 4*n + 1
            endif
            nrm0 = nrm0 - hh(i) * hh(i)
            if (nrm0.lt.zero) nrm0 = zero
            thr = nrm0 * reorth
 40      continue
      endif
!
      do 70 i = 1, ind-1
         fct = distdot(n,vec(1,ind),1,vec(1,i),1)
         hh(i) = fct
         do 50 k = 1, n
            vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 50      continue
         ops = ops + 4 * n + 2
         if (fct*fct.gt.thr) then
            fct = distdot(n,vec(1,ind),1,vec(1,i),1)
            hh(i) = hh(i) + fct
            do 60 k = 1, n
               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
 60         continue
            ops = ops + 4*n + 1
         endif
         nrm0 = nrm0 - hh(i) * hh(i)
         if (nrm0.lt.zero) nrm0 = zero
         thr = nrm0 * reorth
 70   continue
!
!     test the resulting vector
!
      nrm1 = sqrt(distdot(n,vec(1,ind),1,vec(1,ind),1))
      ops = ops + n + n
      hh(ind) = nrm1
      if (nrm1.le.zero) then
         ierr = -3
         return
      endif
!
!     scale the resulting vector
!
      fct = one / nrm1
      do 80 k = 1, n
         vec(k,ind) = vec(k,ind) * fct
 80   continue
      ops = ops + n + 1
!
!     normal return
!
      ierr = 0
      return
!     end surbotine mgsro
      end
