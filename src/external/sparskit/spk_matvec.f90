
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
!    The original file is called BLASSM/matvec.f
!
!--------------------------------------------------------------------------

! revision log :
!
! 29.10.2023    ggu     imported into SHYFEM

!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c
!         Matrix-vector Mulitiplications and Triang. Solves            c
!----------------------------------------------------------------------c
! contents: (as of Nov 18, 1991)                                       c
!----------                                                            c
! 1) Matrix-vector products:                                           c
!---------------------------                                           c
! amux  : A times a vector. Compressed Sparse Row (CSR) format.        c
! amuxms: A times a vector. Modified Compress Sparse Row format.       c
! atmux : Transp(A) times a vector. CSR format.                        c
! atmuxr: Transp(A) times a vector. CSR format. A rectangular.         c
! amuxe : A times a vector. Ellpack/Itpack (ELL) format.               c
! amuxd : A times a vector. Diagonal (DIA) format.                     c
! amuxj : A times a vector. Jagged Diagonal (JAD) format.              c
! vbrmv : Sparse matrix-full vector product, in VBR format             c
!                                                                      c
! 2) Triangular system solutions:                                      c
!-------------------------------                                       c
! lsol  : Unit Lower Triang. solve. Compressed Sparse Row (CSR) format.c
! ldsol : Lower Triang. solve.  Modified Sparse Row (MSR) format.      c
! lsolc : Unit Lower Triang. solve. Comp. Sparse Column (CSC) format.  c
! ldsolc: Lower Triang. solve. Modified Sparse Column (MSC) format.    c
! ldsoll: Lower Triang. solve with level scheduling. MSR format.       c
! usol  : Unit Upper Triang. solve. Compressed Sparse Row (CSR) format.c
! udsol : Upper Triang. solve.  Modified Sparse Row (MSR) format.      c
! usolc : Unit Upper Triang. solve. Comp. Sparse Column (CSC) format.  c
! udsolc: Upper Triang. solve.  Modified Sparse Column (MSC) format.   c
!----------------------------------------------------------------------c
! 1)     M A T R I X    B Y    V E C T O R     P R O D U C T S         c
!----------------------------------------------------------------------c
      subroutine amux (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
!-----------------------------------------------------------------------
!         A times a vector
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in compressed sparse row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
      real*8 t
      integer i, k
!-----------------------------------------------------------------------
      do 100 i = 1,n
!
!     compute the inner product of row i with vector x
! 
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
!
!     store result in y(i) 
!
         y(i) = t
 100  continue
!
      return
!---------end-of-amux---------------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine amuxms (n, x, y, a,ja)
      real*8  x(*), y(*), a(*)
      integer n, ja(*)
!-----------------------------------------------------------------------
!         A times a vector in MSR format
!-----------------------------------------------------------------------
! multiplies a matrix by a vector using the dot product form
! Matrix A is stored in Modified Sparse Row storage.
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,= input matrix in modified compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=Ax
!
!-----------------------------------------------------------------------
! local variables
!
      integer i, k
!-----------------------------------------------------------------------
        do 10 i=1, n
        y(i) = a(i)*x(i)
 10     continue
      do 100 i = 1,n
!
!     compute the inner product of row i with vector x
!
         do 99 k=ja(i), ja(i+1)-1
            y(i) = y(i) + a(k) *x(ja(k))
 99      continue
 100  continue
!
      return
!---------end-of-amuxm--------------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine atmux (n, x, y, a, ja, ia)
      real*8 x(*), y(*), a(*) 
      integer n, ia(*), ja(*)
!-----------------------------------------------------------------------
!         transp( A ) times a vector
!----------------------------------------------------------------------- 
! multiplies the transpose of a matrix by a vector when the original
! matrix is stored in compressed sparse row storage. Can also be
! viewed as the product of a matrix by a vector when the original
! matrix is stored in the compressed sparse column format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=transp(A)*x
!
!-----------------------------------------------------------------------
!     local variables 
!
      integer i, k 
!-----------------------------------------------------------------------
!
!     zero out output vector
! 
      do 1 i=1,n
         y(i) = 0.0
 1    continue
!
! loop over the rows
!
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
!
      return
!-------------end-of-atmux---------------------------------------------- 
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine atmuxr (m, n, x, y, a, ja, ia)
      real*8 x(*), y(*), a(*) 
      integer m, n, ia(*), ja(*)
!-----------------------------------------------------------------------
!         transp( A ) times a vector, A can be rectangular
!----------------------------------------------------------------------- 
! See also atmux.  The essential difference is how the solution vector
! is initially zeroed.  If using this to multiply rectangular CSC 
! matrices by a vector, m number of rows, n is number of columns.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! m     = column dimension of A
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! a, ja,
!    ia = input matrix in compressed sparse row format.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=transp(A)*x
!
!-----------------------------------------------------------------------
!     local variables 
!
      integer i, k 
!-----------------------------------------------------------------------
!
!     zero out output vector
! 
      do 1 i=1,m
         y(i) = 0.0
 1    continue
!
! loop over the rows
!
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1 
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue
!
      return
!-------------end-of-atmuxr--------------------------------------------- 
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine amuxe (n,x,y,na,ncol,a,ja) 
      real*8 x(n), y(n), a(na,*)  
      integer  n, na, ncol, ja(na,*)
!-----------------------------------------------------------------------
!        A times a vector in Ellpack Itpack format (ELL)               
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector when the original matrix is stored 
! in the ellpack-itpack sparse format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! na    = integer. The first dimension of arrays a and ja
!         as declared by the calling program.
! ncol  = integer. The number of active columns in array a.
!         (i.e., the number of generalized diagonals in matrix.)
! a, ja = the real and integer arrays of the itpack format
!         (a(i,k),k=1,ncol contains the elements of row i in matrix
!          ja(i,k),k=1,ncol contains their column numbers) 
!
! on return:
!-----------
! y     = real array of length n, containing the product y=y=A*x
!
!-----------------------------------------------------------------------
! local variables
!
      integer i, j 
!-----------------------------------------------------------------------
      do 1 i=1, n
         y(i) = 0.0 
 1    continue
      do 10 j=1,ncol
         do 25 i = 1,n
            y(i) = y(i)+a(i,j)*x(ja(i,j))
 25      continue
 10   continue
!
      return
!--------end-of-amuxe--------------------------------------------------- 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine amuxd (n,x,y,diag,ndiag,idiag,ioff) 
      integer n, ndiag, idiag, ioff(idiag) 
      real*8 x(n), y(n), diag(ndiag,idiag)
!-----------------------------------------------------------------------
!        A times a vector in Diagonal storage format (DIA) 
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector when the original matrix is stored 
! in the diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n     = row dimension of A
! x     = real array of length equal to the column dimension of
!         the A matrix.
! ndiag  = integer. The first dimension of array adiag as declared in
!         the calling program.
! idiag  = integer. The number of diagonals in the matrix.
! diag   = real array containing the diagonals stored of A.
! idiag  = number of diagonals in matrix.
! diag   = real array of size (ndiag x idiag) containing the diagonals
!          
! ioff   = integer array of length idiag, containing the offsets of the
!   	   diagonals of the matrix:
!          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
!
! on return:
!-----------
! y     = real array of length n, containing the product y=A*x
!
!-----------------------------------------------------------------------
! local variables 
!
      integer j, k, io, i1, i2 
!-----------------------------------------------------------------------
      do 1 j=1, n
         y(j) = 0.0d0
 1    continue
      do 10 j=1, idiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)
         do 9 k=i1, i2
            y(k) = y(k)+diag(k,j)*x(k+io)
 9       continue
 10   continue
! 
      return
!----------end-of-amuxd-------------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine amuxj (n, x, y, jdiag, a, ja, ia)
      integer n, jdiag, ja(*), ia(*)
      real*8 x(n), y(n), a(*)  
!-----------------------------------------------------------------------
!        A times a vector in Jagged-Diagonal storage format (JAD) 
!----------------------------------------------------------------------- 
! multiplies a matrix by a vector when the original matrix is stored 
! in the jagged diagonal storage format.
!-----------------------------------------------------------------------
!
! on entry:
!----------
! n      = row dimension of A
! x      = real array of length equal to the column dimension of
!         the A matrix.
! jdiag  = integer. The number of jadded-diagonals in the data-structure.
! a      = real array containing the jadded diagonals of A stored
!          in succession (in decreasing lengths) 
! j      = integer array containing the colum indices of the 
!          corresponding elements in a.
! ia     = integer array containing the lengths of the  jagged diagonals
!
! on return:
!-----------
! y      = real array of length n, containing the product y=A*x
!
! Note:
!------- 
! Permutation related to the JAD format is not performed.
! this can be done by:
!     call permvec (n,y,y,iperm) 
! after the call to amuxj, where iperm is the permutation produced
! by csrjad.
!-----------------------------------------------------------------------
! local variables 
!
      integer i, ii, k1, len, j 
!-----------------------------------------------------------------------
      do 1 i=1, n
         y(i) = 0.0d0
 1    continue
      do 70 ii=1, jdiag
         k1 = ia(ii)-1
         len = ia(ii+1)-k1-1
         do 60 j=1,len
            y(j)= y(j)+a(k1+j)*x(ja(k1+j)) 
 60      continue
 70   continue
!
      return
!----------end-of-amuxj------------------------------------------------- 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine vbrmv(nr, nc, ia, ja, a, kvstr, kvstc, x, b)
!-----------------------------------------------------------------------
      integer nr, nc, ia(nr+1), ja(*), kvstr(nr+1), kvstc(*)
      real*8  a(*), x(*), b(*)
!-----------------------------------------------------------------------
!     Sparse matrix-full vector product, in VBR format.
!-----------------------------------------------------------------------
!     On entry:
!--------------
!     nr, nc  = number of block rows and columns in matrix A
!     ia,ja,(),a,kvstr,kvstc = matrix A in variable block row format
!     x       = multiplier vector in full format
!
!     On return:
!---------------
!     b = product of matrix A times vector x in full format
!
!     Algorithm:
!---------------
!     Perform multiplication by traversing a in order.
!
!-----------------------------------------------------------------------
!-----local variables
      integer n, i, j, ii, jj, k, istart, istop
      real*8  xjj
!---------------------------------
      n = kvstc(nc+1)-1
      do i = 1, n
         b(i) = 0.d0
      enddo
!---------------------------------
      k = 1
      do i = 1, nr
         istart = kvstr(i)
         istop  = kvstr(i+1)-1
         do j = ia(i), ia(i+1)-1
            do jj = kvstc(ja(j)), kvstc(ja(j)+1)-1
               xjj = x(jj)
               do ii = istart, istop
                  b(ii) = b(ii) + xjj*a(k)
                  k = k + 1
               enddo
            enddo
         enddo
      enddo
!---------------------------------
      return
      end
!-----------------------------------------------------------------------
!----------------------end-of-vbrmv-------------------------------------
!-----------------------------------------------------------------------
!----------------------------------------------------------------------c
! 2)     T R I A N G U L A R    S Y S T E M    S O L U T I O N S       c
!----------------------------------------------------------------------c
      subroutine lsol (n,x,y,al,jal,ial)
      integer n, jal(*),ial(n+1) 
      real*8  x(n), y(n), al(*) 
!-----------------------------------------------------------------------
!   solves    L x = y ; L = lower unit triang. /  CSR format
!----------------------------------------------------------------------- 
! solves a unit lower triangular system by standard (sequential )
! forward elimination - matrix stored in CSR format. 
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real array containg the right side.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in compressed sparse row
!          format. 
!
! On return:
!----------- 
!	x  = The solution of  L x  = y.
!--------------------------------------------------------------------
! local variables 
!
      integer k, j 
      real*8  t
!-----------------------------------------------------------------------
      x(1) = y(1) 
      do 150 k = 2, n
         t = y(k) 
         do 100 j = ial(k), ial(k+1)-1
            t = t-al(j)*x(jal(j))
 100     continue
         x(k) = t 
 150  continue
!
      return
!----------end-of-lsol-------------------------------------------------- 
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine ldsol (n,x,y,al,jal) 
      integer n, jal(*) 
      real*8 x(n), y(n), al(*) 
!----------------------------------------------------------------------- 
!     Solves L x = y    L = triangular. MSR format 
!-----------------------------------------------------------------------
! solves a (non-unit) lower triangular system by standard (sequential) 
! forward elimination - matrix stored in MSR format 
! with diagonal elements already inverted (otherwise do inversion,
! al(1:n) = 1.0/al(1:n),  before calling ldsol).
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real array containg the right hand side.
!
! al,
! jal,   = Lower triangular matrix stored in Modified Sparse Row 
!          format. 
!
! On return:
!----------- 
!	x = The solution of  L x = y .
!--------------------------------------------------------------------
! local variables 
!
      integer k, j 
      real*8 t 
!-----------------------------------------------------------------------
      x(1) = y(1)*al(1) 
      do 150 k = 2, n
         t = y(k) 
         do 100 j = jal(k), jal(k+1)-1
            t = t - al(j)*x(jal(j))
 100     continue
         x(k) = al(k)*t 
 150  continue
      return
!----------end-of-ldsol-------------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine lsolc (n,x,y,al,jal,ial)
      integer n, jal(*),ial(*) 
      real*8  x(n), y(n), al(*) 
!-----------------------------------------------------------------------
!       SOLVES     L x = y ;    where L = unit lower trang. CSC format
!-----------------------------------------------------------------------
! solves a unit lower triangular system by standard (sequential )
! forward elimination - matrix stored in CSC format. 
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real*8 array containg the right side.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in compressed sparse column 
!          format. 
!
! On return:
!----------- 
!	x  = The solution of  L x  = y.
!-----------------------------------------------------------------------
! local variables 
!
      integer k, j
      real*8 t
!-----------------------------------------------------------------------
      do 140 k=1,n
         x(k) = y(k) 
 140  continue
      do 150 k = 1, n-1
         t = x(k) 
         do 100 j = ial(k), ial(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j) 
 100     continue
 150  continue
!
      return
!----------end-of-lsolc------------------------------------------------- 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine ldsolc (n,x,y,al,jal) 
      integer n, jal(*)
      real*8 x(n), y(n), al(*)
!-----------------------------------------------------------------------
!    Solves     L x = y ;    L = nonunit Low. Triang. MSC format 
!----------------------------------------------------------------------- 
! solves a (non-unit) lower triangular system by standard (sequential) 
! forward elimination - matrix stored in Modified Sparse Column format 
! with diagonal elements already inverted (otherwise do inversion,
! al(1:n) = 1.0/al(1:n),  before calling ldsol).
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real array containg the right hand side.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in Modified Sparse Column
!           format.
!
! On return:
!----------- 
!	x = The solution of  L x = y .
!--------------------------------------------------------------------
! local variables
!
      integer k, j
      real*8 t 
!-----------------------------------------------------------------------
      do 140 k=1,n
         x(k) = y(k) 
 140  continue
      do 150 k = 1, n 
         x(k) = x(k)*al(k) 
         t = x(k) 
         do 100 j = jal(k), jal(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j) 
 100     continue
 150  continue
!
      return
!----------end-of-lsolc------------------------------------------------ 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------  
      subroutine ldsoll (n,x,y,al,jal,nlev,lev,ilev) 
      integer n, nlev, jal(*), ilev(nlev+1), lev(n)
      real*8 x(n), y(n), al(*)
!-----------------------------------------------------------------------
!    Solves L x = y    L = triangular. Uses LEVEL SCHEDULING/MSR format 
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real array containg the right hand side.
!
! al,
! jal,   = Lower triangular matrix stored in Modified Sparse Row 
!          format. 
! nlev   = number of levels in matrix
! lev    = integer array of length n, containing the permutation
!          that defines the levels in the level scheduling ordering.
! ilev   = pointer to beginning of levels in lev.
!          the numbers lev(i) to lev(i+1)-1 contain the row numbers
!          that belong to level number i, in the level shcheduling
!          ordering.
!
! On return:
!----------- 
!	x = The solution of  L x = y .
!--------------------------------------------------------------------
      integer ii, jrow, i 
      real*8 t 
!     
!     outer loop goes through the levels. (SEQUENTIAL loop)
!     
      do 150 ii=1, nlev
!     
!     next loop executes within the same level. PARALLEL loop
!     
         do 100 i=ilev(ii), ilev(ii+1)-1 
            jrow = lev(i)
!
! compute inner product of row jrow with x
! 
            t = y(jrow) 
            do 130 k=jal(jrow), jal(jrow+1)-1 
               t = t - al(k)*x(jal(k))
 130        continue
            x(jrow) = t*al(jrow) 
 100     continue
 150  continue
      return
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine usol (n,x,y,au,jau,iau)
      integer n, jau(*),iau(n+1) 
      real*8  x(n), y(n), au(*) 
!----------------------------------------------------------------------- 
!             Solves   U x = y    U = unit upper triangular. 
!-----------------------------------------------------------------------
! solves a unit upper triangular system by standard (sequential )
! backward elimination - matrix stored in CSR format. 
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real array containg the right side.
!
! au,
! jau,
! iau,    = Lower triangular matrix stored in compressed sparse row
!          format. 
!
! On return:
!----------- 
!	x = The solution of  U x = y . 
!-------------------------------------------------------------------- 
! local variables 
!
      integer k, j 
      real*8  t
!-----------------------------------------------------------------------
      x(n) = y(n) 
      do 150 k = n-1,1,-1 
         t = y(k) 
         do 100 j = iau(k), iau(k+1)-1
            t = t - au(j)*x(jau(j))
 100     continue
         x(k) = t 
 150  continue
!
      return
!----------end-of-usol-------------------------------------------------- 
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine udsol (n,x,y,au,jau) 
      integer n, jau(*) 
      real*8  x(n), y(n),au(*) 
!----------------------------------------------------------------------- 
!             Solves   U x = y  ;   U = upper triangular in MSR format
!-----------------------------------------------------------------------
! solves a non-unit upper triangular matrix by standard (sequential )
! backward elimination - matrix stored in MSR format. 
! with diagonal elements already inverted (otherwise do inversion,
! au(1:n) = 1.0/au(1:n),  before calling).
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real array containg the right side.
!
! au,
! jau,    = Lower triangular matrix stored in modified sparse row
!          format. 
!
! On return:
!----------- 
!	x = The solution of  U x = y .
!--------------------------------------------------------------------
! local variables 
!
      integer k, j
      real*8 t
!-----------------------------------------------------------------------
      x(n) = y(n)*au(n)
      do 150 k = n-1,1,-1
         t = y(k) 
         do 100 j = jau(k), jau(k+1)-1
            t = t - au(j)*x(jau(j))
 100     continue
         x(k) = au(k)*t 
 150  continue
!
      return
!----------end-of-udsol-------------------------------------------------
!-----------------------------------------------------------------------
      end
!----------------------------------------------------------------------- 
      subroutine usolc (n,x,y,au,jau,iau)
      real*8  x(*), y(*), au(*) 
      integer n, jau(*),iau(*)
!-----------------------------------------------------------------------
!       SOUVES     U x = y ;    where U = unit upper trang. CSC format
!-----------------------------------------------------------------------
! solves a unit upper triangular system by standard (sequential )
! forward elimination - matrix stored in CSC format. 
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real*8 array containg the right side.
!
! au,
! jau,
! iau,    = Uower triangular matrix stored in compressed sparse column 
!          format. 
!
! On return:
!----------- 
!	x  = The solution of  U x  = y.
!-----------------------------------------------------------------------
! local variables 
!     
      integer k, j
      real*8 t
!-----------------------------------------------------------------------
      do 140 k=1,n
         x(k) = y(k) 
 140  continue
      do 150 k = n,1,-1
         t = x(k) 
         do 100 j = iau(k), iau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j) 
 100     continue
 150  continue
!
      return
!----------end-of-usolc------------------------------------------------- 
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine udsolc (n,x,y,au,jau)   
      integer n, jau(*) 
      real*8 x(n), y(n), au(*)  
!-----------------------------------------------------------------------
!    Solves     U x = y ;    U = nonunit Up. Triang. MSC format 
!----------------------------------------------------------------------- 
! solves a (non-unit) upper triangular system by standard (sequential) 
! forward elimination - matrix stored in Modified Sparse Column format 
! with diagonal elements already inverted (otherwise do inversion,
! auuuul(1:n) = 1.0/au(1:n),  before calling ldsol).
!-----------------------------------------------------------------------
!
! On entry:
!---------- 
! n      = integer. dimension of problem.
! y      = real*8 array containg the right hand side.
!
! au,
! jau,   = Upper triangular matrix stored in Modified Sparse Column
!          format.
!
! On return:
!----------- 
!	x = The solution of  U x = y .
!--------------------------------------------------------------------
! local variables 
! 
      integer k, j
      real*8 t
!----------------------------------------------------------------------- 
      do 140 k=1,n
         x(k) = y(k) 
 140  continue
      do 150 k = n,1,-1
         x(k) = x(k)*au(k) 
         t = x(k) 
         do 100 j = jau(k), jau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j) 
 100     continue
 150  continue
!
      return
!----------end-of-udsolc------------------------------------------------ 
!-----------------------------------------------------------------------
      end
