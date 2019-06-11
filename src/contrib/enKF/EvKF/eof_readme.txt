
*********************************************
This is the quickstart to doing EOF analysis.
*********************************************

Put your data into a matrix so that the rows indicate temporal development and
the columns are variables or spatial data points.  The temporal relationship
between rows is unimportant (ie. doesnt have to be uniform).  Same for the
spatial relationship between columns.
Detrend the columns of the resulting matrix.  Some EOF routines do this for
you, but I prefer to do it separately.
Use singular value decomposition (svd) to break up your data into 3 matrices:
Z = U * D * Vt  
where U and V are orthonormal and D is diagonal.  Then,
EOFs = V
ECs = U * D
covariance matrix = ECst * ECs / (n-1) = D2 / (n-1)
communalities matrix = ECs * ECst 
That is really all there is to it.  The EOFs are really the columns of the
EOFs matrix.  I have included matlab code that performs step 3 above.  See the
references for a more detailed discussion.

After finishing these calculations, you will probably want to reduce the EOFs
and ECs to only those which explain a significant percentage of the overall
variance by just selecting out those columns of the ECs and EOFs.  You then
may or may not wish to rotate the EOFs to increase the physical explainability
of the resulting patterns.  Finally, there are a number of useful ways to
visualize the results of your analysis.  I will not discuss visualization
here.  I will also not discuss EOF analysis of several fields.


*********************************************
DGESVD computes the singular value decomposition (SVD) of a real
*********************************************

 M-by-N matrix A, optionally computing the left and/or right singular
 vectors. The SVD is written

      A = U * SIGMA * transpose(V)

 where SIGMA is an M-by-N matrix which is zero except for its
 min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 are the singular values of A; they are real and non-negative, and
 are returned in descending order.  The first min(m,n) columns of
 U and V are the left and right singular vectors of A.

 Note that the routine returns V**T, not V.

Parameters
[in]	JOBU	
          JOBU is CHARACTER*1
          Specifies options for computing all or part of the matrix U:
          = 'A':  all M columns of U are returned in array U:
          = 'S':  the first min(m,n) columns of U (the left singular
                  vectors) are returned in the array U;
          = 'O':  the first min(m,n) columns of U (the left singular
                  vectors) are overwritten on the array A;
          = 'N':  no columns of U (no left singular vectors) are
                  computed.
[in]	JOBVT	
          JOBVT is CHARACTER*1
          Specifies options for computing all or part of the matrix
          V**T:
          = 'A':  all N rows of V**T are returned in the array VT;
          = 'S':  the first min(m,n) rows of V**T (the right singular
                  vectors) are returned in the array VT;
          = 'O':  the first min(m,n) rows of V**T (the right singular
                  vectors) are overwritten on the array A;
          = 'N':  no rows of V**T (no right singular vectors) are
                  computed.

          JOBVT and JOBU cannot both be 'O'.
[in]	M	
          M is INTEGER
          The number of rows of the input matrix A.  M >= 0.
[in]	N	
          N is INTEGER
          The number of columns of the input matrix A.  N >= 0.
[in,out]	A	
          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the M-by-N matrix A.
          On exit,
          if JOBU = 'O',  A is overwritten with the first min(m,n)
                          columns of U (the left singular vectors,
                          stored columnwise);
          if JOBVT = 'O', A is overwritten with the first min(m,n)
                          rows of V**T (the right singular vectors,
                          stored rowwise);
          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
                          are destroyed.
[in]	LDA	
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
[out]	S	
          S is DOUBLE PRECISION array, dimension (min(M,N))
          The singular values of A, sorted so that S(i) >= S(i+1).
[out]	U	
          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
          if JOBU = 'S', U contains the first min(m,n) columns of U
          (the left singular vectors, stored columnwise);
          if JOBU = 'N' or 'O', U is not referenced.
[in]	LDU	
          LDU is INTEGER
          The leading dimension of the array U.  LDU >= 1; if
          JOBU = 'S' or 'A', LDU >= M.
[out]	VT	
          VT is DOUBLE PRECISION array, dimension (LDVT,N)
          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
          V**T;
          if JOBVT = 'S', VT contains the first min(m,n) rows of
          V**T (the right singular vectors, stored rowwise);
          if JOBVT = 'N' or 'O', VT is not referenced.
[in]	LDVT	
          LDVT is INTEGER
          The leading dimension of the array VT.  LDVT >= 1; if
          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
[out]	WORK	
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
          superdiagonal elements of an upper bidiagonal matrix B
          whose diagonal is in S (not necessarily sorted). B
          satisfies A = U * B * VT, so it has the same singular values
          as A, and singular vectors related by U and VT.
[in]	LWORK	
          LWORK is INTEGER
          The dimension of the array WORK.
          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
             - PATH 1  (M much larger than N, JOBU='N')
             - PATH 1t (N much larger than M, JOBVT='N')
          LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths
          For good performance, LWORK should generally be larger.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]	INFO	
          INFO is INTEGER
          = 0:  successful exit.
          < 0:  if INFO = -i, the i-th argument had an illegal value.
          > 0:  if DBDSQR did not converge, INFO specifies how many
                superdiagonals of an intermediate bidiagonal form B
                did not converge to zero. See the description of WORK
                above for details.

*********************************************
DGEMM  performs one of the matrix-matrix operations
*********************************************

    C := alpha*op( A )*op( B ) + beta*C,

 where  op( X ) is one of

    op( X ) = X   or   op( X ) = X**T,

 alpha and beta are scalars, and A, B and C are matrices, with op( A )
 an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
Parameters
[in]	TRANSA	
          TRANSA is CHARACTER*1
           On entry, TRANSA specifies the form of op( A ) to be used in
           the matrix multiplication as follows:

              TRANSA = 'N' or 'n',  op( A ) = A.

              TRANSA = 'T' or 't',  op( A ) = A**T.

              TRANSA = 'C' or 'c',  op( A ) = A**T.
[in]	TRANSB	
          TRANSB is CHARACTER*1
           On entry, TRANSB specifies the form of op( B ) to be used in
           the matrix multiplication as follows:

              TRANSB = 'N' or 'n',  op( B ) = B.

              TRANSB = 'T' or 't',  op( B ) = B**T.

              TRANSB = 'C' or 'c',  op( B ) = B**T.
[in]	M	
          M is INTEGER
           On entry,  M  specifies  the number  of rows  of the  matrix
           op( A )  and of the  matrix  C.  M  must  be at least  zero.
[in]	N	
          N is INTEGER
           On entry,  N  specifies the number  of columns of the matrix
           op( B ) and the number of columns of the matrix C. N must be
           at least zero.
[in]	K	
          K is INTEGER
           On entry,  K  specifies  the number of columns of the matrix
           op( A ) and the number of rows of the matrix op( B ). K must
           be at least  zero.
[in]	ALPHA	
          ALPHA is DOUBLE PRECISION.
           On entry, ALPHA specifies the scalar alpha.
[in]	A	
          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is
           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
           part of the array  A  must contain the matrix  A,  otherwise
           the leading  k by m  part of the array  A  must contain  the
           matrix A.
[in]	LDA	
          LDA is INTEGER
           On entry, LDA specifies the first dimension of A as declared
           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
           LDA must be at least  max( 1, m ), otherwise  LDA must be at
           least  max( 1, k ).
[in]	B	
          B is DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is
           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
           part of the array  B  must contain the matrix  B,  otherwise
           the leading  n by k  part of the array  B  must contain  the
           matrix B.
[in]	LDB	
          LDB is INTEGER
           On entry, LDB specifies the first dimension of B as declared
           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
           LDB must be at least  max( 1, k ), otherwise  LDB must be at
           least  max( 1, n ).
[in]	BETA	
          BETA is DOUBLE PRECISION.
           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
           supplied as zero then C need not be set on input.
[in,out]	C	
          C is DOUBLE PRECISION array, dimension ( LDC, N )
           Before entry, the leading  m by n  part of the array  C must
           contain the matrix  C,  except when  beta  is zero, in which
           case C need not be set on entry.
           On exit, the array  C  is overwritten by the  m by n  matrix
           ( alpha*op( A )*op( B ) + beta*C ).
[in]	LDC	
          LDC is INTEGER
           On entry, LDC specifies the first dimension of C as declared
           in  the  calling  (sub)  program.   LDC  must  be  at  least
           max( 1, m ).




*********************************************
DGEMV  performs one of the matrix-vector operations
*********************************************

    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,

 where alpha and beta are scalars, x and y are vectors and A is an
 m by n matrix.
Parameters
[in]	TRANS	
          TRANS is CHARACTER*1
           On entry, TRANS specifies the operation to be performed as
           follows:

              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.

              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.

              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
[in]	M	
          M is INTEGER
           On entry, M specifies the number of rows of the matrix A.
           M must be at least zero.
[in]	N	
          N is INTEGER
           On entry, N specifies the number of columns of the matrix A.
           N must be at least zero.
[in]	ALPHA	
          ALPHA is DOUBLE PRECISION.
           On entry, ALPHA specifies the scalar alpha.
[in]	A	
          A is DOUBLE PRECISION array, dimension ( LDA, N )
           Before entry, the leading m by n part of the array A must
           contain the matrix of coefficients.
[in]	LDA	
          LDA is INTEGER
           On entry, LDA specifies the first dimension of A as declared
           in the calling (sub) program. LDA must be at least
           max( 1, m ).
[in]	X	
          X is DOUBLE PRECISION array, dimension at least
           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
           and at least
           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
           Before entry, the incremented array X must contain the
           vector x.
[in]	INCX	
          INCX is INTEGER
           On entry, INCX specifies the increment for the elements of
           X. INCX must not be zero.
[in]	BETA	
          BETA is DOUBLE PRECISION.
           On entry, BETA specifies the scalar beta. When BETA is
           supplied as zero then Y need not be set on input.
[in,out]	Y	
          Y is DOUBLE PRECISION array, dimension at least
           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
           and at least
           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
           Before entry with BETA non-zero, the incremented array Y
           must contain the vector y. On exit, Y is overwritten by the
           updated vector y.
[in]	INCY	
          INCY is INTEGER
           On entry, INCY specifies the increment for the elements of
           Y. INCY must not be zero.


