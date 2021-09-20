      PROGRAM gradient
	    IMPLICIT NONE
      DOUBLE PRECISION A,b, x,x0,r0, eps, d0, alpha, res, beta
	    INTEGER i, j, k, n, kmax
      DIMENSION A(1000,1000),b(1000),x(1000),x0(1000),r0(1000)
      DIMENSION d0(1000)

c	    Variables temporaires de calcul
	    DOUBLE PRECISION Ax0, Ad0, rkScal, rk1Scal, vecTemp, scal
      DIMENSION Ax0(1000), Ad0(1000), vecTemp(1000)

c     taille de la matrice
	    n=3
      eps = 0.1

c     INIT MATRICE
	    A(1,1) = 1
  	  A(2,2) = 2
  	  A(3,3) = 3
  	  A(1,2) = eps
  	  A(1,3) = eps
  	  A(2,1) = eps
  	  A(3,1) = eps
  	  A(2,3) = eps
  	  A(3,2) = eps

c     Solution exacte
  	  x(1) = 1
  	  x(2) = 1
  	  x(3) = 1

c	    b=A*x
  	  call produit(A, x, n, b)

c     Valeur aléatoire de départ
  	  x0(1) = 200
  	  x0(2) = 4
  	  x0(3) = 6


c	    produit A*x0 pour définir r0
  	  call produit(A, x0, n, Ax0)

c	    Definition de r0
  	  DO i=1,3
  	   r0(i) = b(i)-Ax0(i)
  	   d0(i) = r0(i)
  	  ENDDO


  	  res = 1000
  	  eps = 10e-5
      kmax = 1000
      i = 0
  	  DO WHILE ((res.gt.eps).AND.(i.le.kmax))
c       Calcul de alpha et de la norme du résidu
  	    call scalaire(r0, r0, n, rkScal)
        res = sqrt(rkScal)
  	    call produit(A, d0, n, vecTemp)
  	    call scalaire(vecTemp, d0, n, scal)
  	    alpha = rkScal/scal

c       Calcul de x_k+1
  	    DO j = 1,n
  	     x0(j) = x0(j) + alpha*d0(j)
        ENDDO

c       Calcul de r_k+1
  	    call produit(A, d0, n, Ad0)
        DO j = 1,n
         r0(j) = r0(j)-alpha*Ad0(j)
        ENDDO

c       Calcul de beta
        call scalaire(r0, r0, n, rk1Scal)
  	    beta = rk1Scal/rkScal

c       Calcul de la direction d_k+1
  	    DO j = 1,n
  	     d0(j) = r0(j) - beta*d0(j)
        ENDDO
        i=i+1
  	  ENDDO
      WRITE(*,*) "Solution obtenue en", i, "itérations"
      WRITE(*,*) x0(1:3)

      END

      SUBROUTINE conj_grad(A, b, n, x0)
      IMPLICIT NONE
      DOUBLE PRECISION A,b, x,x0,r0, eps, d0, alpha, res, beta
      INTEGER i, j, k, n, kmax
      DIMENSION A(1000,1000),b(1000),x(1000),x0(1000),r0(1000)
      DIMENSION d0(1000)

c	    Variables temporaires de calcul
      DOUBLE PRECISION Ax0, Ad0, rkScal, rk1Scal, vecTemp, scal
      DIMENSION Ax0(1000), Ad0(1000), vecTemp(1000)
c	   produit A*x0 pour définir r0
      call produit(A, x0, n, Ax0)

c	    Definition de r0
      DO i=1,3
        r0(i) = b(i)-Ax0(i)
        d0(i) = r0(i)
      ENDDO


      res = 1000
      eps = 10e-5
      kmax = 1000
      i = 0
      DO WHILE ((res.gt.eps).AND.(i.le.kmax))
c       Calcul de alpha et de la norme du résidu
       call scalaire(r0, r0, n, rkScal)
       res = sqrt(rkScal)
       call produit(A, d0, n, vecTemp)
       call scalaire(vecTemp, d0, n, scal)
       alpha = rkScal/scal

c       Calcul de x_k+1
       DO j = 1,n
        x0(j) = x0(j) + alpha*d0(j)
       ENDDO

c       Calcul de r_k+1
       call produit(A, d0, n, Ad0)
       DO j = 1,n
       r0(j) = r0(j)-alpha*Ad0(j)
       ENDDO

c       Calcul de beta
        call scalaire(r0, r0, n, rk1Scal)
        beta = rk1Scal/rkScal

c       Calcul de la direction d_k+1
        DO j = 1,n
         d0(j) = r0(j) - beta*d0(j)
        ENDDO
        i=i+1
      ENDDO
      WRITE(*,*) "Solution obtenue en", i, "itérations"
      WRITE(*,*) x0(1:3)

      RETURN
      END


c	    PRODUIT MATRICE VECTEUR -> VECTEUR
      SUBROUTINE produit(A, b, n, y)
      IMPLICIT NONE
      DOUBLE PRECISION A, b, y
      INTEGER n,i, j
      DIMENSION A(1000,1000), b(1000), y(1000)
       do i=1,n
        y(i) = 0
        do j=1,n
         y(i) = y(i) + A(i, j)*b(j)
        enddo
      enddo
      RETURN
      END

c     PRODUIT SCALAIRE
      SUBROUTINE scalaire(a, b, n, y)
      DOUBLE PRECISION a, b, y
      DIMENSION a(1000), b(1000)
      INTEGER n, i
      y=0
      DO i=1,n
        y = y + a(i)*b(i)
      ENDDO
      RETURN
      END


c	    AFFICHAGE CONTENU MATRICE de taille nxn
  	  SUBROUTINE printMatrix(A, n)
  	  IMPLICIT NONE
  	  DOUBLE PRECISION A
  	  INTEGER n, i
  	  DIMENSION A(1000, 1000)
  	  DO i=1,n
  	   WRITE(*,*) A(i, 1:n)
  	  ENDDO

  	  RETURN
  	  END
