      INCLUDE "maths_utils.f"

c     Méthode du gradient conjugué pour Ax = b
c     Renvoit le vecteur x0
      SUBROUTINE grad_conj(A, b, n, x0)
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
      DO i=1,n
        r0(i) = b(i)-Ax0(i)
        d0(i) = r0(i)
      ENDDO


      res = 1000
      eps = 10e-5
      kmax = 10000
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
      print*,"résolution en ", i, "itérations"

      RETURN
      END
