      INCLUDE "gradient.f"
      PROGRAM main
      IMPLICIT NONE
      DOUBLE PRECISION A,b, x,x0, eps, h, u
      INTEGER i, j, n
      DIMENSION A(1000,1000),b(1000),x(1000),x0(1000),u(1000)
      DOUBLE PRECISION Dh, K
      DIMENSION Dh(1000, 1000), K(1000,1000)
	  DOUBLE PRECISION f, S, temp, E
	  DIMENSION f(1000)

c     taille de la matrice
      n=31
      h = 1.d0/(n+1)


c      call Identity(n, K)
      call fill_diag(0, -1.0d0/(h*h), n, K)


c	   Remplissage de la matrice Dh
      call fill_diag(0, 4.0d0/(h*h), n, Dh)
      call fill_diag(-1, -1.0d0/(h*h), n, Dh)
      call fill_diag(1, -1.0d0/(h*h), n, Dh)
c	  call printMatrix(Dh, n)

c      Construction de la matrice A
	    call fill_diag_bloc(0, Dh, n, n*n, A)
	    call fill_diag_bloc(-1, K, n, n*n, A)
	    call fill_diag_bloc(1, K, n, n*n, A)
c	  call printMatrix(A, n*n)


	    DO i=1, n*n
	     f(i) = 0.0d0
	    ENDDO
	    DO i=n*(n-1)+1, n*n
	     f(i) = 100.0d0/(h*h)
	    ENDDO

c   Résolution par la méthode du gradient conjugué
	    call grad_conj(A, f, n*n, u)


c	  Affichage des résultats dans la console
c 	  DO i=1, n
c 	   print*, u(n*(n-i)+1:n*(n-i)+n)
c 	  ENDDO

c 	  Enregistrement données dans un fichier texte
 	   open(1, file='heat.dat')
 	   DO i=1, n
 	  	DO j=1,n
 	  	 WRITE(1, *) i*h, j*h, u((j-1)*N+i)
 	  	ENDDO
 	   ENDDO

c     Calcul de la moyenne des solutions
c      DO i=1,n*n
c       S = S + u(i)
c      ENDDO
c      print*, S/(n*n)

c     Affichage de la heat map avec gnuplot
      call system("gnuplot --persist plot2d.plot")

      S = 0.0d0
      DO i=1,n
      	DO j=1,n
      	  call temp_exacte(dble(i)*h, dble(j)*h, temp)
      	  S = S + DABS(u((j-1)*N+i)-temp)**2.0d0
      	ENDDO
      ENDDO
      E = (1.0d0/dble(n))*S
      print*, "Erreur quadratique moyenne =", E

      END

	    SUBROUTINE temp_exacte(x, y, u)
	    IMPLICIT NONE
	    DOUBLE PRECISION x, y, u, S, PI
	    INTEGER i, j, n

	    PI = 4.D0*DATAN(1.D0)
	    S = 0.0d0
	    DO i=0, 100
	  	S = S + (DSIN((2.0d0*dble(i)+1.0d0)*PI*x)*
     &  DSINH((2.0d0*dble(i)+1.0d0)*PI*y))/((2.0d0*dble(i)+1.0d0)*
     &  DSINH((2.0d0*dble(i)+1.0d0)*PI))

	    ENDDO

	    u = 4.0d0*100.0d0*S/PI

	    RETURN
	    END
