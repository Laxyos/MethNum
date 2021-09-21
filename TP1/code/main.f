      INCLUDE "gradient.f"
      PROGRAM main
      IMPLICIT NONE
      DOUBLE PRECISION E, T
      INTEGER i, j, k


c      open(3, file='plaque_slow.dat')
c      DO i=3,31
c       call plaque(i, E, T)
c       WRITE(3,*) i, E, T
c      ENDDO
      call plaque2(5, E, T)

      END




  	  SUBROUTINE plaque2(n, E, T)
  	  IMPLICIT NONE
      DOUBLE PRECISION A,b, h, u
      INTEGER i, j, n
      DIMENSION A(1000,1000),u(1000)
      DOUBLE PRECISION Dh, K
      DIMENSION Dh(1000, 1000), K(1000,1000)
	    DOUBLE PRECISION f, S, temp, E
	    DIMENSION f(1000)


      DOUBLE PRECISION start, finish, T

c     pas du maillage
      h = 1.d0/(n+1)



c      call Identity(n, K)
      call fill_diag(0, -1.0d0/(h*h), n, K)


c	   Remplissage de la matrice Dh
      call fill_diag(0, 4.0d0/(h*h), n, Dh)
      call fill_diag(-1, -1.0d0/(h*h), n, Dh)
      call fill_diag(1, -1.0d0/(h*h), n, Dh)

c      Construction de la matrice A
	  call fill_diag_bloc(0, Dh, n, n*n, A)
	  call fill_diag_bloc(-1, K, n, n*n, A)
	  call fill_diag_bloc(1, K, n, n*n, A)




	  DO i=1, n*n
	   u(i) = 1.0d0
	   f(i) = 0.0d0
	  ENDDO
c     (j-1)*N+i
      DO i=1,n
        print*,n*(n-1)+i
        f(n*(n-1)+i) = 100.0d0/(h*h)
      ENDDO
c	  DO i=n*(n-1)+1, n*n
c	   f(i) = 100.0d0/(h*h)
c	  ENDDO

	  call cpu_time(start)
c   Résolution par la méthode du gradient conjugué
	  call grad_conj2(A, f, n*n, u)
	  call cpu_time(finish)
	  print*, "Execution du gradient conjugué en",
     & (finish-start)*1000.d0,"ms"
      T = (finish-start)*1000.d0


c 	  Enregistrement données dans un fichier texte
 	    open(1, file='heat.dat')
      DO i=0, n+1
       WRITE(1, *) 0.0d0, i*h, 0.0d0
        WRITE(1, *) 1.0d0, i*h, 0.0d0
        WRITE(1, *) i*h, 0.0d0, 0.0d0
        WRITE(1, *) i*h, 1.0d0, 100.0d0
      ENDDO
 	    DO i=1, n
 	  	 DO j=1,n
 	  	  WRITE(1, *) i*h, j*h, u((j-1)*N+i)
 	  	 ENDDO
 	    ENDDO


c     Affichage de la heat map avec gnuplot
      call system("gnuplot --persist plot2d.plot")

c	  Calcul de l'erreur quadratique moyenne
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


c     VERSION LENTE
  	  SUBROUTINE plaque(n, E, T)
  	  IMPLICIT NONE
      DOUBLE PRECISION A, h, u
      INTEGER i, j, n
      DIMENSION A(1000,1000),u(1000)
      DOUBLE PRECISION Dh, K
      DIMENSION Dh(1000, 1000), K(1000,1000)
	  DOUBLE PRECISION f, S, temp, E
	  DIMENSION f(1000)

      DOUBLE PRECISION y1, y2, x1
      DIMENSION y1(1000), y2(1000),x1(1000)

      DOUBLE PRECISION start, finish, T

c     taille du maillage
      h = 1.d0/(n+1)



c      call Identity(n, K)
      call fill_diag(0, -1.0d0/(h*h), n, K)


c	   Remplissage de la matrice Dh
      call fill_diag(0, 4.0d0/(h*h), n, Dh)
      call fill_diag(-1, -1.0d0/(h*h), n, Dh)
      call fill_diag(1, -1.0d0/(h*h), n, Dh)

c      Construction de la matrice A
	  call fill_diag_bloc(0, Dh, n, n*n, A)
	  call fill_diag_bloc(-1, K, n, n*n, A)
	  call fill_diag_bloc(1, K, n, n*n, A)



	  DO i=1, n*n
	   u(i) = 1.0d0
	   f(i) = 0.0d0
	  ENDDO
	  DO i=n*(n-1)+1, n*n
	   f(i) = 100.0d0/(h*h)
	  ENDDO

	  call cpu_time(start)
c   Résolution par la méthode du gradient conjugué
	  call grad_conj(A, f, n*n, u)
	  call cpu_time(finish)
	  print*, "Execution du gradient conjugué en",
     & (finish-start)*1000.d0,"ms"
      T = (finish-start)*1000.d0


c 	  Enregistrement données dans un fichier texte
 	   open(1, file='heat.dat')
 	   DO i=1, n
 	  	DO j=1,n
 	  	 WRITE(1, *) i*h, j*h, u((j-1)*N+i)
 	  	ENDDO
 	   ENDDO

c     Affichage de la heat map avec gnuplot
C      call system("gnuplot --persist plot2d.plot")

c	  Calcul de l'erreur quadratique moyenne
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


c		Calcul analytique de la température au point (x,y)
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
