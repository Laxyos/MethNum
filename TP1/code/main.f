      INCLUDE "gradient.f"
      PROGRAM main
      IMPLICIT NONE
      DOUBLE PRECISION E, T
      INTEGER i, j, k


      open(3, file='data/plaque_slow.dat')
      call plaque2(30, E, T)


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
      DOUBLE PRECISION test

      DOUBLE PRECISION start, finish, T


c     pas du maillage
      h = 1.d0/(n+1)



      call fill_diag(0, -1.0d0/(h*h), n, K)


c	   Remplissage de la matrice Dh
      call fill_diag(0, 4.0d0/(h*h), n, Dh)
      call fill_diag(-1, -1.0d0/(h*h), n, Dh)
      call fill_diag(1, -1.0d0/(h*h), n, Dh)

c      Construction de la matrice A
	    call fill_diag_bloc(0, Dh, n, n*n, A)
	    call fill_diag_bloc(-1, K, n, n*n, A)
	    call fill_diag_bloc(1, K, n, n*n, A)


	    open(20, file='data/matrice.dat')
	    DO i=1,n*n
	      WRITE(20,*) A(i, 1:n*n)
	    ENDDO


	    DO i=1, n*n
	      u(i) = 1.0d0
	      f(i) = 0.0d0
	    ENDDO


c     Conditions aux limites
      DO i=1,n
        f(n*(n-1)+i) = 100.0d0/(h*h)
      ENDDO

	  call cpu_time(start)
c   Résolution par la méthode du gradient conjugué
	  call grad_conj2(A, f, n*n, u)
	  call cpu_time(finish)
	  print*, "Execution du gradient conjugué en",
     & (finish-start)*1000.d0,"ms"
      T = (finish-start)*1000.d0


c 	  Enregistrement données dans un fichier texte
 	  open(1, file='data/heat.dat')
c     Ecriture des valeurs imposées aux limites

	  DO i=1, n
	 	DO j=1,n
	 	  WRITE(1, *) i*h, j*h, u((j-1)*N+i)
	 	ENDDO
	  ENDDO

c     Affichage de la heat map avec gnuplot
	  call system('gnuplot plot2d.plot')




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


  	  SUBROUTINE plaque(n, E, T)
  	  IMPLICIT NONE
      DOUBLE PRECISION A,b, h, u
      INTEGER i, j, n
      DIMENSION A(1000,1000),u(1000)
      DOUBLE PRECISION Dh, K
      DIMENSION Dh(1000, 1000), K(1000,1000)
	    DOUBLE PRECISION f, S, temp, E
	    DIMENSION f(1000)
      DOUBLE PRECISION test

      DOUBLE PRECISION start, finish, T


c     pas du maillage
      h = 1.d0/(n+1)



      call fill_diag(0, -1.0d0/(h*h), n, K)


c	   Remplissage de la matrice Dh
      call fill_diag(0, 4.0d0/(h*h), n, Dh)
      call fill_diag(-1, -1.0d0/(h*h), n, Dh)
      call fill_diag(1, -1.0d0/(h*h), n, Dh)

c      Construction de la matrice A
	    call fill_diag_bloc(0, Dh, n, n*n, A)
	    call fill_diag_bloc(-1, K, n, n*n, A)
	    call fill_diag_bloc(1, K, n, n*n, A)


	    open(20, file='data/matrice.dat')
	    DO i=1,n*n
	      WRITE(20,*) A(i, 1:n*n)
	    ENDDO


	    DO i=1, n*n
	      u(i) = 1.0d0
	      f(i) = 0.0d0
	    ENDDO


c     Conditions aux limites
      DO i=1,n
        f(n*(n-1)+i) = 100.0d0/(h*h)
      ENDDO

	  call cpu_time(start)
c   Résolution par la méthode du gradient conjugué
	  call grad_conj(A, f, n*n, u)
	  call cpu_time(finish)
	  print*, "Execution du gradient conjugué en",
     & (finish-start)*1000.d0,"ms"
      T = (finish-start)*1000.d0


c 	  Enregistrement données dans un fichier texte
 	  open(1, file='data/heat.dat')
c     Ecriture des valeurs imposées aux limites
c	  DO i=1, n
c		WRITE(1, *) 0.0d0, i*h, 0.0d0
c		WRITE(1, *) 1.0d0, i*h, 0.0d0
c		WRITE(1, *) i*h, 0.0d0, 0.0d0
c		WRITE(1, *) i*h, 1.0d0, 100.0d0
c	  ENDDO
c     Coins de la plaque
c	  WRITE(1,*) 0.0d0, 1.0d0, 100.0d0
c	  WRITE(1,*) 1.0d0, 1.0d0, 100.0d0
c	  WRITE(1,*) 0.0d0, 0.0d0, 0.0d0
c	  WRITE(1,*) 1.0d0, 0.0d0, 0.0d0

	  DO i=1, n
	 	DO j=1,n
	 	  WRITE(1, *) i*h, j*h, u((j-1)*N+i)
	 	ENDDO
	  ENDDO

c     Affichage de la heat map avec gnuplot
	  call system('gnuplot plot2d.plot')




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
