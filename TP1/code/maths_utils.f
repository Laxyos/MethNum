    
      SUBROUTINE diag(D, n, m, A)
      IMPLICIT NONE
      INTEGER i, j, n, m, k
      DOUBLE PRECISION D, A
      DIMENSION D(1000, 1000), A(1000, 1000)

      DO k=1,m
        DO i = 1,n
          DO j=1,n
            A((k-1)*n+i, (k-1)*n+j) = D(i, j)
          ENDDO
        ENDDO

      ENDDO
      RETURN
      END

c     Cette fonction permet de remplir les diagonales de A par des matrices blocs
c     pos : position de la diagonale (0 diag principale, 1 sur-diag, -1 sous-diag)
c     D matrice bloc de taille nxn à placer dans A
      SUBROUTINE fill_diag_bloc(pos, D, n, sizeA, A)
      IMPLICIT NONE
      INTEGER i, j, n, sizeA, k, pos
      DOUBLE PRECISION D, A, nbre
      DIMENSION D(1000, 1000), A(1000, 1000)

      nbre = dble(sizeA)/dble(n)
      DO k=1,int(nbre)-abs(pos)
        DO i = 1,n
          DO j=1,n
            IF(pos.eq.0) THEN
              A((k-1)*n+i, (k-1)*n+j) = D(i, j)
            ELSE IF (pos.gt.0) THEN
              A((k-1)*n+i, (k-1+pos)*n+j) = D(i, j)
            ELSE
              A((k-1+abs(pos))*n+i, (k-1)*n+j) = D(i, j)
            ENDIF
          ENDDO
        ENDDO

      ENDDO
      RETURN
      END


c	  Permet de remplir la diagonale numero D sur une matrice de taille sizeA
      SUBROUTINE fill_diag(pos, d, sizeA, A)
      IMPLICIT NONE
      INTEGER i, j, n, sizeA, k, pos
      DOUBLE PRECISION d, A
      DIMENSION  A(1000, 1000)

      DO k=1, sizeA-abs(pos)

        IF(pos.eq.0) THEN
          A(k, k) = d
        ELSE IF(pos.gt.0) THEN
          A(k, k+pos) = d
        ELSE
          A(k+abs(pos), k) = d
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE top_diag(D, n, m, A)
      IMPLICIT NONE
      INTEGER i, j, n, m, k
      DOUBLE PRECISION D, A
      DIMENSION D(1000, 1000), A(1000, 1000)

      DO k=1, m-1
        DO i = 1,n
          DO j=1,n
            A((k-1)*n+i, k*n+j) = D(i, j)
         ENDDO
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE bot_diag(D, n, m, A)
      IMPLICIT NONE
      INTEGER i, j, n, m, k
      DOUBLE PRECISION D, A
      DIMENSION D(1000, 1000), A(1000, 1000)

      DO k=1, m-1
        DO i = 1,n
          DO j=1,n
            A(k*n+i, (k-1)*n+j) = D(i, j)
         ENDDO
        ENDDO
      ENDDO
      RETURN
      END





c     RENVOIE UNE MATRICE IDENTITÉ DE TAILLE NxN
      SUBROUTINE identity(n, Id)
      IMPLICIT NONE
      INTEGER n, i, j
      DOUBLE PRECISION Id
      DIMENSION Id(1000, 1000)
      DO i=1,n+1
        DO j=1,n+1
          IF (i.eq.j) THEN
            Id(i, j) = 1
          ELSE
            Id(i, j) = 0
          ENDIF
        ENDDO
      ENDDO
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
