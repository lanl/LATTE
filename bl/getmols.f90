PROGRAM CFGIR
  
  IMPLICIT NONE

  INTEGER :: I, J, K, STEP, COUNT, II, JJ, KK,  NOMOL, MYATS
  INTEGER, ALLOCATABLE :: NEBS(:,:,:), NNEB(:)
  INTEGER :: STACKINDEX, ATOMINDEX1, MOLECULEINDEX=1
  INTEGER, PARAMETER :: NATS = 336
  INTEGER, PARAMETER :: NC = 16*7, NH = 16*5, NN = 16*3, NO = 16*6
  INTEGER, ALLOCATABLE :: NEIGHBORINDEXSTACK(:), ATOMINDEXSTACK(:)
  INTEGER, ALLOCATABLE :: MOLID(:)
  REAL :: X(3,NATS)
  REAL ::  F(3,NATS), Y2
  REAL :: V(3,NATS), Q(NATS), BOX(3,3), R2, RIJ(3), MAGR, MYQ
  CHARACTER(LEN=50) :: Y, FLNM
  CHARACTER(LEN=2) :: NSPEC(4)
  CHARACTER(LEN=2) :: SPEC(NATS)

  
  OPEN(UNIT=11, STATUS="OLD", FILE="dump2atomeye.494500.cfg")
  
  READ(11,*) Y, Y, Y, Y, I
  READ(11,*) Y, Y, Y2, Y
  READ(11,*) Y, Y, Y, BOX(1,1), Y
  READ(11,*) Y, Y, Y, BOX(1,2), Y
  READ(11,*) Y, Y, Y, BOX(1,3), Y
  READ(11,*) Y, Y, Y, BOX(2,1), Y
  READ(11,*) Y, Y, Y, BOX(2,2), Y
  READ(11,*) Y, Y, Y, BOX(2,3), Y
  READ(11,*) Y, Y, Y, BOX(3,1), Y
  READ(11,*) Y, Y, Y, BOX(3,2), Y
  READ(11,*) Y, Y, Y, BOX(3,3), Y
  READ(11,*) Y, Y, I
  READ(11,*) Y, Y, Y
  READ(11,*) Y, Y, Y
  READ(11,*) Y, Y, Y
  READ(11,*) Y, Y, Y
  READ(11,*) Y, Y, Y
  READ(11,*) Y, Y, Y
  
  READ(11,*) Y2
  READ(11,*) Y
  
  COUNT = 0
  
  DO I = 1, NC
     
     COUNT = COUNT + 1
     READ(11,*) X(1,COUNT), X(2,COUNT), X(3,COUNT), &
          V(1,COUNT), V(2,COUNT), V(3,COUNT), &
          F(1,COUNT), F(2,COUNT), F(3,COUNT), &
          Y2, Y2, Q(COUNT)
     
     SPEC(COUNT) = "C"
     
  ENDDO
  
  READ(11,*) Y2
  READ(11,*) Y
  
  DO I = 1, NH
     
     COUNT = COUNT + 1
     READ(11,*) X(1,COUNT), X(2,COUNT), X(3,COUNT), &
          V(1,COUNT), V(2,COUNT), V(3,COUNT), &
          F(1,COUNT), F(2,COUNT), F(3,COUNT), &
          Y2, Y2, Q(COUNT)
     
     SPEC(COUNT) = "H"
     
  ENDDO
  
  READ(11,*) Y2
  READ(11,*) Y
  
  DO I = 1, NN
     
     COUNT = COUNT + 1
     READ(11,*) X(1,COUNT), X(2,COUNT), X(3,COUNT), &
          V(1,COUNT), V(2,COUNT), V(3,COUNT), &
          F(1,COUNT), F(2,COUNT), F(3,COUNT), &
          Y2, Y2, Q(COUNT)
     
     SPEC(COUNT) = "N"
     
  ENDDO
     
  READ(11,*) Y2
  READ(11,*) Y
  
  DO I = 1, NO
     
     COUNT = COUNT + 1
     READ(11,*) X(1,COUNT), X(2,COUNT), X(3,COUNT), &
          V(1,COUNT), V(2,COUNT), V(3,COUNT), &
          F(1,COUNT), F(2,COUNT), F(3,COUNT), &
          Y2, Y2, Q(COUNT)
     
     SPEC(COUNT) = "O"
     
  ENDDO
  
  CLOSE(11)
  
  
  DO I = 1, NATS
     X(1,I) = X(1,I)*BOX(1,1)
     X(2,I) = X(2,I)*BOX(2,2)
     X(3,I) = X(3,I)*BOX(3,3)
  ENDDO
  
  ALLOCATE(NEBS(NATS,NATS,4), NNEB(NATS), MOLID(NATS))
  
  NEBS = 0
  NNEB = 0
  
  DO I = 1, NATS
     DO J = 1, NATS
        
        DO II = -1, 1
           DO JJ = -1, 1
              DO KK = -1, 1

                 RIJ(1) = X(1,J) + REAL(II)*BOX(1,1) + &
                      REAL(JJ)*BOX(2,1) + REAL(KK)*BOX(3,1) - X(1,I)
                 
                 RIJ(2) = X(2,J) + REAL(II)*BOX(1,2) + &
                      REAL(JJ)*BOX(2,2) + REAL(KK)*BOX(3,2) - X(2,I)
                 
                 RIJ(3) = X(3,J) + REAL(II)*BOX(1,3) + &
                      REAL(JJ)*BOX(2,3) + REAL(KK)*BOX(3,3) - X(3,I)

                 R2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

                 MAGR = SQRT(R2)

                 IF (SPEC(I) .EQ. "C" .AND. SPEC(J) .EQ. "C" &
                         .AND. MAGR .LT. 1.6) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "C" .AND. SPEC(J) .EQ. "H" &
                         .AND. MAGR .LT. 1.2) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "C" .AND. SPEC(J) .EQ. "N" &
                         .AND. MAGR .LT. 1.6) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "C" .AND. SPEC(J) .EQ. "O" &
                         .AND. MAGR .LT. 1.4) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J


                    ELSEIF (SPEC(I) .EQ. "O" .AND. SPEC(J) .EQ. "C" &
                         .AND. MAGR .LT. 1.4) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "O" .AND. SPEC(J) .EQ. "H" &
                         .AND. MAGR .LT. 1.2) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J


                    ELSEIF (SPEC(I) .EQ. "O" .AND. SPEC(J) .EQ. "N" &
                         .AND. MAGR .LT. 1.5) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "O" .AND. SPEC(J) .EQ. "O" &
                         .AND. MAGR .LT. 1.4) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "N" .AND. SPEC(J) .EQ. "C" &
                         .AND. MAGR .LT. 1.6) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "N" .AND. SPEC(J) .EQ. "H" &
                         .AND. MAGR .LT. 1.2) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "N" .AND. SPEC(J) .EQ. "N" &
                         .AND. MAGR .LT. 1.5) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "N" .AND. SPEC(J) .EQ. "O" &
                         .AND. MAGR .LT. 1.5) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "H" .AND. SPEC(J) .EQ. "C" &
                         .AND. MAGR .LT. 1.2) THEN
                       
                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J
                       
                    ELSEIF (SPEC(I) .EQ. "H" .AND. SPEC(J) .EQ. "N" &
                         .AND. MAGR .LT. 1.2) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "H" .AND. SPEC(J) .EQ. "O" &
                         .AND. MAGR .LT. 1.2) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ELSEIF (SPEC(I) .EQ. "H" .AND. SPEC(J) .EQ. "H" &
                         .AND. MAGR .LT. 0.85) THEN

                       NNEB(I) = NNEB(I) + 1
                       NEBS(I,NNEB(I), 1) = J

                    ENDIF

!                 IF (R2 .LT. 1.55**2) THEN
!                    NNEB(I) = NNEB(I) + 1
!!                    NEBS(I, NNEB(I), 1) = J
!                    NEBS(I, NNEB(I), 2) = II
!                    NEBS(I, NNEB(I), 3) = JJ
!                    NEBS(I, NNEB(I), 4) = KK
!                 ENDIF

              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ALLOCATE(ATOMINDEXSTACK(NATS+1))
  ALLOCATE(NEIGHBORINDEXSTACK(NATS+1))

  DO ATOMINDEX1 = 1, NATS
     MOLID(ATOMINDEX1) = 0
  END DO
  
  !                                                                               
  ! WE NEED TO TRAVERSE THE NEIGHBOR TREE, SO WE'LL USE A STACK                   
  ! OF INDICES INTO THE NEIGHBOR LISTS, AND SKIP OVER BRANCHES                    
  ! WE'VE ALREADY DONE                                                            
  !                                                                               

  STACKINDEX = 1
  ATOMINDEXSTACK(1) = 0

  DO

     !                                                                            
     ! SET UP ATOMINDEXSTACK                                                      
     !                                                                            

     IF (STACKINDEX .GT. 1) THEN

        ATOMINDEXSTACK(STACKINDEX) = &
             NEBS(ATOMINDEXSTACK(STACKINDEX-1), &
             NEIGHBORINDEXSTACK(STACKINDEX), 1)

     END IF

     !                                                                            
     ! DEAL WITH CURRENT ATOM                                                     
     !                                                                            

     IF (ATOMINDEXSTACK(STACKINDEX) .GT. 0) THEN

        IF (MOLID(ATOMINDEXSTACK(STACKINDEX)) .EQ. 0) THEN

           MOLID(ATOMINDEXSTACK(STACKINDEX)) = MOLECULEINDEX

           STACKINDEX = STACKINDEX+1

           NEIGHBORINDEXSTACK(STACKINDEX) = 0

        END IF
     END IF

     !                                                                            
     ! INCREMENT THE CURRENT ATOM                                                 
     !                                                                            

     IF (STACKINDEX .GT. 1) THEN

        NEIGHBORINDEXSTACK(STACKINDEX) = NEIGHBORINDEXSTACK(STACKINDEX) + 1

     ELSE IF (STACKINDEX .EQ. 1) THEN

        ATOMINDEXSTACK(STACKINDEX) = ATOMINDEXSTACK(STACKINDEX) + 1

     END IF

     !                                                                            
     ! CHECK BOUNDS                                                               
     !                                                                            

     IF (STACKINDEX .GT. 1) THEN

        IF (NEIGHBORINDEXSTACK(STACKINDEX) > &
             NNEB(ATOMINDEXSTACK(STACKINDEX-1))) THEN

           STACKINDEX = STACKINDEX - 1

           IF (STACKINDEX .EQ. 1) THEN

              MOLECULEINDEX = MOLECULEINDEX + 1

           END IF

        END IF

     ELSE IF (STACKINDEX .EQ. 1) THEN

        IF (ATOMINDEXSTACK(STACKINDEX) > NATS) THEN

           EXIT

        END IF

     END IF
  END DO

  DEALLOCATE(ATOMINDEXSTACK, NEIGHBORINDEXSTACK)

  ! MJC put in this bit to get the total number of molecules                      
  !                                                                               


  NOMOL = MAXVAL(MOLID)
  

  print*, NOMOL

  DO II = 1, NOMOL
     
     IF (II .LT. 10) THEN
        WRITE(FLNM, '("mymols.",I1,".xyz")') II
     ELSEIF (II .GE. 10 .AND. II .LT. 100) THEN
        WRITE(FLNM, '("mymols.",I2,".xyz")') II
     ENDIF

     OPEN(UNIT=22, STATUS="UNKNOWN", FILE=FLNM)
     
     MYQ = 0.0
     MYATS = 0

     DO J = 1, NATS
        IF (MOLID(J) .EQ. II) THEN
           MYQ = MYQ + Q(J)
           MYATS = MYATS + 1
        ENDIF
     ENDDO


     WRITE(22,*) MYATS
     WRITE(22,*) MYQ
     
     DO J = 1, NATS
        IF (MOLID(J) .EQ. II) THEN
           WRITE(22,*) SPEC(J), X(1,J), X(2,J), X(3,J) 
        ENDIF
     ENDDO

     CLOSE(22)
     
  ENDDO

        
!  DO J = 1, NATS
!     IF (MOLID(J) .EQ. 1) THEN
!        PRINT*, SPEC(J), X(1,J), X(2,J), X(3,J)
!     ENDIF
!  ENDDO
     


END PROGRAM CFGIR
