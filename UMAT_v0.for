      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2)
	  
	  INTEGER, PARAMETER :: DP = selected_real_kind(p=15, r=307)
	  PARAMETER (maxN = 256, maxL = 7)
	  
	  REAL(DP) INP_MIN(6), INP_MAX(6), OUT_MIN(6), OUT_MAX(6)
	  REAL W(maxL+1, maxN, maxN+1)
	  INTEGER Struc(maxL+2),idx(2),L
C#***********************************************************************************#
	  DATA L /4/, Struc /6, 128, 256, 256, 128, 6/
C#***********************************************************************************#
	  COMMON /surrogate/ INP_MIN, INP_MAX, OUT_MIN, OUT_MAX, W, Struc, L


	  IF (LOP.EQ.0) THEN 
		
		OPEN(15,file='D:\FE-dNN_example\inputs_min.txt')
		READ(15,*) INP_MIN
		CLOSE(15)
		
		OPEN(16,file='D:\FE-dNN_example\inputs_max.txt')
		READ(16,*) INP_MAX
		CLOSE(16)
		
		OPEN(17,file='D:\FE-dNN_example\outputs_min.txt')
		READ(17,*) OUT_MIN
		CLOSE(17)
		
		OPEN(18,file='D:\FE-dNN_example\outputs_max.txt')
		READ(18,*) OUT_MAX
		CLOSE(18)
		
		
		idx =(/0, 0/)
		OPEN(20,file='D:\FE-dNN_example\model_nh.txt')

		DO i = 1,L+1

			idx(1) = idx(2)
			idx(2) = idx(2)+Struc(i+1)

			DO j = idx(1)+1,idx(2)
				READ(20,*) W(i,j-idx(1),:Struc(i)+1)
			ENDDO

		ENDDO
		CLOSE(20)
		
	  END IF


      RETURN
      END


C ################## UMAT ##############

      SUBROUTINE UMAT(STRESS ,STATEV ,DDSDDE, SSE, SPD,
     1 SCD,RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,
     2 DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI ,NSHR,
     3 NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1 , NOEL , NPT , LAYER , KSPT , KSTEP ,
     5 KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
	  DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

	  	  
	  DIMENSION EELAS(6), EELASP(3), BBAR(6), BBARP(3), BBARN(3,3), DISTGR(3,3)
      Parameter (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0, FOUR= 4.0D0, SIX= 6.0D0)
	  
	  INTEGER, PARAMETER :: DP = selected_real_kind(p=15, r=307)
	  INTEGER M1, M2, N1, N2, K1, K2
	  REAL(DP) R(3,3), U(3,3), E(3,3), CG(3,3), PK2(3,3), PK1(3,3), SIGMA(3,3), J, CIJKL(3,3,3,3), DELTA(3,3),TANGENT(6,6),
     1 TEM(3,3),FTINV(3,3), DPTS(6), Y(6), GRAD(6,6), TERM1(3,3,3,3), TERM2(3,3,3,3), TERM3(3,3,3,3),CSE(3,3,3,3)
	 

C#***********************************************************************************#	  
	  
	  DO K1 = 1, 3
		DO K2 = 1, 3
			IF (K1.EQ.K2) THEN
				DELTA(K1,K2) = 1.0D0
			ELSE
				DELTA(K1,K2) = 0.0D0
			ENDIF
	    END DO
	  END DO
C	  CALL POLAR(DFGRD1, U , R)
	  
C#***********************************************************************************#		  
	  
	  CALL DETERMINANT(DFGRD1, J)	

	  CG = MATMUL(TRANSPOSE(DFGRD1), DFGRD1)
	  E =  (CG - DELTA)/ 2.0D0
	  
	  
	  DPTS(1) = E(1,1)
	  DPTS(2) = E(2,1)
	  DPTS(3) = E(3,1)
	  DPTS(4) = E(2,2)
	  DPTS(5) = E(3,2)
	  DPTS(6) = E(3,3)	
	  
	  CALL NN(DPTS,GRAD)
		
	  PK2(1,1) = DPTS(1)
	  PK2(1,2) = DPTS(2)
	  PK2(1,3) = DPTS(3)
	  PK2(2,2) = DPTS(4)
	  PK2(2,3) = DPTS(5)
	  PK2(3,3) = DPTS(6) 
		
	  
	  DO K1 =1 , 3
		DO K2 =1 , K1-1
			PK2(K1, K2) = PK2(K2, K1)
		END DO
	  END DO

	  
	  PK1 = MATMUL(PK2,TRANSPOSE(DFGRD1))
	  
	  SIGMA = MATMUL(DFGRD1,PK1)/J
	  
C	  Print*, 'SIGMA:', SIGMA(1,1),SIGMA(2,2),SIGMA(3,3),SIGMA(1,2),SIGMA(1,3),SIGMA(2,3)
	  
	  STRESS(1) = SIGMA(1,1)
	  STRESS(2) = SIGMA(2,2)
	  STRESS(3) = SIGMA(3,3)
	  STRESS(4) = SIGMA(1,2)
	  STRESS(5) = SIGMA(1,3)
	  STRESS(6) = SIGMA(2,3)
	  
	  
	  GRAD(:,2) = GRAD(:,2)/2.0D0
	  GRAD(:,3) = GRAD(:,3)/2.0D0
	  GRAD(:,5) = GRAD(:,5)/2.0D0
	  
	  GRAD = (GRAD + TRANSPOSE(GRAD)) / 2.0D0
	  
      DO K1 = 1, 6
        DO K2 = 1, 6
			CALL get_mn_index(K1, M1, M2)
			CALL get_mn_index(K2, N1, N2)
			CSE(M1,M2,N1,N2) = GRAD(K1,K2)
		END DO
	  END DO
	  
	  DO K1 = 1, 3
		DO K2 = K1, 3
			IF (K1.NE.K2) THEN
				CSE(:,:,K2,K1) = CSE(:,:,K1,K2)
			ENDIF
	    END DO
	  END DO
	  
	  DO K1 = 1, 3
		DO K2 = K1, 3
			IF (K1.NE.K2) THEN
				CSE(K2,K1,:,:) = CSE(K1,K2,:,:)
			ENDIF
	    END DO
	  END DO
	  
C#***********************************************************************************#	  
	  TERM1 = 0.0
	  DO M1 = 1, 3
        DO M2 = 1, 3
            DO N1 = 1, 3
                DO N2 = 1, 3
                    TERM1(M1,M2,N1,N2) = DELTA(M1,N1) * PK2(N2,M2)
                END DO
            END DO
        END DO
	  END DO
	  
	  TERM2 = 0.0
	  DO M1 = 1, 3
        DO M2 = 1, 3
            DO N1 = 1, 3
                DO N2 = 1, 3
                    DO K1 = 1, 3
						TERM2(M1,M2,N1,N2) = TERM2(M1,M2,N1,N2) + DFGRD1(K1,M1) * CSE(K1,M2,N1,N2)
					END DO
                END DO
            END DO
        END DO
	  END DO

	  TERM3 = 0.0
	  DO M1 = 1, 3
        DO M2 = 1, 3
            DO N1 = 1, 3
                DO N2 = 1, 3
                    DO K1 = 1, 3
						TERM3(M1,M2,N1,N2) = TERM3(M1,M2,N1,N2) + DFGRD1(K1,N1) * TERM2(M1,M2,K1,N2)
					END DO
                END DO
            END DO
        END DO
	  END DO  

	  CIJKL = TERM1 + TERM3
	  
C#***********************************************************************************#	  
	  CALL INVERSE(DFGRD1,FTINV)
	  FTINV = TRANSPOSE(FTINV)
	  
	  TERM1 = 0.0
	  DO M1 = 1, 3
        DO M2 = 1, 3
            DO N1 = 1, 3
                DO N2 = 1, 3
                    TERM1(M1,M2,N1,N2) = -FTINV(N1,N2) * SIGMA(M1,M2)
                END DO
            END DO
        END DO
	  END DO
	  
	  TERM2 = 0.0
	  DO M1 = 1, 3
        DO M2 = 1, 3
            DO N1 = 1, 3
                DO N2 = 1, 3
					DO K1 = 1, 3
						TERM2(M1,M2,N1,N2) = TERM2(M1,M2,N1,N2) + CIJKL(M1,K1,N1,N2) * DFGRD1(K1,M2) / J
					END DO
				END DO
            END DO
        END DO
	  END DO
	  	  
	  TERM3 = 0.0
	  DO M1 = 1, 3
        DO M2 = 1, 3
            DO N1 = 1, 3
                DO N2 = 1, 3			
					TERM3(M1,M2,N1,N2) = PK1(N2,M1) * DELTA(M2,N1) / J
				END DO
            END DO
        END DO
	  END DO
	  	  
	  CIJKL = TERM1 + TERM2 + TERM3
	  
	  TERM1 = 0.0
	  DO M1 = 1, 3
        DO M2 = 1, 3
            DO N1 = 1, 3
                DO N2 = 1, 3
                    TERM1(M1,M2,N1,N2) = SIGMA(M1,M2) * DELTA(N1,N2)
                END DO
            END DO
        END DO
	  END DO
	  
	  TERM2 = 0.0
	  DO M1 = 1, 3
        DO M2 = 1, 3
            DO N1 = 1, 3
                DO N2 = 1, 3
                    DO K1 = 1, 3
						TERM2(M1,M2,N1,N2) = TERM2(M1,M2,N1,N2) + CIJKL(M1,M2,N1,K1) * DFGRD1(K1,N2)
					END DO
                END DO
            END DO
        END DO
	  END DO
	  
	  CIJKL = TERM1 + TERM2
	  
C#***********************************************************************************#	 	  
	  TANGENT(1,1) = CIJKL(1,1,1,1)
	  TANGENT(1,2) = CIJKL(1,1,2,2)
	  TANGENT(1,3) = CIJKL(1,1,3,3)
	  TANGENT(1,4) = (CIJKL(1,1,1,2)+ CIJKL(1,1,2,1))/2.0D0
	  TANGENT(1,5) = (CIJKL(1,1,1,3)+ CIJKL(1,1,3,1))/2.0D0
	  TANGENT(1,6) = (CIJKL(1,1,2,3)+ CIJKL(1,1,3,2))/2.0D0
	  
	  TANGENT(2,2) = CIJKL(2,2,2,2)
	  TANGENT(2,3) = CIJKL(2,2,3,3)
	  TANGENT(2,4) = (CIJKL(2,2,1,2)+ CIJKL(2,2,2,1))/2.0D0
	  TANGENT(2,5) = (CIJKL(2,2,1,3)+ CIJKL(2,2,3,1))/2.0D0
	  TANGENT(2,6) = (CIJKL(2,2,2,3)+ CIJKL(2,2,3,2))/2.0D0
	  
	  TANGENT(3,3) = CIJKL(3,3,3,3)
	  TANGENT(3,4) = (CIJKL(3,3,1,2)+ CIJKL(3,3,2,1))/2.0D0
	  TANGENT(3,5) = (CIJKL(3,3,1,3)+ CIJKL(3,3,3,1))/2.0D0
	  TANGENT(3,6) = (CIJKL(3,3,2,3)+ CIJKL(3,3,3,2))/2.0D0
	  
	  TANGENT(4,4) = (CIJKL(1,2,1,2)+ CIJKL(1,2,2,1))/2.0D0
	  TANGENT(4,5) = (CIJKL(1,2,1,3)+ CIJKL(1,2,3,1))/2.0D0
	  TANGENT(4,6) = (CIJKL(1,2,2,3)+ CIJKL(1,2,3,2))/2.0D0
	  
	  TANGENT(5,5) = (CIJKL(1,3,1,3)+ CIJKL(1,3,3,1))/2.0D0
	  TANGENT(5,6) = (CIJKL(1,3,2,3)+ CIJKL(1,3,3,2))/2.0D0
	  
	  TANGENT(6,6) = (CIJKL(2,3,2,3)+ CIJKL(2,3,3,2))/2.0D0


	  DO K1 = 1, 6
		DO K2 = 1, 6
			IF (K1.NE.K2) THEN
				TANGENT(K2, K1) = TANGENT(K1, K2)
			ENDIF
	    END DO
	  END DO
	  
C	  Print*, 'TANGENT:', TANGENT
		

		
	  DDSDDE = TANGENT 
	

	  IF (NOEL.EQ.111) THEN
		IF (KINC.EQ.1) THEN
			IF (NPT.EQ.3) THEN
		Print*, 'KINC :', KINC 
		Print*, 'NPT:', NPT
		Print*, 'DFGRD1:', DFGRD1
		Print*, 'STRESS:', STRESS
		Print*, 'DDSDDE:', DDSDDE
		Print*, '----------------------------'
		END IF
			END IF
	  END IF

	  
      RETURN
      END

C ################## right polar decomposition##############
	  SUBROUTINE POLAR(F, U ,R)
	  
	  INCLUDE 'ABA_PARAM.INC'
	  
	  DIMENSION F(3,3)
	  INTEGER, PARAMETER :: DP = selected_real_kind(p=15, r=307) 
	  REAL(DP) C(3,3), V(3,3), RT(3,3), Q(3,3), UPPER3(3,3), R(3,3), U(3,3), INVU(3,3), det, max_off_diag, theta
	  INTEGER IND(2)
	  PARAMETER (tol=1.0e-6, max_iter=1000, pi  = 4 * ATAN(1.0_8))
		
	  DO K1 = 1, 3
		DO K2 = 1, 3
			IF (K1.EQ.K2) THEN
				Q(K1,K2) = 1.0D0
			ELSE
				Q(K1,K2) = 0.0D0
			ENDIF
	    END DO
	  END DO

	  C = MATMUL(TRANSPOSE(F),F)
	  
	  DO K = 1, max_iter
	  
		UPPER3 = 0.0D0
		UPPER3(3,1) = C(3,1)
	    UPPER3(3,2) = C(3,2)
	    UPPER3(2,1) = C(2,1)
		max_off_diag = MAXVAL(ABS(UPPER3))
		
		IF (max_off_diag.LT.tol) EXIT

		
		IND = MAXLOC(ABS(UPPER3))
		
		IF (C(IND(1),IND(1)).EQ.C(IND(2),IND(2))) THEN 
			theta = pi/4 
        ELSE
			theta = 0.5D0 * ATAN(2.0D0 * C(IND(1),IND(2)) / (C(IND(2),IND(2)) - C(IND(1),IND(1))))
        END IF 
		
		
		
		DO K1 = 1, 3
			DO K2 = 1, 3
				IF (K1.EQ.K2) THEN
					V(K1,K2) = 1.0D0
				ELSE
					V(K1,K2) = 0.0D0
				ENDIF
			END DO
		END DO
		V(IND(2),IND(2)) = DCOS(theta)
		V(IND(1),IND(1)) = DCOS(theta)
		V(IND(1),IND(2)) = -DSIN(theta)
		V(IND(2),IND(1)) = DSIN(theta)
		

		C = MATMUL(C, TRANSPOSE(V))
		C = MATMUL(V, C)
		Q = MATMUL(V, Q)
		
	  ENDDO
	  
	  C(1,1) = SQRT(C(1,1))
	  C(2,2) = SQRT(C(2,2))
	  C(3,3) = SQRT(C(3,3))
	  U = MATMUL(C, Q)
	  U = MATMUL(TRANSPOSE(Q), U)
	  
C     	  
	  CALL INVERSE(U, INVU)
	  
	  R = MATMUL(F,INVU) 
	  
	  RETURN
      END
C ################## INVERSE ##############
	  SUBROUTINE INVERSE(M, INVM)
	  
	  INCLUDE 'ABA_PARAM.INC'
	  
	  INTEGER, PARAMETER :: DP = selected_real_kind(p=15, r=307) 
	  REAL(DP) M(3,3), INVM(3,3), J

	  CALL DETERMINANT(M, J)
	  
	  INVM(1, 1) = M(2, 2) * M(3, 3) - M(2, 3) * M(3, 2)
	  INVM(1, 2) = -(M(2, 1) * M(3, 3) - M(2, 3) * M(3, 1))
      INVM(1, 3) = M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1)
      INVM(2, 1) = -(M(1, 2) * M(3, 3) - M(1, 3) * M(3, 2))
      INVM(2, 2) = M(1, 1) * M(3, 3) - M(1, 3) * M(3, 1)
      INVM(2, 3) = -(M(1, 1) * M(3, 2) - M(1, 2) * M(3, 1))
      INVM(3, 1) = M(1, 2) * M(2, 3) - M(1, 3) * M(2, 2)
      INVM(3, 2) = -(M(1, 1) * M(2, 3) - M(1, 3) * M(2, 1))
      INVM(3, 3) = M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)
	  
	  INVM = INVM/J
	  
	  RETURN
      END 
	  
C ################## DETERMINANT( ##############	  
	  SUBROUTINE DETERMINANT(M, J)
	  
	  INCLUDE 'ABA_PARAM.INC'
	  
	  INTEGER, PARAMETER :: DP = selected_real_kind(p=15, r=307) 
	  REAL(DP) M(3,3), J
	  PARAMETER (tol=1.0e-6, ZERO=0.0D0)
	  
	  
	  J = M(1, 1) * (M(2, 2) * M(3, 3) - M(2, 3) * M(3, 2))
     1 - M(1, 2) * (M(2, 1) * M(3, 3) - M(2, 3) * M(3, 1))
     2 + M(1, 3) * (M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1))

	  IF (J.LT.tol) THEN		
		J = ZERO
		PRINT *, 'Exiting program due to error.(DET=0)'
		RETURN
	  END IF
	  
	  RETURN
      END         

C ################## Surrogate model ##############
	  SUBROUTINE NN(DPTS,GRAD)
	  
	  INCLUDE 'ABA_PARAM.INC'
	  
	  INTEGER, PARAMETER :: DP = selected_real_kind(p=15, r=307)
	  PARAMETER (maxN = 256, maxL = 7)
	  
	  REAL(DP) DPTS(6), GRAD(6,6), dN(6,6), INP_MIN(6), INP_MAX(6), OUT_MIN(6), OUT_MAX(6)	
      INTEGER Struc(maxL+2), L
	  REAL W(maxL+1, maxN, maxN+1), Inter_x(maxL+1, maxN, 1), Inter_dx(maxL+1, maxN, maxN), x(maxN, 1), dx(maxN, maxN),  da(maxN, 1)

	  
C#***********************************************************************************#
	  COMMON /surrogate/ INP_MIN, INP_MAX, OUT_MIN, OUT_MAX, W, Struc, L
	  

	  DPTS = (DPTS- INP_MIN) / (INP_MAX- INP_MIN)
	  
      x(:Struc(1),1) = DPTS

	  DO i = 1, L+1
		x(:Struc(i+1),1) = matmul(W(i,:Struc(i+1),:Struc(i)),x(:Struc(i),1)) + W(i,:Struc(i+1),Struc(i)+1)
		Inter_x(i, :Struc(i+1), 1) = x(:Struc(i+1),1)
C		PRINT*, Inter_x(i, :Struc(i+1), 1)
		IF (i.NE.L+1) THEN
			CALL Sigmoid(x(:Struc(i+1),1), Struc(i+1))
C			PRINT*, 'a:', x(:Struc(i+1),1)
		ELSE
C			PRINT*, 'OUTPUT:', x(:Struc(i+1),1)
		END IF
 	  ENDDO
	  
	  DPTS = x(:Struc(i+1),1)
	  
	  DPTS = DPTS* (OUT_MAX - OUT_MIN) + OUT_MIN


C 	  PRINT*, '******************'
	  
 	  DO i = 1, L+1
 		da(:Struc(i+1), 1) = Inter_x(i, :Struc(i+1), 1)
 		CALL d_Sigmoid(da(:Struc(i+1), 1), Struc(i+1))

 		DO j = 1 , Struc(i)
 			Inter_dx(i, :Struc(i+1), j) = W(i,:Struc(i+1), j)*da(:Struc(i+1), 1)
 		ENDDO

 	  ENDDO

	  DO i = 1, Struc(1)
		DO j = 1, Struc(1)
			IF(i.EQ.j) THEN
				dx(i, j) = 1
			ELSE
				dx(i, j) = 0
			END IF
		ENDDO
	  ENDDO

 	  DO i = 1, L
 		dx(:Struc(i+1), :Struc(1)) =  matmul(Inter_dx(i, :Struc(i+1), :Struc(i)), dx(:Struc(i), :Struc(1)))

	  ENDDO
	  dx(:Struc(L+1+1), :Struc(1)) =  matmul(W(L+1, :Struc(L+1+1), :Struc(L+1)), dx(:Struc(L+1), :Struc(1)))
	  
C#***********************************************************************************#	
	  GRAD = dx(:Struc(L+1+1), :Struc(1))
	  
	  DO i = 1, 6
		DO j = 1, 6
			IF(i.EQ.j) THEN
				dN(i, j) = 1 / (INP_MAX(i)- INP_MIN(i))
			ELSE
				dN(i, j) = 0
			END IF
		ENDDO
	  ENDDO
  
	  GRAD = matmul(GRAD, dN)
	  
	  DO i = 1, 6
		DO j = 1, 6
			IF(i.EQ.j) THEN
				dN(i, j) = OUT_MAX(i) - OUT_MIN(i)
			ELSE
				dN(i, j) = 0
			END IF
		ENDDO
	  ENDDO
	  
	  GRAD = matmul(dN, GRAD)
	  

	  RETURN
      END 

C ################## Activation func ##############	 
      SUBROUTINE ReLu(x, n)

      INCLUDE 'aba_param.inc'
      
	  INTEGER n 
	  REAL x(n,1)
	  
	  DO i = 1, n
		IF (x(i,1).GT.0) THEN
			x(i,1) = x(i,1)
		ELSE
		 	x(i,1) = 0
		END IF
	  ENDDO
	  
      RETURN
      END
C ################## Activation func ##############	 
      SUBROUTINE d_ReLu(x, n)

      INCLUDE 'aba_param.inc'
      
	  INTEGER n 
	  REAL x(n,1)
	  
	  DO i = 1, n
		IF (x(i,1).GT.0) THEN
			x(i,1) = 1
		ELSE
		 	x(i,1) = 0
		END IF
	  ENDDO
	  
      RETURN
      END
C ################## Activation func ##############		  
	  SUBROUTINE Sigmoid(x, n)

      INCLUDE 'aba_param.inc'
      
      INTEGER n
      REAL x(n,1)

      DO i = 1, n
         x(i,1) = 1.0 / (1.0 + EXP(-x(i,1)))
      ENDDO
      
      RETURN
      END
C ################## Activation func ##############	
      SUBROUTINE d_Sigmoid(x, n)

      INCLUDE 'aba_param.inc'
      
      INTEGER n
      REAL x(n,1)
      REAL sigmoid_val
      
      DO i = 1, n
         sigmoid_val = 1.0 / (1.0 + EXP(-x(i,1))) 
         x(i,1) = sigmoid_val * (1.0 - sigmoid_val)
      ENDDO
      
      RETURN
      END
C ################## matrix_to_tensor ##############  	  
      SUBROUTINE get_mn_index(i, m, n)
      
	  INCLUDE 'aba_param.inc'
      
	  INTEGER i, m, n

      SELECT CASE (i)
		CASE (1)
			m = 1
			n = 1
		CASE (2)
			m = 1
			n = 2
		CASE (3)
			m = 1
			n = 3
		CASE (4)
			m = 2
			n = 2
		CASE (5)
			m = 2
			n = 3
		CASE (6)
			m = 3
			n = 3
      END SELECT
	
	  RETURN
      END