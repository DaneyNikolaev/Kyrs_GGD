Program Main
  include 'omp_lib.h'

  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt', ResFile='Residuals.dat', ProfileFile='Profile.dat'
  character MeshFile*30,ctmp        ! name of file with computational mesh
  integer, parameter:: IO = 12 ! input-output unit
  integer :: i,j,NI,NJ, NT, IP, JP, Scheme
  real,allocatable,dimension(:,:):: X,Y,CellVolume ! scalar arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector !geometry vector arrays
  real,allocatable,dimension(:,:):: P, T, Ro, M, U, V
  real :: CFL, Pin, Uin, Vin, Tin, gamma, Cp, Rm
  
  
!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  READ(IO,*) Pin     ! read key of method for calculating gradient
  READ(IO,*) Uin    ! scheme for scalar
  READ(IO,*) Vin        ! velocity scale
  READ(IO,*) Tin        ! length scale
  READ(IO,*) gamma      ! Reynolds number
  READ(IO,*) Cp        ! Prandtl number
  READ(IO,*) Rm       ! CFL
  READ(IO,*) CFL       ! VNM
  READ(IO,*) NT     ! number of iterations
  READ(IO,*) Scheme
  CLOSE(IO)
!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ
  WRITE(*,*) 'Pin = ',Pin
  WRITE(*,*) 'Uin = ',Uin,'Vin = ', Vin
  WRITE(*,*) 'Tin = ',Tin,'gamma = ', gamma,'Cp = ', Cp,'Rm = ', Rm
  WRITE(*,*) 'CFL = ',CFL,'NT = ', NT
  WRITE(*,*) 'Scheme = ',Scheme
  
  

!=== ALLOCATE ALL ARRAYS ===
  !WRITE(*,*) 'Allocate arrays'
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(T(0:NI,0:NJ))
  allocate(U(0:NI,0:NJ))
  allocate(V(0:NI,0:NJ))
  allocate(Ro(0:NI,0:NJ))
  allocate(M(0:NI,0:NJ))


!===  READ GRID ===
  !WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J),rtmp,I=1,NI),J=1,NJ)
  CLOSE(IO)



!=== CALCULATE METRIC ===
  !WRITE(*,*) 'Calculate metric'
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)






  !Ro = 0.0
  !M = 0.0
  
	Call B_Euler(NI, NJ, IFaceVector, JFaceVector, CellVolume, CellCenter, IFaceCenter, JFaceCenter, &
P, U, V, T, gamma, Cp, Rm, CFL, NT, Pin, Uin, Vin, Tin, Scheme, X, Y)
  Ro(:,:) = P(:,:)/(Rm*T(:,:))
  M(:,:) = sqrt((U(:,:)*U(:,:)+V(:,:)*V(:,:))/(gamma*Rm*T(:,:)))
  
  
  
  
  
  !=== OUTPUT FIELDS ===
  !WRITE(*,*) 'Output fields to file: ', OutputFile
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,U,V,T,Ro,M)
  Close(IO)
END PROGRAM Main



SUBROUTINE B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 
  REAL X(NI,NJ),Y(NI,NJ),&                               ! input: nodes coordinates
       CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   !output: cell centers and volumes
       IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& !        face centers and vectors for I-faces
       JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)   !        face centers and vectors for J-faces
  REAL r(2)

  !=== FACE CENTERS AND FACE VECTORS ===
  ! I-DIRECTION
  DO J = 1,NJ-1
    DO I = 1,NI
      r(1) = X(I,J+1) - X(I,J)  ! r = vector from one node to another
      r(2) = Y(I,J+1) - Y(I,J)
      IFaceVector(I,J,1) = r(2) ! IFaceVector = r rotated on 90 degree
      IFaceVector(I,J,2) =-r(1) ! IFaceVector directed to increasing I-index
      IFaceCenter(I,J,1) = 0.5*(X(i,j)+x(i,j+1))
      IFaceCenter(I,J,2) = 0.5*(Y(i,j)+Y(i,j+1))
    ENDDO
  ENDDO

  ! J-DIRECTION
  DO J = 1,NJ
    DO I = 1,NI-1
      r(1) = X(I+1,J) - X(I,J)  ! r = vector from one node to another
      r(2) = Y(I+1,J) - Y(I,J)
      JFaceVector(I,J,1) =-r(2) ! JFaceVector = r rotated on -90 degree
      JFaceVector(I,J,2) = r(1) ! JFaceVector directed to increasing J-index 
      JFaceCenter(I,J,1) = 0.5*(X(i,j)+x(i+1,j))
      JFaceCenter(I,J,2) = 0.5*(Y(i,j)+Y(i+1,j))
    ENDDO
  ENDDO


 !=== CELL VOLUMES ===
  DO J = 1,NJ-1
    DO I = 1,NI-1
      r(1)=X(I+1,J+1) - X(I,J)
      r(2)=Y(I+1,J+1) - Y(I,J)
      CellVolume(I,J) = 0.5*DOT_PRODUCT(IFaceVector(I,J,:),r)& ! sum surfaces of two triangles
                      + 0.5*DOT_PRODUCT(JFaceVector(I,J,:),r)
    ENDDO
  ENDDO


  !=== CELL CENTERS ===
  ! FOR INNER CELLS: CENTER OF CONTOUR (sum of FaceCenter*FaceLength/Perimeter)
  DO J = 1,NJ-1
    DO  I = 1,NI-1
      CellCenter(I,J,:) = ( IFaceCenter(I  ,J,:)*Norm2(IFaceVector(I  ,J,:))+&
                            IFaceCenter(I+1,J,:)*Norm2(IFaceVector(I+1,J,:))+&
                            JFaceCenter(I,J  ,:)*Norm2(JFaceVector(I,J  ,:))+&
                            JFaceCenter(I,J+1,:)*Norm2(JFaceVector(I,J+1,:)) )&
                         /( Norm2(IFaceVector(I,J,:))+Norm2(IFaceVector(I+1,J,:))+&
                            Norm2(JFaceVector(I,J,:))+Norm2(JFaceVector(I,J+1,:)) )
    ENDDO
  ENDDO

  ! FOR DUMMY CELLS ON BOUNDARIES: CELL CENTER = FACE CENTER
  ! I-BOUNDARIES -----------------------------------------------------
  DO NBOUND = 1,2
    IF (NBOUND.EQ.1) THEN
      IBOUND =  1; IOUT =  0
    ELSE 
      IBOUND = NI; IOUT =  NI
    ENDIF
    DO J = 1,NJ-1
      CellCenter(IOUT,J,:) = IFaceCenter(IBOUND,J,:)
    ENDDO
  ENDDO

  ! J-BOUNDARIES -----------------------------------------------------
  DO NBOUND = 1,2
    IF (NBOUND.EQ.1) THEN
      JBOUND = 1;  JOUT =  0
    ELSE 
      JBOUND = NJ; JOUT =  NJ
    ENDIF
    DO  I = 1,NI-1
      CellCenter(I,JOUT,:) = JFaceCenter(I,JBOUND,:) 
    ENDDO
  ENDDO

END SUBROUTINE




Subroutine B_Euler(NI, NJ, IFaceVector, JFaceVector, CellVolume, CellCenter, IFaceCenter, JFaceCenter,&
P, U, V, T, gamma, Cp, Rm, CFL, NT, Pin, Uin, Vin, Tin, Scheme, X, Y)

Integer Scheme, NI, NJ, IP, JP
Real,Dimension(0:NI,0:NJ):: P, U, V, T
Real IFaceVector(NI, NJ-1, 2), JFaceVector(NI-1, NJ, 2), CellVolume(NI-1, NJ-1),&
CellCenter(0:NI,0:NJ,2), IFaceCenter( NI,NJ-1,2), JFaceCenter( NI-1,NJ,2), X(NI,NJ),Y(NI,NJ)
Real gamma, Cp, Rm, dt, CFL, Pin, Uin, Vin, Tin, Roin
Real, dimension(0:NI,0:NJ) :: ro, roU, roV, roE
Real, dimension(0:NI,0:NJ) :: ro1, roU1, roV1, roE1
Real, dimension(0:NI,0:NJ) :: Res1, Res2, Res3, Res4
Real Flux(4), qL(4), qR(4), SF(2)
Real VN, Vf(2), SS, VNL
Real Pr, Hl, mod_Sf
Real q_im1(4), q_im2(4), q_ip1(4), q_ip2(4)
Real d(2), dn(2)
integer ::IP2, JP2, P_X_1, P_X_2, P_Y_1, P_Y_2

  P = Pin
  U = Uin
  V = Vin
  T = Tin
  Ro_in = Pin/Rm/Tin
  
  DO I=1, NI-1
	DO J=1, NJ-1
		ro(I,J) = P(I,J)/(Rm*T(I,J))
		roU(I,J) = ro(I,J)*U(I,J)
		roV(I,J) = ro(I,J)*V(I,J)
		roE(I,J) =(U(I,J)*U(I,J) + V(I,J)*V(I,J))/2.0*ro(I,J) + P(I,J)/(gamma - 1)
	end do
  end do

	call C_Location (27.0, 3.0, NI, NJ, X, Y, CellVolume, IP2, JP2)
	IP = IP2
	JP = JP2
	OPEN(11, FILE = 'MonitorFile1.plt')
	OPEN(12,FILE = 'Residuals1.dat')
	OPEN(10, FILE = 'ProfileFile1.plt')
	write(11, *) 'Variables = "Iterations","U"'
	
  do k =1, NT
  
  Res1(:,:) = 0.0
  Res2(:,:) = 0.0
  Res3(:,:) = 0.0
  Res4(:,:) = 0.0 
  
!I-DIRECTION
  do J = 1, NJ-1
	do I = 2, NI-1
		
		IF (I.EQ.2) THEN
			d(1) = 2.0
			dn(1) = 0.0
		ELSE
			d(1) = norm2(CellCenter(I-1,J,:) - CellCenter(I-2,J,:))
			dn(1) = norm2(IFaceCenter(I,J,:) - CellCenter(I-1,J,:))
		ENDIF
		
		IF (I.EQ.NI-1) THEN
			d(2) = 2.0
			dn(2) = 0.0
		ELSE
			d(2) = norm2(CellCenter(I+1,J,:) - CellCenter(I,J,:))
			dn(2) = norm2(IFaceCenter(I,J,:) - CellCenter(I,J,:))
		ENDIF
		!write(*,*) I, J, dn(1)/d(1), dn(2)/d(2)
		!write(*,*) JFaceCenter(1,3,:)		
				
		IF (I.EQ.2) THEN
			q_im2(1) = P(I-1,J)
			q_im2(2) = U(I-1,J)
			q_im2(3) = V(I-1,J)
			q_im2(4) = T(I-1,J)
			! q_im2(1) = 2*P(I-1,J) - P(I,J)
			! q_im2(2) = 2*U(I-1,J) - U(I,J)
			! q_im2(3) = 2*V(I-1,J) - V(I,J)
			! q_im2(4) = 2*T(I-1,J) - T(I,J)
		ELSE
			q_im2(1) = P(I-2,J)
			q_im2(2) = U(I-2,J)
			q_im2(3) = V(I-2,J)
			q_im2(4) = T(I-2,J)
		ENDIF
		
		q_im1(1) = P(I-1,J)
		q_im1(2) = U(I-1,J)
		q_im1(3) = V(I-1,J)
		q_im1(4) = T(I-1,J)		
		
		q_ip1(1) = P(I,J)
		q_ip1(2) = U(I,J)
		q_ip1(3) = V(I,J)
		q_ip1(4) = T(I,J)
		
		IF (I.EQ.NI-1) THEN
			q_ip2(1) = P(I,J)
			q_ip2(2) = U(I,J)
			q_ip2(3) = V(I,J)
			q_ip2(4) = T(I,J)
			! q_ip2(1) = 2*P(I,J) - P(I-1,J)
			! q_ip2(2) = 2*U(I,J) - U(I-1,J)
			! q_ip2(3) = 2*V(I,J) - V(I-1,J)
			! q_ip2(4) = 2*T(I,J) - T(I-1,J)
		ELSE
			q_ip2(1) = P(I+1,J)
			q_ip2(2) = U(I+1,J)
			q_ip2(3) = V(I+1,J)
			q_ip2(4) = T(I+1,J)
		ENDIF
		
		call CALC_Variables(q_im1, q_im2, q_ip1, q_ip2, d, dn, qL, qR)
		
		SF(:) = IFaceVector(I,J,:)
		
		call CALC_FLUX(SF, qR, qL, gamma,Cp, Rm, FLUX, Scheme)
		
		Res1(I, J) = Res1(I, J) - FLUX(1)
        Res2(I, J) = Res2(I, J) - FLUX(2)
	    Res3(I, J) = Res3(I, J) - FLUX(3)
		Res4(I, J) = Res4(I, J) - FLUX(4)
		
		Res1(I-1, J) = Res1(I-1, J) + FLUX(1)
        Res2(I-1, J) = Res2(I-1, J) + FLUX(2)
	    Res3(I-1, J) = Res3(I-1, J) + FLUX(3)
		Res4(I-1, J) = Res4(I-1, J) + FLUX(4)
	end do
  end do
  
!J-DIRECTION
    do J = 2, NJ-1
	 do I = 1, NI-1
	 
		IF (J.EQ.2) THEN
			d(1) = 2.0
			dn(1) = 0.0
		ELSE
			d(1) = norm2(CellCenter(I,J-1,:) - CellCenter(I,J-2,:))
			dn(1) = norm2(JFaceCenter(I,J,:) - CellCenter(I,J-1,:))
		ENDIF
		
		IF (J.EQ.NJ-1) THEN
			d(2) = 2.0
			dn(2) = 0.0
		ELSE
			d(2) = norm2(CellCenter(I,J+1,:) - CellCenter(I,J,:))
			dn(2) = norm2(JFaceCenter(I,J,:) - CellCenter(I,J,:))
		ENDIF
		!write(*,*) I, J, dn(1)/d(1), dn(2)/d(2)
		
		IF (J.EQ.2) THEN
			q_im2(1) = P(I,J-1)
			q_im2(2) = U(I,J-1)
			q_im2(3) = V(I,J-1)
			q_im2(4) = T(I,J-1)
			! q_im2(1) = 2*P(I,J-1) - P(I,J)
			! q_im2(2) = 2*U(I,J-1) - U(I,J)
			! q_im2(3) = 2*V(I,J-1) - V(I,J)
			! q_im2(4) = 2*T(I,J-1) - T(I,J)
		ELSE
			q_im2(1) = P(I,J-2)
			q_im2(2) = U(I,J-2)
			q_im2(3) = V(I,J-2)
			q_im2(4) = T(I,J-2)
		ENDIF
		
		q_im1(1) = P(I,J-1)
		q_im1(2) = U(I,J-1)
		q_im1(3) = V(I,J-1)
		q_im1(4) = T(I,J-1)
		
		
		q_ip1(1) = P(I,J)
		q_ip1(2) = U(I,J)
		q_ip1(3) = V(I,J)
		q_ip1(4) = T(I,J)
		
		IF (J.EQ.NJ-1) THEN
			 q_ip2(1) = P(I,J)
			 q_ip2(2) = U(I,J)
			 q_ip2(3) = V(I,J)
			 q_ip2(4) = T(I,J)
			! q_ip2(1) = 2*P(I,J) - P(I,J-1)
			! q_ip2(2) = 2*U(I,J) - U(I,J-1)
			! q_ip2(3) = 2*V(I,J) - V(I,J-1)
			! q_ip2(4) = 2*T(I,J) - T(I,J-1)
		ELSE
			q_ip2(1) = P(I,J+1)
			q_ip2(2) = U(I,J+1)
			q_ip2(3) = V(I,J+1)
			q_ip2(4) = T(I,J+1)
		ENDIF
		
		call CALC_Variables(q_im1, q_im2, q_ip1, q_ip2, d, dn, qL, qR)
		
		SF(:) = JFaceVector(I,J,:)
		
		call CALC_FLUX(SF, qR, qL, gamma,Cp, Rm, FLUX, Scheme)
		
		Res1(I, J) = Res1(I, J) - FLUX(1)
        Res2(I, J) = Res2(I, J) - FLUX(2)
	    Res3(I, J) = Res3(I, J) - FLUX(3)
		Res4(I, J) = Res4(I, J) - FLUX(4)
		
		Res1(I, J-1) = Res1(I, J-1) + FLUX(1)
        Res2(I, J-1) = Res2(I, J-1) + FLUX(2)
	    Res3(I, J-1) = Res3(I, J-1) + FLUX(3)
		Res4(I, J-1) = Res4(I, J-1) + FLUX(4)
	end do
  end do


	do J = 1, NJ-1
	
!==================Inlet=========================
	SF(:) = IFaceVector(1,J,:)
	mod_Sf = sqrt(SF(1)**2 + SF(2)**2)
	
	FLUX(1) = Ro_in*Uin*mod_Sf
	FLUX(2) = Ro_in*Uin**2*mod_Sf + Sf(1)*Pin
	FLUX(3) = Sf(2)*Pin
	FLUX(4) = Ro_in*Uin*(Cp*Tin + Uin**2/2.0)*mod_Sf
	
	Res1(1, J) = Res1(1, J) - FLUX(1)
    Res2(1, J) = Res2(1, J) - FLUX(2)
	Res3(1, J) = Res3(1, J) - FLUX(3)
	Res4(1, J) = Res4(1, J) - FLUX(4)

!==================Outlet=========================
	SF(:) = IFaceVector(NI,J,:)

	Vn = U(NI-1,J)*Sf(1) + V(NI-1,J)*Sf(2)
	
	FLUX(1) = ro(NI-1,J)*Vn
	FLUX(2) = ro(NI-1,J)*U(NI-1,J)*Vn + Sf(1)*P(NI-1,J)
	FLUX(3) = ro(NI-1,J)*V(NI-1,J)*Vn + Sf(2)*P(NI-1,J)
	FLUX(4) = ro(NI-1,J)*Vn*(Cp*T(NI-1,J) + (U(NI-1,J)*U(NI-1,J) + V(NI-1,J)*V(NI-1,J))/2.0)
		
	Res1(NI-1, J) = Res1(NI-1, J) + FLUX(1)
    Res2(NI-1, J) = Res2(NI-1, J) + FLUX(2)
	Res3(NI-1, J) = Res3(NI-1, J) + FLUX(3)
	Res4(NI-1, J) = Res4(NI-1, J) + FLUX(4)
	end do
	
	
	
	do I = 1, NI-1
	
!==================Wall=========================
	SF(:) = JFaceVector(I,1,:)
	
	FLUX(1) = 0.0
	FLUX(2) = SF(1)*P(I,1)
	FLUX(3) = SF(2)*P(I,1)
	FLUX(4) = 0.0
	
	Res1(I, 1) = Res1(I, 1) - FLUX(1)
    Res2(I, 1) = Res2(I, 1) - FLUX(2)
	Res3(I, 1) = Res3(I, 1) - FLUX(3)
	Res4(I, 1) = Res4(I, 1) - FLUX(4)

!==================Outlet=========================
	SF(:) = JFaceVector(I,NJ,:)
	
	! FLUX(1) = 0.0
	! FLUX(2) = SF(1)*P(I,NJ-1)
	! FLUX(3) = SF(2)*P(I,NJ-1)
	! FLUX(4) = 0.0
	
	Vn = U(I,NJ-1)*Sf(1) + V(I,NJ-1)*Sf(2)
	
	FLUX(1) = ro(I,NJ-1)*Vn
	FLUX(2) = ro(I,NJ-1)*U(I,NJ-1)*Vn + Sf(1)*P(I,NJ-1)
	FLUX(3) = ro(I,NJ-1)*V(I,NJ-1)*Vn + Sf(2)*P(I,NJ-1)
	FLUX(4) = ro(I,NJ-1)*Vn*(Cp*T(I,NJ-1) + (U(I,NJ-1)*U(I,NJ-1) + V(I,NJ-1)*V(I,NJ-1))/2.0)
		
	Res1(I,NJ-1) = Res1(I,NJ-1) + FLUX(1)
    Res2(I,NJ-1) = Res2(I,NJ-1) + FLUX(2)
	Res3(I,NJ-1) = Res3(I,NJ-1) + FLUX(3)
	Res4(I,NJ-1) = Res4(I,NJ-1) + FLUX(4)
	end do
	



  DO I=1, NI-1
	DO J=1, NJ-1
	
	dt = CFL*sqrt(CellVolume(I, J))/(sqrt(gamma*Rm*Tin) + max(U(I,J),V(I,J)))
	 
	ro1(I,J) = ro(I,J) - dt*Res1(I,J)/CellVolume(I, J)
	roU1(I,J) = roU(I,J) - dt*Res2(I,J)/CellVolume(I, J)
	roV1(I,J) = roV(I,J) - dt*Res3(I,J)/CellVolume(I, J)
	roE1(I,J) =roE(I,J) - dt*Res4(I,J)/CellVolume(I, J)
	
	end do
  end do
  
	ro(1:NI-1,1:NJ-1) = ro1(1:NI-1,1:NJ-1)
	roU(1:NI-1,1:NJ-1) = roU1(1:NI-1,1:NJ-1)
	roV(1:NI-1,1:NJ-1) = roV1(1:NI-1,1:NJ-1)
	roE(1:NI-1,1:NJ-1) =roE1(1:NI-1,1:NJ-1)
	
	
	DO I=1, NI-1
		DO J=1, NJ-1
		
			U(I,J) = roU(I,J)/ro(I,J)
			V(I,J) = roV(I,J)/ro(I,J)
			P(I,J) = (gamma-1)*(roE(I,J) - ro(I,J)*(U(I,J)*U(I,J) + V(I,J)*V(I,J))/2.0)
			T(I,J) = P(I,J)/(ro(I,J)*Rm)

		end do
	end do
  
	write(*,*) k, maxval(Res1(1:NI-1,1:NJ-1)),maxval(Res2(1:NI-1,1:NJ-1)),maxval(Res3(1:NI-1,1:NJ-1)),maxval(Res4(1:NI-1,1:NJ-1))
	write(12,*) k, maxval(Res1(1:NI-1,1:NJ-1)),maxval(Res2(1:NI-1,1:NJ-1)),maxval(Res3(1:NI-1,1:NJ-1)),maxval(Res4(1:NI-1,1:NJ-1))
	
	write(11, *) k, U(IP,JP)
  end do
  
	CLOSE(12) 
	
	call C_Location (35.0, 8.5, NI, NJ, X, Y, CellVolume, IP2, JP2)
	P_X_1 = IP2
	P_Y_1 = JP2
	call C_Location (43.0, 15.0, NI, NJ, X, Y, CellVolume, IP2, JP2)
	P_X_2 = IP2
	P_Y_2 = JP2
	write(10, *) 'Variables = "Y_P1", "U_P1", "T_P1", "Y_P2", "U_P2", "T_P2"'
		DO J=1, NJ-1
			write(10, *) Y(P_X_1, J), U(P_X_1, J), T(P_X_1,J), Y(P_X_2, J), U(P_X_2, J), T(P_X_2,J)
		end do

End Subroutine



Subroutine CALC_FLUX(SF, qR, qL, gamma,Cp, Rm, FLUX, Scheme)

Integer Scheme
Real gamma, Cp, Rm
Real Flux(4), qL(4), qR(4), SF(2)
Real nf(2), modSF, Vn, Vn_F, Ro_L, Ro_R, H_L, H_R, E_L, E_R
Real cL, cR, sL, sR, S_cs, p_csL, p_csR
Real Mach_L, Mach_R, p_plus, p_minus, p_F, M_F, M_minus, M_plus, FLUX_C(4), roc, Koef, w_cs(4)
Real u_roy, H_roy, c_roy, alpha, betta, c_F

modSF = sqrt(SF(1)**2 + SF(2)**2)
nf(1) = SF(1)/modSF
nf(2) = SF(2)/modSF

!Vn_R = qR(2)*Sf(1) + qR(3)*Sf(2)

Vn_L = qL(2)*nf(1) + qL(3)*nf(2)
Vn_R = qR(2)*nf(1) + qR(3)*nf(2)
Vn_F = Vn_L + Vn_R

Ro_L = qL(1)/(Rm*qL(4))
Ro_R = qR(1)/(Rm*qR(4))
H_L = (qL(2)**2 + qL(3)**2)/2.0 + Cp*qL(4)
H_R = (qR(2)**2 + qR(3)**2)/2.0 + Cp*qR(4)

SELECT CASE(Scheme)

	CASE(1)!Схема, в котрой скорости берутся против потока, а давление по потоку
	
		IF (Vn_F.GT.0) THEN
		
			FLUX(1) = Ro_L*Vn_L*modSF
			FLUX(2) = (Ro_L*qL(2)*Vn_L + nf(1)*qR(1))*modSF
			FLUX(3) = (Ro_L*qL(3)*Vn_L + nf(2)*qR(1))*modSF
			FLUX(4) = Ro_L*Vn_L*H_L*modSF
		ELSE
	  
			FLUX(1) = Ro_R*Vn_R*modSF
			FLUX(2) = (Ro_R*qR(2)*Vn_R + nf(1)*qL(1))*modSF
			FLUX(3) = (Ro_R*qR(3)*Vn_R + nf(2)*qL(1))*modSF
			FLUX(4) = Ro_R*Vn_R*H_R*modSF
		ENDIF
	
	
	
	
	CASE(2)!HLLC
		
		cL = sqrt(Cp*(gamma - 1)*qL(4))
		cR = sqrt(Cp*(gamma - 1)*qR(4))
		
		!sL = minval([Vn_L - cL,Vn_R - cR])
		!sR = maxval([Vn_L + cL,Vn_R + cR])
		
		u_roy = (sqrt(Ro_L)*Vn_L + sqrt(Ro_R)*Vn_R)/(sqrt(Ro_L) + sqrt(Ro_R))
		H_roy = (sqrt(Ro_L)*H_L + sqrt(Ro_R)*H_R)/(sqrt(Ro_L) + sqrt(Ro_R))
		c_roy = sqrt((gamma - 1)*(H_roy - 0.5*u_roy**2.0))
		
		sL = minval([Vn_L - cL, u_roy - c_roy])
		sR = maxval([Vn_R + cR, u_roy + c_roy])

		
		S_cs = ((qR(1) - qL(1)) + Ro_L*Vn_L*(sL - Vn_L) - Ro_R*Vn_R*(sR - Vn_R))/(Ro_L*(sL - Vn_L) - Ro_R*(sR - Vn_R))
		
		IF (sL.GE.0) THEN
			
			FLUX(1) = Ro_L*Vn_L*modSF
			FLUX(2) = (Ro_L*qL(2)*Vn_L + nf(1)*qL(1))*modSF
			FLUX(3) = (Ro_L*qL(3)*Vn_L + nf(2)*qL(1))*modSF
			FLUX(4) = Ro_L*Vn_L*H_L*modSF
			
		ELSE IF (sL.LE.0.AND.S_cs.GE.0) THEN
		
			E_L = (Cp - Rm)*qL(4) + (qL(2)**2 + qL(3)**2)/2.0
			p_csL = qL(1) + Ro_L*(sL - Vn_L)*(S_cs - Vn_L)
	  
			FLUX(1) = (S_cs*(sL*Ro_L - Ro_L*Vn_L) + sL*p_csl*0.0)/(sL - S_cs)*modSF
			FLUX(2) = (S_cs*(sL*Ro_L*qL(2) - (Ro_L*qL(2)*Vn_L + nf(1)*qL(1))) + sL*p_csl*nF(1))/(sL - S_cs)*modSF
			FLUX(3) = (S_cs*(sL*Ro_L*qL(3) - (Ro_L*qL(3)*Vn_L + nf(2)*qL(1))) + sL*p_csl*nF(2))/(sL - S_cs)*modSF
			FLUX(4) = (S_cs*(sL*Ro_L*E_L - Ro_L*Vn_L*H_L) + sL*p_csl*S_cs)/(sL - S_cs)*modSF
			
		ELSE IF (S_cs.LE.0.AND.sR.GE.0) THEN
		
			E_R = (Cp - Rm)*qR(4) + (qR(2)**2 + qR(3)**2)/2.0
			p_csR = qR(1) + Ro_R*(sR - Vn_R)*(S_cs - Vn_R)
	  
			FLUX(1) = (S_cs*(sR*Ro_R - Ro_R*Vn_R) + sR*p_csR*0.0)/(sR - S_cs)*modSF
			FLUX(2) = (S_cs*(sR*Ro_R*qR(2) - (Ro_R*qR(2)*Vn_R + nf(1)*qR(1))) + sR*p_csR*nF(1))/(sR - S_cs)*modSF
			FLUX(3) = (S_cs*(sR*Ro_R*qR(3) - (Ro_R*qR(3)*Vn_R + nf(2)*qR(1))) + sR*p_csR*nF(2))/(sR - S_cs)*modSF
			FLUX(4) = (S_cs*(sR*Ro_R*E_R - Ro_R*Vn_R*H_R) + sR*p_csR*S_cs)/(sR - S_cs)*modSF
		
		ELSE IF (sR.LE.0) THEN	
		
			FLUX(1) = Ro_R*Vn_R*modSF
			FLUX(2) = (Ro_R*qR(2)*Vn_R + nf(1)*qR(1))*modSF
			FLUX(3) = (Ro_R*qR(3)*Vn_R + nf(2)*qR(1))*modSF
			FLUX(4) = Ro_R*Vn_R*H_R*modSF
		
		ENDIF
	


	
	CASE(3)!AUSM
	
		Mach_L = Vn_L/sqrt(gamma*Rm*qL(4))
		Mach_R = Vn_R/sqrt(gamma*Rm*qR(4))
		
		IF (abs(Mach_L).GT.1.0) THEN
			p_plus = 1.0/2.0*(Mach_L + abs(Mach_L))/Mach_L
		ELSE
			p_plus = 1.0/4.0*(Mach_L + 1.0)**2.0*(2.0 - Mach_L)
		ENDIF
		
		IF (abs(Mach_R).GT.1.0) THEN
			p_minus = 1.0/2.0*(Mach_R - abs(Mach_R))/Mach_R
		ELSE 
			p_minus = 1.0/4.0*(Mach_R - 1.0)**2.0*(2.0 + Mach_R)
		ENDIF
		
		p_F = p_plus*qL(1) + p_minus*qR(1)
		
		IF (abs(Mach_L).GT.1.0) THEN
			M_plus = 0.5*(Mach_L + abs(Mach_L))
		ELSE 
			M_plus = 0.25*(Mach_L + 1.0)**2! + 0.125*(Mach_L**2 - 1.0)**2
		ENDIF
		
		IF (abs(Mach_R).GT.1.0) THEN
			M_minus = 0.5*(Mach_R - abs(Mach_R))
		ELSE 
			M_minus = -0.25*(Mach_R - 1.0)**2! - 0.125*(Mach_R**2 - 1.0)**2
		ENDIF
		
		M_F = M_plus + M_minus
		
		IF (M_F.GE.0.0) THEN
			
			roc = qL(1)/(Rm*qL(4))*sqrt(gamma*Rm*qL(4))
			
			FLUX_C(1) = M_F*roc
			FLUX_C(2) = M_F*roc*qL(2)
			FLUX_C(3) = M_F*roc*qL(3)
			FLUX_C(4) = M_F*roc*(Cp*qL(4) + (qL(2)**2 + qL(3)**2)/2.0)
		ELSE 
			
			roc = qR(1)/(Rm*qR(4))*sqrt(gamma*Rm*qR(4))
			
			FLUX_C(1) = M_F*roc
			FLUX_C(2) = M_F*roc*qR(2)
			FLUX_C(3) = M_F*roc*qR(3)
			FLUX_C(4) = M_F*roc*(Cp*qR(4) + (qR(2)**2 + qR(3)**2)/2.0)
		ENDIF
		
			FLUX(1) = FLUX_C(1)*modSF
			FLUX(2) = (FLUX_C(2) + p_F*nF(1))*modSF
			FLUX(3) = (FLUX_C(3) + p_F*nF(2))*modSF
			FLUX(4) = FLUX_C(4)*modSF
	


	
	case(4)!AUSM+
	
		alpha = 3.0/16.0
		betta = 1.0/8.0
	
		c_F = sqrt(sqrt(gamma*Rm*qL(4))*sqrt(gamma*Rm*qR(4)))
		
		Mach_L = Vn_L/c_F
		Mach_R = Vn_R/c_F
		
		IF (abs(Mach_L).GE.1.0) THEN
			p_plus = 0.5*(Mach_L + abs(Mach_L))/Mach_L
		ELSE 
			p_plus = 0.25*(Mach_L + 1.0)**2*(2.0 - Mach_L) + alpha*Mach_L*(Mach_L**2 - 1.0)**2
		ENDIF
		
		IF (abs(Mach_R).GE.1.0) THEN
			p_minus = 0.5*(Mach_R - abs(Mach_R))/Mach_R
		ELSE 
			p_minus = 0.25*(Mach_R - 1.0)**2*(2.0 + Mach_R) - alpha*Mach_R*(Mach_R**2 - 1.0)**2
		ENDIF
		
		p_F = p_plus*qL(1) + p_minus*qR(1)
		
		IF (abs(Mach_L).GT.1.0) THEN
			M_plus = 0.5*(Mach_L + abs(Mach_L))
		ELSE 
			M_plus = 0.25*(Mach_L + 1.0)**2 + betta*(Mach_L**2 - 1.0)**2
		ENDIF
		
		IF (abs(Mach_R).GT.1.0) THEN
			M_minus = 0.5*(Mach_R - abs(Mach_R))
		ELSE 
			M_minus = -0.25*(Mach_R - 1.0)**2 - betta*(Mach_R**2 - 1.0)**2
		ENDIF
		
		M_F = M_plus + M_minus
		
		IF (M_F.GE.0.0) THEN
			
			roc = qL(1)/(Rm*qL(4))*c_F
			
			FLUX_C(1) = M_F*roc
			FLUX_C(2) = M_F*roc*qL(2)
			FLUX_C(3) = M_F*roc*qL(3)
			FLUX_C(4) = M_F*roc*(Cp*qL(4) + (qL(2)**2 + qL(3)**2)/2.0)
		ELSE 
			
			roc = qR(1)/(Rm*qR(4))*c_F
			
			FLUX_C(1) = M_F*roc
			FLUX_C(2) = M_F*roc*qR(2)
			FLUX_C(3) = M_F*roc*qR(3)
			FLUX_C(4) = M_F*roc*(Cp*qR(4) + (qR(2)**2 + qR(3)**2)/2.0)
		ENDIF
		
			FLUX(1) = FLUX_C(1)*modSF
			FLUX(2) = (FLUX_C(2) + p_F*nF(1))*modSF
			FLUX(3) = (FLUX_C(3) + p_F*nF(2))*modSF
			FLUX(4) = FLUX_C(4)*modSF
END SELECT
	
End Subroutine


Subroutine CALC_Variables(q_im1, q_im2, q_ip1, q_ip2, d, dn, qL, qR)

Real q_im1(4), q_im2(4), q_ip1(4), q_ip2(4)
Real qL(4), qR(4)
Real r_L(4),r_R(4), TVD_L(4), TVD_R(4)
Real d(2), dn(2)

DO I=1,4


	IF (abs(q_im1(I) - q_im2(I)).LE.1e-10) THEN
		r_L(I) = 1.0
	ELSE
		r_L(I) = (q_ip1(I) - q_im1(I))/(q_im1(I) - q_im2(I))
	ENDIF
	
	IF (r_L(I).LT.0.0) THEN
		TVD_L(I) = 0.0
		qL(I) = q_im1(I)
	ELSE
		TVD_L(I) = (r_L(I)**2 + r_L(I))/(r_L(I)**2 + 1.0)
		!TVD_L(I) = min(1.0, r_L(I))
		qL(I) = q_im1(I) + dn(1)/d(1)*TVD_L(I)*(q_im1(I) - q_im2(I))
		!qL(I) = q_im1(I) + 0.5*TVD_L(I)*(q_im1(I) - q_im2(I))
	ENDIF
	
	
	IF (abs(q_ip2(I) - q_ip1(I)).LE.1e-10) THEN
		r_R(I) = 1.0
	ELSE
		r_R(I) = (q_ip1(I) - q_im1(I))/(q_ip2(I) - q_ip1(I))
	ENDIF
	
	IF (r_R(I).LT.0.0) THEN
		TVD_R(I) = 0.0
		qR(I) = q_ip1(I)
	ELSE
		TVD_R(I) = (r_R(I)**2 + r_R(I))/(r_R(I)**2 + 1.0)
		!TVD_R(I) = min(1.0, r_R(I))
		qR(I) = q_ip1(I) + dn(2)/d(2)*TVD_R(I)*(q_ip1(I) - q_ip2(I))
		!qR(I) = q_ip1(I) - 0.5*TVD_R(I)*(q_ip2(I) - q_ip1(I))
	ENDIF
	
ENDDO
End Subroutine

Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,U,V,T,Ro,M)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P,U,V,T,Ro,M

  !Write(IO,*) 'VARIABLES = "X","Y","U","V","T","Ro","M"'!,"P"'
  Write(IO,*) 'VARIABLES = "X","Y","U","V","T","Ro","M","P"'
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F20.8)') X(1:NI,1:NJ)
  Write(IO,'(100F20.8)') Y(1:NI,1:NJ)
  Write(IO,'(100F20.8)') U(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') V(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') T(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') Ro(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') M(1:NI-1,1:NJ-1)
  Write(IO,'(100F20.8)') P(1:NI-1,1:NJ-1)
End Subroutine



Subroutine C_Location (x_m,y_m,NI,NJ,X,Y,CellVolume,Ip,Jp)

  REAL  X(NI,NJ),Y(NI,NJ),&							! input: nodes coordinates
        CellVolume(NI-1,NJ-1)  						! cell volumes
  REAL :: x_m, y_m
  integer NI, NJ, IP, JP
  real :: eps
  real :: p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y
  real :: T1, T2, T3, T4

IP = -1
JP = -1
EPS = 1.0e-6

do j = 1, NJ-1
	do i = 1, NI-1

p1x = X(I,J)
p1y = Y(I,J)

p2x = X(I,J+1)
p2y = Y(I,J+1)

p3x = X(I+1,J+1)
p3y = Y(I+1,J+1)

p4x = X(I+1,J)
p4y = Y(I+1,J)

T1 = TS(p1x, p1y, p2x, p2y, x_m, y_m)
T2 = TS(p2x, p2y, p3x, p3y, x_m, y_m)
T3 = TS(p3x, p3y, p4x, p4y, x_m, y_m)
T4 = TS(p4x, p4y, p1x, p1y, x_m, y_m)

S = T1 + T2 + T3 + T4

	If (abs(S - CellVolume(i,j)) .LE. eps) then
		IP = I
		JP = J	
	endif	
	
	end do
end do

End Subroutine



REAL FUNCTION TS(x1,y1,x2,y2,x3,y3)
real x1,y1,x2,y2,x3,y3
real P

A = SQRT((x1-x2)**2.0 + (y1-y2)**2.0)
B = SQRT((x1-x3)**2.0 + (y1-y3)**2.0)
C = SQRT((x2-x3)**2.0 + (y2-y3)**2.0)
		
P = (A+B+C)/2.0
	
TS = SQRT(P*(P-A)*(P-B)*(P-C))

end function