Program Main

  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt', &
  SolutionFile='res_cav.plt',InputParticle='input_particle.txt',&
  OutputParticle='particle.plt',OutputForce='Force.plt' ! names of input and output files
  character MeshFile*30,ctmp        ! name of file with computational mesh
  integer, parameter:: IO = 12 ! input-output unit
  integer :: IGrad,scheme,iter,niter,i,j,NI,NJ,cavity, Nt, m, IP, JP, IP1, JP1, St
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume ! scalar arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector !geometry vector arrays
  real,allocatable,dimension(:,:,:):: GradP
  real,allocatable,dimension(:,:,:):: V                              
  real,allocatable,dimension(:,:):: RotV, Vel_F, Vel_P
  real,allocatable,dimension(:,:,:):: GradV_x, GradV_y, G_GP_x, G_GP_y               !flos velocity and temperature gradient
  real,allocatable,dimension(:,:):: Vmod, T

  real :: Reyn, Pr, Vs, Ls, CFL, VNM, a_diff,time1,time2, rtmp
  real :: ro1, ro2, mu, dp, x0, y0, u0, v0, w0, dt, Vref, Lref, Sk, fb_x, fb_y
  real :: x_m, y_m, u_m, v_m, w_m, x_m1, y_m1, u_m1, v_m1, w_m1, Fx, Fy, Mz, Vx_p, Vy_p, Vrot, Xcros, Ycros,GradP_x, GradP_y
  real :: FDx, FDy, FAx, FAy, FMx, FMy, FBx, FBy, C_FB, pi, mp
!===  READ INPUT FILE ===
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  READ(IO,*) IGrad     ! read key of method for calculating gradient
  READ(IO,*) scheme    ! scheme for scalar
  READ(IO,*) Vs        ! velocity scale
  READ(IO,*) Ls        ! length scale
  READ(IO,*) Reyn      ! Reynolds number
  READ(IO,*) Pr        ! Prandtl number
  READ(IO,*) CFL       ! CFL
  READ(IO,*) VNM       ! VNM
  READ(IO,*) niter     ! number of iterations
  READ(IO,*) cavity    ! 1 - cavity, 0 - operators
  CLOSE(IO)
!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ
  
  

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for I-faces
  allocate(GradP(0:NI,0:NJ,2)) ! Pressure gradient
  allocate(RotV(0:NI,0:NJ)) ! Velocity rotor
  allocate(GradV_x(0:NI,0:NJ,2))
  allocate(GradV_y(0:NI,0:NJ,2))
  allocate(G_GP_x(0:NI,0:NJ,2))
  allocate(G_GP_y(0:NI,0:NJ,2))
  allocate(Vmod(0:NI,0:NJ))    ! delta tau
  allocate(V(0:NI,0:NJ,2)) ! Velocity field
  allocate(T(0:NI,0:NJ))


!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  !READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
  READ(IO,*) ((X(I,J),Y(I,J),rtmp,I=1,NI),J=1,NJ)
  CLOSE(IO)



!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)



!=== READ PARTICLE INFO ===
  OPEN(IO, FILE = InputParticle)
  WRITE(*,*) 'Read particle info from file:', InputParticle
  READ(IO,*) ro1, ro2
  READ(IO,*) mu
  READ(IO,*) dp
  READ(IO,*) x0, y0
  READ(IO,*) u0, v0
  READ(IO,*) w0
  READ(IO,*) dt
  READ(IO,*) Nt
  CLOSE(IO)

  print *, ro1, ro2, mu1, dp, x0, y0, u0, v0, w0, dt, Nt



  allocate(Vel_F(2,1:Nt))
  Vel_F = 0.0
  allocate(Vel_P(2,1:Nt))
  Vel_P = 0.0


!=== READ SOLUTION ===
  OPEN(IO, FILE = SolutionFile)
  READ(IO, *) ctmp
  READ(IO, *) ctmp
  WRITE(*,*) 'Read solution data from file:', SolutionFile
  READ(IO,*) ((rtmp,rtmp,V(I,J,1),V(I,J,2),Vmod(I,J),P(I,J),T(I,J),rtmp,rtmp,I=0,NI),J=0,NJ)
  CLOSE(IO)

!V = 0.0
GradV_x = 0.0
GradV_y = 0.0
! P = 0.0
  ! DO  J = 0,NJ
    ! DO  I = 0,NI
      ! P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(i,j,2))
    ! ENDDO
  ! ENDDO

!=== CALCULATE GRADIENTS ===
  WRITE(*,*) 'Calculate derivatives'
  Call B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,IGrad,P,GradP)
  Call B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,IGrad,GradP(:,:,1),G_GP_x)
  Call B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,IGrad,GradP(:,:,2),G_GP_y)
  
  Call B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,IGrad,V(:,:,1),GradV_x)
  Call B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,IGrad,V(:,:,2),GradV_y)



!=== CALCULATE ROTOR ===
  Call B_CalcRot(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,V,RotV)


!GradP = 0.0
!RotV = 0.0
!=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', OutputFile
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,GradP,V,RotV,GradV_x,GradV_y,Vmod)
  Close(IO)

  Vref = 1.0!Vs
  Lref = 1.0!Ls

  Sk = ro2*dp**2*Vref/(18.0*mu*Lref)
  WRITE(*,*) 'Particle: ', 'd=', dp, 'ro=', ro2
  WRITE(*,*) 'Stokes number Sk= ', Sk
  WRITE(*,*) 'Iteration process: ', 'NT=', Nt, 'dt=', dt

  x_m = x0
  y_m = y0
  u_m = u0
  v_m = v0
  w_m = w0
  
  f_bx = 0.0
  f_by = 0.0
  FDx = 0.0
  FDy = 0.0
  FAx = 0.0
  FAy = 0.0
  FMx  = 0.0
  FMy = 0.0
  FBx = 0.0
  FBy = 0.0
  
  OPEN(IO, FILE = OutputParticle)
  write(IO, *) 'Variables = "it","t","X","Y","u","v","w"'
  OPEN(IO+1, FILE = OutputForce)
  write(IO+1,*) 'Variables = "it","t","x","y","FDx","FDy","FAx","FAy","FMx","FMy","FBx","FBy","Fx","Fy","Mz"'

  
  IP = -1
  JP = -1
  call C_Location (x_m, y_m, NI, NJ, X, Y, CellVolume, IP, JP)
  
  	pi = atan(1.0)*4.0
	mp = ro2*pi*dp**3.0/6.0
  
  	write(*,*) 'Kd=',mp
	!write(*,*) 'Km=',pi*ro1*dp**3/mp/8*1000
	!write(*,*) 'Ka=',pi*dp**3/mp/6

  do m = 1,Nt

	
    Vx_p = V(IP,JP,1) + GradV_x(IP,JP,1)*(x_m - X(IP,JP)) + GradV_x(IP,JP,2)*(y_m - Y(IP,JP))
	Vy_p = V(IP,JP,2) + GradV_y(IP,JP,1)*(x_m - X(IP,JP)) + GradV_y(IP,JP,2)*(y_m - Y(IP,JP))
	Vel_F(1,m) = Vx_p
	Vel_F(2,m) = Vy_p
	
	Vrot = RotV(IP,JP)
	
	Vel_P(1,m) = u_m
	Vel_P(2,m) = v_m
	
	GradP_x = GradP(IP,JP,1) + G_GP_x(IP,JP,1)*(x_m - X(IP,JP)) + G_GP_x(IP,JP,2)*(y_m - Y(IP,JP))
	GradP_y = GradP(IP,JP,2) + G_GP_y(IP,JP,1)*(x_m - X(IP,JP)) + G_GP_y(IP,JP,2)*(y_m - Y(IP,JP))
	!write(*,*) 'GradP_x=',GradP_x
	
	Fx = 0.0
    Fy = 0.0
    Mz = 0.0
	
	C_FB = 6.0*(dp/2.0)**2.0*sqrt(pi*ro1*mu)/mp
	
	if (m < Nt) then
	f_bx = f_bx + C_FB*(Vx_p - u_m - Vel_F(1,m-1) + Vel_P(1,m-1))/sqrt(dt*(Nt - m))
	f_by = f_by + C_FB*(Vy_p - v_m - Vel_F(2,m-1) + Vel_P(2,m-1))/sqrt(dt*(Nt - m))
  end if
	!write(*,*) 'f_bx=',f_bx
	
	call C_Force (ro1, ro2, mu, dp, u_m, v_m, w_m, Vx_p, Vy_p, Vrot, GradP_x, GradP_y, Fx, Fy, Mz, f_bx, f_by, &
	FDx, FDy, FAx, FAy, FMx, FMy, FBx, FBy)
	
    !write(*,*) 'it=',m,'x=',x_m,'y=',y_m,'IP=',IP,'JP=',JP,'u=',u_m,'v=',v_m,'Fx=',Fx,'Fy=',Fy,'Mz=',Mz
	
	write(IO, *) m, m*dt, x_m, y_m, u_m, v_m, w_m
	write(IO+1,*) m, m*dt, x_m, y_m, FDx, FDy, FAx, FAy, FMx, FMy, FBx, FBy, Fx, Fy, Mz
	
	
    x_m1 = x_m + dt*u_m
    y_m1 = y_m + dt*v_m
    u_m1 = u_m + dt*Fx
    v_m1 = v_m + dt*Fy
    w_m1 = w_m + dt*Mz
	call C_Location (x_m1, y_m1, NI, NJ, X, Y, CellVolume, IP1, JP1)
	
	IP = IP1
	JP = JP1
	!write(*,*) 'it=',JP
	
	if(IP1.eq.-1.and.JP1.eq.-1) then
		call C_Boundary(X,Y, IFaceVector, JFaceVector, NI, NJ, x_m, y_m,u_m, v_m, w_m,x_m1, y_m1,u_m1, v_m1, w_m1, Xcros, Ycros,St)
		write(*,*) 'it=',m,'x=',Xcros,'y=',Ycros,'IP=',IP,'JP=',JP,'u=',u_m,'v=',v_m,'Fx=',Fx,'Fy=',Fy,'Mz=',Mz
		write(IO, *) m,m*dt,Xcros,Ycros,u_m,v_m,w_m
		write(IO+1,*) m, m*dt, x_m, y_m, FDx, FDy, FAx, FAy, FMx, FMy, FBx, FBy, Fx, Fy, Mz
		!exit
		call C_Location (x_m1, y_m1, NI, NJ, X, Y, CellVolume, IP2, JP2)
		IP = IP2
		JP = JP2
		if(IP2.eq.-1.and.JP2.eq.-1) then
			write(*,*) 'Particle unexpected left computational domain'
			exit
		end if
	end if

    x_m = x_m1
    y_m = y_m1
    u_m = u_m1
    v_m = v_m1
    w_m = w_m1
  end do
  
END PROGRAM Main


real function Pressure(X,Y) result(Pres)
  Pres = 5*sin(5*y)+5*cos(5*x)
End Function


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




Subroutine B_CalcGradient(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
  IGrad,P,GradP)
  REAL  X(NI,NJ),Y(NI,NJ),&                               ! input: nodes coordinates
        CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   ! cell centers and volumes
        IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! face centers and vectors for I-faces
        JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! face centers and vectors for J-faces
        P(0:NI,0:NJ),&                                    ! Pressure
        GradP_Exact(0:NI,0:NJ,2), GradP_Error(0:NI,0:NJ,2), & !exact value and error for pressure gradient
        GradP(0:NI,0:NJ,2),GradP_tmp(0:NI,0:NJ,2)         ! output: pressure gradient
  integer IGrad     !Gradient method
  real Sf(4,2),rf(4,2),Pf(4),A_ls(2,2),b_ls(2),ri(4,2),GPE(2),rE(2)
  integer NeighCell(4,2),GrG_iter
  real d,dn,det_A, weight, p_dum
  integer i,j,inc,jnc,iface,icell

      select case (IGrad)

      !Green-Gauss
      case(1,2)
        do i = 1,NI-1
          do j = 1,NJ-1

            Sf(1,:) = -IFaceVector(i,j,:)
            Sf(2,:) =  IFaceVector(i+1,j,:)
            Sf(3,:) = -JFaceVector(i,j,:)
            Sf(4,:) =  JFaceVector(i,j+1,:)

            rf(1,:) = IFaceCenter(i,j,:)
            rf(2,:) = IFaceCenter(i+1,j,:)
            rf(3,:) = JFaceCenter(i,j,:)
            rf(4,:) = JFaceCenter(i,j+1,:)

            NeighCell(1,:) = [i-1,j]
            NeighCell(2,:) = [i+1,j]
            NeighCell(3,:) = [i,j-1]
            NeighCell(4,:) = [i,j+1]

            GradP(i,j,:) = 0

            do iface = 1,4
              inc = NeighCell(iface,1)
              jnc = NeighCell(iface,2)
              d = norm2(rf(iface,:) - CellCenter(i,j,:))
              dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:))
              Pf(iface) = RLinearInterp(d,dn,P(i,j),P(inc,jnc))
              ! if (dn < 1e-7) then
                ! p_dum = Pressure(2*rf(iface,1) - CellCenter(i,j,1),2*rf(iface,2) - CellCenter(i,j,2))
                ! Pf(iface) = 0.5*(p_dum + P(i,j))
              ! end if
              GradP(i,j,:) = GradP(i,j,:) + Pf(iface)*Sf(iface,:)
            end do

            GradP(i,j,:) = GradP(i,j,:)/CellVolume(i,j)
            !call Calc_GradP_Exact(CellCenter(I,J,1),CellCenter(i,j,2),GradP_Exact(I,J,:))

          end do
        end do
        !GradP_Error = ABS((GradP_Exact-GradP)/GradP_Exact)
        !write(*,*) '1', maxval(GradP_Error(1:NI-1,1:NJ-1,:))

      !Green-Gauss with iterations
      if (IGrad == 2) then

        GrG_iter = 9

        do k = 1,GrG_iter

        do i = 1,NI-1
          do j = 1,NJ-1

            GradP_tmp(i,j,:) = 0

            Sf(1,:) = -IFaceVector(i,j,:)
            Sf(2,:) =  IFaceVector(i+1,j,:)
            Sf(3,:) = -JFaceVector(i,j,:)
            Sf(4,:) =  JFaceVector(i,j+1,:)

            rf(1,:) = IFaceCenter(i,j,:)
            rf(2,:) = IFaceCenter(i+1,j,:)
            rf(3,:) = JFaceCenter(i,j,:)
            rf(4,:) = JFaceCenter(i,j+1,:)

            NeighCell(1,:) = [i-1,j]
            NeighCell(2,:) = [i+1,j]
            NeighCell(3,:) = [i,j-1]
            NeighCell(4,:) = [i,j+1]

            do iface = 1,4
              inc = NeighCell(iface,1)
              jnc = NeighCell(iface,2)
              d = norm2(rf(iface,:) - CellCenter(i,j,:))
              dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:))
              !Pf(iface) = RLinearInterp(d,dn,P(i,j),P(inc,jnc))
              PE = RLinearInterp(d,dn,P(i,j),P(inc,jnc))
              !if (dn < 1e-7) then
              !  p_dum = Pressure(2*rf(iface,1) - CellCenter(i,j,1),2*rf(iface,2) - CellCenter(i,j,2))
              !  PE = 0.5*(p_dum + P(i,j))
              !end if

              rE(1) = RLinearInterp(d,dn,CellCenter(i,j,1),CellCenter(inc,jnc,1))
              rE(2) = RLinearInterp(d,dn,CellCenter(i,j,2),CellCenter(inc,jnc,2))

              GPE(1) = RLinearInterp(d,dn,GradP(i,j,1),GradP(inc,jnc,1))
              GPE(2) = RLinearInterp(d,dn,GradP(i,j,2),GradP(inc,jnc,2))

              Pf(iface) = PE + dot_product(GPE(:),rf(iface,:) - rE(:))

              GradP_tmp(i,j,:) = GradP_tmp(i,j,:) + Pf(iface)*Sf(iface,:)
            end do

            GradP_tmp(i,j,:) = GradP_tmp(i,j,:)/CellVolume(i,j)

            GradP(i,j,:) = GradP_tmp(i,j,:)

          end do
        end do

        GradP_Error = ABS((GradP_Exact-GradP)/GradP_Exact)
        !write(*,*) k+1, maxval(GradP_Error(1:NI-1,1:NJ-1,:))

        end do

      end if

      !Least squares
      case(3)
        do j = 1,NJ-1
          do i = 1,NI-1

            ri(1,:) = CellCenter(i-1,j,:) - CellCenter(i,j,:)
            ri(2,:) = CellCenter(i+1,j,:) - CellCenter(i,j,:)
            ri(3,:) = CellCenter(i,j-1,:) - CellCenter(i,j,:)
            ri(4,:) = CellCenter(i,j+1,:) - CellCenter(i,j,:)

            NeighCell(1,:) = [i-1,j]
            NeighCell(2,:) = [i+1,j]
            NeighCell(3,:) = [i,j-1]
            NeighCell(4,:) = [i,j+1]

            GradP(i,j,:) = 0
            A_ls(:,:) = 0
            b_ls(:) = 0

            do icell = 1,4
              weight = 1.0/norm2(ri(icell,:))
              inc = NeighCell(icell,1)
              jnc = NeighCell(icell,2)
              A_ls(1,1) = A_ls(1,1) + ri(icell,1)*ri(icell,1)*weight*weight
              A_ls(1,2) = A_ls(1,2) + ri(icell,1)*ri(icell,2)*weight*weight
              A_ls(2,1) = A_ls(2,1) + ri(icell,2)*ri(icell,1)*weight*weight
              A_ls(2,2) = A_ls(2,2) + ri(icell,2)*ri(icell,2)*weight*weight
              b_ls(1) = b_ls(1) + ri(icell,1)*(P(inc,jnc) - P(i,j))*weight*weight
              b_ls(2) = b_ls(2) + ri(icell,2)*(P(inc,jnc) - P(i,j))*weight*weight
            end do

            det_A = A_ls(1,1)*A_ls(2,2) - A_ls(1,2)*A_ls(2,1)

            GradP(i,j,1) = (A_ls(2,2)*b_ls(1) - A_ls(1,2)*b_ls(2))/det_A
            GradP(i,j,2) = (-A_ls(2,1)*b_ls(1) + A_ls(1,1)*b_ls(2))/det_A


          end do
        end do
      end select

End Subroutine



Subroutine B_CalcRot(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
  V,RotV)
  REAL  X(NI,NJ),Y(NI,NJ),&                               ! input: nodes coordinates
        CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   ! cell centers and volumes
        IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& ! face centers and vectors for I-faces
        JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! face centers and vectors for J-faces
        V(0:NI,0:NJ,2),RotV(0:NI,0:NJ)                    ! Velocity and RotV arrays
  real Sf(4,2),rf(4,2),Vf(2),GPE(2),rE(2)
  integer NeighCell(4,2)
  real d,dn
  integer i,j,inc,jnc,iface,icell

        do i = 1,NI-1
          do j = 1,NJ-1

            Sf(1,:) = -IFaceVector(i,j,:)
            Sf(2,:) =  IFaceVector(i+1,j,:)
            Sf(3,:) = -JFaceVector(i,j,:)
            Sf(4,:) =  JFaceVector(i,j+1,:)

            rf(1,:) = IFaceCenter(i,j,:)
            rf(2,:) = IFaceCenter(i+1,j,:)
            rf(3,:) = JFaceCenter(i,j,:)
            rf(4,:) = JFaceCenter(i,j+1,:)

            NeighCell(1,:) = [i-1,j]
            NeighCell(2,:) = [i+1,j]
            NeighCell(3,:) = [i,j-1]
            NeighCell(4,:) = [i,j+1]

            RotV(i,j) = 0

            do iface = 1,4
              inc = NeighCell(iface,1)
              jnc = NeighCell(iface,2)
              d = norm2(rf(iface,:) - CellCenter(i,j,:))
              dn = norm2(rf(iface,:) - CellCenter(inc,jnc,:))
              Vf(1) = RLinearInterp(d,dn,V(i,j,1),V(inc,jnc,1))
              Vf(2) = RLinearInterp(d,dn,V(i,j,2),V(inc,jnc,2))

              RotV(i,j) = RotV(i,j) + (Sf(iface,1)*Vf(2) - Sf(iface,2)*Vf(1))

            end do

            RotV(i,j) = RotV(i,j)/CellVolume(i,j)

          end do
        end do

End Subroutine



real function RLinearInterp(d1,d2,p1,p2)
  real :: d1, d2, p1, p2
  RLinearInterp = (d1*p2 + d2*p1)/(d1+d2)
End Function



REAL FUNCTION TS(x1,y1,x2,y2,x3,y3)
real x1,y1,x2,y2,x3,y3
real P

A = SQRT((x1-x2)**2.0 + (y1-y2)**2.0)
B = SQRT((x1-x3)**2.0 + (y1-y3)**2.0)
C = SQRT((x2-x3)**2.0 + (y2-y3)**2.0)
		
P = (A+B+C)/2.0
	
TS = SQRT(P*(P-A)*(P-B)*(P-C))

end function



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



Subroutine C_Force (ro1, ro2, mu, dp, u_m, v_m, w_m, Vx_p, Vy_p, Vrot, GradP_x, GradP_y, Fx, Fy, Mz, f_bx, f_by, &
	FDx, FDy, FAx, FAy, FMx, FMy, FBx, FBy)
	
  real :: ro1, ro2, mu, dp, u_m, v_m, w_m, Vx_p, Vy_p, Vrot, GradP_x, GradP_y, f_bx, f_by
  real :: Fx, Fy, Mz, FDx, FDy, FAx, FAy, FMx, FMy, FBx, FBy
  real :: pi, mp, Vr,  Rep, Cd, Fd(2), Fa(2), Fm(2), Cm
  real :: Vrm, Rem, Cl, MJp
  
  pi = atan(1.0)*4.0
  mp = ro2*pi*dp**3.0/6.0
  Vr = sqrt((u_m-Vx_p)**2.0 + (v_m-Vy_p)**2.0)
  Rep = dp*ro1*Vr/mu
  Cd = 1.0 + 0.179*sqrt(Rep)+0.013*Rep
  
  Vrm = 0.5*Vrot-w_m
  Rem = dp**2.0*ro1*abs(Vrm)/mu
!write(*,*) 'Rem=',Rep
  if (Rem < 32.0) then
	Cl = 64.0/Rem
    else
    Cl = 64.0/Rem*(0.091*sqrt(Rem) + 0.0017*Rem)
  end if
  MJp = 0.4*dp**2.0*0.25*mp
  
  Cm = 0.75
  
  Fd(1) = 3.0*mu*pi*dp*Cd*(Vx_p-u_m)/mp!*0.0
  Fd(2) = 3.0*mu*pi*dp*Cd*(Vy_p-v_m)/mp!*0.0
  !write(*,*) 'Rem=',Rep
  !write(*,*) 'Rem=',Fd(2)
  Fa(1) = -4.0/3.0*pi*(dp/2.0)**3.0*GradP_x/mp!*0.0
  Fa(2) = -4.0/3.0*pi*(dp/2.0)**3.0*GradP_y/mp!*0.0
  Fm(1) = -4.0/3.0*pi*Cm*(dp/2.0)**3.0*ro1*(Vy_p-v_m)*Vrm/mp!*0.0
  Fm(2) = 4.0/3.0*pi*Cm*(dp/2.0)**3.0*ro1*(Vx_p-u_m)*Vrm/mp!*0.0
  
  
  FDx = Fd(1)
  FDy = Fd(2)
  FAx = Fa(1)
  FAy = Fa(2)
  FMx = Fm(1)
  FMy = Fm(2)
  FBx = f_bx!*0.0
  FBy = f_by!*0.0
  
  Fx = FDx + FAx + FMx + FBx
  Fy = FDy + FAy + FMy + FBy
  
  Mz = 1.0/64.0*pi*Cl*dp**5.0*ro1*abs(Vrm)*Vrm/MJp
End Subroutine



Subroutine C_Boundary(X,Y, IFaceVector, JFaceVector, NI, NJ, x_m, y_m,u_m, v_m, w_m,x_m1, y_m1,u_m1, v_m1, w_m1, Xcros, Ycros,St)

REAL  X(NI,NJ),Y(NI,NJ),IFaceVector( NI,NJ-1,2),JFaceVector( NI-1,NJ,2)
REAL :: x_m, y_m, u_m, v_m, w_m, x_m1, y_m1,u_m1, v_m1, w_m1, Xcros, Ycros, D, NV(2)
integer NI, NJ, IP, JP, Icr, St


! J=1  - wall, St=1 
! J=NJ - wall, St=2 
! I=1  - wall, St=3 
! I=NI - wall, St=4 
  
  
St = 0
Icr = 0

DO I=1, NI-1
	
	CALL Cross_edges(x_m, y_m,x_m1, y_m1, X(I,1),Y(I,1), X(I+1,1),Y(I+1,1),Xcros, Ycros, Icr)
	IF(Icr.eq.1) then
		St=1
		D = sqrt(JFaceVector( I,1,1)**2+JFaceVector(I,1,2)**2)
		NV(1) = -JFaceVector( I,1,1)/D
		NV(2) = -JFaceVector( I,1,2)/D
		CALL Boundary_Condition(NV, x_m, y_m, x_m1, y_m1, Xcros, Ycros, u_m, v_m, w_m,u_m1, v_m1, w_m1)
		write(*,*) 'Particle cross boundary 1', St, 'in coordinates', Xcros, Ycros
		return
	endif
	
	CALL Cross_edges(x_m, y_m,x_m1, y_m1, X(I,NJ),Y(I,NJ), X(I+1,NJ),Y(I+1,NJ),Xcros, Ycros, Icr)
	IF(Icr.eq.1) then
		St=2
		D = sqrt(JFaceVector( I,NJ,1)**2+JFaceVector(I,NJ,2)**2)
		NV(1) = -JFaceVector( I,NJ,1)/D
		NV(2) = -JFaceVector( I,NJ,2)/D
		CALL Boundary_Condition(NV, x_m, y_m, x_m1, y_m1, Xcros, Ycros, u_m, v_m, w_m,u_m1, v_m1, w_m1)
		write(*,*) 'Particle cross boundary 2', St, 'in coordinates', Xcros, Ycros
		return
	endif
end do

DO J=1, NJ-1
	
	CALL Cross_edges(x_m, y_m,x_m1, y_m1, X(1,J),Y(1,J), X(1,J+1),Y(1,J+1),Xcros, Ycros, Icr)
	IF(Icr.eq.1) then
		St=3
		D = sqrt(IFaceVector( 1,J,1)**2+IFaceVector(1,J,2)**2)
		NV(1) = -IFaceVector( 1,J,1)/D
		NV(2) = -IFaceVector( 1,J,2)/D
		CALL Boundary_Condition(NV, x_m, y_m, x_m1, y_m1, Xcros, Ycros, u_m, v_m, w_m,u_m1, v_m1, w_m1)
		write(*,*) 'Particle cross boundary 3', St, 'in coordinates', Xcros, Ycros
		return
	endif
	
	CALL Cross_edges(x_m, y_m,x_m1, y_m1, X(NI,J),Y(NI,J), X(NI,J+1),Y(NI,J+1),Xcros, Ycros, Icr)
	IF(Icr.eq.1) then
		St=4
		D = sqrt(IFaceVector( NI,J,1)**2+IFaceVector(NI,J,2)**2)
		NV(1) = -IFaceVector( NI,J,1)/D
		NV(2) = -IFaceVector( NI,J,2)/D
		CALL Boundary_Condition(NV, x_m, y_m, x_m1, y_m1, Xcros, Ycros, u_m, v_m, w_m,u_m1, v_m1, w_m1)
		write(*,*) 'Particle cross boundary 4', St, 'in coordinates', Xcros, Ycros
		return
	endif
end do

End Subroutine



Subroutine Cross_edges(x, y, x1, y1, Cx, Cy, Dx, Dy, Xcros, Ycros, Icr)

REAL :: x, y, x1, y1, Cx, Cy, Dx, Dy, Xcros, Ycros
REAL :: t, rx, ry, sx, sy, u
integer Icr

Icr = 0


rx = x1-x
ry = y1-y

sx = Dx-Cx
sy = Dy-Cy

t = ((Cx-x)*sy-(Cy-y)*sx)/(rx*sy-ry*sx)
u = ((Cx-x)*ry-(Cy-y)*rx)/(rx*sy-ry*sx)

!Icr = 1
Xcros = x + t*rx
Ycros = y + t*ry
		
if(t.le.1.and.t.ge.0) then
	if (u.le.1.and.u.ge.0) then
		Icr = 1
		Xcros = x + t*rx
		Ycros = y + t*ry
	end if
end if

End Subroutine



Subroutine Boundary_Condition(NV, x_m, y_m, x_m1, y_m1, Xcros, Ycros, u_m, v_m, w_m,u_m1, v_m1, w_m1)

REAL ::u_m, v_m, w_m, x_m, y_m, x_m1, y_m1,u_m1, v_m1, w_m1, Xcros, Ycros, D, NV(2), TV(2), Kn, Kt, Vn, Vt, delT, A, r, B
integer NI, NJ, IP, JP, Icr, St
!write(*,*) 'NV1',NV(1),'NV2',NV(2)

Kt = 0.9
Kn = 0.8

TV(1) = NV(2)
TV(2) = - NV(1)

!write(*,*) 'NV(1)=',NV(1), 'NV(2)=',NV(2), 'TV(1)=',TV(1), 'TV(2)=',TV(2)

Vn = u_m*NV(1) + v_m*NV(2)
Vt = u_m*TV(1) + v_m*TV(2)

r = sqrt((x_m1 - x_m)**2 + (y_m1 - y_m)**2)

!==========ГУ 0=================================================

! u_m1 = Kt*Vt*TV(1) - Kn*Vn*NV(1)
! v_m1 = Kt*Vt*TV(2) - Kn*Vn*NV(2)
! w_m1 = Kt*(u_m1*TV(1) + v_m1*TV(2))/r

!==========ГУ 3=================================================

A = 1.0/7.0*((5.0+2.0*Kt)*Vt+2*r*w_m*(1-Kt))
B = -Kn*Vn
u_m1 = A*NV(2) - B*TV(2)
write(*,*) 'u_m1=',u_m1
v_m1 = B*TV(1) - A*NV(1)
w_m1 = w_m-5.0/7.0*(1.0 - Kt)*(w_m - (u_m1*TV(1) + v_m1*TV(2))/r)

!===============================================================

delT = sqrt((x_m1-Xcros)**2+(y_m1-Ycros)**2)/sqrt(u_m**2+v_m**2)
x_m1 = Xcros + u_m1*delT
y_m1 = Ycros + v_m1*delT



End Subroutine



Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,GradP,V,RotV,GradV_x,GradV_y,Vmod)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P,DivV,LaplP,RotV,Vmod
  Real,Dimension(0:NI,0:NJ,2)::GradP,V,GradV_x,GradV_y

  Write(IO,*) 'VARIABLES = "X", "Y", "P","GradP_X","GradP_Y","V_x","V_y","Vmod","RotV"'
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ)
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') Vmod(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') RotV(1:NI-1,1:NJ-1)
End Subroutine