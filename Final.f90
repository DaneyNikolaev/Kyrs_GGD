Program Main
include 'omp_lib.h'

  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt' ! names of input and output files
  character MeshFile*30        ! name of file with computational mesh
  integer, parameter:: IO = 12 ! input-output unit
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume ! scalar arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! vector arrays
  real,allocatable,dimension(:,:,:):: GradP, GradP_Exact, GradP_Error 
  real,allocatable,dimension(:,:):: DivV, DivV_Exact, DivV_Error
  real,allocatable,dimension(:,:,:):: V  
  real,allocatable,dimension(:,:) :: LapP, RLapP_Exact, LapP_Error
  real,allocatable,dimension(:,:):: RotV, RotV_Exact, RotV_Error
  real,allocatable,dimension(:,:,:):: VFlos, GradT
  real,allocatable,dimension(:,:):: T, ResT, TFlos, T_Error, dtau
  Integer :: niter, iter, IGrad, scheme, i,j,NI, NJ, var
  real :: Rei, Pr, Vs, L, CFL, VNM, Koef_Temp,time1,time2


!===============  READ INPUT FILE ===============
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
 READ(IO,*) MeshFile  ! read name of file with computational mesh
 READ(IO,*) niter  !число итераций
 READ(IO,*) scheme  !схема для расчета конвективных потоков
 READ(IO,*) IGrad  !метод для градиента
 READ(IO,*) var  !вариант задачи (== 1 если течение в каверне)
 READ(IO,*) Vs  !масштаб скорости
 READ(IO,*) L  !масштаб длины
 READ(IO,*) Rei  !число Рейнольдса
 READ(IO,*) Pr  !число Прандтля
 READ(IO,*) CFL  !число Куранта
 READ(IO,*) VNM  !число фон Неймана
  CLOSE(IO)
!===============  READ INPUT FILE ===============



!===============   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===============
  WRITE(*,*) 'Read nodes number from file: ', MeshFile !чтение размеров сетки из файла с ней
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ
!===============   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===============



!=============== ALLOCATE ALL ARRAYS ===============
  WRITE(*,*) 'Allocate arrays'
  allocate(X(NI,NJ)) ! x-координаты узлов сетки
  allocate(Y(NI,NJ)) ! y-координаты узлов сетки
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(CellVolume(NI-1,NJ-1))   ! объем ячеек
  allocate(CellCenter(0:NI,0:NJ,2)) ! центры ячеек
  allocate(IFaceCenter( NI,NJ-1,2)) ! центры граней для I-граней
  allocate(IFaceVector( NI,NJ-1,2)) ! норм. векторы для I-граней
  allocate(JFaceCenter( NI-1,NJ,2)) ! центры граней для J-граней
  allocate(JFaceVector( NI-1,NJ,2)) ! норм. векторы для J-граней
  allocate(GradP(0:NI,0:NJ,2))		! градиент давления
  allocate(GradP_Exact(0:NI,0:NJ,2))! аналит градиент давления
  allocate(GradP_Error(0:NI,0:NJ,2))! ошибка град давления
  allocate(DivV(0:NI,0:NJ))			! дивергенция
  allocate(DivV_Exact(0:NI,0:NJ))   ! аналит дивергенци
  allocate(DivV_Error(0:NI,0:NJ))   ! ошибка дивергенции
  allocate(V(0:NI,0:NJ,2))			! вектор скорости
  allocate(LapP(0:NI,0:NJ))			! лаплас
  allocate(RLapP_Exact(0:NI,0:NJ))	! аналит лаплас
  allocate(LapP_Error(0:NI,0:NJ))	! ошибка лаплас
  allocate(RotV(0:NI,0:NJ))			! ротор
  allocate(RotV_Exact(0:NI,0:NJ))	! аналит ротор
  allocate(RotV_Error(0:NI,0:NJ))	! ошибка ротор
  allocate(GradT(0:NI,0:NJ,2))		! градиент температуры
  allocate(VFlos(0:NI,0:NJ,2))		! поле скорости из Flos
  allocate(T(0:NI,0:NJ))			! температура
  allocate(ResT(0:NI,0:NJ))			! невязка
  allocate(TFlos(0:NI,0:NJ))		! температура из Flos
  allocate(T_Error(0:NI,0:NJ))		! ошибка температуры
  allocate(dtau(0:NI,0:NJ))			! шаг по посевдовремени
!=============== ALLOCATE ALL ARRAYS ===============



!===============  READ GRID ===============
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
  CLOSE(IO)
!===============  READ GRID =============== 
  
  
  
  
!=============== Задача о каверне ===============
if (var == 1) then
!==  READ FLOS VELOCITY FIELD ==
  OPEN(IO,FILE = 'VelocityField.txt')
  !OPEN(IO,FILE = 'v.txt')
  READ(IO,*) ((VFlos(I,J,1),VFlos(I,J,2),I=1,NI),J=1,NJ)
  CLOSE(IO)

!==  READ FLOS TEMPERATURE FIELD ==
  OPEN(IO,FILE = 'Temperature.txt')
  READ(IO,*) ((TFlos(I,J),I=1,NI),J=1,NJ)
  CLOSE(IO)
end if
!=============== Задача о каверне ===============




!=============== CALCULATE METRIC ===============
  WRITE(*,*) 'Calculate metric'
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) !вызов расчета геом.характеристик
!=============== CALCULATE METRIC ===============



!=============== INITIATE FIELDS ===============
WRITE(*,*) 'Initiate fields'
DO  J = 0,NJ
    DO  I = 0,NI
		P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(i,j,2))
		Call CalcGradPressureExact(GradP_Exact(I,J,:), CellCenter(I,J,1), CellCenter(I,J,2))  ! расчет точных значений
		Call VELOCITY(CellCenter(I,J,1), CellCenter(I,J,2),V(I,J,:))
		DivV_Exact(I,j)=DivVelocityExact(CellCenter(I,J,1),CellCenter(i,j,2))
		RLapP_Exact(I,j)=RLapPExact(CellCenter(I,J,1),CellCenter(i,j,2))
		RotV_Exact(I,J) = RotVelocityExact(CellCenter(I,J,1),CellCenter(i,j,2))
	ENDDO
ENDDO
!=============== INITIATE FIELDS ===============



!=============== CALCULATE DIF OPERATORS ===============

gradP=0

!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'CALCULATE GRADIENT'
  Call B_CalcGradient(NI, NJ, X, Y, GradP,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector, p, IGrad)
  GradP_Error = ABS((GradP_Exact-GradP)/GradP_Exact)
  write(*,*) maxval(GradP_Error(1:NI-1,1:NJ-1,:))

!=== CALCULATE DIVERGENCE ===
WRITE(*,*) 'CALCULATE DIVERGENCE'
  Call B_CalcDIVERGANCE(NI, NJ, X, Y, DivV, CellCenter,CellVolume,IFaceCenter,IFaceVector,&
    JFaceCenter,JFaceVector, V,p, GradP, scheme)
  DivV_Error = ABS((DivV_Exact-DivV)/DivV_Exact)
  write(*,*) maxval(DivV_Error(1:NI-1,1:NJ-1))

!=== CALCULATE LAPLACIAN ===
WRITE(*,*) 'CALCULATE LAPLACIAN'
  Call B_CalcLaplacian(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
  P,GradP,LapP)
  LapP_Error = ABS((RLapP_Exact-LapP)/RLapP_Exact)
  write(*,*) maxval(LapP_Error(1:NI-1,1:NJ-1))

!=== CALCULATE ROTOR ===
WRITE(*,*) 'CALCULATE ROTOR'
  Call B_CalcRotor(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
  V,RotV)
  RotV_Error = ABS((RotV_Exact-RotV)/RotV_Exact)
  write(*,*) maxval(RotV_Error(1:NI-1,1:NJ-1))
!=============== CALCULATE DIF OPERATORS ===============




!=============== CALCULATE RESIDUAL T ===============

  T(:,:) = 1.0  !инициализация поля
  T(:,NJ) = 2.0	!ГУ на крышке	
  dtau = 0.01	!инициализация поля
  Koef_Temp = Vs*L/(Rei*Pr)  !коэффициент температуропроводности
  
  OPEN(IO,FILE = 'ResT.plt')
  write(IO,*) 'Variables = "iterations", "Res T"'
  !GradT = 0
!$ time1 = OMP_GET_WTIME()
!$OMP parallel
  do iter = 1,niter
    Call B_CalcGradient(NI, NJ, X, Y, GradT,CellCenter,CellVolume,IFaceCenter,&
    IFaceVector,JFaceCenter,JFaceVector, T, IGrad) !вычисление градиента температуры
    Call B_CalcResidualT(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
  Rei,Pr,VFlos,T,GradT,ResT,scheme,CFL,VNM,Koef_Temp,dtau) !вычисление невязки температуры
!$OMP single
    write(*,*) iter, maxval(abs(ResT(1:NI-1,1:NJ-1)))
    write(IO,*) iter, maxval(abs(ResT(1:NI-1,1:NJ-1)))
!$OMP end single
!$OMP DO private(i,j)
!Метод установления. Явная схема Эйлера
    do j = 1, NJ-1
      do i = 1, NI-1
        T(i,j) = T(i,j) - ResT(i,j)*dtau(i,j)
      end do
    end do
!$OMP END DO
  end do
!$OMP end parallel
!$ time2 = OMP_GET_WTIME()
  print *, 'time = ', time2-time1
  close(IO)
  !write(*,*) T(:,NJ)
if (var == 1) then
  T_Error = ABS((TFlos-T)/TFlos)
end if
write(*,*) 'TEMPERATURE error', maxval(T_Error(1:NI-1,1:NJ-1))
!=============== CALCULATE RESIDUAL T ===============



!=============== OUTPUT FIELDS ===============
  WRITE(*,*) 'Output fields to file: ', OutputFile
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO, NI, NJ, X, Y, P, GradP, GradP_Error,V,DivV,DivV_Exact,&
    DivV_Error,LapP,RLapP_Exact,LapP_Error,RotV,RotV_Exact,RotV_Error,&
    VFlos,T,GradT,ResT,TFlos,T_Error)
  Close(IO)
!=============== OUTPUT FIELDS ===============
END PROGRAM Main



SUBROUTINE B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector)
  REAL X(NI,NJ),Y(NI,NJ),&                               
       CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   
       IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& 
       JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)
  REAL r(2)

  !=== FACE CENTERS AND FACE VECTORS ===
  ! I-DIRECTION
  DO J = 1,NJ-1
    DO I = 1,NI
      r(1) = X(I,J+1) - X(I,J)  ! вектор от одного узла к другому
      r(2) = Y(I,J+1) - Y(I,J)
      IFaceVector(I,J,1) = r(2) ! IFaceVector = r повернутый на 90 градусов
      IFaceVector(I,J,2) =-r(1) ! IFaceVector в направлении роста I
      IFaceCenter(I,J,1) = 0.5*(X(i,j)+x(i,j+1)) !центры для I-граней
      IFaceCenter(I,J,2) = 0.5*(Y(i,j)+Y(i,j+1))
    ENDDO
  ENDDO

  ! J-DIRECTION
  DO J = 1,NJ
    DO I = 1,NI-1
      r(1) = X(I+1,J) - X(I,J)  ! вектор от одного узла к другому
      r(2) = Y(I+1,J) - Y(I,J)
      JFaceVector(I,J,1) =-r(2) ! JFaceVector = r повернутый на -90 градусов
      JFaceVector(I,J,2) = r(1) ! JFaceVector в направлении роста J
      JFaceCenter(I,J,1) = 0.5*(X(i,j)+x(i+1,j))
      JFaceCenter(I,J,2) = 0.5*(Y(i,j)+Y(i+1,j))
    ENDDO
  ENDDO


 !=== CELL VOLUMES ===
  DO J = 1,NJ-1
    DO I = 1,NI-1
      r(1)=X(I+1,J+1) - X(I,J)
      r(2)=Y(I+1,J+1) - Y(I,J)
      CellVolume(I,J) = 0.5*DOT_PRODUCT(IFaceVector(I,J,:),r)&
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



!=============== FUNCTIONS ===============
Function Pressure(X,Y) !Скалярная функция давления
  !Pres = (x+y)**2
  Pres = x**3+y**3
  !Pres = (x+y)**3
End Function

Subroutine VELOCITY(x,y,V) !Вектор скорости
Real V(2)
  V(1) = y**2
  V(2) = -x**1
End  Subroutine

Function rLinearInterp(d1,d2,p1,p2) !линейная интерполяция
    !d1 - расстояние до p2, d2 - расстояния до p1
	rLinearInterp=(d1*p2+d2*p1)/(d1+d2)
End Function

Subroutine CalcGradPressureExact(GP_E, x, y) !Аналитическое значение градиента давления
Real GP_E(2)
  GP_E(1) = 1
  GP_E(2) = 1
  !GP_E(1) = 2*x
  !GP_E(2) = 2*y
  !GP_E(1) = 3*(x+y)**2
  !GP_E(2) = 3*(x+y)**2
End  Subroutine

Function DivVelocityExact(X,Y) !Аналитическое значение дивергенции
	DivVelocityExact = 2
	!DivVelocityExact = 4*x**3+2*x*y**2+2*y*x**2+4*y**3
End Function

Function RLapPExact(X,Y) !Аналитическое значение лапласиана
RLapPExact = 4
End Function

real function RotVelocityExact(X,Y) !Аналитическое значение ротора
RotVelocityExact = -(1+2*y)
End Function

!=============== FUNCTIONS ===============



Subroutine B_CalcGradient(NI, NJ, X, Y, GradP,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector, p, IGrad)

 real,dimension(0:NI,0:NJ,2):: GradP, GradP_tmp, GradP_Exact, GradP_Error
 real,dimension(0:NI,0:NJ):: p
 INTEGER :: NeighCell(4,2), GGI, IGrad, i,j , inc, jnc,iface
 REAL SF(4,2), rF(4,2), Pf(4)
 REAL CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   ! центры и объемы ячеек
       IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& !центры и нормали граней для I-граней
       JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! центры и нормали граней для J-граней
       X(NI,NJ),Y(NI,NJ) ! координаты узлов сетки
 REAL p_dum, d, dn
 REAL rE(2), GPE(2), Pe(4)


 select case (IGrad)

      !метод Грина-Гаусса
      case(1)
!$OMP DO private(i,j,Sf,rf,NeighCell,iface,IN,JN,d,dn,Pf)

     DO I=1, NI-1
     DO J=1, NJ-1

SF(1,:)=-IFaceVector(I,J,:) 
SF(2,:)=IFaceVector(I+1,J,:)
SF(3,:)=-JFaceVector(I,J,:)
SF(4,:)=JFaceVector(I,J+1,:)

rF(1,:)=IFaceCenter(I,J,:) 
rF(2,:)=IFaceCenter(I+1,J,:)
rF(3,:)=JFaceCenter(I,J,:)
rF(4,:)=JFaceCenter(I,J+1,:)


NeighCell(1,:)=[i-1,j]
NeighCell(2,:)=[I+1,j]
NeighCell(3,:)=[i,j-1]
NeighCell(4,:)=[i,j+1]

GradP(I,J,:)=0

	DO IFACE=1,4 
	IN=NeighCell(iface,1)
	JN=NeighCell(iface,2)

	d=Norm2(rF(iface,:)-CellCenter(I,j,:))  !расстояние от центра ячейки до центра грани
	dn=Norm2(rF(iface,:)-CellCenter(IN,JN,:)) !расстояние до центра соседней ячейки

	Pf(iface)=rLinearInterp(d,dn,P(I,J),P(IN,JN))

if (dn<1e-7) then !значение в центре заграничной ячейки
p_dum=Pressure(2*rf(iface,1)-CellCenter(i,j,1),2*rf(iface,2)-CellCenter(i,j,2))
Pf(iface)=0.5*(p_dum+P(i,j))
end if

GradP(I,J,:)=GradP(I,J,:)+Pf(iface)*SF(iface,:)
enddo

GradP(I,J,:)=GradP(I,J,:)/CellVolume(I,J)

END DO
END DO
!$OMP END DO

!метод Грина-Гуасса с итерациями
IF (IGrad==2) then

    GGI=5 !итерации самого метода
    do k=1, GGI
!$OMP DO private(i,j,Sf,rf,NeighCell,iface,IN,JN,d,dn,PE,rE,GPE,Pf)
    DO I=1, NI-1
    DO J=1, NJ-1

        GradP_tmp(i,j,:)=0

SF(1,:)=-IFaceVector(I,J,:) 
SF(2,:)=IFaceVector(I+1,J,:)
SF(3,:)=-JFaceVector(I,J,:)
SF(4,:)=JFaceVector(I,J+1,:)

rF(1,:)=IFaceCenter(I,J,:) 
rF(2,:)=IFaceCenter(I+1,J,:)
rF(3,:)=JFaceCenter(I,J,:)
rF(4,:)=JFaceCenter(I,J+1,:)


NeighCell(1,:)=[i-1,j]
NeighCell(2,:)=[I+1,j]
NeighCell(3,:)=[i,j-1]
NeighCell(4,:)=[i,j+1]

    DO IFACE=1,4
	IN=NeighCell(iface,1)
	JN=NeighCell(iface,2)

	d=Norm2(rF(iface,:)-CellCenter(I,j,:)) !расстояние от центра ячейки до центра грани
	dn=Norm2(rF(iface,:)-CellCenter(IN,JN,:)) !расстояние до центра соседней ячейки

    Pe(iface)=rLinearInterp(d,dn,P(I,J),P(IN,JN))

	rE(1)=rLinearInterp(d,dn,CellCenter(I,J,1),CellCenter(IN,JN,1))
	rE(2)=rLinearInterp(d,dn,CellCenter(I,J,2),CellCenter(IN,JN,2))

	GPE(1)=rLinearInterp(d,dn,GradP(I,J,1),GradP(IN,JN,1))
	GPE(2)=rLinearInterp(d,dn,GradP(I,J,2),GradP(IN,JN,2))

    Pf(iface)=Pe(iface)+DOT_PRODUCT(GPE(:),rF(iface,:)-rE(:)) ! с поправкой для скошенных ячеек

	GradP_tmp(I,J,:)=GradP_tmp(I,J,:)+Pf(iface)*SF(iface,:)
	ENDDO

GradP_tmp(I,J,:)=GradP_tmp(I,J,:)/CellVolume(I,J)

GradP(I,J,:)=GradP_tmp(I,J,:)

    END DO
    END DO
!$OMP END DO

enddo
endif

end select
End Subroutine



Subroutine B_CalcDIVERGANCE(NI, NJ, X, Y, DivV, CellCenter,CellVolume,IFaceCenter,IFaceVector,&
    JFaceCenter,JFaceVector, V,p, GradP, scheme)
 real,dimension(0:NI,0:NJ,2):: V, GradP
 real,dimension(0:NI,0:NJ):: DivV, p
 INTEGER :: NeighCell(4,2)
 REAL SF(4,2), rF(4,2), Pf(4)
 REAL X(NI,NJ),Y(NI,NJ),& ! координаты узлов сетки
 CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   ! центры и объемы ячеек
       IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& !центры и нормали граней для I-граней
       JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)   ! центры и нормали граней для J-граней
REAL rE(2), GPE(2), Pe(4), Vf(2)
integer scheme 

    DO I=1, NI-1
      DO J=1, NJ-1

        SF(1,:)=-IFaceVector(I,J,:)
        SF(2,:)=IFaceVector(I+1,J,:)
        SF(3,:)=-JFaceVector(I,J,:)
        SF(4,:)=JFaceVector(I,J+1,:)

        rF(1,:)=IFaceCenter(I,J,:)
        rF(2,:)=IFaceCenter(I+1,J,:)
        rF(3,:)=JFaceCenter(I,J,:)
        rF(4,:)=JFaceCenter(I,J+1,:)

        NeighCell(1,:)=[i-1,j]
        NeighCell(2,:)=[I+1,j]
        NeighCell(3,:)=[i,j-1]
        NeighCell(4,:)=[i,j+1]

    DivV(I,J)=0

	DO IFACE=1,4
	IN=NeighCell(iface,1)
	JN=NeighCell(iface,2)

	d=Norm2(rF(iface,:)-CellCenter(I,j,:))!расстояние от центра ячейки до центра грани
	dn=Norm2(rF(iface,:)-CellCenter(IN,JN,:)) !расстояние до центра соседней ячейки


    Vf(1)=rLinearInterp(d,dn,V(I,J,1),V(IN,JN,1))
    Vf(2)=rLinearInterp(d,dn,V(I,J,2),V(IN,JN,2))

select case(scheme)
     case(0)
     Pf(iface)=1

     case(1) ! FOU
       If (DOT_PRODUCT(Vf,SF(iface,:)).gt.0.0) then !Знак расхода на грани
         Pf(iface)=P(I,j)
       else
         Pf(iface)=P(In,jn)
		 if (dn.lt.1e-7) then !для приграничных ячеек
			Pf(iface)=Pf(iface)+(P(in,jn)-P(I,j))
		end if
       end if

          case(2) !central
          Pf(iface)=rLinearInterp(d,dn,P(I,J),P(IN,JN))

          case(3) !SOU
            if (dot_product(Vf(:),Sf(iface,:)) > 0.0) then !Знак расхода на грани
                Pf(iface) = P(i,j) + dot_product(GradP(i,j,:),rf(iface,:) - CellCenter(i,j,:))
            else
                Pf(iface) = P(in,jn) + dot_product(GradP(in,jn,:),rf(iface,:) - CellCenter(in,jn,:))
                if (dn < 1e-7) then !для приграничных ячеек
                    Pf(iface) = 2*P(in,jn) - P(i,j) -&
                    4*P(in,jn) + 4*P(i,j) +&
                    3*dot_product(GradP(i,j,:),rf(iface,:)-CellCenter(i,j,:))
                    !Pf(iface) = P(i,j) + dot_product(GradP(i,j,:),rf(iface,:) - CellCenter(i,j,:))
                end if
            end if
          end select
		DivV(I,j)=DivV(I,j)+DOT_PRODUCT(Pf(iface)*Vf(:),Sf(iface,:))
	ENDDO
	DivV(I,j)=DivV(I,j)/CellVolume(I,J)
	end do
end do

End Subroutine



Subroutine B_CalcLaplacian(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
  P,GradP,LapP)

 real,dimension(0:NI,0:NJ,2):: GradP
 real,dimension(0:NI,0:NJ):: p, LapP
 INTEGER :: NeighCell(4,2)
 REAL SF(4,2), rF(4,2), Pf(4)
 REAL X(NI,NJ),Y(NI,NJ),&  ! координаты узлов сетки
       CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   !центры и объемы ячеек
       IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& !центры и нормали граней для I-граней
       JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)   ! центры и нормали граней для J-граней
 REAL rE(2), GPE(2), Pe(4), NF(2), RNC(2)
 REAL d, dn, DNC, dpdn

LapP=0 

DO I=1, NI-1
DO J=1, NJ-1

SF(1,:)=-IFaceVector(I,J,:) 
SF(2,:)=IFaceVector(I+1,J,:)
SF(3,:)=-JFaceVector(I,J,:)
SF(4,:)=JFaceVector(I,J+1,:)

rF(1,:)=IFaceCenter(I,J,:) 
rF(2,:)=IFaceCenter(I+1,J,:)
rF(3,:)=JFaceCenter(I,J,:)
rF(4,:)=JFaceCenter(I,J+1,:)


NeighCell(1,:)=[i-1,j]
NeighCell(2,:)=[I+1,j]
NeighCell(3,:)=[i,j-1]
NeighCell(4,:)=[i,j+1]

	DO IFACE=1,4 
	IN=NeighCell(iface,1)
	JN=NeighCell(iface,2)

	d=Norm2(rF(iface,:)-CellCenter(I,j,:)) !расстояние от центра ячейки до центра грани
	dn=Norm2(rF(iface,:)-CellCenter(IN,JN,:))  !расстояние до центра соседней ячейки

	Pf(iface)=rLinearInterp(d,dn,P(I,J),P(IN,JN))

    DNC=Norm2(CellCenter(IN,JN,:)-CellCenter(I,J,:)) !расстояние между соседними ячейками
    dpdn=(P(IN,JN)-P(I,J))/DNC !производная по направлению к соседней ячейке
	RNC=(CellCenter(IN,JN,:)-CellCenter(I,J,:))/DNC !вектор, задающий направление от одного центра к другому

	rE(1)=rLinearInterp(d,dn,CellCenter(I,J,1),CellCenter(IN,JN,1))
	rE(2)=rLinearInterp(d,dn,CellCenter(I,J,2),CellCenter(IN,JN,2))

	GPE(1)=rLinearInterp(d,dn,GradP(I,J,1),GradP(IN,JN,1))
	GPE(2)=rLinearInterp(d,dn,GradP(I,J,2),GradP(IN,JN,2))

    NF(:)=SF(iface,:)/Norm2(SF(iface,:))

    if (dn.le.1e-7) then !для приграничных ячеек

	   dpdn_c=DOT_PRODUCT(GradP(i,j,:),NF(:))
	   !dpdn=dpdn+(dpdn-dpdn_c)   !первый порядок
	   dpdn=5./3.*dpdn-2./3.*dpdn_c   !второй порядок
       GPE(:)=GradP(i,j,:)	   !на границе берем значение из центра ячейки
	endif

    dpdn=dpdn+DOT_PRODUCT(NF(:)-RNC(:),GPE(:))
	LapP(I,J)=LapP(I,J)+dpdn*Norm2(SF(iface,:))
	ENDDO

LapP(I,J)=LapP(I,J)/CellVolume(I,J)

END DO
END DO


End Subroutine



Subroutine B_CalcRotor(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
  V,RotV)
 real,dimension(0:NI,0:NJ,2):: V
 INTEGER :: NeighCell(4,2)
 REAL SF(4,2), rF(4,2)
 REAL X(NI,NJ),Y(NI,NJ),& ! координаты узлов сетки
 CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   !центры и объемы ячеек
       IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& !центры и нормали граней для I-граней
       JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),&   ! центры и нормали граней для J-граней
       RotV(0:NI,0:NJ)
 REAL rE(2), GPE(2), Vf(2)
 REAL  d, dn


    DO I=1, NI-1
      DO J=1, NJ-1

        SF(1,:)=-IFaceVector(I,J,:)
        SF(2,:)=IFaceVector(I+1,J,:)
        SF(3,:)=-JFaceVector(I,J,:)
        SF(4,:)=JFaceVector(I,J+1,:)

        rF(1,:)=IFaceCenter(I,J,:)
        rF(2,:)=IFaceCenter(I+1,J,:)
        rF(3,:)=JFaceCenter(I,J,:)
        rF(4,:)=JFaceCenter(I,J+1,:)

        NeighCell(1,:)=[i-1,j]
        NeighCell(2,:)=[I+1,j]
        NeighCell(3,:)=[i,j-1]
        NeighCell(4,:)=[i,j+1]

    RotV(I,J)=0

	DO IFACE=1,4
	IN=NeighCell(iface,1)
	JN=NeighCell(iface,2)

	d=Norm2(rF(iface,:)-CellCenter(I,j,:))!расстояние от центра ячейки до центра грани
	dn=Norm2(rF(iface,:)-CellCenter(IN,JN,:))  !расстояние до центра соседней ячейки

    Vf(1)=rLinearInterp(d,dn,V(I,J,1),V(IN,JN,1))
    Vf(2)=rLinearInterp(d,dn,V(I,J,2),V(IN,JN,2))

              RotV(i,j) = RotV(i,j) + (Sf(iface,1)*Vf(2) - Sf(iface,2)*Vf(1))

            end do

            RotV(i,j) = RotV(i,j)/CellVolume(i,j)

          end do
        end do

End Subroutine



Subroutine B_CalcResidualT(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
  Rei,Pr,VFlos,T,GradT,ResT,scheme,CFL,VNM,Koef_Temp,dtau)
 INTEGER :: NeighCell(4,2), scheme, i,j,iface, in,jn
  real Sf(4,2),rf(4,2),GTE(2),rE(2),Tf,rnc(2),Nf(2),Vf(2), VTf(2)
 REAL CellCenter(0:NI,0:NJ,2),CellVolume(NI-1,NJ-1),&   !центры и объемы ячеек
       IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2),& !центры и нормали граней для I-граней
       JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),& ! центры и нормали граней для J-граней
       X(NI,NJ),Y(NI,NJ),& ! координаты узлов сетки
       T(0:NI,0:NJ),GradT(0:NI,0:NJ,2),&
       VFlos(0:NI,0:NJ,2),&
       ResT(0:NI,0:NJ),dtau(0:NI,0:NJ)
 REAL d, dn, DNC,dTdn,dTdn_c,Rei,Pr,CFL,VNM,Koef_Temp,dtau_c,dtau_d
	

!$OMP DO private(i,j,Sf,rf,NeighCell,dtau_c,dtau_d,iface,IN,JN,d,dn,Vf,Tf,VTf,DNC,dTdn,dTdn_c,RNC,rE,GTE,Nf)

DO j=1, NJ-1
      DO i=1, NJ-1
		ResT(i,j) = 0.0
        

        SF(1,:)=-IFaceVector(I,J,:)
        SF(2,:)=IFaceVector(I+1,J,:)
        SF(3,:)=-JFaceVector(I,J,:)
        SF(4,:)=JFaceVector(I,J+1,:)

        rF(1,:)=IFaceCenter(I,J,:)
        rF(2,:)=IFaceCenter(I+1,J,:)
        rF(3,:)=JFaceCenter(I,J,:)
        rF(4,:)=JFaceCenter(I,J+1,:)

        NeighCell(1,:)=[i-1,j]
        NeighCell(2,:)=[I+1,j]
        NeighCell(3,:)=[i,j-1]
        NeighCell(4,:)=[i,j+1]

            dtau_c = 0.0 !конвективный шаг
            dtau_d = 0.0 !диффузионный шаг
            dtau(i,j) = 0.0

            DO IFACE=1,4
	        IN=NeighCell(iface,1)
	        JN=NeighCell(iface,2)

	        d=Norm2(rF(iface,:)-CellCenter(I,j,:)) !расстояние от центра ячейки до центра грани
	        dn=Norm2(rF(iface,:)-CellCenter(IN,JN,:)) !расстояние до центра соседней ячейки

	        Vf(1)=rLinearInterp(d,dn,VFlos(I,J,1),VFlos(IN,JN,1))
	        Vf(2)=rLinearInterp(d,dn,VFlos(I,J,2),VFlos(IN,JN,2))

            Tf = RLinearInterp(d,dn,T(i,j),T(in,jn)) !central

              select case(scheme)
                case(0)  !нулевой конквективный поток
                  Vf(:) = 0.0
				  
				case (2) ! схема QUICK
					if (dot_product(Vf(:),Sf(iface,:)) > 0.0) then !Знак расхода на грани
						Tf = T(i,j) + dot_product(GradT(i,j,:),rf(iface,:) - CellCenter(i,j,:))
					else
						Tf = T(in,jn) + dot_product(GradT(in,jn,:),rf(iface,:) - CellCenter(in,jn,:))
						if (dn < 1e-7) then !для приграничных ячеек
							Tf = 2*T(in,jn) - T(i,j) -&
							4*T(in,jn) + 4*T(i,j) +&
							3*dot_product(GradT(i,j,:),rf(iface,:)-CellCenter(i,j,:))
						end if
					end if
					Tf = 0.05*(Tf + RLinearInterp(d,dn,T(i,j),T(in,jn)))
              end select
			  
              VTf(:) = Vf(:)*Tf !итоговое ковективное слагаемое
			  
	!диффузионный поток
    DNC=Norm2(CellCenter(IN,JN,:)-CellCenter(I,J,:))
    dTdn=(T(IN,JN)-T(I,J))/DNC
	RNC=(CellCenter(IN,JN,:)-CellCenter(I,J,:))/DNC

	rE(1)=rLinearInterp(d,dn,CellCenter(I,J,1),CellCenter(IN,JN,1))
	rE(2)=rLinearInterp(d,dn,CellCenter(I,J,2),CellCenter(IN,JN,2))

	GTE(1)=rLinearInterp(d,dn,GradT(I,J,1),GradT(IN,JN,1))
	GTE(2)=rLinearInterp(d,dn,GradT(I,J,2),GradT(IN,JN,2))

    NF(:)=SF(iface,:)/Norm2(SF(iface,:))

   if (dn.le.1e-7) then !для заграничных ячеек

	   dTdn_c=DOT_PRODUCT(GradT(i,j,:),NF(:))
	   dTdn=5./3.*dTdn-2./3.*dTdn_c   !второй порядок
       GTE(:)=GradT(i,j,:)	   !на границе берем значение из центра ячейки
	endif

    dTdn=dTdn+DOT_PRODUCT(NF(:)-RNC(:),GTE(:))

              if ((IN == 0) .or. (IN == NI)) then  !ГУ - адиабатические левая и правая стенки
                dTdn = 0.0
                T(IN,JN) = T(i,j) + DNC*3./5.*(2./3.*dot_product(GradT(i,j,:),Nf(:)) -&
                dot_product(GradT(i,j,:),Nf(:) - RNC(:)))
              end if

              ResT(i,j) = ResT(i,j) + dot_product(VTf(:),Sf(iface,:)) - dTdn*norm2(Sf(iface,:))*Koef_Temp !невязка на м3

              dtau_c = dtau_c + abs(dot_product(Vf(:),Sf(iface,:))) !конвективный шаг
              dtau_d = dtau_d + Koef_Temp/norm2(Sf(iface,:)) !диффузионный шаг
            end do

            ResT(i,j) = ResT(i,j)/CellVolume(i,j)
			!шаг по псевдовремени
            if (dtau_c > 1e-7) then
              dtau(i,j) = dtau(i,j) + dtau_c/(CFL*CellVolume(i,j))
            end if

            if (dtau_d > 1e-7) then
              dtau(i,j) = dtau(i,j) + 2*dtau_d/VNM
            end if

            dtau(i,j) = 1.0/dtau(i,j)

          end do
        end do
!$OMP END DO

End Subroutine



Subroutine B_OutputFields(IO, NI, NJ, X, Y, P, GradP, GradP_Error,V,DivV,DivV_Exact,&
    DivV_Error,LapP,RLapP_Exact,LapP_Error,RotV,RotV_Exact,RotV_Error,&
    VFlos,T,GradT,ResT,TFlos,T_Error)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ):: P, DivV, DivV_Exact, DivV_Error,LapP,&
  RLapP_Exact,LapP_Error, RotV, RotV_Exact, RotV_Error,T,ResT,TFlos,T_Error
  real,dimension(0:NI,0:NJ,2):: GradP, V,GradT,VFlos
  real,dimension(0:NI,0:NJ,2):: GradP_Error


  Write(IO,*) 'VARIABLES = "X", "Y", "P", "GradP_X", "GradP_Y", "GradP_X_Error", "GradP_Y_Error","V_x","V_y",&
   "DivV","DivV_Exact","DivV_Error","LapP","LapP_Exact","LapP_Error","RotV","RotV_Exact","RotV_Error",&
   "VFlos_x","VFlos_y","T","GradT_X", "GradT_Y","ResT","TFlos","T_Error"'
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ)
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,1) !GradP_X
  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,2) !GradP_Y
  Write(IO,'(100F14.7)') GradP_Error(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradP_Error(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1) !V_X
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2) !V_Y
  Write(IO,'(100F14.7)') DivV(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') DivV_Exact(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') DivV_Error(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') LapP(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') RLapP_Exact(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') LapP_Error(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') RotV(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') RotV_Exact(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') RotV_Error(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') VFlos(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') VFlos(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') T(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') GradT(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradT(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') ResT(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') TFlos(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') T_Error(1:NI-1,1:NJ-1)


End Subroutine