MODULE common_data
    save
    integer, parameter          :: inputRe = 1000
    integer, parameter          :: N = 32
    real*8, parameter           :: Len=1.d0, u0=1.d0, nu=1.d0/real(inputRe), CFL=0.325d0, t_tilde=3.2d0
    real*8, dimension(N+2)      :: u, uOld, x
    real*8, dimension(N)        :: Hj, Fj, expVj
    real*8                      :: dt, dx, pi, Re, t
    real*8, dimension(N)        :: a, b, c, d
    real*8, dimension(N/2+1)    :: EnergySpectrum, DissipationSpectrum
END MODULE

MODULE common_fft
    save
    integer, parameter          :: nn = 32
END MODULE

PROGRAM MAIN
    USE common_data
    implicit none
    integer                     :: i, j, it, nsteps
    integer, parameter          :: Ntime = 1
    real*8                      :: Tend
    real*8                      :: ke, epsilon, PP, etSAV, dkedt
    integer                     :: timeStart, timeEnd
    real*8                      :: cc
    real*8                      :: ke_fft, epsilon_fft
    logical                     :: convergent_judge
    character(len=100)          :: fileName
    
    write(*,*) &
    '==================================================================='
    write(*,*)'Scription:'
    write(*,*)'The code is used to calculate the non-linear Burgers Equation!'
    write(*,*)'Methedology:'
    write(*,*)'Finite-Difference Method'
    write(*,*)'Spatial Scheme ---> Second Central Difference'
    write(*,*)'Time Scheme:'
    write(*,*)'Nonlinear term ---> Adams-Bashforth Scheme'
    write(*,*)'Forcing term   ---> Adams-Bashforth Scheme'
    write(*,*)'Viscous term   ---> Crank-Nicholson Scheme'
    write(*,*)'Matrix Solver:'
    write(*,*)'tri-dag and cyclic'
    write(*,*)'Numerical Integral Method:'
    write(*,*)'Trapezoidal rule'
    write(*,*)' '
    write(*,*)'CopyRight@Weiyuan Liang'
    write(*,*) &
    '==================================================================='

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                            Mesh setting
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dx = Len/real(N)
    dt = CFL*dx/u0
    Re = u0*Len/nu
    Tend = t_tilde*Len/u0

    write(*,*) &
    '==================================================================='
    write(*,*) 'Mesh Initializing !!!'
    write(*,*) 'N                ---> ', N
    write(*,*) 'dx               ---> ', dx
    write(*,*) 'dt               ---> ', dt
    write(*,*) 'Tend             ---> ', Tend
    write(*,*) 'Reynolds numbers ---> ', Re
    write(*,*) 'CFL              ---> ', CFL
    write(*,*) &
    '==================================================================='

    do i = 1, N+2
        x(i) = real(i-1)*dx
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                            file operation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! open(55, file='uSolution.dat', status='replace', action='write')
    ! write(fileName, '("spatialAveragedRe", I5.5, "N", I4.4, ".dat")') int(Re), N
    ! open(56, file=fileName, status='replace', action='write')
    ! write(fileName, '("uSteadySolutionRe", I5.5, "N", I4.4, ".dat")') int(Re), N
    ! open(57, file=fileName, status='replace', action='write')
    ! write(fileName, '("spatialAveragedInkSpaceRe", I5.5, "N", I4.4, ".dat")') int(Re), N
    ! open(58, file=fileName, status='replace', action='write')
    ! write(fileName, '("spectrumRe", I5.5, "N", I4.4, ".dat")') int(Re), N
    ! open(59, file=fileName, status='replace', action='write')
    open(60, file='output.dat', status='unknown', position='append', action='write')

    pi = 4.0d0*atan(1.0d0)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                           Initialization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) 'Initializing Velocity Field!'
    ! Initial condition
    u = 0.d0
    uOld = u

    etSAV = 0.d0
    ke=0.d0

    ! initialize the tri-matrix coefficient
    cc = real(N)*CFL/2.d0/Re
    a = -cc
    b = 1.d0 + 2.d0*cc
    c = -cc

    !
    convergent_judge = .False.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                           Time Loop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    t = 0.0d0
    nsteps = Tend/dt
    write(*,*) 'Time loop start!'
    call system_clock(timeStart)
    do it = 1, nsteps

        ! Next time step
        CALL uNewcal()

        ! Update
        uOld = u
        t = t + dt

        ! Result output
        CALL averageCal(ke, epsilon, PP, dkedt, etSAV)
        CALL calInFourierSpace(ke_fft, epsilon_fft)
        if ((mod(it, Ntime) .eq. 0) .or. (it .eq. 1)) then
            ! do j = 1, N+1
            !     write(55, 100)t, x(j), u(j)
            ! end do
            ! write(56, 101)t, ke, epsilon, PP, (PP - epsilon), dkedt
            ! write(58, 102)t, ke_fft, epsilon_fft
            ! do j = 1, N/2+1
            !     write(59, 103)t, real(j-1), EnergySpectrum(j), DissipationSpectrum(j) 
            ! end do

            if ((mod(it, 100*Ntime) .eq. 0).or. (it .eq. 1)) then
                write(*,*)'Current time=', t, 'steps=', it
                write(*,*)'ke=', ke, 'epsilon=', epsilon, 'P=', PP, 'P-epsilon', (PP - epsilon)
                write(*,*)'***<u>=', etSAV, '***d<ke>/dt', dkedt
                write(*,*)'***Calculated by FFT:'
                write(*,*)'ke=', ke_fft, 'epsilon=', epsilon_fft, 'P-epsilon', (PP - epsilon_fft)
                write(*,*) &
                          '==================================================================='
            end if
        end if

        ! judge whether converge
        if ((abs(PP-epsilon_fft) .LT. 2e-3) .and. (.not.convergent_judge)) then
            if (t .GT. 1) then
                write(*,*)'This case is converged !!!'
                convergent_judge = .True.
                write(60, 104)int(Re), N
            end if
        elseif ((abs(PP-epsilon_fft) .GT. 2e-3) .and. it .eq. nsteps) then
            print *, 'Finally, this case is not converged !!!'
        end if

        ! Output the steady solution
        ! if (it .eq. nsteps) then
        !     do j = 1, N+1
        !         write(57, 100)t, x(j), u(j)
        !     end do
        ! end if


    end do
    write(*,*) &
    '==================================================================='
    write(*,*) 'Time loop end!'
    call system_clock(timeEnd)
    write(*,*)'Run time = ', (timeEnd-timeStart)/1000

    100     format(2x, 3f16.12)
    101     format(2x, 6f16.12)
    102     format(2x, 3f16.12)
    103     format(2x, 4f16.12)
    104     format(2x, 2I8)

    ! close(56)
    ! close(57)
    ! close(58)
    ! close(59)
    ! close(60)

END PROGRAM MAIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          Define a some subrutine for calculating
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Matrix solver
!--------------------------------------------------------------------- 
subroutine solveMatrix(a, b, c, alpha, beta, r, u_new, n)
    implicit none
    integer, intent(in)             :: n
    real*8, intent(in)              :: a(n), b(n), c(n), alpha, beta, r(n)
    real*8, intent(inout)           :: u_new(*)

    call cyclic(a, b, c, alpha, beta, r, u_new, n)

end subroutine solveMatrix
!--------------------------------------------------------------------- 

! Hcal
!--------------------------------------------------------------------- 
subroutine Hjcal()
    USE common_data
    implicit none
    integer                             :: i, j
    real*8                              :: term1, term2

    j = 1
    do i = 2, N+1
        term1 = (u(i+1)**2.d0)/2.d0 - (u(i-1)**2.d0)/2.d0
        term2 = (uOld(i+1)**2.d0)/2.d0 - (uOld(i-1)**2.d0)/2.d0
        Hj(j) = -3.d0/2.d0*(dt/2.d0/dx)*term1 + 1.d0/2.d0*(dt/2.d0/dx)*term2
        j = j + 1
    end do
    
end subroutine Hjcal
!--------------------------------------------------------------------- 

! ExpVcal
!--------------------------------------------------------------------- 
subroutine ExpVcal()
    USE common_data
    implicit none
    integer                             :: i, j
    real*8                              :: cc

    cc = real(N)*CFL/2.d0/Re
    j = 1
    do i = 2, N+1
        expVj(j) = cc*(u(i+1) - 2.d0*u(i) + u(i-1))
        j = j + 1
    end do

end subroutine ExpVcal
!--------------------------------------------------------------------- 

! Fjcal
!--------------------------------------------------------------------- 
subroutine Fjcal()
    USE common_data
    implicit none
    integer                             :: j
    logical                             :: unsteady_switch
    real*8                              :: phi1, phi2

    unsteady_switch = .False.
    if (unsteady_switch) then
        phi1 = pi*dsin(3.d0*pi*u0*t/Len)
        phi2 = pi*dsin(5.d0*pi*u0*t/Len)
    else
        phi1 = 0.d0
        phi2 = 0.d0
    end if

    do j = 1, N
        Fj(j) = dt*u0*u0/Len*(dsin(2.d0*pi*x(j+1)/Len + phi1) + dcos(4.d0*pi*x(j+1)/Len + phi2))
    end do

end subroutine Fjcal
!--------------------------------------------------------------------- 

! uNewcal
!--------------------------------------------------------------------- 
subroutine uNewcal()
    USE common_data
    implicit none
    integer                             :: i, j
    real*8, dimension(0:N+1)            :: u_new

    Call Hjcal()
    Call ExpVcal()
    Call Fjcal()

    ! Construct b
    do i = 1, N
        d(i) = u(i+1) + Hj(i) + expVj(i) + Fj(i)
    end do

    ! Calculate the U at next time
    Call solveMatrix(a, b, c, c(N), c(N), d, u(2:N+1), N)

    ! Periodical boundary condition
    u(1) = u(N+1)
    u(N+2) = u(2)

end subroutine uNewcal
!--------------------------------------------------------------------- 

! Spatial averaged terms calculation
!--------------------------------------------------------------------- 
subroutine averageCal(ke, epsilon, PP, dkedt, etSAV)
    USE common_data
    implicit none
    real*8, intent(inout)         :: ke, epsilon, PP, dkedt, etSAV
    integer                       :: j
    real*8, dimension(N+1)        :: tempke, tempepsilon, tempPP
    real*8                        :: keOld

    ! The numerical integral are using the Trapezoidal rule

    keOld = ke

    ! Averaged kinetic energy
    tempke = 0.d0
    do j = 1, N+1
        tempke(j) = u(j)**2.d0 / 2.d0 
    end do
    ke = (tempke(1) + tempke(N+1))/2.d0*dx/Len
    do j = 2, N
        ke = ke + tempke(j)*dx/Len
    end do

    ! Averaged viscous dissipation rate
    tempepsilon = 0.d0
    tempepsilon(1) = ((u(2) - u(N))/2.d0/dx)**2.d0
    do j = 2, N+1
        tempepsilon(j) = ((u(j+1) - u(j-1))/2.d0/dx)**2.d0
    end do
    epsilon = (tempepsilon(1) + tempepsilon(N+1))/2.d0*dx/Len*nu
    do j = 2, N
        epsilon = epsilon + tempepsilon(j)*dx/Len*nu
    end do

    ! Averaged energy input by the forcing term
    tempPP = 0.d0
    do j = 1, N+1
        tempPP(j) = u0*u0/Len*(dsin(2.d0*pi*x(j)/Len) + dcos(4.d0*pi*x(j)/Len)) * u(j)
    end do
    PP = (tempPP(1) + tempPP(N+1))/2.d0*dx/Len
    do j = 2, N
        PP = PP + tempPP(j)*dx/Len
    end do

    ! Every time step Spartial Averaged velocity
    etSAV = (u(1) + u(N+1))/2.d0*dx/Len
    do j = 2, N
        etSAV = etSAV + u(j)*dx/Len
    end do

    dkedt = (ke - keOld)/dt

end subroutine averageCal
!--------------------------------------------------------------------- 

! tridag-cyclic matrix solver
!--------------------------------------------------------------------- 
subroutine tridag(a,b,c,r,u,n)
    integer j,n
    real*8, dimension(n) ::  a,b,c,r,u,gam
    real*8               ::  alpha,beta,fact,gamma,bet
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo

end subroutine tridag

subroutine cyclic(a,b,c,alpha,beta,r,u_new,n)
    integer i,n
    real*8, dimension(n) ::  a,b,c,r,u_new,bb,u,z
    real*8               ::  alpha,beta,fact,gamma
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do i=2,n-1
        bb(i)=b(i)
      enddo
      call tridag(a,bb,c,r,u_new,n)
      u(1)=gamma
      u(n)=alpha
      do i=2,n-1
        u(i)=0.
      enddo
      call tridag(a,bb,c,u,z,n)
      fact=(u_new(1)+beta*u_new(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do i=1,n
        u_new(i)=u_new(i)-fact*z(i)
      enddo

end subroutine cyclic
!--------------------------------------------------------------------- 


! 1-D Fast Fourier Transform
!--------------------------------------------------------------------- 
subroutine fft(cur)
    USE common_fft
    real, intent(inout)             :: cur(nn+2)

    integer                         :: INCX, JUMPX
    INTEGER, dimension(13)          :: IFAXX
	REAL, dimension(3*nn/2+1)       :: TRIGSX
    REAL, dimension(nn+1)           :: WORKX

    ! Define common blocks
	COMMON /TRIG/  TRIGSX
	COMMON /FAX/   IFAXX

    ! Arrays used for FFTs
    INCX = 1
    JUMPX = nn + 2


    ! Initialize arrays for FFT
    CALL FFTFAX(nn,IFAXX,TRIGSX)

    ! Fast Fourier Transform
    CALL FFT991(cur(1),WORKX,TRIGSX,IFAXX,INCX,JUMPX,nn,1,-1)

    do i = 1, nn+2
        if (abs(cur(i)) .LT. .00003 ) then
            cur(i) = 0
        end if
    end do


end subroutine fft
!--------------------------------------------------------------------- 

! 1-D Inverse Fast Fourier Transform
!--------------------------------------------------------------------- 
subroutine ifft(cur)
    USE common_fft
    real, intent(inout)             :: cur(nn+2)

    integer                         :: INCX, JUMPX
    INTEGER, dimension(13)          :: IFAXX
	REAL, dimension(3*nn/2+1)       :: TRIGSX
    REAL, dimension(nn+1)           :: WORKX

    integer                         :: j

    ! Define common blocks
	COMMON /TRIG/  TRIGSX
	COMMON /FAX/   IFAXX

    ! Arrays used for FFTs
    INCX = 1
    JUMPX = nn + 2

    ! Initialize arrays for FFT
    CALL FFTFAX(nn,IFAXX,TRIGSX)

    ! Fast Fourier Transform
    CALL FFT991(cur(1),WORKX,TRIGSX,IFAXX,INCX,JUMPX,nn,1,+1)

    do i = 1, nn+2
        if (abs(cur(i)) .LT. .00003 ) then
            cur(i) = 0
        end if
    end do

    
end subroutine ifft
!--------------------------------------------------------------------- 

! Sepctrum
!--------------------------------------------------------------------- 
subroutine calSpectrum()
    USE common_data
    implicit none
    real, dimension(N+2)              :: ukArray
    integer                           :: j, i
    real*8                            :: coeff
    real*8, dimension(N+2)            :: doubleukArray(N+2) 
    real*8, dimension(N/2+1)          :: realPart(N/2+1), imagePart(N/2+1)

    ! Convert to single precision for fft
    ukArray(1:N+2) = u(1:N+2)
    
    ! velocity component in Fourier space
    CALL fft(ukArray)
    doubleukArray(1:N+2) = ukArray(1:N+2)

    i = 1
    do j = 1, N+2, 2
        realPart(i) = doubleukArray(j)
        imagePart(i) = doubleukArray(j+1)
        i = i + 1
    end do

    ! Energy spectrum
    do j = 1, N/2+1
        EnergySpectrum(j) = realPart(j)**2.d0 + imagePart(j)**2.d0
    end do

    ! Dissipation spectrum
    coeff = 8.d0*nu*pi*pi/Len/Len
    do j = 1, N/2+1
        DissipationSpectrum(j) = coeff*real(j-1)*real(j-1) & 
                                 *(realPart(j)**2.d0 + imagePart(j)**2.d0)
    end do
    
end subroutine calSpectrum
!--------------------------------------------------------------------- 


! Calculate the kenimatic energy and dissipation rate based on spectrum
!--------------------------------------------------------------------- 
subroutine calInFourierSpace(ke_fft, epsilon_fft)
    USE common_data
    implicit none
    real*8, intent(inout)               :: ke_fft, epsilon_fft

    CALL calSpectrum()

    ke_fft = sum(EnergySpectrum)
    epsilon_fft = sum(DissipationSpectrum)
end subroutine calInFourierSpace
!--------------------------------------------------------------------- 