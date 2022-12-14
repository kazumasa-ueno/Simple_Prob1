! solve the discrete Poisson equation with Dirichlet boundary conditions in the 2D unit square using Multigrid methods
!
! equations:
!   -(phi_xx + phi_yy) = rho / e0
!	  phi[boundary] = 0
!
! smoother: GS-RB relaxation
! restriction: full weighting
! prolongation: bilinear interpolation
! standard coarsening: h_{k+1} = h_k / 2
! size of coarsest grid: h_0 = 1/2

program main
	use,intrinsic :: iso_fortran_env
	implicit none

	integer(int32) :: time_begin_c,time_end_c, CountPerSec, CountMax

	integer, parameter :: k = 8 !一番細かい格子の深さ
	integer, parameter :: gamma = 1
	integer, parameter :: nu1 = 1, nu2 = 1 !smoothing のステップ数
	integer, parameter :: maxnt = 12 !mgcycの回数
	integer, parameter :: N = 2**k+1 !横方向格子数
	integer, parameter :: Nc = 2**(k-1)+1 !1段階分粗い格子の横方向格子数
	real(8), parameter :: X = 1.d0 !計算領域サイズ
	real(8), parameter :: hf = X / (N-1) !格子間隔
	real(8), parameter :: e0 = 8.85d-12 !真空の誘電率
	integer :: ios, nt
	real(8) :: Prev(N,N), tmp(N**2) !前のステップの値を格納
	
	real(8) :: phi(N,N) !計算する電位
	real(8) :: rho(N,N) !電荷密度
	real(8) :: f(N,N) !右辺の値
	real(8) :: L(-1:1,-1:1) !L on fine grid
	real(8) :: Ifc(-1:1,-1:1) !Restriction operator
	real(8) :: Icf(-1:1,-1:1) !Interpolation operation

	call system_clock(time_begin_c, CountPerSec, CountMax)

	!Initialize parameters
	call initialize(phi, rho, f, L, Ifc, Icf, N, e0, hf)

	!main loop
	do nt = 1, maxnt
		Prev(:,:) = phi(:,:)
		call MGCYC(k, gamma, phi, L, f, nu1, nu2, X, Ifc, Icf, N, Nc)
		tmp(:) = reshape(Prev(:,:) - phi(:,:),(/N**2/))
		write(*,*) (dot_product(tmp(:),tmp(:)))**0.5
	end do

  call system_clock(time_end_c)
	print *,time_begin_c,time_end_c, CountPerSec,CountMax
  print *,real(time_end_c - time_begin_c)/CountPerSec,"sec"

	open(unit=10, file="./output/output.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/output.txt"
	
	write(10,*) phi(:,:)

	stop
contains

	subroutine initialize(phi, rho, f, L, Ifc, Icf, N, e0, hf)
		implicit none

		real(8), intent(inout) :: phi(:,:), rho(:,:), f(:,:), L(-1:1,-1:1), Ifc(-1:1,-1:1), Icf(-1:1,-1:1)
		real(8), intent(in) :: e0, hf
		integer, intent(in) :: N
		integer :: i, j

		phi(:,:) = 0.d0
		rho(:,:) = 0.d0
		
		do j = 1, N
			do i = 1, N
				if (((N/2-i)**2+(N/2-j)**2)*hf**2 < 0.05**2) then
					rho(i,j) = 1.d-8
				end if
			end do
		end do

		f(:,:) =  rho(:,:) / e0

		L( :, :) = 0.d0
		L( 0, 0) = 4.d0
		L(-1, 0) = -1.d0
		L( 0,-1) = -1.d0
		L( 1, 0) = -1.d0
		L( 0, 1) = -1.d0

		Ifc( :, :) = 1.d0
		Ifc( 0, 0) = 4.d0
		Ifc(-1, 0) = 2.d0
		Ifc( 0,-1) = 2.d0
		Ifc( 1, 0) = 2.d0
		Ifc( 0, 1) = 2.d0
		Ifc(:,:) = Ifc / 16.d0

		Icf( :, :) = 1.d0
		Icf( 0, 0) = 4.d0
		Icf(-1, 0) = 2.d0
		Icf( 0,-1) = 2.d0
		Icf( 1, 0) = 2.d0
		Icf( 0, 1) = 2.d0
		Icf(:,:) = Icf / 4.d0

	end subroutine initialize

	subroutine smooth(phi, f, hf, nu)
		implicit none

		real(8), intent(in) :: f(:,:), hf
		integer, intent(in) :: nu
		real(8), intent(inout) :: phi(:,:)

		integer :: i, j, nt, N, b(2)

		b = shape(phi)
		N = b(1)

		do nt = 1, nu
			do j = 2, N-1
				do i = 2, N-1
					if(mod(i+j,2) == 0) then
						phi(i,j) = (hf**2 * f(i,j) + phi(i-1,j) + phi(i+1,j) + phi(i,j-1) + phi(i,j+1)) / 4.d0
					end if
				end do
			end do
			do j = 2, N-1
				do i = 2, N-1
					if(mod(i+j,2) /= 0) then
						phi(i,j) = (hf**2 * f(i,j) + phi(i-1,j) + phi(i+1,j) + phi(i,j-1) + phi(i,j+1)) / 4.d0
					end if
				end do
			end do
		end do

	end subroutine smooth

	subroutine stencil(M, u, N)
		implicit none

		real(8), intent(in) :: M(-1:1,-1:1)
		real(8), intent(inout) :: u(:,:)
		integer, intent(in) :: N

		integer :: i, j
		real(8) :: tmp(N,N)

		tmp(:,:) = u(:,:)
		do j = 2, N-1
			do i = 2, N-1
				tmp(i,j) = sum(M(:,:) * u(i-1:i+1,j-1:j+1))
			end do
		end do
		u(:,:) = tmp(:,:)

	end subroutine stencil

	subroutine Restriction(M, uf, uc, N, Nc)
		implicit none
		
		real(8), intent(in) :: M(-1:1,-1:1), uf(:,:)
		real(8), intent(inout) :: uc(:,:)
		integer, intent(in) :: N, Nc

		integer :: i, j
		real(8) :: tmp(Nc,Nc)

		tmp(:,:) = uc(:,:)
		do j = 2, Nc-1
			do i = 2, Nc-1
				tmp(i,j) = sum(M(:,:) * uf((i-1)*2:(i-1)*2+2,(j-1)*2:(j-1)*2+2))
			end do
		end do
		uc(:,:) = tmp(:,:)
		
	end subroutine Restriction

	subroutine Interpolation(M, uf, uc, N, Nc)
		implicit none

		real(8), intent(in) :: M(-1:1,-1:1), uc(:,:)
		real(8), intent(inout) :: uf(:,:)
		integer, intent(in) :: N, Nc

		integer :: i, j
		real(8) :: tmp(N,N)

		tmp(:,:) = uf(:,:)
		do j = 2, N-1
			do i = 2, N-1
				if( (mod(i,2) /= 0) .and. (mod(j,2) /= 0) ) then
					tmp(i,j) = M(0,0) * uc((i-1)/2+1,(j-1)/2+1)
				elseif ( (mod(i,2) /= 0 ) .and. (mod(j,2) == 0) ) then
					tmp(i,j) = M(0,1) * uc((i-1)/2+1,j/2+1) + M(0,-1) * uc((i-1)/2+1,j/2)
				elseif ( (mod(i,2) == 0 ) .and. (mod(j,2) /= 0) ) then
					tmp(i,j) = M(1,0) * uc(i/2+1,(j-1)/2+1) + M(-1,0) * uc(i/2,(j-1)/2+1)
				else
					tmp(i,j) = M(1,1) * uc(i/2+1,j/2+1) + M(1,-1) * uc(i/2+1,j/2) &
										+ M(-1,1) * uc(i/2,j/2+1) + M(-1,-1) * uc(i/2,j/2)
				end if
			end do
		end do
		uf(:,:) = tmp(:,:)

	end subroutine Interpolation

	recursive subroutine MGCYC(k,gamma,u,L,f,nu1,nu2, X,Ifc,Icf,N,Nc)
		implicit none

		integer, intent(in) :: k, gamma, nu1, nu2, N, Nc
		real(8), intent(in) :: L(-1:1,-1:1), f(:,:), X, Ifc(-1:1,-1:1), Icf(-1:1,-1:1)
		real(8), intent(inout) :: u(:,:)

		integer :: Ncc !coarser grid number
		integer :: nt
		real(8) :: df(N,N), dc(Nc,Nc), vf(N,N), vc(Nc,Nc), tmp(N,N) !defects and errors etc.
		real(8) :: hf, hc !grid size of fine and coarse grid


		hf = X/(N-1)
		hc = X/(Nc-1)

		!Presmoothing
		call smooth(u, f, hf, nu1)

		!Coarse grid correction
		!Compute the defect
		tmp(:,:) = u(:,:)
		call stencil(L/hf**2, tmp, N)
		df(:,:) = f(:,:) - tmp(:,:)

		!Restrict the defect
		dc(:,:) = 0.d0
		call Restriction(Ifc, df, dc, N, Nc)

		!Compute an approximate solution v of the defect equation on k-1
		if(k==2) then
			vc(:,:) = 0.d0
			vc(2,2) = dc(2,2)*hc**2 / 4.d0
		else
			do nt = 1, gamma
				Ncc = 2**(k-2)+1	
				vc(:,:) = 0.d0
				call MGCYC(k-1,gamma,vc,L,dc,nu1,nu2,X,Ifc,Icf,Nc,Ncc)
			end do
		end if

		!Interpolate the correction
		vf(:,:) = 0.d0
		call Interpolation(Icf, vf, vc, N, Nc)

		!Compute the corrected approximation on k
		u(:,:) = u(:,:) + vf(:,:)

		!Postsmoothing
		call smooth(u, f, hf, nu2)

	end subroutine MGCYC

end program main