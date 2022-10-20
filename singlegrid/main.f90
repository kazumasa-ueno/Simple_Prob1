program main
	use,intrinsic :: iso_fortran_env
	implicit none

	integer(int32) :: time_begin_c,time_end_c, CountPerSec, CountMax

	integer, parameter :: N = 129 !横方向格子数
	real(8), parameter :: X = 1.d0 !計算領域サイズ
	real(8), parameter :: hf = X / N !細かい格子間隔
	real(8), parameter :: e0 = 8.85d-12 !真空の誘電率
	real(8), parameter :: Conv = 1.d-8 !収束判定用
	integer, parameter :: nu = 100000 !smoothing のステップ数上限
	integer :: ios
	
	real(8) :: phi(N,N) !計算する電位
	real(8) :: rho(N,N) !電荷密度
	real(8) :: f(N,N) !右辺の値
	real(8) :: L(-1:1,-1:1) !L on fine grid

	call system_clock(time_begin_c, CountPerSec, CountMax)

	call initialize(phi, rho, f, L, N, e0, hf)

	call smooth(phi, f, L, hf, Conv, nu, N)

  call system_clock(time_end_c)
	print *,time_begin_c,time_end_c, CountPerSec,CountMax
  print *,real(time_end_c - time_begin_c)/CountPerSec,"sec"

	open(unit=10, file="./output/output.txt", iostat=ios, status="replace", action="write")
	if ( ios /= 0 ) stop "Error opening file ./output/output.txt"
	
	write(10,*) phi(:,:)

	stop
contains

	subroutine initialize(phi, rho, f, L, N, e0, hf)
		implicit none

		real(8), intent(inout) :: phi(:,:), rho(:,:), f(:,:), L(-1:1,-1:1)
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

		f(:,:) = rho(:,:) / e0

		L( :, :) = 0.d0
		L( 0, 0) = 4.d0
		L(-1, 0) = -1.d0
		L( 0,-1) = -1.d0
		L( 1, 0) = -1.d0
		L( 0, 1) = -1.d0
		L(:,:) = L / hf**2

	end subroutine initialize

	subroutine smooth(phi, f, L, hf, Conv, nu, N)
		implicit none

		real(8), intent(in) :: f(:,:), L(-1:1,-1:1), hf, Conv
		integer, intent(in) :: nu, N
		real(8), intent(inout) :: phi(:,:)

		integer :: i, j, nt
		real(8) :: norm
		real(8) :: Prev(N,N) !前のループの値
		real(8) :: tmp(N**2)

		do nt = 1, nu
			Prev(:,:) = phi(:,:)
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
			tmp(:) = reshape(Prev(:,:) - phi(:,:),(/N**2/))
			norm = (dot_product(tmp(:),tmp(:)))**0.5
			if(norm < Conv) then 
				exit
			end if
		end do

	end subroutine smooth


end program main