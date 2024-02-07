!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Read 2D data on a domain and resample using the axis_input.nc file
!>
program resample
  use xios
  use mpi

  implicit none

  integer :: comm = -1
  integer :: rank = -1
  integer :: npar = 0

  call initialise()
  call simulate()
  call finalise()
contains

  subroutine initialise()

    type(xios_date) :: origin
    type(xios_date) :: start
    type(xios_duration) :: tstep
    integer :: mpi_error
    integer :: lenx
    integer :: lenrx
    integer :: leny
    integer :: lenry
    integer :: isegment
    integer :: jsegment
    integer :: iseg_st
    integer :: jseg_st
    integer :: iseg_n
    integer :: jseg_n
    integer :: irseg_st
    integer :: jrseg_st
    integer :: irseg_n
    integer :: jrseg_n

    ! Arbitrary datetime setup, required for XIOS but unused
    origin = xios_date(2022, 2, 2, 12, 0, 0)
    start = xios_date(2022, 12, 13, 12, 0, 0)
    tstep = xios_hour

    ! Initialise MPI and XIOS
    call MPI_INIT(mpi_error)

    call xios_initialize('client', return_comm=comm)

    call MPI_Comm_rank(comm, rank, mpi_error)
    call MPI_Comm_size(comm, npar, mpi_error)

    ! use the axis_check context to obtain sizing information on all arrays
    ! for use in defining the main context interpretation
    call xios_context_initialize('axis_check', comm)
    call xios_set_time_origin(origin)
    call xios_set_start_date(start)
    call xios_set_timestep(tstep)

    call xios_close_context_definition()

    call xios_get_axis_attr('x', n_glo=lenx)
    call xios_get_axis_attr('x_resample', n_glo=lenrx)
    call xios_get_axis_attr('y', n_glo=leny)
    call xios_get_axis_attr('y_resample', n_glo=lenry)

    ! print *, 'x, y', lenx, ', ', leny
    ! print *, 'rx, ry', lenrx, ', ', lenry

    ! initialize the main context for interacting with the data.
    call xios_context_initialize('main', comm)

    call xios_set_time_origin(origin)
    call xios_set_start_date(start)
    call xios_set_timestep(tstep)

    ! analyse domain sizes to assign sub-domains to ranks
    iseg_st = rank * (lenx / npar)
    jseg_st = rank * (leny / npar)

    if( npar == rank + 1 ) then
      iseg_n = lenx - iseg_st
      jseg_n = leny - jseg_st
    else
      iseg_n = ((rank + 1) * (lenx / npar)) - iseg_st 
      jseg_n = ((rank + 1) * (leny / npar)) - jseg_st 
    end if

    ! analyse domain sizes to assign sub-domains to ranks
    irseg_st = rank * (lenrx / npar)
    jrseg_st = rank * (lenry / npar)

    if( npar == rank + 1 ) then
      irseg_n = lenrx - irseg_st
      jrseg_n = lenry - jrseg_st
    else
      irseg_n = ((rank + 1) * (lenrx / npar)) - irseg_st 
      jrseg_n = ((rank + 1) * (lenry / npar)) - jrseg_st 
    end if

    print*, 'rank:', rank, 'i', iseg_st, '+', iseg_n, '| j', jseg_st, '+', jseg_n
    print*, 'rank:', rank, 'ir', irseg_st, '+', irseg_n, '| jr', jrseg_st, '+', jrseg_n

    call xios_set_domain_attr("original_domain", ni=iseg_n,&
                              & nj=jseg_n, ibegin=iseg_st, jbegin=jseg_st)
    call xios_set_domain_attr("resampled_domain", ni=irseg_n,&
                              & nj=jrseg_n, ibegin=irseg_st, jbegin=jrseg_st)

    call xios_close_context_definition()

  end subroutine initialise

  subroutine finalise()

    integer :: mpi_error

    ! Finalise all XIOS contexts and MPI
    call xios_set_current_context('axis_check')
    call xios_context_finalize()
    call xios_set_current_context('main')
    call xios_context_finalize()
    call MPI_Comm_free(comm, mpi_error)
    call xios_finalize()
    call MPI_Finalize(mpi_error)

  end subroutine finalise

  subroutine simulate()

    type(xios_date) :: current
    integer :: ts
    integer :: lenx
    integer :: lenrx
    integer :: leny
    integer :: lenry

    ! Allocatable arrays, size is taken from input file
    double precision, dimension (:,:), allocatable :: inodata
    double precision, dimension (:,:), allocatable :: inedata

    ! call xios_get_domain_attr('original_domain', ni_glo=lenx)
    ! call xios_get_domain_attr('original_domain', nj_glo=leny)
    ! call xios_get_domain_attr('resampled_domain', ni_glo=lenrx)
    ! call xios_get_domain_attr('resampled_domain', nj_glo=lenry)

    call xios_get_domain_attr('original_domain', ni=lenx)
    call xios_get_domain_attr('original_domain', nj=leny)
    call xios_get_domain_attr('resampled_domain', ni=lenrx)
    call xios_get_domain_attr('resampled_domain', nj=lenry)

    allocate ( inodata(leny, lenx) )
    allocate ( inedata(lenry, lenrx) )

    ! Load data from the input file
    call xios_recv_field('odatain', inodata)
    call xios_recv_field('edatain', inedata)

    do ts=1, 1
      call xios_update_calendar(ts)
      call xios_get_current_date(current)
      ! Send (copy) the original data to the output file.
      call xios_send_field('odata', inodata)
      ! Send (copy) the expected data to the output file.
      call xios_send_field('edata', inedata)
    enddo

    deallocate (inodata)
    deallocate (inedata)

  end subroutine simulate

end program resample
