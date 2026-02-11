module signal_handler
  use iso_c_binding
  implicit none

  logical(c_bool), volatile :: checkpoint_flag
  logical(c_bool) :: writing_output

  interface
    function c_signal(signum, handler) bind(C, name="signal")
      import :: c_int, c_funptr
      integer(c_int), value :: signum
      type(c_funptr), value :: handler
      type(c_funptr) :: c_signal
    end function c_signal
  end interface

contains

  subroutine handle_signal(slurm_sig) bind(C)
    use iso_c_binding
    integer(c_int), value :: slurm_sig

    ! Only set flag, safe inside signal handler
    checkpoint_flag = .true.
  end subroutine handle_signal

end module signal_handler
