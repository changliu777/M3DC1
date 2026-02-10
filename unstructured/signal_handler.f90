module signal_handler
  use iso_c_binding
  implicit none

  logical(c_bool), volatile :: checkpoint_flag = .false.
  logical(c_bool) :: writing_output = .false.  ! lock to avoid double write

contains

  subroutine handle_signal(sig) bind(C)
    use iso_c_binding
    integer(c_int), value :: sig

    ! Set flag
    checkpoint_flag = .true.
  end subroutine handle_signal

end module signal_handler

