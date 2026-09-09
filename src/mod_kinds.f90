module mod_kinds
  ! Central precision definitions used by every module.
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: i4 = kind(1)
end module mod_kinds
