module palm_test
    use topography_mod, only: init_topography

    implicit none

    private

    public :: say_hello
contains
    subroutine say_hello
        call init_topography
    end subroutine say_hello
end module palm_test
