module io_ascii

    use fabm_driver
    use fabm_config_types
    use fabm_yaml,yaml_parse=>parse,yaml_error_length=>error_length
    use fabm_types, only: rk

    implicit none
    private
    
    type brom_par ! Type for sparse matrix
        real(rk) :: realvalue
        character(len=64) :: initname
        type(brom_par), pointer :: next
    end type brom_par
    
    type (brom_par), pointer :: first
    
    public find_index, porting_initial_state_variables_data, &
        init_brom_par, get_brom_par, close_brom_par, &
        saving_state_variables_data 
contains
    
    subroutine init_brom_par()
    !from brom.yaml initialization
    
        class (type_node), pointer       :: node
        character(len=yaml_error_length) :: yaml_error
        integer                          :: unit_eff
        character(len=256)               :: path_eff
             
        unit_eff = 1
        path_eff = 'brom.yaml'
        ! Parse YAML file.
        node => yaml_parse(trim(path_eff),unit_eff,yaml_error) 
        if (yaml_error/='') call fatal_error('fabm_create_model_from_yaml_file',trim(yaml_error))
        if (.not.associated(node)) call fatal_error('fabm_create_model_from_yaml_file', &
             'No configuration information found in '//trim(path_eff)//'.')
        select type (node)
             class is (type_dictionary)
                call tree(node)
             class is (type_node)
                call fatal_error('brom', trim(path_eff)//' must contain a dictionary &
                   &at the root (non-indented) level, not a single value. Are you missing a trailing colon?')
        end select
    
    end subroutine init_brom_par
    
    subroutine tree(mapping)
    
        class (type_dictionary), intent(in) :: mapping
        class (type_node), pointer          :: node
        type (type_key_value_pair), pointer :: pair
        character(len=64)                   :: instancename
        
        node => mapping%get('instances')
        if (.not.associated(node)) &
            call fatal_error('create_model_tree_from_dictionary', 'No "instances" dictionary found at root level.')
        select type (node)
            class is (type_dictionary)
                pair => node%first
            class is (type_node)
                nullify(pair)
                call fatal_error('create_model_tree_from_dictionary',trim(node%path)// &
                   ' must be a dictionary with (model name : information) pairs, not a single value.')
        end select
            
        ! only brom acceptance, iterate throw other models with no action
        do while (associated(pair))
            instancename = trim(pair%key)
            if (instancename == "brom") then
                select type (dict=>pair%value)
                    class is (type_dictionary)
                        call from_tree(instancename, dict)
                    class is (type_null)
                        call fatal_error('create_model_tree_from_dictionary','Configuration information for model "'// &
                            trim(instancename)//'" must be a dictionary, not a single value or nothing')
                    class is (type_node)
                        call fatal_error('create_model_tree_from_dictionary','Configuration information for model "'// &
                            trim(instancename)//'" must be a dictionary, not a single value.')
                    end select
            end if
            pair => pair%next
        end do
    
    end subroutine tree
    
    subroutine from_tree(instancename, node)
    
        character(len=*),       intent(in)          :: instancename
        class (type_dictionary),intent(in)          :: node
        class (type_dictionary),pointer             :: childmap
        type (type_error),pointer                   :: config_error
        
        childmap => node%get_dictionary('initialization',required=.false.,error=config_error)
        if (associated(config_error)) call fatal_error('create_model_from_dictionary',config_error%message)
        if (associated(childmap)) call parse_initialization(childmap)
      
    end subroutine from_tree
    
    subroutine parse_initialization(node)
    
        class (type_dictionary),intent(in)  :: node
        type (type_key_value_pair), pointer :: pair
        type (brom_par), pointer :: current
        
        logical                             :: success
        real(rk)                            :: realvalue
        character(len=64)                   :: initname
            
        ! Transfer user-specified initial state to the model.
        nullify(first)
        pair => node%first
        do while (associated(pair))
            select type (value=>pair%value)
                class is (type_scalar)
                    initname = trim(pair%key)
                    realvalue = value%to_real(default=real(0,real_kind),success=success)
                    if (.not.success) call fatal_error('parse_initialization', &
                        trim(value%path)//': "'//trim(value%string)//'" is not a real number.')
                class is (type_null)
                    call fatal_error('parse_initialization',&
                        trim(value%path)//' must be set to a real number, not to null.')
                class is (type_dictionary)
                    call fatal_error('parse_initialization',&
                        trim(value%path)//' must be set to a real number, not to a dictionary.')
                end select
            allocate(current)
            current = brom_par(realvalue, initname, first)
            first => current
            pair => pair%next
        end do

    end subroutine parse_initialization
    
    function get_brom_par(initname) result (realvalue)
    
    real(rk)                            :: realvalue
    character(len=*), intent(in)        :: initname
    type (brom_par), pointer            :: current
    
        current => first
        do
            if (.not.associated(current)) exit
            if (trim(initname) == trim(current%initname)) then
                realvalue = current%realvalue
                return
            end if
            current => current%next
        end do
        
        print *, "Check brom.yaml or name of the input variable please"
        stop
    
    end function get_brom_par

    subroutine close_brom_par()
        
        nullify(first)

    end subroutine close_brom_par
    
    function find_index(names, name) result(index)
    
        character(len=*), intent(in) :: names(:)
        character(len=*), intent(in) :: name
        integer :: index
        
        do index = 1, size(names)
            if (names(index) == name) return
        end do
        index = -1
    
    end function
    
    subroutine porting_initial_state_variables_data(lev_max, par_name, cc)
    
        integer, intent(in)                     :: lev_max
        character(len = *), intent(in)          :: par_name(:)
        real(rk), dimension(:, :), intent(out)  :: cc
        integer                                 :: d_year_start, julianday_start, lev_max_start, par_max_start
        integer, allocatable                    :: column2state(:)
        character(20000)                        :: labels
        character(200)                          :: comments
        integer                                 :: i, j, istart, istop, foo
        real(rk)                                :: value
        
        open(9, file = 'input.dat')
        read(9, '( 1x,i4,1x,i5,1x,i5,1x,i5 )') d_year_start, julianday_start, lev_max_start, par_max_start
        if (lev_max_start /= lev_max) stop 'number of levels in start-up file does not match number of levels in model'
        allocate(column2state(par_max_start))

        !"labels" is a space-separated string with state variable names.
        !find indiviudal names and look up their state variable index.
        read(9, '( 1x, a )') labels
        istart = 1
        do i = 1, par_max_start
            istop = istart - 1 + index(labels(istart:), ' ')
            column2state(i) = find_index(par_name, labels(istart:istop-1))
            istart = istop + 1
        end do
        read (9, *) comments
        do i = 1, lev_max_start
            read(9, '( i5 )', advance = 'no') foo
            read(9, '( 1x,f6.2 )', advance = 'no') value
            do j = 1, par_max_start
                read(9, '( 1x, f15.9 )', advance = 'no') value
                if (column2state(j) /= -1) then
                    if (value /= 0.) then
                        cc(i, column2state(j)) = value
                    else
                        cc(i, column2state(j)) = 0.0001
                    end if
                end if
            end do
            read(9, *)
        end do
        deallocate(column2state)
        close (9)
    
    end subroutine porting_initial_state_variables_data
    
    subroutine saving_state_variables_data(dYear, julianday, LevMax, ParMax, ParName, z, Cc)
    
        integer, intent(in)                         :: dYear, julianday, LevMax, ParMax
        real(rk), dimension(:, :), intent(in)       :: Cc
        character(len = *), intent(in)              :: ParName(:)
        real(rk), intent(in)                        :: z(:)
        
        integer :: ip, k
        
        open(10,FILE = 'output.dat') 

        ! First line with year, julian day, number of levels, number of state variables
        write(10,'( i4,1x,i5,1x,i5,1x,i5 )') dYear, julianday, LevMax, ParMax
        do ip=1, ParMax
            write(10,'(1x,a)',advance='NO') trim(ParName(ip))
        end do
        write (10,*)
        
        write(10,'(3h n ,6hDepth )',advance='NO')
        do ip=1, ParMax
        !       write(10,'(1x,a)',advance='NO') trim(ParName(ip)(len_trim(ParName(ip))-3:))
            write(10,'(1x,a)',advance='NO') trim(ParName(ip)(15:))
        end do
        write (10,*)

        ! Subsequent lines: one for each depth level
        do k=1, LevMax
            write(10,'( i5 )',advance='NO') k
            write(10,'(1x,f6.2 )',advance='NO') z(k)
            do ip=1,ParMax
                write(10,'( 1x,f15.9 )',advance='NO') Cc(k,ip)
            end do
            write (10,*)
        end do
        close (10)
    
    end subroutine saving_state_variables_data
    
end module io_ascii

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
