module xxh3_hashmap_module
    use, intrinsic :: iso_c_binding
    use iso_fortran_env, only: int32, int64
    implicit none
    private
    
    public :: hashmap_type, hashmap_create, hashmap_destroy, hashmap_get, hashmap_put
    
    ! C interface for XXH3 hashing
    interface
        function xxh3_hash_c(key, length) bind(C, name="xxh3_hash_c")
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: key
            integer(c_int), value :: length
            integer(c_int64_t) :: xxh3_hash_c
        end function xxh3_hash_c
    end interface

    ! Node type for separate chaining
    type :: hashmap_node_type
        character(len=:), allocatable :: key
        integer(int32) :: value
        type(hashmap_node_type), pointer :: next => null()
    end type hashmap_node_type
    
    ! Parameters
    integer, parameter :: DEFAULT_KEY_LENGTH = 256
    real, parameter :: MAX_LOAD_FACTOR = 0.75
    logical, parameter :: DEBUG = .false.
    
    ! Hash table structure
    type :: hashmap_type
        integer(int32) :: size = 0
        integer(int32) :: count = 0
        type(hashmap_node_type), pointer :: buckets(:) => null()
    end type hashmap_type

contains

! Find the next power of two greater than or equal to n
function next_power_of_two(n) result(power)
    integer(int32), intent(in) :: n
    integer(int32) :: power
    
    power = 1
    do while (power < n)
        power = ishft(power, 1)
    end do
end function next_power_of_two

! Create the hashmap
subroutine hashmap_create(map, initial_size)
    type(hashmap_type), intent(out) :: map
    integer(int32), intent(in), optional :: initial_size
    
    integer(int32) :: table_size, i
    
    ! Calculate table size (power of two)
    if (present(initial_size)) then
        table_size = next_power_of_two(int(initial_size / MAX_LOAD_FACTOR, int32))
    else
        table_size = 16  ! Default size
    end if
    
    ! Ensure minimum size
    table_size = max(table_size, 16)
    
    ! Allocate buckets
    allocate(map%buckets(table_size))
    
    ! Initialize buckets to null
    do i = 1, table_size
        nullify(map%buckets(i)%next)
    end do
    
    map%size = table_size
    map%count = 0
    
    if (DEBUG) print *, "Hashmap created with size:", table_size
end subroutine hashmap_create

! Destroy the hashmap
subroutine hashmap_destroy(map)
    type(hashmap_type), intent(inout) :: map
    
    integer(int32) :: i
    type(hashmap_node_type), pointer :: current, next
    
    if (associated(map%buckets)) then
        do i = 1, map%size
            current => map%buckets(i)%next
            do while (associated(current))
                next => current%next
                deallocate(current)
                current => next
            end do
        end do
        deallocate(map%buckets)
    end if
    
    map%size = 0
    map%count = 0
end subroutine hashmap_destroy

! Compute XXH3 hash of a string
function xxh3_hash_fortran(key, table_size) result(hash_idx)
    character(len=*), intent(in) :: key
    integer(int32), intent(in) :: table_size
    integer(int32) :: hash_idx
    
    integer(int64) :: hash_val
    integer :: key_len
    character(len=:), allocatable, target :: trimmed_key
    
    trimmed_key = normalize_gene_id(key)
    key_len = len(trimmed_key)
    
    hash_val = xxh3_hash_c(c_loc(trimmed_key), key_len)
    hash_idx = int(iand(hash_val, int(table_size - 1, int64)) + 1, int32)
end function xxh3_hash_fortran

! Insert a key-value pair
subroutine hashmap_put(map, key, value)
    type(hashmap_type), intent(inout) :: map
    character(len=*), intent(in) :: key
    integer(int32), intent(in) :: value
    
    integer(int32) :: hash_idx
    type(hashmap_node_type), pointer :: new_node, current
    character(len=:), allocatable :: normalized_key
    
    ! Normalize key
    normalized_key = normalize_gene_id(key)
    
    if (DEBUG) print *, "PUT: ", trim(normalized_key), " -> ", value
    
    ! Check if we need to resize
    if (real(map%count) / real(map%size) > MAX_LOAD_FACTOR) then
        call resize_hashmap(map)
    end if
    
    ! Compute hash index
    hash_idx = xxh3_hash_fortran(normalized_key, map%size)
    
    if (DEBUG) print *, "  Hash index:", hash_idx
    
    ! Check if key already exists
    current => map%buckets(hash_idx)%next
    do while (associated(current))
        if (current%key == normalized_key) then
            ! Key exists - update value
            current%value = value
            if (DEBUG) print *, "  Updated existing key"
            return
        end if
        current => current%next
    end do
    
    ! Key doesn't exist - create new node
    allocate(new_node)
    new_node%key = normalized_key
    new_node%value = value
    new_node%next => map%buckets(hash_idx)%next
    map%buckets(hash_idx)%next => new_node
    map%count = map%count + 1
    
    if (DEBUG) print *, "  Added new key, count:", map%count
end subroutine hashmap_put

! Lookup a key
function hashmap_get(map, key) result(value)
    type(hashmap_type), intent(in) :: map
    character(len=*), intent(in) :: key
    integer(int32) :: value
    
    integer(int32) :: hash_idx
    type(hashmap_node_type), pointer :: current
    character(len=:), allocatable :: normalized_key
    
    value = 0  ! Default: not found
    
    ! Normalize key
    normalized_key = normalize_gene_id(key)
    
    if (DEBUG) print *, "GET: ", trim(normalized_key)
    
    ! Compute hash index
    hash_idx = xxh3_hash_fortran(normalized_key, map%size)
    
    if (DEBUG) print *, "  Hash index:", hash_idx
    
    ! Search for the key in the bucket
    current => map%buckets(hash_idx)%next
    do while (associated(current))
        if (current%key == normalized_key) then
            value = current%value
            if (DEBUG) print *, "  FOUND: ", value
            return
        end if
        current => current%next
    end do
    
    if (DEBUG) print *, "  NOT FOUND"
end function hashmap_get

! Resize the hashmap when load factor is too high
subroutine resize_hashmap(map)
    type(hashmap_type), intent(inout) :: map
    
    type(hashmap_type) :: new_map
    integer(int32) :: i, new_size
    type(hashmap_node_type), pointer :: current, next
    
    if (DEBUG) print *, "Resizing hashmap from ", map%size, " to ", map%size * 2
    
    ! Create a new larger map
    new_size = map%size * 2
    call hashmap_create(new_map, new_size)
    
    ! Reinsert all elements
    do i = 1, map%size
        current => map%buckets(i)%next
        do while (associated(current))
            call hashmap_put(new_map, current%key, current%value)
            next => current%next
            deallocate(current)
            current => next
        end do
    end do
    
    ! Replace the old map with the new one
    deallocate(map%buckets)
    map%size = new_map%size
    map%count = new_map%count
    map%buckets => new_map%buckets
    
    ! Nullify the new_map's buckets to prevent deallocation
    nullify(new_map%buckets)
end subroutine resize_hashmap

! Normalize gene ID by trimming spaces and converting to uppercase
function normalize_gene_id(gene_id) result(normalized)
    character(len=*), intent(in) :: gene_id
    character(len=:), allocatable :: normalized
    
    integer :: i, j
    character(len=len(gene_id)) :: buffer
    
    buffer = trim(adjustl(gene_id))
    
    ! Remove any internal whitespace
    j = 1
    do i = 1, len_trim(buffer)
        if (buffer(i:i) /= ' ' .and. buffer(i:i) /= char(9)) then
            buffer(j:j) = buffer(i:i)
            j = j + 1
        end if
    end do
    
    ! Convert to uppercase
    buffer(1:j-1) = to_upper(buffer(1:j-1))
    normalized = trim(buffer(1:j-1))
end function normalize_gene_id

! Convert string to uppercase
function to_upper(str) result(upper)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: upper
    
    integer :: i, ic
    
    do i = 1, len(str)
        ic = iachar(str(i:i))
        if (ic >= iachar('a') .and. ic <= iachar('z')) then
            upper(i:i) = achar(ic - 32)
        else
            upper(i:i) = str(i:i)
        end if
    end do
end function to_upper

end module xxh3_hashmap_module