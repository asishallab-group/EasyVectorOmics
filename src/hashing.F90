module hashmap_module
    use iso_fortran_env, only: int32, int64
    implicit none
    private
    
    public :: hashmap_type, hashmap_create, hashmap_destroy, hashmap_get, hashmap_put
    
    type :: hashmap_node_type
        character(len=:), allocatable :: key
        integer(int32) :: value
        type(hashmap_node_type), pointer :: next => null()
    end type hashmap_node_type
    
    type :: hashmap_type
        type(hashmap_node_type), pointer :: buckets(:) => null()
        integer(int32) :: size = 0
    end type hashmap_type

contains

! DJB2 Hash Algorithmus
integer function djb2_hash(key, table_size) result(hash)
    character(len=*), intent(in) :: key
    integer(int32), intent(in) :: table_size
    integer(int32) :: i
    integer(int64) :: hash_val
    character(len=1) :: c
    
    hash_val = 257
    do i = 1, len_trim(key)
        c = key(i:i)
        hash_val = (ishft(hash_val, 5) + hash_val) + iachar(c)  ! hash * 33 + c
    end do
    
    hash = mod(abs(hash_val), table_size) + 1
end function djb2_hash

! Hashmap erstellen
subroutine hashmap_create(map, size)
    type(hashmap_type), intent(out) :: map
    integer(int32), intent(in) :: size
    integer(int32) :: i
    
    allocate(map%buckets(size))
    do i = 1, size
        nullify(map%buckets(i)%next)  ! Korrekte Initialisierung
    end do
    map%size = size
end subroutine hashmap_create

! Hashmap zerstören
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
end subroutine hashmap_destroy

! Wert in Hashmap speichern
subroutine hashmap_put(map, key, value)
    type(hashmap_type), intent(inout) :: map
    character(len=*), intent(in) :: key
    integer(int32), intent(in) :: value
    
    integer(int32) :: hash_idx
    type(hashmap_node_type), pointer :: new_node, current
    character(len=:), allocatable :: trimmed_key
    
    trimmed_key = normalize_gene_id(key)
    hash_idx = djb2_hash(trimmed_key, map%size)
    
    ! Prüfe ob Key bereits existiert
    current => map%buckets(hash_idx)%next
    do while (associated(current))
        if (current%key == trimmed_key) then
            current%value = value
            return
        end if
        current => current%next
    end do
    
    ! Neuer Node erstellen
    allocate(new_node)
    new_node%key = trimmed_key
    new_node%value = value
    new_node%next => map%buckets(hash_idx)%next
    map%buckets(hash_idx)%next => new_node
end subroutine hashmap_put

! Wert aus Hashmap lesen
integer function hashmap_get(map, key) result(value)
    type(hashmap_type), intent(in) :: map
    character(len=*), intent(in) :: key
    
    integer(int32) :: hash_idx
    type(hashmap_node_type), pointer :: current
    character(len=:), allocatable :: trimmed_key
    
    value = 0  ! 0 = nicht gefunden
    trimmed_key = trim(adjustl(key))
    hash_idx = djb2_hash(trimmed_key, map%size)
    
    current => map%buckets(hash_idx)%next
    do while (associated(current))
        if (current%key == trimmed_key) then
            value = current%value
            return
        end if
        current => current%next
    end do
end function hashmap_get

function normalize_gene_id(gene_id) result(normalized)
    character(len=*), intent(in) :: gene_id
    character(len=:), allocatable :: normalized
    integer :: i, j
    character(len=len(gene_id)) :: buffer

    buffer = trim(adjustl(gene_id))
    j = 1
    do i = 1, len_trim(buffer)
        if (buffer(i:i) /= ' ' .and. buffer(i:i) /= char(9)) then
            buffer(j:j) = buffer(i:i)
            j = j + 1
        end if
    end do
    normalized = trim(buffer(1:j-1))
end function normalize_gene_id

end module hashmap_module