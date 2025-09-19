!! TEST PROGRAMM FOR XXH3 NOT ACTUAL TEST CASES !!
program test_hashmap
    use xxh3_hashmap_module
    implicit none
    
    type(hashmap_type) :: map
    integer :: i, value
    
    ! Create hashmap
    call hashmap_create(map, 1)
    
    ! Insert some values
    call hashmap_put(map, "NP_001000000.1", 42)
    call hashmap_put(map, "Gene2", 123)
    call hashmap_put(map, "Gene3", 789)
    call hashmap_put(map, "Gene3", 123)
    call hashmap_put(map, "Gene4", 789)
    call hashmap_put(map, "Gene5", 12)
    call hashmap_put(map, "Gene6", 3)
    call hashmap_put(map, "Gene8", 789)
    call hashmap_put(map, "Gene99", 4)
    call hashmap_put(map, "Gene1", 5)
    call hashmap_put(map, "Gene23", 9)
    call hashmap_put(map, "Gene0", 6)
    call hashmap_put(map, "Gene142", 12)
    
    ! Retrieve values
    value = hashmap_get(map, "NP_001000000.1")
    print *, "Gene1:", value
    
    value = hashmap_get(map, "Gene2")
    print *, "Gene2:", value
    
    value = hashmap_get(map, "Gene3")
    print *, "Gene3:", value
    
    ! Test non-existent key
    value = hashmap_get(map, "Gene4")
    print *, "Gene4:", value
    
    ! Clean up
    call hashmap_destroy(map)
end program test_hashmap