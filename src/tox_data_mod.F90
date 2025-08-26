!> Main module that aggregates all tox_data functionality
module tox_data_mod
    use tox_data_type
    ! use tox_data_readers
    !! under development
    ! use accessors
    !! under development
    implicit none

    public :: tox_data_t
    public :: create_tox_data, destroy_tox_data
    public :: read_gene_ids_file, read_gene_families_file, read_expression_tables
    public :: compute_family_centroids, compute_shift_vectors
    public :: validate_tox_data

end module tox_data_mod