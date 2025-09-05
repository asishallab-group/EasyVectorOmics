gfortran -O3 -Wall -std=f2008 \
    src/error_handling.F90\
    src/hashing.F90\
    src/array_utils.F90\
    src/tox_gene_centroids.f90\
    src/serialize_int.F90\
    src/serialize_real.F90\
    src/serialize_char.F90\
    src/deserialize_int.F90\
    src/deserialize_real.F90\
    src/deserialize_char.F90\
    src/tox_data_read_write.F90\
    src/tox_data_tools.F90\
    src/tox_data_validation.F90\
    test/test_expression_readers.f90 -o test_expression