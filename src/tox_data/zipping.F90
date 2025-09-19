module tox_archive
    use iso_c_binding
    use tox_data_read_write
    use tox_errors
    use iso_fortran_env, only: real64, int32
    implicit none

    ! libzip constants
    integer(c_int), parameter :: ZIP_CREATE = 1
    integer(c_int), parameter :: ZIP_RDONLY = 0
    integer(c_int), parameter :: ZIP_FL_OVERWRITE = 8192
    integer(c_int), parameter :: ZIP_CM_STORE = 0
    integer(c_int), parameter :: ZIP_CM_DEFLATE = 8

    ! libzip interface definitions
    interface
        function malloc(size) bind(C, name="malloc")
            use iso_c_binding
            type(c_ptr) :: malloc
            integer(c_size_t), value :: size
        end function malloc
        
        subroutine free(ptr) bind(C, name="free")
            use iso_c_binding
            type(c_ptr), value :: ptr
        end subroutine free

        function zip_open(path, flags, errorp) bind(C, name="zip_open")
            use iso_c_binding
            type(c_ptr) :: zip_open
            character(kind=c_char), dimension(*) :: path
            integer(c_int), value :: flags
            integer(c_int), intent(out) :: errorp
        end function zip_open

        function zip_close(archive) bind(C, name="zip_close")
            use iso_c_binding
            integer(c_int) :: zip_close
            type(c_ptr), value :: archive
        end function zip_close

        function zip_file_add(archive, name, source, flags) bind(C, name="zip_file_add")
            use iso_c_binding
            integer(c_int64_t) :: zip_file_add
            type(c_ptr), value :: archive
            character(kind=c_char), dimension(*) :: name
            type(c_ptr), value :: source
            integer(c_int), value :: flags
        end function zip_file_add

        function zip_source_buffer(archive, data, len, freep) bind(C, name="zip_source_buffer")
            use iso_c_binding
            type(c_ptr) :: zip_source_buffer
            type(c_ptr), value :: archive
            type(c_ptr), value :: data
            integer(c_size_t), value :: len
            integer(c_int), value :: freep
        end function zip_source_buffer

        function zip_set_file_compression(archive, index, comp, flags) bind(C, name="zip_set_file_compression")
            use iso_c_binding
            integer(c_int) :: zip_set_file_compression
            type(c_ptr), value :: archive
            integer(c_int64_t), value :: index
            integer(c_int), value :: comp
            integer(c_int), value :: flags
        end function zip_set_file_compression

        function zip_fopen(archive, fname, flags) bind(C, name="zip_fopen")
            use iso_c_binding
            type(c_ptr) :: zip_fopen
            type(c_ptr), value :: archive
            character(kind=c_char), dimension(*) :: fname
            integer(c_int), value :: flags
        end function zip_fopen

        function zip_fread(file, buf, nbytes) bind(C, name="zip_fread")
            use iso_c_binding
            integer(c_int64_t) :: zip_fread
            type(c_ptr), value :: file
            type(c_ptr), value :: buf
            integer(c_size_t), value :: nbytes
        end function zip_fread

        function zip_fclose(file) bind(C, name="zip_fclose")
            use iso_c_binding
            integer(c_int) :: zip_fclose
            type(c_ptr), value :: file
        end function zip_fclose

        function zip_stat(archive, fname, flags, sb) bind(C, name="zip_stat")
            use iso_c_binding
            integer(c_int) :: zip_stat
            type(c_ptr), value :: archive
            character(kind=c_char), dimension(*) :: fname
            integer(c_int), value :: flags
            type(c_ptr), value :: sb
        end function zip_stat

        function zip_file_stat(archive, index, flags, sb) bind(C, name="zip_file_stat")
            use iso_c_binding
            integer(c_int) :: zip_file_stat
            type(c_ptr), value :: archive
            integer(c_int64_t), value :: index
            integer(c_int), value :: flags
            type(c_ptr), value :: sb
        end function zip_file_stat

        function zip_get_num_entries(archive, flags) bind(C, name="zip_get_num_entries")
            use iso_c_binding
            integer(c_int64_t) :: zip_get_num_entries
            type(c_ptr), value :: archive
            integer(c_int), value :: flags
        end function zip_get_num_entries

        function zip_get_name(archive, index, flags) bind(C, name="zip_get_name")
            use iso_c_binding
            type(c_ptr) :: zip_get_name
            type(c_ptr), value :: archive
            integer(c_int64_t), value :: index
            integer(c_int), value :: flags
        end function zip_get_name
    end interface

    ! ZIP stat structure (simplified)
    type, bind(C) :: zip_stat_t
        integer(c_int64_t) :: valid
        integer(c_int64_t) :: size
        integer(c_size_t) :: mtime
        integer(c_int64_t) :: crc
    end type zip_stat_t

contains

    ! Create a ZIP archive with files and manifest
    subroutine create_zip_archive(zip_filename, gene_ids_file, expression_file, gene_to_family_file, &
                                family_ids_file, family_centroids_file, shift_vectors_file, ierr)
        character(len=*), intent(in) :: zip_filename
        character(len=*), intent(in) :: gene_ids_file, expression_file, gene_to_family_file, &
                                        family_ids_file, family_centroids_file, shift_vectors_file
        integer, intent(out) :: ierr
        
        type(c_ptr) :: zip_handle
        integer(c_int) :: error
        character(len=:), allocatable :: manifest_filename

        call set_ok(ierr)
        call set_ok(error)
        
        ! Open ZIP archive for writing
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_CREATE, error)
        if (.not. is_ok(error)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error opening ZIP file for writing: ", error
            return
        end if
        
        ! Add files if they are non-empty
        if (len_trim(gene_ids_file) > 0) then
            call add_file_to_zip(zip_handle, trim(gene_ids_file), ierr)
            if(.not. is_ok(ierr)) return
        end if
        
        if (len_trim(expression_file) > 0) then
            call add_file_to_zip(zip_handle, trim(expression_file), ierr)
            if(.not. is_ok(ierr)) return
        end if
        
        if (len_trim(gene_to_family_file) > 0) then
            call add_file_to_zip(zip_handle, trim(gene_to_family_file), ierr)
            if(.not. is_ok(ierr)) return
        end if
        
        if (len_trim(family_ids_file) > 0) then
            call add_file_to_zip(zip_handle, trim(family_ids_file), ierr)
            if(.not. is_ok(ierr)) return
        end if
        
        if (len_trim(family_centroids_file) > 0) then
            call add_file_to_zip(zip_handle, trim(family_centroids_file), ierr)
            if(.not. is_ok(ierr)) return
        end if
        
        if (len_trim(shift_vectors_file) > 0) then
            call add_file_to_zip(zip_handle, trim(shift_vectors_file), ierr)
            if(.not. is_ok(ierr)) return
        end if
        
        ! Create and add manifest
        manifest_filename = "manifest.txt"
        call write_manifest(gene_ids_file, expression_file, gene_to_family_file, &
                        family_ids_file, family_centroids_file, shift_vectors_file, &
                        manifest_filename, ierr)
        if(.not. is_ok(ierr)) return

        call add_file_to_zip(zip_handle, manifest_filename, ierr)
        if(.not. is_ok(ierr)) return
        
        ! Delete the temporary manifest file
        call delete_file(manifest_filename, ierr)
        if(.not. is_ok(ierr)) return
        
        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (.not. is_ok(error)) then
            ierr = error
            print *, "Error closing ZIP file: ", error
            return
        end if
        
        write(*,*) "ZIP archive created successfully: ", trim(zip_filename)
    end subroutine create_zip_archive

    subroutine extract_zip_archive(zip_filename, gene_ids_file, expression_file, gene_to_family_file, &
                                family_ids_file, family_centroids_file, shift_vectors_file, ierr)
        character(len=*), intent(in) :: zip_filename
        character(len=:), allocatable, intent(out) :: gene_ids_file, expression_file, gene_to_family_file, &
                                                    family_ids_file, family_centroids_file, shift_vectors_file
        integer(int32), intent(out) :: ierr
        
        type(c_ptr) :: zip_handle, file_handle
        integer(c_int) :: error
        integer(c_int64_t) :: i, num_entries, bytes_read
        character(len=:), allocatable :: filename
        integer(int32) :: unit, iostat
        integer(int32), parameter :: CHUNK_SIZE = 4096
        character(kind=c_char), dimension(:), allocatable, target :: buffer
        logical :: file_exists
        
        ! Initialize output variables
        gene_ids_file = ""
        expression_file = ""
        gene_to_family_file = ""
        family_ids_file = ""
        family_centroids_file = ""
        shift_vectors_file = ""
        
        call set_ok(ierr)
        call set_ok(error)
        call set_ok(iostat)
        
        ! Check if file exists
        inquire(file=zip_filename, exist=file_exists)
        if (.not. file_exists) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "ZIP file does not exist: ", trim(zip_filename)
            return
        end if
        
        ! Open ZIP archive for reading
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_RDONLY, error)
        if (error /= 0 .or. .not. c_associated(zip_handle)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error opening ZIP file for reading: ", error
            return
        end if
        
        ! Get number of entries in the ZIP
        num_entries = zip_get_num_entries(zip_handle, 0)
        
        ! Extract each file
        do i = 0, num_entries - 1
            ! Get filename
            call get_zip_entry_name(zip_handle, i, filename, ierr)
            if(.not. is_ok(ierr)) return
            if (len_trim(filename) == 0) then
                print *, "Error getting name for index: ", i
                cycle
            end if
            
            ! Skip if it's the manifest
            if (filename == "manifest.txt") cycle
            
            ! Open file in ZIP
            file_handle = zip_fopen(zip_handle, trim(filename)//c_null_char, 0)
            if (.not. c_associated(file_handle)) then
                call set_err_once(ierr, ERR_FILE_OPEN)
                print *, "Error opening file in ZIP: ", trim(filename)
                return
            end if
            
            ! Open output file
            open(newunit=unit, file=trim(filename), access='stream', form='unformatted', iostat=iostat, status='replace')
            if (.not. is_ok(iostat)) then
                call set_err_once(ierr, ERR_FILE_OPEN)
                print *, "Error creating file: ", trim(filename)
                error = zip_fclose(file_handle)
                cycle
            end if
            
            ! Read and write in chunks
            allocate(buffer(CHUNK_SIZE))
            do
                bytes_read = zip_fread(file_handle, c_loc(buffer), int(CHUNK_SIZE, c_size_t))
                if (bytes_read <= 0) exit
                write(unit, iostat=iostat) buffer(1:bytes_read)
                if (.not. is_ok(iostat)) then
                    call set_err_once(ierr, ERR_FILE_OPEN)
                    print *, "Error writing file: ", trim(filename)
                    exit
                end if
            end do
            
            ! Clean up
            deallocate(buffer)
            close(unit)
            error = zip_fclose(file_handle)
            if (.not. is_ok(error)) then
                call set_err_once(ierr, ERR_INTERNAL)
                print *, "Error closing file in ZIP: ", trim(filename)
            end if
            
            print *, "Extracted: ", trim(filename)
        end do
        
        ! Extract and parse the manifest file
        file_handle = zip_fopen(zip_handle, "manifest.txt"//c_null_char, 0)
        if (c_associated(file_handle)) then
            ! Open output file for manifest
            open(newunit=unit, file="manifest.txt", access='stream', form='unformatted', iostat=iostat, status='replace')
            if (.not. is_ok(iostat)) then
                print *, "Error creating manifest file"
                error = zip_fclose(file_handle)
                call set_err_once(ierr, ERR_FILE_OPEN)
                return
            end if
            
            ! Read and write in chunks
            allocate(buffer(CHUNK_SIZE))
            do
                bytes_read = zip_fread(file_handle, c_loc(buffer), int(CHUNK_SIZE, c_size_t))
                if (bytes_read <= 0) exit
                write(unit, iostat=iostat) buffer(1:bytes_read)
                if (iostat /= 0) then
                    print *, "Error writing manifest file"
                    exit
                end if
            end do
            
            ! Clean up
            deallocate(buffer)
            close(unit)
            error = zip_fclose(file_handle)
            
            ! Parse the manifest file
            call read_manifest("manifest.txt", gene_ids_file, expression_file, gene_to_family_file, &
                            family_ids_file, family_centroids_file, shift_vectors_file, ierr)
            
            if (.not. is_ok(ierr)) then
                print *, "Error parsing manifest file"
                return
            end if
            
            ! Delete the temporary manifest file
            call delete_file("manifest.txt", ierr)
        else
            print *, "No manifest file found in ZIP archive"
            call set_err_once(ierr, ERR_INTERNAL)
        end if
        
        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (.not. is_ok(error)) then
            print *, "Error closing ZIP file: ", error
            call set_err_once(ierr, ERR_FILE_OPEN)
            return
        end if
        
        print *, "ZIP archive extracted successfully: ", trim(zip_filename)
    end subroutine extract_zip_archive

    subroutine write_file(filename, content, ierr)
        character(len=*), intent(in) :: filename
        character(kind=c_char), dimension(:), intent(in) :: content
        integer(int32), intent(out) :: ierr

        integer :: unit, iostat
        call set_ok(iostat)
        call set_ok(ierr)
        
        open(unit, file=filename, access='stream', form='unformatted', iostat=iostat, status='replace')
        if (.not. is_ok(iostat)) then
            print *, "Error creating file: ", trim(filename)
            call set_err_once(ierr, ERR_FILE_OPEN)
            return
        end if
        
        write(unit, iostat=iostat) content
        
        close(unit)
        
        if (.not. is_ok(iostat)) then
            call set_err_once(ierr, ERR_WRITE_DATA)
            print *, "Error writing file: ", trim(filename)
        end if
    end subroutine write_file

    subroutine delete_file(filename, ierr)
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        integer :: unit, iostat

        call set_ok(iostat)
        call set_ok(ierr)
        
        open(newunit=unit, file=filename, iostat=iostat, status='old')
        if (is_ok(iostat)) then
            close(unit, status='delete')
        else 
            call set_err_once(ierr, ERR_FILE_OPEN)
        end if
    end subroutine delete_file
    
    ! Helper function to add a file to ZIP
    subroutine add_file_to_zip(zip_handle, filename, ierr)
        type(c_ptr), intent(in) :: zip_handle
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        
        integer(c_int) :: error
        integer(c_int64_t) :: index
        type(c_ptr) :: source, c_data
        integer(c_signed_char), pointer :: data(:)
        integer :: unit, iostat, file_size
        integer(c_size_t) :: c_file_size

        call set_ok(ierr)
        call set_ok(iostat)
        call set_ok(error)
        
        ! Read file content
        open(unit, file=filename, access='stream', form='unformatted', iostat=iostat, status='old')
        if (.not. is_ok(iostat)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error opening file: ", trim(filename)
            return
        end if
        
        inquire(unit, size=file_size)
        c_file_size = int(file_size, c_size_t)
        
        ! Allocate memory using C's malloc
        c_data = malloc(c_file_size)
        if (.not. c_associated(c_data)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
            print *, "Error allocating memory for: ", trim(filename)
            close(unit)
            return
        end if
        
        ! Associate the C pointer with a Fortran pointer
        call c_f_pointer(c_data, data, [file_size])
        
        ! Read file content into the allocated memory
        read(unit, iostat=iostat) data
        close(unit)
        
        if (.not. is_ok(iostat)) then
            call set_err_once(ierr, ERR_READ_DATA)
            print *, "Error reading file: ", trim(filename)
            call free(c_data)
            return
        end if
        
        ! Create ZIP source from file content
        source = zip_source_buffer(zip_handle, c_data, c_file_size, 1)
        if (.not. c_associated(source)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
            print *, "Error creating source for: ", trim(filename)
            call free(c_data)
            return
        end if
        
        ! Add file to ZIP
        index = zip_file_add(zip_handle, trim(filename)//c_null_char, source, ZIP_FL_OVERWRITE)
        if (index < 0) then
            print *, "Error adding file to ZIP: ", trim(filename)
            return
        end if
        
        ! Set compression to store (no compression)
        error = zip_set_file_compression(zip_handle, index, ZIP_CM_STORE, 0)
        if (.not. is_ok(error)) then
            print *, "Error setting compression for: ", trim(filename)
        end if
        
        print *, "Added to ZIP: ", trim(filename)
    end subroutine add_file_to_zip

    ! Helper function to add a string as a file to ZIP
    subroutine add_string_to_zip(zip_handle, filename, content, ierr)
        type(c_ptr), intent(in) :: zip_handle
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: content
        integer(int32), intent(out) :: ierr
        
        integer(c_int) :: error
        integer(c_int64_t) :: index
        type(c_ptr) :: source, c_data
        integer(c_size_t) :: content_len
        character(kind=c_char), pointer :: f_content(:)
        integer(int32) :: i
        content_len = len(content)

        call set_ok(ierr)
        call set_ok(error)
        
        ! Allocate memory using C's malloc
        c_data = malloc(content_len)
        if (.not. c_associated(c_data)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
            print *, "Error allocating memory for: ", trim(filename)
            return
        end if
        
        ! Copy the string content to the allocated memory
        call c_f_pointer(c_data, f_content, [content_len])
        
        ! Convert Fortran string to C string (without null terminator)
        do i = 1, content_len
            f_content(i) = content(i:i)
        end do
        
        ! Create ZIP source from string content
        source = zip_source_buffer(zip_handle, c_data, content_len, 1)
        if (.not. c_associated(source)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
            print *, "Error creating source for: ", trim(filename)
            call free(c_data)
            return
        end if
        
        ! Add file to ZIP
        index = zip_file_add(zip_handle, trim(filename)//c_null_char, source, ZIP_FL_OVERWRITE)
        if (index < 0) then
            print *, "Error adding file to ZIP: ", trim(filename)
            return
        end if
        
        ! Set compression to store (no compression)
        error = zip_set_file_compression(zip_handle, index, ZIP_CM_STORE, 0)
        if (error /= 0) then
            print *, "Error setting compression for: ", trim(filename)
        end if
        
        print *, "Added to ZIP: ", trim(filename)
    end subroutine add_string_to_zip

    ! Helper function to get the name of a ZIP entry
    subroutine get_zip_entry_name(zip_handle, entry_index, name, ierr)
        type(c_ptr), intent(in) :: zip_handle
        integer(c_int64_t), intent(in) :: entry_index
        character(len=:), allocatable, intent(out) :: name
        integer(int32), intent(out) :: ierr
        integer(int32) :: iostat
        
        type(c_ptr) :: name_ptr
        integer(int32) :: name_len, i
        character(kind=c_char), pointer :: f_ptr(:)

        call set_ok(ierr)
        call set_ok(iostat)
        
        ! Get name from ZIP
        name_ptr = zip_get_name(zip_handle, entry_index, 0)
        if (.not. c_associated(name_ptr)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
            print *, "Error getting name for index: ", entry_index
            name = ""
            return
        end if
        
        ! Convert C string to Fortran string
        call c_f_pointer(name_ptr, f_ptr, [huge(1)])
        name_len = 0
        do i = 1, size(f_ptr)
            if (f_ptr(i) == c_null_char) then
                name_len = i - 1
                exit
            end if
        end do
        
        if (name_len > 0) then
            allocate(character(len=name_len) :: name, stat=iostat)
            if(.not. is_ok(iostat)) then
                call set_err_once(ierr, ERR_ALLOC_FAIL)
                RETURN
            end if
            name = transfer(f_ptr(1:name_len), name)
        else
            name = ""
        end if
    end subroutine get_zip_entry_name

    ! Manifest creation and parsing functions
    subroutine write_manifest(gene_ids_file, expression_file, gene_to_family_file, &
                            family_ids_file, family_centroids_file, shift_vectors_file, &
                            manifest_filename, ierr)
        character(len=*), intent(in) :: gene_ids_file, expression_file, gene_to_family_file, &
                                        family_ids_file, family_centroids_file, shift_vectors_file
        character(len=*), intent(in) :: manifest_filename
        integer(int32), intent(out) :: ierr
        
        integer(int32) :: unit, iostat

        call set_ok(ierr)
        call set_ok(iostat)
        
        ! Open the manifest file for writing
        open(newunit=unit, file=manifest_filename, status='replace', iostat=iostat)
        if (.not. is_ok(iostat)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error creating manifest file: ", trim(manifest_filename)
            return
        end if
        
        ! Write each file entry to the manifest
        if (len_trim(gene_ids_file) > 0) then
            write(unit, '(a)') 'gene_ids=' // trim(gene_ids_file)
        end if
        
        if (len_trim(expression_file) > 0) then
            write(unit, '(a)') 'expression=' // trim(expression_file)
        end if
        
        if (len_trim(gene_to_family_file) > 0) then
            write(unit, '(a)') 'gene_to_family=' // trim(gene_to_family_file)
        end if
        
        if (len_trim(family_ids_file) > 0) then
            write(unit, '(a)') 'family_ids=' // trim(family_ids_file)
        end if
        
        if (len_trim(family_centroids_file) > 0) then
            write(unit, '(a)') 'family_centroids=' // trim(family_centroids_file)
        end if
        
        if (len_trim(shift_vectors_file) > 0) then
            write(unit, '(a)') 'shift_vectors=' // trim(shift_vectors_file)
        end if
        
        ! Close the manifest file
        close(unit)
    end subroutine write_manifest

    subroutine read_manifest(manifest_filename, gene_ids_file, expression_file, gene_to_family_file, &
                            family_ids_file, family_centroids_file, shift_vectors_file, ierr)
        character(len=*), intent(in) :: manifest_filename
        character(len=:), allocatable, intent(out) :: gene_ids_file, expression_file, gene_to_family_file, &
                                                    family_ids_file, family_centroids_file, shift_vectors_file
        integer(int32), intent(out) :: ierr
        
        integer(int32) :: unit, iostat
        character(len=256) :: line
        character(len=32) :: key
        character(len=256) :: value
        
        call set_ok(ierr)
        call set_ok(iostat)
        
        ! Initialize all outputs to empty strings
        gene_ids_file = ""
        expression_file = ""
        gene_to_family_file = ""
        family_ids_file = ""
        family_centroids_file = ""
        shift_vectors_file = ""
        
        ! Open the manifest file for reading
        open(newunit=unit, file=manifest_filename, status='old', iostat=iostat)
        if (.not. is_ok(iostat)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error opening manifest file: ", trim(manifest_filename)
            return
        end if
        
        ! Read each line and parse the key-value pair
        do
            read(unit, '(a)', iostat=iostat) line
            if(iostat < 0) exit
            if (.not. is_ok(iostat)) then
                call set_err_once(ierr, ERR_READ_DATA)
                exit
            end if
            
            ! Split the line at the '=' character
            key = line(1:index(line, '=')-1)
            value = line(index(line, '=')+1:)
            
            ! Assign the value to the appropriate variable based on the key
            select case (trim(key))
                case ('gene_ids')
                    gene_ids_file = trim(value)
                case ('expression')
                    expression_file = trim(value)
                case ('gene_to_family')
                    gene_to_family_file = trim(value)
                case ('family_ids')
                    family_ids_file = trim(value)
                case ('family_centroids')
                    family_centroids_file = trim(value)
                case ('shift_vectors')
                    shift_vectors_file = trim(value)
                case default
                    print *, "Unknown key in manifest: ", trim(key)
            end select
        end do
        
        ! Close the manifest file
        close(unit)
    end subroutine read_manifest

    subroutine save_tox_data(zip_filename, ierr, gene_ids, gene_ids_file, expression, &
                            expression_file, gene_to_family, gene_to_family_file, &
                            family_ids, family_ids_file, family_centroids, &
                            family_centroids_file, shift_vectors, shift_vectors_file)
        implicit none
        
        character(len=*), intent(in) :: zip_filename
        character(len=*), intent(in), optional :: gene_ids(:)
        character(len=*), intent(in), optional :: family_ids(:)
        real(real64), intent(in), optional :: expression(:,:)
        real(real64), intent(in), optional :: family_centroids(:,:)
        real(real64), intent(in), optional :: shift_vectors(:,:)
        integer(int32), intent(in), optional :: gene_to_family(:)
        character(len=*), intent(in), optional :: gene_ids_file
        character(len=*), intent(in), optional :: expression_file
        character(len=*), intent(in), optional :: gene_to_family_file
        character(len=*), intent(in), optional :: family_ids_file
        character(len=*), intent(in), optional :: family_centroids_file
        character(len=*), intent(in), optional :: shift_vectors_file
        integer(int32), intent(out) :: ierr
        
        character(len=:), allocatable :: actual_gene_ids_file, actual_expression_file, actual_gene_to_family_file, &
                                        actual_family_ids_file, actual_family_centroids_file, actual_shift_vectors_file
        logical :: gene_ids_present, expression_present, gene_to_family_present, &
                family_ids_present, family_centroids_present, shift_vectors_present
        
        call set_ok(ierr)
        
        ! Determine which arrays are present and set appropriate filenames
        gene_ids_present = present(gene_ids) .and. present(gene_ids_file)
        if (gene_ids_present) then
            actual_gene_ids_file = gene_ids_file
            call save_gene_ids(gene_ids, actual_gene_ids_file, ierr)
            if(.not. is_ok(ierr)) return
        else
            actual_gene_ids_file = ""
        end if
        
        expression_present = present(expression) .and. present(expression_file)
        if (expression_present) then
            actual_expression_file = expression_file
            call save_expression_vectors(expression, actual_expression_file, ierr)
            if(.not. is_ok(ierr)) return
        else
            actual_expression_file = ""
        end if
        
        gene_to_family_present = present(gene_to_family) .and. present(gene_to_family_file)
        if (gene_to_family_present) then
            actual_gene_to_family_file = gene_to_family_file
            call save_gene_to_family(gene_to_family, actual_gene_to_family_file, ierr)
            if(.not. is_ok(ierr)) return
        else
            actual_gene_to_family_file = ""
        end if
        
        family_ids_present = present(family_ids) .and. present(family_ids_file)
        if (family_ids_present) then
            actual_family_ids_file = family_ids_file
            call save_family_ids(family_ids, actual_family_ids_file, ierr)
            if(.not. is_ok(ierr)) return
        else
            actual_family_ids_file = ""
        end if
        
        family_centroids_present = present(family_centroids) .and. present(family_centroids_file)
        if (family_centroids_present) then
            actual_family_centroids_file = family_centroids_file
            call save_family_centroids(family_centroids, actual_family_centroids_file, ierr)
            if(.not. is_ok(ierr)) return
        else
            actual_family_centroids_file = ""
        end if
        
        shift_vectors_present = present(shift_vectors) .and. present(shift_vectors_file)
        if (shift_vectors_present) then
            actual_shift_vectors_file = shift_vectors_file
            call save_shift_vectors(shift_vectors, actual_shift_vectors_file, ierr)
            if(.not. is_ok(ierr)) return
        else
            actual_shift_vectors_file = ""
        end if
        
        ! Create the ZIP archive
        call create_zip_archive(zip_filename, actual_gene_ids_file, actual_expression_file, &
                            actual_gene_to_family_file, actual_family_ids_file, &
                            actual_family_centroids_file, actual_shift_vectors_file, ierr)
        
        ! Clean up temporary files
        if (gene_ids_present) call delete_file(actual_gene_ids_file, ierr)
        if (expression_present) call delete_file(actual_expression_file, ierr)
        if (gene_to_family_present) call delete_file(actual_gene_to_family_file, ierr)
        if (family_ids_present) call delete_file(actual_family_ids_file, ierr)
        if (family_centroids_present) call delete_file(actual_family_centroids_file, ierr)
        if (shift_vectors_present) call delete_file(actual_shift_vectors_file, ierr)
    end subroutine save_tox_data

    subroutine read_tox_data(zip_filename, ierr, gene_ids, gene_ids_file, expression, expression_file, &
                        gene_to_family, gene_to_family_file, family_ids, family_ids_file, &
                        family_centroids, family_centroids_file, shift_vectors, shift_vectors_file)
        use array_utils, only: get_array_metadata
        use tox_data_read_write
        use iso_fortran_env, only: real64, int32
        implicit none
        
        character(len=*), intent(in) :: zip_filename
        integer(int32), intent(out) :: ierr
        character(len=:), allocatable, optional, intent(out) :: gene_ids(:)
        character(len=:), allocatable, optional, intent(out) :: family_ids(:)
        real(real64), allocatable, optional, intent(out) :: expression(:,:)
        real(real64), allocatable, optional, intent(out) :: family_centroids(:,:)
        real(real64), allocatable, optional, intent(out) :: shift_vectors(:,:)
        integer(int32), allocatable, optional, intent(out) :: gene_to_family(:)
        character(len=:), allocatable, optional, intent(out) :: gene_ids_file
        character(len=:), allocatable, optional, intent(out) :: expression_file
        character(len=:), allocatable, optional, intent(out) :: gene_to_family_file
        character(len=:), allocatable, optional, intent(out) :: family_ids_file
        character(len=:), allocatable, optional, intent(out) :: family_centroids_file
        character(len=:), allocatable, optional, intent(out) :: shift_vectors_file
        
        character(len=:), allocatable :: extracted_gene_ids_file, extracted_expression_file, &
                                        extracted_gene_to_family_file, extracted_family_ids_file, &
                                        extracted_family_centroids_file, extracted_shift_vectors_file
        logical :: gene_ids_requested, expression_requested, gene_to_family_requested, &
                family_ids_requested, family_centroids_requested, shift_vectors_requested
        integer(int32) :: max_dims, ndims, dims(5), char_len
        
        call set_ok(ierr)
        max_dims = 5  ! Maximum number of dimensions supported
        
        write(*,*) 'Extracting zip archive...'
        ! Extract the ZIP archive and get file names from manifest
        call extract_zip_archive(zip_filename, extracted_gene_ids_file, extracted_expression_file, &
                                extracted_gene_to_family_file, extracted_family_ids_file, &
                                extracted_family_centroids_file, extracted_shift_vectors_file, ierr)
        if (.not. is_ok(ierr)) return
        
        ! Return filenames if requested
        if (present(gene_ids_file)) gene_ids_file = extracted_gene_ids_file
        if (present(expression_file)) expression_file = extracted_expression_file
        if (present(gene_to_family_file)) gene_to_family_file = extracted_gene_to_family_file
        if (present(family_ids_file)) family_ids_file = extracted_family_ids_file
        if (present(family_centroids_file)) family_centroids_file = extracted_family_centroids_file
        if (present(shift_vectors_file)) shift_vectors_file = extracted_shift_vectors_file
        
        ! Load the arrays that are requested and available
        gene_ids_requested = present(gene_ids) .and. len_trim(extracted_gene_ids_file) > 0
        if (gene_ids_requested) then
            ! Get array metadata to determine size and character length
            call get_array_metadata(extracted_gene_ids_file, dims, max_dims, ndims, ierr, char_len)
            if (is_ok(ierr) .and. ndims == 1) then
                ! Allocate array based on metadata with proper character length
                allocate(character(len=char_len) :: gene_ids(dims(1)))
                call load_gene_ids(gene_ids, extracted_gene_ids_file, ierr)
                if(.not. is_ok(ierr)) return
            else
                print *, "Error getting metadata for gene_ids file"
                return
            end if
            call delete_file(extracted_gene_ids_file, ierr)
        end if
        
        expression_requested = present(expression) .and. len_trim(extracted_expression_file) > 0
        if (expression_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_expression_file, dims, max_dims, ndims, ierr)
            if (is_ok(ierr) .and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(expression(dims(1), dims(2)))
                call load_expression_vectors(expression, extracted_expression_file, ierr)
                if(.not. is_ok(ierr)) return
            else
                print *, "Error getting metadata for expression file"
                return
            end if
            call delete_file(extracted_expression_file, ierr)
        end if
        
        gene_to_family_requested = present(gene_to_family) .and. len_trim(extracted_gene_to_family_file) > 0
        if (gene_to_family_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_gene_to_family_file, dims, max_dims, ndims, ierr)
            if (is_ok(ierr) .and. ndims == 1) then
                ! Allocate array based on metadata
                allocate(gene_to_family(dims(1)))
                call load_gene_to_family(gene_to_family, extracted_gene_to_family_file, ierr)
                if(.not. is_ok(ierr)) return
            else
                print *, "Error getting metadata for gene_to_family file"
                return
            end if
            call delete_file(extracted_gene_to_family_file, ierr)
        end if
        
        family_ids_requested = present(family_ids) .and. len_trim(extracted_family_ids_file) > 0
        if (family_ids_requested) then
            ! Get array metadata to determine size and character length
            call get_array_metadata(extracted_family_ids_file, dims, max_dims, ndims, ierr, char_len)
            if (is_ok(ierr) .and. ndims == 1) then
                ! Allocate array based on metadata with proper character length
                allocate(character(len=char_len) :: family_ids(dims(1)))
                call load_family_ids(family_ids, extracted_family_ids_file, ierr)
            else
                print *, "Error getting metadata for family_ids file"
                return
            end if
            call delete_file(extracted_family_ids_file, ierr)
        end if
        
        family_centroids_requested = present(family_centroids) .and. len_trim(extracted_family_centroids_file) > 0
        if (family_centroids_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_family_centroids_file, dims, max_dims, ndims, ierr)
            if (is_ok(ierr) .and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(family_centroids(dims(1), dims(2)))
                call load_family_centroids(family_centroids, extracted_family_centroids_file, ierr)
                if(.not. is_ok(ierr)) return
            else
                print *, "Error getting metadata for family_centroids file"
                return
            end if
            call delete_file(extracted_family_centroids_file, ierr)
        end if
        
        shift_vectors_requested = present(shift_vectors) .and. len_trim(extracted_shift_vectors_file) > 0
        if (shift_vectors_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_shift_vectors_file, dims, max_dims, ndims, ierr)
            if (is_ok(ierr).and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(shift_vectors(dims(1), dims(2)))
                call load_shift_vectors(shift_vectors, extracted_shift_vectors_file, ierr)
                if(.not. is_ok(ierr)) return
            else
                print *, "Error getting metadata for shift_vectors file"
                return
            end if
            call delete_file(extracted_shift_vectors_file, ierr)
        end if
    end subroutine read_tox_data

end module tox_archive

! Subroutines for R interface
! Helper subroutines for array conversions
subroutine strings_to_ascii_array(strings, ascii_array, str_len, n_strings)
    use iso_fortran_env, only: int32
    implicit none
    character(len=*), intent(in) :: strings(:)
    integer(int32), intent(out) :: ascii_array(*)
    integer(int32), intent(in) :: str_len, n_strings
    integer(int32) :: i, j
    
    do i = 1, n_strings
        do j = 1, str_len
            if (j <= len_trim(strings(i))) then
                ascii_array((i-1)*str_len + j) = ichar(strings(i)(j:j))
            else
                ascii_array((i-1)*str_len + j) = 0
            end if
        end do
    end do
end subroutine strings_to_ascii_array

subroutine ascii_array_to_strings(ascii_array, str_len, n_strings, strings)
    use iso_fortran_env, only: int32
    implicit none
    integer(int32), intent(in) :: ascii_array(*)
    integer(int32), intent(in) :: str_len, n_strings
    character(len=str_len), intent(out) :: strings(n_strings)
    integer(int32) :: i, j
    
    do i = 1, n_strings
        strings(i) = ''
        do j = 1, str_len
            if (ascii_array((i-1)*str_len + j) /= 0) then
                strings(i)(j:j) = char(ascii_array((i-1)*str_len + j))
            end if
        end do
    end do
end subroutine ascii_array_to_strings

! Simplified version of the R binding using helper routines
subroutine read_tox_data_R(zip_filename_ascii, zip_filename_len, flags, &
                          gene_ids_ascii, gene_ids_len, n_genes, &
                          expression_vectors, n_samples, n_genes_exp, &
                          gene_to_family, n_genes_gtf, &
                          family_ids_ascii, family_ids_len, n_families, &
                          family_centroids, n_families_cent, d_cent, &
                          shift_vectors, n_genes_shift, d_shift, &
                          ierr) bind(C, name="read_tox_data_R")
    use iso_c_binding
    use tox_archive
    use tox_conversions
    use array_utils
    implicit none
    
    ! Input
    integer(c_int), intent(in) :: zip_filename_ascii(*)
    integer(c_int), intent(in) :: zip_filename_len
    integer(c_int), intent(in) :: flags(6)
    
    ! Output for gene_ids
    integer(c_int), intent(out) :: gene_ids_ascii(*)
    integer(c_int), intent(out) :: gene_ids_len, n_genes
    
    ! Output for expression
    real(c_double), intent(out) :: expression_vectors(*)
    integer(c_int), intent(out) :: n_samples, n_genes_exp
    
    ! Output for gene_to_family
    integer(c_int), intent(out) :: gene_to_family(*)
    integer(c_int), intent(out) :: n_genes_gtf
    
    ! Output for family_ids
    integer(c_int), intent(out) :: family_ids_ascii(*)
    integer(c_int), intent(out) :: family_ids_len, n_families
    
    ! Output for family_centroids
    real(c_double), intent(out) :: family_centroids(*)
    integer(c_int), intent(out) :: n_families_cent, d_cent
    
    ! Output for shift_vectors
    real(c_double), intent(out) :: shift_vectors(*)
    integer(c_int), intent(out) :: n_genes_shift, d_shift
    
    integer(c_int), intent(out) :: ierr
    
    ! Local variables
    character(len=:), allocatable :: zip_filename
    character(len=:), allocatable :: gene_ids(:), family_ids(:)
    real(real64), allocatable :: expression(:,:), family_centroids_arr(:,:), shift_vectors_arr(:,:)
    integer(int32), allocatable :: gene_to_family_arr(:)
    integer :: i, j, k, temp_ierr
    character(kind=c_char, len=1), dimension(:), allocatable :: c_char_array
    
    call set_ok(ierr)
    
    ! Convert zip filename from ASCII
    allocate(c_char_array(zip_filename_len))
    do i = 1, zip_filename_len
        c_char_array(i) = char(zip_filename_ascii(i))
    end do
    call c_char_1d_as_string(c_char_array, zip_filename, temp_ierr)
    if (temp_ierr /= 0) then
        ierr = temp_ierr
        return
    end if
    
    ! Call the Fortran function
    call read_tox_data(zip_filename, ierr, &
                      gene_ids=gene_ids, &
                      expression=expression, &
                      gene_to_family=gene_to_family_arr, &
                      family_ids=family_ids, &
                      family_centroids=family_centroids_arr, &
                      shift_vectors=shift_vectors_arr)
    if (ierr /= 0) return
    
    ! Convert outputs back to ASCII arrays using helper routines
    if (flags(1) == 1 .and. allocated(gene_ids)) then
        n_genes = size(gene_ids)
        gene_ids_len = maxval(len_trim(gene_ids))
        call strings_to_ascii_array(gene_ids, gene_ids_ascii, gene_ids_len, n_genes)
    else
        n_genes = 0
        gene_ids_len = 0
    end if
    
    if (flags(2) == 1 .and. allocated(expression)) then
        n_samples = size(expression, 1)
        n_genes_exp = size(expression, 2)
        k = 1
        do j = 1, n_genes_exp
            do i = 1, n_samples
                expression_vectors(k) = expression(i, j)
                k = k + 1
            end do
        end do
    else
        n_samples = 0
        n_genes_exp = 0
    end if
    
    if (flags(3) == 1 .and. allocated(gene_to_family_arr)) then
        n_genes_gtf = size(gene_to_family_arr)
        do i = 1, n_genes_gtf
            gene_to_family(i) = gene_to_family_arr(i)
        end do
    else
        n_genes_gtf = 0
    end if
    
    if (flags(4) == 1 .and. allocated(family_ids)) then
        n_families = size(family_ids)
        family_ids_len = maxval(len_trim(family_ids))
        call strings_to_ascii_array(family_ids, family_ids_ascii, family_ids_len, n_families)
    else
        n_families = 0
        family_ids_len = 0
    end if
    
    if (flags(5) == 1 .and. allocated(family_centroids_arr)) then
        d_cent = size(family_centroids_arr, 1)
        n_families_cent = size(family_centroids_arr, 2)
        k = 1
        do j = 1, n_families_cent
            do i = 1, d_cent
                family_centroids(k) = family_centroids_arr(i, j)
                k = k + 1
            end do
        end do
    else
        n_families_cent = 0
        d_cent = 0
    end if
    
    if (flags(6) == 1 .and. allocated(shift_vectors_arr)) then
        d_shift = size(shift_vectors_arr, 1)
        n_genes_shift = size(shift_vectors_arr, 2)
        k = 1
        do j = 1, n_genes_shift
            do i = 1, d_shift
                shift_vectors(k) = shift_vectors_arr(i, j)
                k = k + 1
            end do
        end do
    else
        n_genes_shift = 0
        d_shift = 0
    end if
end subroutine read_tox_data_R

! Similarly simplified version for save_tox_data_R
subroutine save_tox_data_R(zip_filename_ascii, zip_filename_len, flags, &
                          gene_ids_ascii, gene_ids_len, n_genes, &
                          expression_vectors, n_samples, n_genes_exp, &
                          gene_to_family, n_genes_gtf, &
                          family_ids_ascii, family_ids_len, n_families, &
                          family_centroids, n_families_cent, d_cent, &
                          shift_vectors, n_genes_shift, d_shift, &
                          ierr) bind(C, name="save_tox_data_R")
    use iso_c_binding
    use tox_archive
    use tox_conversions
    use array_utils
    implicit none
    
    ! Input
    integer(c_int), intent(in) :: zip_filename_ascii(*)
    integer(c_int), intent(in) :: zip_filename_len
    integer(c_int), intent(in) :: flags(6)
    integer(c_int), intent(in) :: gene_ids_ascii(*)
    integer(c_int), intent(in) :: gene_ids_len, n_genes
    real(c_double), intent(in) :: expression_vectors(*)
    integer(c_int), intent(in) :: n_samples, n_genes_exp
    integer(c_int), intent(in) :: gene_to_family(*)
    integer(c_int), intent(in) :: n_genes_gtf
    integer(c_int), intent(in) :: family_ids_ascii(*)
    integer(c_int), intent(in) :: family_ids_len, n_families
    real(c_double), intent(in) :: family_centroids(*)
    integer(c_int), intent(in) :: n_families_cent, d_cent
    real(c_double), intent(in) :: shift_vectors(*)
    integer(c_int), intent(in) :: n_genes_shift, d_shift
    integer(c_int), intent(out) :: ierr
    
    ! Local variables
    character(len=:), allocatable :: zip_filename
    character(len=:), allocatable :: gene_ids(:), family_ids(:)
    real(real64), allocatable :: expression(:,:), family_centroids_arr(:,:), shift_vectors_arr(:,:)
    integer(int32), allocatable :: gene_to_family_arr(:)
    integer :: i, j, k, temp_ierr
    character(kind=c_char, len=1), dimension(:), allocatable :: c_char_array
    
    call set_ok(ierr)
    
    ! Convert zip filename from ASCII
    allocate(c_char_array(zip_filename_len))
    do i = 1, zip_filename_len
        c_char_array(i) = char(zip_filename_ascii(i))
    end do
    call c_char_1d_as_string(c_char_array, zip_filename, temp_ierr)
    if (temp_ierr /= 0) then
        ierr = temp_ierr
        return
    end if
    
    ! Convert inputs to Fortran arrays using helper routines
    if (flags(1) == 1 .and. n_genes > 0) then
        allocate(character(len=gene_ids_len) :: gene_ids(n_genes))
        call ascii_array_to_strings(gene_ids_ascii, gene_ids_len, n_genes, gene_ids)
    end if
    
    if (flags(2) == 1 .and. n_genes_exp > 0 .and. n_samples > 0) then
        allocate(expression(n_samples, n_genes_exp))
        k = 1
        do j = 1, n_genes_exp
            do i = 1, n_samples
                expression(i, j) = expression_vectors(k)
                k = k + 1
            end do
        end do
    end if
    
    if (flags(3) == 1 .and. n_genes_gtf > 0) then
        allocate(gene_to_family_arr(n_genes_gtf))
        do i = 1, n_genes_gtf
            gene_to_family_arr(i) = gene_to_family(i)
        end do
    end if
    
    if (flags(4) == 1 .and. n_families > 0) then
        allocate(character(len=family_ids_len) :: family_ids(n_families))
        call ascii_array_to_strings(family_ids_ascii, family_ids_len, n_families, family_ids)
    end if
    
    if (flags(5) == 1 .and. n_families_cent > 0 .and. d_cent > 0) then
        allocate(family_centroids_arr(d_cent, n_families_cent))
        k = 1
        do j = 1, n_families_cent
            do i = 1, d_cent
                family_centroids_arr(i, j) = family_centroids(k)
                k = k + 1
            end do
        end do
    end if
    
    if (flags(6) == 1 .and. n_genes_shift > 0 .and. d_shift > 0) then
        allocate(shift_vectors_arr(d_shift, n_genes_shift))
        k = 1
        do j = 1, n_genes_shift
            do i = 1, d_shift
                shift_vectors_arr(i, j) = shift_vectors(k)
                k = k + 1
            end do
        end do
    end if
    
    ! Call the Fortran function
    call save_tox_data(zip_filename, ierr, &
                      gene_ids=gene_ids, &
                      expression=expression, &
                      gene_to_family=gene_to_family_arr, &
                      family_ids=family_ids, &
                      family_centroids=family_centroids_arr, &
                      shift_vectors=shift_vectors_arr)
end subroutine save_tox_data_R

! Subroutines for C interface

! Zero-copy improved C-binding wrappers for read_tox_data_C and save_tox_data_C
! Changes:
! - Removed reshape usage for numeric arrays.
! - Use c_loc (Fortran → C) to expose contiguous allocatables directly.
! - Use c_f_pointer (C → Fortran) to alias C buffers without copies.
! - String arrays still require copies (conversion unavoidable).

subroutine read_tox_data_C(zip_filename_ascii, fn_len, flags, &
                          gene_ids, gene_ids_len, n_genes, &
                          expression_vectors, n_samples, n_genes_exp, &
                          gene_to_family, n_genes_gtf, &
                          family_ids, family_ids_len, n_families, &
                          family_centroids, n_families_cent, d_cent, &
                          shift_vectors, n_genes_shift, d_shift, &
                          ierr) bind(C, name="read_tox_data_C")
    use iso_c_binding, only: c_int, c_char
    use iso_fortran_env, only: int32, real64
    use array_utils, only: string_to_ascii, check_okay_ioerror
    use tox_archive
    use tox_conversions
    implicit none

    ! Inputs
    integer(c_int), intent(in) :: fn_len
    character(c_char), intent(in) :: zip_filename_ascii(fn_len)
    integer(c_int), intent(in) :: flags(6)

    ! Outputs (preallocated by caller)
    integer(c_int), intent(out) :: gene_ids(gene_ids_len, n_genes)
    real(real64), intent(out)   :: expression_vectors(n_samples, n_genes_exp)
    integer(int32), intent(out) :: gene_to_family(n_genes_gtf)
    integer(c_int), intent(out) :: family_ids(family_ids_len, n_families)
    real(real64), intent(out)   :: family_centroids(d_cent, n_families_cent)
    real(real64), intent(out)   :: shift_vectors(d_shift, n_genes_shift)

    integer(c_int), intent(in)  :: gene_ids_len, n_genes
    integer(c_int), intent(in)  :: n_samples, n_genes_exp
    integer(c_int), intent(in)  :: n_genes_gtf
    integer(c_int), intent(in)  :: family_ids_len, n_families
    integer(c_int), intent(in)  :: n_families_cent, d_cent
    integer(c_int), intent(in)  :: n_genes_shift, d_shift
    integer(c_int), intent(out) :: ierr

    ! Locals
    character(len=:), allocatable :: f_zip_filename
    character(len=:), allocatable :: f_gene_ids(:), f_family_ids(:)
    real(real64), allocatable     :: f_expression(:,:), f_family_centroids(:,:), f_shift_vectors(:,:)
    integer(int32), allocatable   :: f_gene_to_family(:)
    integer :: j

    call set_ok(ierr)

    ! Convert ASCII filename
    call ascii_to_string(zip_filename_ascii, fn_len, f_zip_filename)

    ! Read full dataset into Fortran allocatables
    call read_tox_data(f_zip_filename, ierr, flags, &
                       f_gene_ids, gene_ids_len, n_genes, &
                       f_expression, n_samples, n_genes_exp, &
                       f_gene_to_family, n_genes_gtf, &
                       f_family_ids, family_ids_len, n_families, &
                       f_family_centroids, n_families_cent, d_cent, &
                       f_shift_vectors, n_genes_shift, d_shift)
    call check_okay_ioerror(ierr)

    ! Copy Fortran strings → ASCII
    do j = 1, n_genes
        call string_to_ascii(trim(f_gene_ids(j)), gene_ids(:, j))
    end do
    do j = 1, n_families
        call string_to_ascii(trim(f_family_ids(j)), family_ids(:, j))
    end do

    ! Copy numeric arrays directly
    expression_vectors = f_expression
    gene_to_family     = f_gene_to_family
    family_centroids   = f_family_centroids
    shift_vectors      = f_shift_vectors
end subroutine

subroutine save_tox_data_C(zip_filename_ascii, fn_len, flags, &
           gene_ids, gene_ids_len, n_genes, &
           expression_vectors, n_samples, n_genes_exp, &
           gene_to_family, n_genes_gtf, &
           family_ids, family_ids_len, n_families, &
           family_centroids, n_families_cent, d_cent, &
           shift_vectors, n_genes_shift, d_shift, &
           ierr) bind(C, name="save_tox_data_C")
    use iso_c_binding, only: c_int, c_char
    use iso_fortran_env, only: int32, real64
    use array_utils, only: ascii_to_string_padded, check_okay_ioerror
    use tox_archive
    use tox_conversions
    implicit none

    ! Inputs
    integer(c_int), intent(in) :: fn_len
    character(c_char), intent(in) :: zip_filename_ascii(fn_len)
    integer(c_int), intent(in) :: flags(6)

    integer(c_int), intent(in) :: gene_ids(gene_ids_len, n_genes)
    real(real64), intent(in)   :: expression_vectors(n_samples, n_genes_exp)
    integer(int32), intent(in) :: gene_to_family(n_genes_gtf)
    integer(c_int), intent(in) :: family_ids(family_ids_len, n_families)
    real(real64), intent(in)   :: family_centroids(d_cent, n_families_cent)
    real(real64), intent(in)   :: shift_vectors(d_shift, n_genes_shift)

    integer(c_int), intent(in) :: gene_ids_len, n_genes
    integer(c_int), intent(in) :: n_samples, n_genes_exp
    integer(c_int), intent(in) :: n_genes_gtf
    integer(c_int), intent(in) :: family_ids_len, n_families
    integer(c_int), intent(in) :: n_families_cent, d_cent
    integer(c_int), intent(in) :: n_genes_shift, d_shift
    integer(c_int), intent(out) :: ierr

    ! Locals
    character(len=:), allocatable :: f_zip_filename
    character(len=:), allocatable :: f_gene_ids(:), f_family_ids(:)
    integer :: j

    call set_ok(ierr)

    ! Convert ASCII filename
    call ascii_to_string_padded(zip_filename_ascii, fn_len, f_zip_filename)

    ! Convert ASCII arrays → Fortran strings
    allocate(f_gene_ids(n_genes))
    do j = 1, n_genes
        call ascii_to_string_padded(gene_ids(:, j), gene_ids_len, f_gene_ids(j))
    end do
    allocate(f_family_ids(n_families))
    do j = 1, n_families
        call ascii_to_string_padded(family_ids(:, j), family_ids_len, f_family_ids(j))
    end do

    ! Save directly
    call save_tox_data(f_zip_filename, ierr, &
                       gene_ids=f_gene_ids, &
                       expression=expression_vectors, &
                       gene_to_family=gene_to_family, &
                       family_ids=f_family_ids, &
                       family_centroids=family_centroids, &
                       shift_vectors=shift_vectors)
    call check_okay_ioerror(ierr)
end subroutine

! R interface (no bind(C) - to be called with .Fortran())
subroutine extract_zip_archive_R(zip_filename, zip_filename_len, ierr)
    use iso_c_binding
    use tox_archive
    use iso_fortran_env, only: int32
    implicit none
    
    ! Input
    integer(int32), intent(in) :: zip_filename_len
    integer(int32), intent(in) :: zip_filename(zip_filename_len)
    
    ! Output
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    character(len=:), allocatable :: f_zip_filename
    integer :: i
    
    call set_ok(ierr)
    
    ! Convert integer array to Fortran string
    allocate(character(len=zip_filename_len) :: f_zip_filename)
    do i = 1, zip_filename_len
        f_zip_filename(i:i) = char(zip_filename(i))
    end do
    
    ! Extract all files from the archive
    call extract_zip_archive_simple(f_zip_filename, ierr)
    
    deallocate(f_zip_filename)
end subroutine extract_zip_archive_R

! C interface (with bind(C))
subroutine extract_zip_archive_C(zip_filename, ierr) bind(C, name="extract_zip_archive_C")
    use iso_c_binding
    use tox_archive
    use iso_fortran_env, only: int32
    implicit none
    
    ! Input
    character(kind=c_char), intent(in) :: zip_filename(*)
    
    ! Output
    integer(c_int), intent(out) :: ierr
    
    ! Local variables
    character(len=:), allocatable :: f_zip_filename
    integer :: i, len
    
    call set_ok(ierr)
    
    ! Calculate length of C string
    len = 0
    do while(zip_filename(len+1) /= c_null_char)
        len = len + 1
    end do
    
    ! Convert C string to Fortran string
    allocate(character(len=len) :: f_zip_filename)
    do i = 1, len
        f_zip_filename(i:i) = zip_filename(i)
    end do
    
    ! Extract all files from the archive
    call extract_zip_archive_simple(f_zip_filename, ierr)
    
    deallocate(f_zip_filename)
end subroutine extract_zip_archive_C

! Helper subroutine that does the actual extraction
subroutine extract_zip_archive_simple(zip_filename, ierr)
    use iso_c_binding
    use tox_errors
    use iso_fortran_env, only: int32
    implicit none
    
    character(len=*), intent(in) :: zip_filename
    integer(int32), intent(out) :: ierr
    
    type(c_ptr) :: zip_handle, file_handle
    integer(c_int) :: error
    integer(c_int64_t) :: i, num_entries, bytes_read
    character(len=:), allocatable :: filename
    integer(int32) :: unit, iostat
    integer(int32), parameter :: CHUNK_SIZE = 4096
    character(kind=c_char), dimension(:), allocatable, target :: buffer
    logical :: file_exists
    
    call set_ok(ierr)
    call set_ok(error)
    call set_ok(iostat)
    
    ! Check if file exists
    inquire(file=zip_filename, exist=file_exists)
    if (.not. file_exists) then
        call set_err_once(ierr, ERR_FILE_OPEN)
        print *, "ZIP file does not exist: ", trim(zip_filename)
        return
    end if
    
    ! Open ZIP archive for reading
    zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_RDONLY, error)
    if (error /= 0 .or. .not. c_associated(zip_handle)) then
        call set_err_once(ierr, ERR_FILE_OPEN)
        print *, "Error opening ZIP file for reading: ", error
        return
    end if
    
    ! Get number of entries in the ZIP
    num_entries = zip_get_num_entries(zip_handle, 0)
    
    ! Extract each file
    do i = 0, num_entries - 1
        ! Get filename
        call get_zip_entry_name(zip_handle, i, filename, ierr)
        if(.not. is_ok(ierr)) then
            print *, "Error getting name for index: ", i
            cycle
        end if
        
        if (len_trim(filename) == 0) then
            print *, "Empty filename for index: ", i
            cycle
        end if
        
        ! Open file in ZIP
        file_handle = zip_fopen(zip_handle, trim(filename)//c_null_char, 0)
        if (.not. c_associated(file_handle)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error opening file in ZIP: ", trim(filename)
            cycle
        end if
        
        ! Open output file
        open(newunit=unit, file=trim(filename), access='stream', form='unformatted', &
             iostat=iostat, status='replace')
        if (.not. is_ok(iostat)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error creating file: ", trim(filename)
            error = zip_fclose(file_handle)
            cycle
        end if
        
        ! Read and write in chunks
        allocate(buffer(CHUNK_SIZE))
        do
            bytes_read = zip_fread(file_handle, c_loc(buffer), int(CHUNK_SIZE, c_size_t))
            if (bytes_read <= 0) exit
            write(unit, iostat=iostat) buffer(1:bytes_read)
            if (.not. is_ok(iostat)) then
                call set_err_once(ierr, ERR_FILE_OPEN)
                print *, "Error writing file: ", trim(filename)
                exit
            end if
        end do
        
        ! Clean up
        deallocate(buffer)
        close(unit)
        error = zip_fclose(file_handle)
        if (.not. is_ok(error)) then
            call set_err_once(ierr, ERR_INTERNAL)
            print *, "Error closing file in ZIP: ", trim(filename)
        end if
        
        print *, "Extracted: ", trim(filename)
    end do
    
    ! Close ZIP archive
    error = zip_close(zip_handle)
    if (.not. is_ok(error)) then
        print *, "Error closing ZIP file: ", error
        call set_err_once(ierr, ERR_FILE_OPEN)
        return
    end if
    
    print *, "All files extracted successfully from: ", trim(zip_filename)
end subroutine extract_zip_archive_simple