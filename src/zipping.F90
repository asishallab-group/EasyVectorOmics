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