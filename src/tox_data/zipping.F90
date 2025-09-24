module tox_archive
    use iso_c_binding
    use tox_data_read_write
    use tox_errors
    use iso_fortran_env, only: real64, int32
    implicit none

    ! libzip constants
    integer(c_int), parameter :: ZIP_CREATE = 1
    integer(c_int), parameter :: ZIP_EXCL = 3
    integer(c_int), parameter :: ZIP_RDONLY = 0
    integer(c_int), parameter :: ZIP_FL_OVERWRITE = 8192
    integer(c_int), parameter :: ZIP_CM_STORE = 0
    integer(c_int), parameter :: ZIP_CM_DEFLATE = 8
    integer(c_int), parameter :: ZIP_ER_EXISTS = -11

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

        subroutine zip_source_free(source) bind(C, name="zip_source_free")
            use iso_c_binding
            type(c_ptr), value :: source
        end subroutine zip_source_free
        
        function c_strlen(str) bind(C, name="strlen")
            use iso_c_binding
            integer(c_size_t) :: c_strlen
            type(c_ptr), value :: str
        end function c_strlen
    end interface

    ! ZIP stat structure (simplified)
    type, bind(C) :: zip_stat_t
        integer(c_int64_t) :: valid
        integer(c_int64_t) :: size
        integer(c_size_t) :: mtime
        integer(c_int64_t) :: crc
    end type zip_stat_t

contains

    !> Creates a zip archive if it does not exist already
    subroutine create_zip_archive(zip_filename, gene_ids_file, expression_file, gene_to_family_file, &
                                family_ids_file, family_centroids_file, shift_vectors_file, ierr)
        character(len=*), intent(in) :: zip_filename
            !! Name of the zip file to create
        character(len=*), intent(in) :: gene_ids_file
            !! Name of the file containing the gene_ids
        character(len=*), intent(in) :: expression_file
            !! Name of the file containing the expression vectors
        character(len=*), intent(in) :: gene_to_family_file
            !! Name of the file containing the gene to family mapping
        character(len=*), intent(in) :: family_ids_file
            !! Name of the file containing the family ids
        character(len=*), intent(in) :: family_centroids_file
            !! Name of the file containing the centroids 
        character(len=*), intent(in) :: shift_vectors_file
            !! Name of the file containing the shift vectors
        integer, intent(out) :: ierr
            !! Error code
        
        type(c_ptr) :: zip_handle
        integer(c_int) :: error
        character(len=:), allocatable :: manifest_filename

        call set_ok(ierr)
        call set_ok(error)
        
        ! Open ZIP archive
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_EXCL, error)
        if (.not. is_ok(error)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            if (error == 10) print *, "Error opening ZIP file for writing: File already exists"
            return
        end if
        
        ! Add all data files
        call add_files_to_zip(zip_handle, gene_ids_file, expression_file, gene_to_family_file, &
                            family_ids_file, family_centroids_file, shift_vectors_file, ierr)
        if (.not. is_ok(ierr)) then
            error = zip_close(zip_handle)
            return
        end if
        
        ! Add manifest
        manifest_filename = "manifest.txt"
        call write_manifest(gene_ids_file, expression_file, gene_to_family_file, &
                        family_ids_file, family_centroids_file, shift_vectors_file, &
                        manifest_filename, ierr)
        if (is_ok(ierr)) then
            call add_data_to_zip(zip_handle, manifest_filename, manifest_filename, 1, ierr)
            if (is_ok(ierr)) call delete_file(manifest_filename, ierr)
        end if
        
        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (.not. is_ok(error)) then
            call set_err_once(ierr, error)
            print *, "Error closing ZIP file: ", error
        else if (is_ok(ierr)) then
            print *, "ZIP archive created successfully: ", trim(zip_filename)
        end if
        
    contains
        subroutine add_files_to_zip(zip_handle, gid_file, exp_file, g2f_file, fam_file, cent_file, shift_file, ierr)
            type(c_ptr), intent(in) :: zip_handle
                !! open zip file
            character(len=*), intent(in) :: gid_file
                !! gene_ids file
            character(len=*), intent(in) :: exp_file 
                !! expression file
            character(len=*), intent(in) :: g2f_file 
                !! gene to family file
            character(len=*), intent(in) :: fam_file 
                !! family file
            character(len=*), intent(in) :: cent_file 
                !! centroids file
            character(len=*), intent(in) :: shift_file
                !! shift vectors file
            integer(int32), intent(out) :: ierr 
                !! Error code
            
            call set_ok(ierr)
            
            if (len_trim(gid_file) > 0) then
                call add_data_to_zip(zip_handle, gid_file, gid_file, 1, ierr)
                if (.not. is_ok(ierr)) return
            end if
            
            if (len_trim(exp_file) > 0) then
                call add_data_to_zip(zip_handle, exp_file, exp_file, 1, ierr)
                if (.not. is_ok(ierr)) return
            end if
            
            if (len_trim(g2f_file) > 0) then
                call add_data_to_zip(zip_handle, g2f_file, g2f_file, 1, ierr)
                if (.not. is_ok(ierr)) return
            end if
            
            if (len_trim(fam_file) > 0) then
                call add_data_to_zip(zip_handle, fam_file, fam_file, 1, ierr)
                if (.not. is_ok(ierr)) return
            end if
            
            if (len_trim(cent_file) > 0) then
                call add_data_to_zip(zip_handle, cent_file, cent_file, 1, ierr)
                if (.not. is_ok(ierr)) return
            end if
            
            if (len_trim(shift_file) > 0) then
                call add_data_to_zip(zip_handle, shift_file, shift_file, 1, ierr)
                if (.not. is_ok(ierr)) return
            end if
        end subroutine add_files_to_zip
    end subroutine create_zip_archive

    !> Extract a zip archive. Reads manifest and fills arrays accordingly.
    subroutine extract_zip_archive(zip_filename, gene_ids_file, expression_file, gene_to_family_file, &
                                family_ids_file, family_centroids_file, shift_vectors_file, ierr)
        character(len=*), intent(in) :: zip_filename
            !! Zip file to read
        character(len=:), allocatable, intent(out) :: gene_ids_file
            !! Gene IDs filename
        character(len=:), allocatable, intent(out) :: expression_file
            !! Expression vectors filename
        character(len=:), allocatable, intent(out) :: gene_to_family_file   
            !! gene to family filename
        character(len=:), allocatable, intent(out) :: family_ids_file
            !! Family ids filename
        character(len=:), allocatable, intent(out) :: family_centroids_file
            !! Family centroids filename
        character(len=:), allocatable, intent(out) :: shift_vectors_file
            !! Shift vectors filename
        integer(int32), intent(out) :: ierr
            !! Error code
        
        type(c_ptr) :: zip_handle
        integer(c_int) :: error
        integer(c_int64_t) :: i, num_entries
        character(len=:), allocatable :: filename
        logical :: file_exists
        
        ! Initialize outputs
        gene_ids_file = ""
        expression_file = ""
        gene_to_family_file = ""
        family_ids_file = ""
        family_centroids_file = ""
        shift_vectors_file = ""
        
        call set_ok(ierr)
        call set_ok(error)
        
        ! Check if file exists
        inquire(file=zip_filename, exist=file_exists)
        if (.not. file_exists) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "ZIP file does not exist: ", trim(zip_filename)
            return
        end if
        
        ! Open ZIP archive
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_RDONLY, error)
        if (error /= 0 .or. .not. c_associated(zip_handle)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error opening ZIP file for reading: ", error
            return
        end if
        
        ! Extract all files
        num_entries = zip_get_num_entries(zip_handle, 0)
        do i = 0, num_entries - 1
            call get_zip_entry_name(zip_handle, i, filename, ierr)
            if (.not. is_ok(ierr)) then
                error = zip_close(zip_handle)
                return 
            end if 
            
            if (filename == "manifest.txt") cycle  ! Handle manifest separately
            
            call extract_file_from_zip(zip_handle, filename, ierr)
            if (.not. is_ok(ierr)) then
                error = zip_close(zip_handle)
                RETURN
            end if
        end do
        
        ! Extract and parse manifest
        call extract_and_parse_manifest(zip_handle, gene_ids_file, expression_file, &
                                    gene_to_family_file, family_ids_file, &
                                    family_centroids_file, shift_vectors_file, ierr)
        if (.not. is_ok(ierr)) then
            error = zip_close(zip_handle)
            return 
        end if
        
        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (.not. is_ok(error)) then
            print *, "Error closing ZIP file: ", error
            call set_err_once(ierr, ERR_FILE_CLOSE)
            return
        end if
        
        print *, "ZIP archive extracted successfully: ", trim(zip_filename)
    end subroutine extract_zip_archive

    !> Delete a file from the disk
    subroutine delete_file(filename, ierr)
        character(len=*), intent(in) :: filename
            !! File to delete
        integer(int32), intent(out) :: ierr
            !! Error code
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

    !> Helper function to get the name of a ZIP entry
    subroutine get_zip_entry_name(zip_handle, entry_index, entry_name, ierr)
        type(c_ptr), intent(in) :: zip_handle
            !! Zip file connection
        integer(c_int64_t), intent(in) :: entry_index
            !! Index of the entry
        character(len=:), allocatable, intent(out) :: entry_name
            !! Name of the entry
        integer(int32), intent(out) :: ierr
            !! Error code
        
        type(c_ptr) :: name_ptr
        integer(int32) :: iostat, name_len, i
        character(kind=c_char), pointer :: f_ptr(:)
        integer, parameter :: MAX_NAME_LENGTH = 4096  ! Reasonable maximum

        call set_ok(ierr)
        call set_ok(iostat)
        
        ! Get name from ZIP
        name_ptr = zip_get_name(zip_handle, entry_index, 0)
        if (.not. c_associated(name_ptr)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
            print *, "Error getting name for index: ", entry_index
            entry_name = ""
            return
        end if
        
        ! Convert C string to Fortran string
        call c_f_pointer(name_ptr, f_ptr, [MAX_NAME_LENGTH])
        name_len = 0
        
        ! Find null terminator with bounds checking
        do i = 1, MAX_NAME_LENGTH
            if (f_ptr(i) == c_null_char) then
                name_len = i - 1
                exit
            end if
        end do
        
        ! Check if we exceeded maximum length
        if (i > MAX_NAME_LENGTH) then
            call set_err_once(ierr, ERR_STRING_TOO_LONG)
            print *, "ZIP entry name too long at index: ", entry_index
            entry_name = ""
            return
        end if
        
        if (name_len > 0) then
            allocate(character(len=name_len) :: entry_name, stat=iostat)
            if(.not. is_ok(iostat)) then
                call set_err_once(ierr, ERR_ALLOC_FAIL)
                entry_name = ""
                return
            end if
            
            ! Copy the string safely
            do i = 1, name_len
                entry_name(i:i) = f_ptr(i)
            end do
        else
            entry_name = ""
            print *, "Warning: Empty name for ZIP entry index: ", entry_index
        end if
    end subroutine get_zip_entry_name

    ! Unified subroutine to extract a file from ZIP archive
    subroutine extract_file_from_zip(zip_handle, filename, ierr)
        type(c_ptr), intent(in) :: zip_handle
            !! Zip connection
        character(len=*), intent(in) :: filename
            !! Name of the file to extract
        integer(int32), intent(out) :: ierr
            !! Error code
        
        type(c_ptr) :: file_handle
        integer(c_int) :: error
        integer(c_int64_t) :: bytes_read
        integer(int32) :: unit, iostat
        integer(int32), parameter :: CHUNK_SIZE = 4096
        character(kind=c_char), dimension(:), allocatable, target :: buffer
        
        call set_ok(ierr)
        
        ! Open file in ZIP
        file_handle = zip_fopen(zip_handle, trim(filename)//c_null_char, 0)
        if (.not. c_associated(file_handle)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error opening file in ZIP: ", trim(filename)
            return
        end if
        
        ! Open output file
        open(newunit=unit, file=trim(filename), access='stream', form='unformatted', &
            iostat=iostat, status='replace', action='write')
        if (.not. is_ok(iostat)) then
            call set_err_once(ierr, ERR_FILE_OPEN)
            print *, "Error creating file: ", trim(filename)
            error = zip_fclose(file_handle)
            return
        end if
        
        ! Read and write in chunks
        allocate(buffer(CHUNK_SIZE), stat=iostat)
        if (.not. is_ok(iostat)) then
            call set_err_once(ierr, ERR_ALLOC_FAIL)
            print *, "Error allocating buffer for: ", trim(filename)
            close(unit, status='delete')
            error = zip_fclose(file_handle)
            return
        end if
        
        do
            bytes_read = zip_fread(file_handle, c_loc(buffer), int(CHUNK_SIZE, c_size_t))
            if (bytes_read <= 0) exit
            if (bytes_read > 0) then
                write(unit, iostat=iostat) buffer(1:bytes_read)
                if (.not. is_ok(iostat)) then
                    call set_err_once(ierr, ERR_WRITE_DATA)
                    print *, "Error writing file: ", trim(filename)
                    exit
                end if
            end if
        end do
        
        ! Clean up
        if (allocated(buffer)) deallocate(buffer)
        close(unit)
        error = zip_fclose(file_handle)
        
        if (.not. is_ok(error)) then
            call set_err_once(ierr, ERR_FILE_CLOSE)
            print *, "Error closing file in ZIP: ", trim(filename)
        end if
        
        if (is_ok(ierr)) then
            print *, "Extracted: ", trim(filename)
        end if
    end subroutine extract_file_from_zip

    ! Unified subroutine to add data to ZIP (handles both files and strings)
    subroutine add_data_to_zip(zip_handle, filename, data_source, data_type, ierr)
        type(c_ptr), intent(in) :: zip_handle
            !! Zip connection
        character(len=*), intent(in) :: filename
            !! Filename to add
        character(len=*), intent(in) :: data_source 
            !! File path or string content
        integer, intent(in) :: data_type 
            !! 1 = file, 2 = string
        integer(int32), intent(out) :: ierr
            !! Error code
        
        integer(c_int) :: error
        integer(c_int64_t) :: index
        type(c_ptr) :: source, c_data
        integer(c_signed_char), pointer :: file_data(:)
        character(kind=c_char), pointer :: string_data(:)
        integer :: unit, iostat, file_size
        integer(c_size_t) :: data_len
        integer(int32) :: i
        
        call set_ok(ierr)
        
        if (len_trim(data_source) == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        select case (data_type)
        case (1) ! File data
            ! Read file content
            open(unit, file=data_source, access='stream', form='unformatted', iostat=iostat, status='old')
            if (.not. is_ok(iostat)) then
                call set_err_once(ierr, ERR_FILE_OPEN)
                print *, "Error opening file: ", trim(data_source)
                return
            end if
            
            inquire(unit, size=file_size)
            data_len = int(file_size, c_size_t)
            
            if (file_size == 0) then
                close(unit)
                ! Add empty file
                call add_empty_file_to_zip(zip_handle, filename, ierr)
                return
            end if
            
            ! Allocate memory and read file
            c_data = malloc(data_len)
            if (.not. c_associated(c_data)) then
                call set_err_once(ierr, ERR_POINTER_NULL)
                close(unit)
                return
            end if
            
            call c_f_pointer(c_data, file_data, [file_size])
            read(unit, iostat=iostat) file_data
            close(unit)
            
            if (.not. is_ok(iostat)) then
                call set_err_once(ierr, ERR_READ_DATA)
                call free(c_data)
                return
            end if
            
            source = zip_source_buffer(zip_handle, c_data, data_len, 1)
            
        case (2) ! String data
            
            data_len = len(data_source)
            c_data = malloc(data_len)
            if (.not. c_associated(c_data)) then
                call set_err_once(ierr, ERR_POINTER_NULL)
                return
            end if
            
            call c_f_pointer(c_data, string_data, [data_len])
            do i = 1, data_len
                string_data(i) = data_source(i:i)
            end do
            
            source = zip_source_buffer(zip_handle, c_data, data_len, 1)
            
        case default
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end select
        
        if (.not. c_associated(source)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
            if (data_type == 1 .or. data_type == 2) call free(c_data)
            return
        end if
        
        ! Add file to ZIP
        index = zip_file_add(zip_handle, trim(filename)//c_null_char, source, ZIP_FL_OVERWRITE)
        if (index < 0) then
            call set_err_once(ierr, ERR_FILE_ADD)
            call zip_source_free(source)
            return
        end if
        
        ! Set compression to store (no compression)
        error = zip_set_file_compression(zip_handle, index, ZIP_CM_STORE, 0)
        if (.not. is_ok(error)) then
            print *, "Warning: Error setting compression for: ", trim(filename)
        end if
        
        print *, "Added to ZIP: ", trim(filename)
        
    contains
        subroutine add_empty_file_to_zip(zip_handle, filename, ierr)
            type(c_ptr), intent(in) :: zip_handle
                !! Zip connection
            character(len=*), intent(in) :: filename
                !! Name of the file to add
            integer(int32), intent(out) :: ierr
                !! Error code
            
            type(c_ptr) :: source
            integer(c_int64_t) :: index
            integer(c_int) :: error
            
            call set_ok(ierr)
            
            source = zip_source_buffer(zip_handle, c_null_ptr, 0_c_size_t, 0)
            if (.not. c_associated(source)) then
                call set_err_once(ierr, ERR_POINTER_NULL)
                return
            end if
            
            index = zip_file_add(zip_handle, trim(filename)//c_null_char, source, ZIP_FL_OVERWRITE)
            if (index < 0) then
                call set_err_once(ierr, ERR_FILE_ADD)
                call zip_source_free(source)
                return
            end if
            
            print *, "Added empty file to ZIP: ", trim(filename)
        end subroutine add_empty_file_to_zip
    end subroutine add_data_to_zip

    ! >Manifest creation
    subroutine write_manifest(gene_ids_file, expression_file, gene_to_family_file, &
                            family_ids_file, family_centroids_file, shift_vectors_file, &
                            manifest_filename, ierr)
        character(len=*), intent(in) :: gene_ids_file
            !! Name of the gene ids file
        character(len=*), intent(in) :: expression_file
            !! Name of the expression file
        character(len=*), intent(in) :: gene_to_family_file
            !! Name of the gene to family mapping file
        character(len=*), intent(in) :: family_ids_file 
            !! Name of the family ids file
        character(len=*), intent(in) :: family_centroids_file
            !! Name of the family centroids file
        character(len=*), intent(in) :: shift_vectors_file
            !! Name of the shift vectors file
        character(len=*), intent(in) :: manifest_filename
            !! Name of the manifest (should be manifest.txt)
        integer(int32), intent(out) :: ierr
            !! Error code
        
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

    !> Read an extracted manifest
    subroutine read_manifest(manifest_filename, gene_ids_file, expression_file, gene_to_family_file, &
                            family_ids_file, family_centroids_file, shift_vectors_file, ierr)
        character(len=*), intent(in) :: manifest_filename
            !! Filename of the manifest (should be manifest.txt)
        character(len=:), allocatable, intent(out) :: gene_ids_file
            !! gene ids filename
        character(len=:), allocatable, intent(out) :: expression_file
            !! expression vectors filename
        character(len=:), allocatable, intent(out) :: gene_to_family_file
            !! gene to family filename
        character(len=:), allocatable, intent(out) :: family_ids_file
            !! family ids filename
        character(len=:), allocatable, intent(out) :: family_centroids_file
            !! family centroids filename
        character(len=:), allocatable, intent(out) :: shift_vectors_file
            !! shift vectors filename
        integer(int32), intent(out) :: ierr
            !! Error code
        
        integer(int32) :: unit, iostat
        character(len=512) :: line
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

    ! Module-level subroutine for manifest extraction and parsing
    subroutine extract_and_parse_manifest(zip_handle, gene_ids_file, expression_file, gene_to_family_file, &
                                        family_ids_file, family_centroids_file, shift_vectors_file, &
                                        ierr)
        type(c_ptr), intent(in) :: zip_handle
            !! Zip file connection
        character(len=:), allocatable, intent(out) :: gene_ids_file
            !! Gene ids filename
        character(len=:), allocatable, intent(out) :: expression_file
            !! expression vectors filename
        character(len=:), allocatable, intent(out) :: gene_to_family_file
            !! gene to family mapping filename
        character(len=:), allocatable, intent(out) :: family_ids_file
            !! family ids filename
        character(len=:), allocatable, intent(out) :: family_centroids_file
            !! family centroids filename
        character(len=:), allocatable, intent(out) :: shift_vectors_file
            !! Shift vectors filename
        integer(int32), intent(out) :: ierr
            !! Error code
        
        type(c_ptr) :: file_handle
        integer(c_int) :: error
        integer(c_int64_t) :: bytes_read
        integer(int32) :: unit, iostat
        integer(int32), parameter :: CHUNK_SIZE = 4096
        character(kind=c_char), dimension(:), allocatable, target :: buffer
        
        ! Extract and parse the manifest file
        file_handle = zip_fopen(zip_handle, "manifest.txt"//c_null_char, 0)
        if (c_associated(file_handle)) then
            ! Open output file for manifest
            open(newunit=unit, file="manifest.txt", access='stream', form='unformatted', &
                iostat=iostat, status='replace', action='write')
            if (.not. is_ok(iostat)) then
                print *, "Error creating manifest file"
                error = zip_fclose(file_handle)
                call set_err_once(ierr, ERR_FILE_OPEN)
                return
            end if
            
            ! Read and write in chunks
            allocate(buffer(CHUNK_SIZE), stat=iostat)
            if (is_ok(iostat)) then
                do
                    bytes_read = zip_fread(file_handle, c_loc(buffer), int(CHUNK_SIZE, c_size_t))
                    if (bytes_read <= 0) exit
                    write(unit, iostat=iostat) buffer(1:bytes_read)
                    if (.not. is_ok(iostat)) then
                        print *, "Error writing manifest file"
                        call set_err(ierr, ERR_WRITE_DATA)
                        exit
                    end if
                end do
            else
                call set_err(ierr, ERR_ALLOC_FAIL)
                error = zip_close(zip_handle)
                return
            end if
            
            ! Clean up manifest extraction
            if (allocated(buffer)) deallocate(buffer)
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
            if (.not. is_ok(ierr)) then
                print *, "Warning: Could not delete temporary manifest file"
                call set_ok(ierr)  ! Not critical
            end if
        else
            print *, "No manifest file found in ZIP archive"
            call set_err_once(ierr, ERR_MISSING_MANIFEST)
        end if
    end subroutine extract_and_parse_manifest

    subroutine save_tox_data(zip_filename, ierr, gene_ids, gene_ids_file, expression, &
                            expression_file, gene_to_family, gene_to_family_file, &
                            family_ids, family_ids_file, family_centroids, &
                            family_centroids_file, shift_vectors, shift_vectors_file)
        implicit none
        
        character(len=*), intent(in) :: zip_filename
            !! Name of the zipfile
        character(len=*), intent(in), optional :: gene_ids(:)
            !! Gene ids array
        character(len=*), intent(in), optional :: family_ids(:)
            !! family ids array
        real(real64), intent(in), optional :: expression(:,:)
            !! Expression vectors array
        real(real64), intent(in), optional :: family_centroids(:,:)
            !! family centroids array
        real(real64), intent(in), optional :: shift_vectors(:,:)
            !! Shift vectors array
        integer(int32), intent(in), optional :: gene_to_family(:)
            !! gene to family mapping array
        character(len=*), intent(in), optional :: gene_ids_file
            !! gene ids filename
        character(len=*), intent(in), optional :: expression_file
            !! expression vectors filename
        character(len=*), intent(in), optional :: gene_to_family_file
            !! Gene to family mapping filename
        character(len=*), intent(in), optional :: family_ids_file
            !! Family ids filename
        character(len=*), intent(in), optional :: family_centroids_file
            !! Family centroids filename
        character(len=*), intent(in), optional :: shift_vectors_file
            !! shift vectors filename
        integer(int32), intent(out) :: ierr
            !! Error code
        
        character(len=:), allocatable :: actual_gene_ids_file, actual_expression_file, actual_gene_to_family_file, &
                                        actual_family_ids_file, actual_family_centroids_file, actual_shift_vectors_file
        logical :: gene_ids_present, expression_present, gene_to_family_present, &
                family_ids_present, family_centroids_present, shift_vectors_present
        integer(int32) :: temp_ierr
        
        call set_ok(ierr)
        call set_ok(temp_ierr)
        
        ! Determine which arrays are present
        gene_ids_present = present(gene_ids) .and. present(gene_ids_file)
        expression_present = present(expression) .and. present(expression_file)
        gene_to_family_present = present(gene_to_family) .and. present(gene_to_family_file)
        family_ids_present = present(family_ids) .and. present(family_ids_file)
        family_centroids_present = present(family_centroids) .and. present(family_centroids_file)
        shift_vectors_present = present(shift_vectors) .and. present(shift_vectors_file)
        
        ! Save data files
        if (is_ok(ierr) .and. gene_ids_present) then
            actual_gene_ids_file = gene_ids_file
            call save_gene_ids(gene_ids, actual_gene_ids_file, ierr)
        else
            actual_gene_ids_file = ""
        end if
        
        if (is_ok(ierr) .and. expression_present) then
            actual_expression_file = expression_file
            call save_expression_vectors(expression, actual_expression_file, ierr)
        else
            actual_expression_file = ""
        end if
        
        if (is_ok(ierr) .and. gene_to_family_present) then
            actual_gene_to_family_file = gene_to_family_file
            call save_gene_to_family(gene_to_family, actual_gene_to_family_file, ierr)
        else
            actual_gene_to_family_file = ""
        end if
        
        if (is_ok(ierr) .and. family_ids_present) then
            actual_family_ids_file = family_ids_file
            call save_family_ids(family_ids, actual_family_ids_file, ierr)
        else
            actual_family_ids_file = ""
        end if
        
        if (is_ok(ierr) .and. family_centroids_present) then
            actual_family_centroids_file = family_centroids_file
            call save_family_centroids(family_centroids, actual_family_centroids_file, ierr)
        else
            actual_family_centroids_file = ""
        end if
        
        if (is_ok(ierr) .and. shift_vectors_present) then
            actual_shift_vectors_file = shift_vectors_file
            call save_shift_vectors(shift_vectors, actual_shift_vectors_file, ierr)
        else
            actual_shift_vectors_file = ""
        end if
        
        ! Create the ZIP archive if all saves were successful
        if (is_ok(ierr)) then
            call create_zip_archive(zip_filename, actual_gene_ids_file, actual_expression_file, &
                                actual_gene_to_family_file, actual_family_ids_file, &
                                actual_family_centroids_file, actual_shift_vectors_file, ierr)
        end if
        
        ! Clean up temporary files, preserving the original error
        call cleanup_temporary_files(gene_ids_present, actual_gene_ids_file, "Gene IDs")
        call cleanup_temporary_files(expression_present, actual_expression_file, "Expression")
        call cleanup_temporary_files(gene_to_family_present, actual_gene_to_family_file, "Gene to family mapping")
        call cleanup_temporary_files(family_ids_present, actual_family_ids_file, "Family IDs")
        call cleanup_temporary_files(family_centroids_present, actual_family_centroids_file, "Centroids")
        call cleanup_temporary_files(shift_vectors_present, actual_shift_vectors_file, "Shift vectors")
    contains
        subroutine cleanup_temporary_files(file_present, filename, description)
            logical, intent(in) :: file_present
            character(len=*), intent(in) :: filename
            character(len=*), intent(in) :: description
            
            if (file_present .and. len_trim(filename) > 0) then
                call delete_file(filename, temp_ierr)
                if(is_err(temp_ierr)) then
                    write(*,*) 'Warning: ', trim(description), ' file could not be removed: ', trim(filename)
                end if
            end if
        end subroutine cleanup_temporary_files
    end subroutine save_tox_data

    subroutine read_tox_data(zip_filename, ierr, gene_ids, gene_ids_file, expression, expression_file, &
                        gene_to_family, gene_to_family_file, family_ids, family_ids_file, &
                        family_centroids, family_centroids_file, shift_vectors, shift_vectors_file)
        use array_utils, only: get_array_metadata
        use tox_data_read_write
        use iso_fortran_env, only: real64, int32
        implicit none
        
        character(len=*), intent(in) :: zip_filename
            !! Zip file to read from
        integer(int32), intent(out) :: ierr
            !! Error code
        character(len=:), allocatable, optional, intent(out) :: gene_ids(:)
            !! Gene ids array
        character(len=:), allocatable, optional, intent(out) :: family_ids(:)
            !! Family ids array
        real(real64), allocatable, optional, intent(out) :: expression(:,:)
            !! expression vectors array
        real(real64), allocatable, optional, intent(out) :: family_centroids(:,:)
            !! family centroids array
        real(real64), allocatable, optional, intent(out) :: shift_vectors(:,:)
            !! shift vectors array
        integer(int32), allocatable, optional, intent(out) :: gene_to_family(:) 
            !! gene to family mapping array
        character(len=:), allocatable, optional, intent(out) :: gene_ids_file
            !! gene ids filename
        character(len=:), allocatable, optional, intent(out) :: expression_file
            !! expression vectors filename
        character(len=:), allocatable, optional, intent(out) :: gene_to_family_file
            !! gene to family mapping filename
        character(len=:), allocatable, optional, intent(out) :: family_ids_file
            !! family ids filename
        character(len=:), allocatable, optional, intent(out) :: family_centroids_file   
            !! family centroids filename
        character(len=:), allocatable, optional, intent(out) :: shift_vectors_file
            !! shift vectors filename
        
        character(len=:), allocatable :: extracted_gene_ids_file
        character(len=:), allocatable :: extracted_expression_file
        character(len=:), allocatable :: extracted_gene_to_family_file
        character(len=:), allocatable :: extracted_family_ids_file 
        character(len=:), allocatable :: extracted_family_centroids_file
        character(len=:), allocatable :: extracted_shift_vectors_file

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
        end if

        if(len_trim(extracted_gene_ids_file) > 0) call delete_file(extracted_gene_ids_file, ierr)
        if(len_trim(extracted_expression_file) > 0) call delete_file(extracted_expression_file, ierr)
        if(len_trim(extracted_gene_to_family_file) > 0) call delete_file(extracted_gene_to_family_file, ierr)
        if(len_trim(extracted_family_ids_file) > 0) call delete_file(extracted_family_ids_file, ierr)
        if(len_trim(extracted_family_centroids_file) > 0) call delete_file(extracted_family_centroids_file, ierr)
        if(len_trim(extracted_shift_vectors_file) > 0) call delete_file(extracted_shift_vectors_file, ierr)
    end subroutine read_tox_data

end module tox_archive