module tox_archive
    use iso_c_binding
    use tox_data_read_write
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
        character(len=:), allocatable, target :: manifest_content
        
        ! Open ZIP archive for writing
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_CREATE, error)
        if (error /= 0) then
            ierr = error
            print *, "Error opening ZIP file for writing: ", error
            return
        end if
        
        ! Add files if they are non-empty
        if (len_trim(gene_ids_file) > 0) then
            call add_file_to_zip(zip_handle, trim(gene_ids_file))
        end if
        
        if (len_trim(expression_file) > 0) then
            call add_file_to_zip(zip_handle, trim(expression_file))
        end if
        
        if (len_trim(gene_to_family_file) > 0) then
            call add_file_to_zip(zip_handle, trim(gene_to_family_file))
        end if
        
        if (len_trim(family_ids_file) > 0) then
            call add_file_to_zip(zip_handle, trim(family_ids_file))
        end if
        
        if (len_trim(family_centroids_file) > 0) then
            call add_file_to_zip(zip_handle, trim(family_centroids_file))
        end if
        
        if (len_trim(shift_vectors_file) > 0) then
            call add_file_to_zip(zip_handle, trim(shift_vectors_file))
        end if
        
        ! Create and add manifest
        call write_manifest(gene_ids_file, expression_file, gene_to_family_file, &
                        family_ids_file, family_centroids_file, shift_vectors_file, &
                        manifest_content)
        call add_string_to_zip(zip_handle, "manifest.txt", manifest_content)
        
        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (error /= 0) then
            ierr = error
            print *, "Error closing ZIP file: ", error
            return
        end if
        
        ierr = 0
        print *, "ZIP archive created successfully: ", trim(zip_filename)
    end subroutine create_zip_archive

    subroutine extract_zip_archive(zip_filename, gene_ids_file, expression_file, gene_to_family_file, &
                                family_ids_file, family_centroids_file, shift_vectors_file, ierr)
        character(len=*), intent(in) :: zip_filename
        character(len=:), allocatable, intent(out) :: gene_ids_file, expression_file, gene_to_family_file, &
                                                    family_ids_file, family_centroids_file, shift_vectors_file
        integer, intent(out) :: ierr
        
        type(c_ptr) :: zip_handle, file_handle
        integer(c_int) :: error
        integer(c_int64_t) :: i, num_entries, bytes_read
        type(zip_stat_t), target :: sb
        character(len=:), allocatable :: filename, manifest_content
        character(kind=c_char), dimension(:), allocatable, target :: buffer
        
        ! Initialize output variables
        gene_ids_file = ""
        expression_file = ""
        gene_to_family_file = ""
        family_ids_file = ""
        family_centroids_file = ""
        shift_vectors_file = ""
        ierr = 0
        
        ! Open ZIP archive for reading
        zip_handle = zip_open(trim(zip_filename)//c_null_char, ZIP_RDONLY, error)
        if (error /= 0) then
            ierr = error
            print *, "Error opening ZIP file for reading: ", error
            return
        end if
        
        ! Get number of entries in the ZIP
        num_entries = zip_get_num_entries(zip_handle, 0)
        
        ! Extract each file
        do i = 0, num_entries - 1
            ! Get filename
            call get_zip_entry_name(zip_handle, i, filename)
            
            ! Skip if it's the manifest (we'll handle it separately)
            if (filename == "manifest.txt") cycle
            
            ! Get file info
            error = zip_stat(zip_handle, trim(filename)//c_null_char, 0, c_loc(sb))
            if (error /= 0) then
                print *, "Error getting file info for: ", trim(filename)
                cycle
            end if
            
            ! Allocate buffer for file content
            allocate(buffer(sb%size))
            
            ! Open file in ZIP
            file_handle = zip_fopen(zip_handle, trim(filename)//c_null_char, 0)
            if (.not. c_associated(file_handle)) then
                print *, "Error opening file in ZIP: ", trim(filename)
                deallocate(buffer)
                cycle
            end if
            
            ! Read file content
            bytes_read = zip_fread(file_handle, c_loc(buffer), int(sb%size, c_size_t))
            if (bytes_read /= sb%size) then
                print *, "Error reading file from ZIP: ", trim(filename)
                deallocate(buffer)
                error = zip_fclose(file_handle)
                cycle
            end if
            
            ! Write to disk
            call write_file(trim(filename), buffer)
            
            ! Clean up
            deallocate(buffer)
            error = zip_fclose(file_handle)
            if (error /= 0) then
                print *, "Error closing file in ZIP: ", trim(filename)
            end if
            
            print *, "Extracted: ", trim(filename)
        end do
        
        ! Read and parse the manifest file
        error = zip_stat(zip_handle, "manifest.txt"//c_null_char, 0, c_loc(sb))
        if (error == 0) then
            ! Allocate buffer for manifest content
            allocate(buffer(sb%size))
            
            ! Open manifest file in ZIP
            file_handle = zip_fopen(zip_handle, "manifest.txt"//c_null_char, 0)
            if (c_associated(file_handle)) then
                ! Read manifest content
                bytes_read = zip_fread(file_handle, c_loc(buffer), int(sb%size, c_size_t))
                if (bytes_read == sb%size) then
                    ! Convert to string
                    allocate(character(len=sb%size) :: manifest_content)
                    manifest_content = transfer(buffer, manifest_content)
                    
                    ! Parse manifest to get the filenames
                    call read_manifest(manifest_content, gene_ids_file, expression_file, gene_to_family_file, &
                                    family_ids_file, family_centroids_file, shift_vectors_file, ierr)
                    
                    if (ierr /= 0) then
                        print *, "Error parsing manifest file"
                    end if
                    
                    deallocate(manifest_content)
                end if
                
                error = zip_fclose(file_handle)
            end if
            
            deallocate(buffer)
            
            ! Extract the manifest file to disk
            error = zip_stat(zip_handle, "manifest.txt"//c_null_char, 0, c_loc(sb))
            if (error == 0) then
                allocate(buffer(sb%size))
                file_handle = zip_fopen(zip_handle, "manifest.txt"//c_null_char, 0)
                if (c_associated(file_handle)) then
                    bytes_read = zip_fread(file_handle, c_loc(buffer), int(sb%size, c_size_t))
                    if (bytes_read == sb%size) then
                        call write_file("manifest.txt", buffer)
                    end if
                    error = zip_fclose(file_handle)
                end if
                deallocate(buffer)
            end if
        else
            print *, "No manifest file found in ZIP archive"
            ierr = -2
        end if
        
        ! Close ZIP archive
        error = zip_close(zip_handle)
        if (error /= 0) then
            print *, "Error closing ZIP file: ", error
            ierr = error
            return
        end if
        
        print *, "ZIP archive extracted successfully: ", trim(zip_filename)
    end subroutine extract_zip_archive

    subroutine write_file(filename, content)
        character(len=*), intent(in) :: filename
        character(kind=c_char), dimension(:), intent(in) :: content
        
        integer :: unit, iostat
        
        open(unit, file=filename, access='stream', form='unformatted', iostat=iostat, status='replace')
        if (iostat /= 0) then
            print *, "Error creating file: ", trim(filename)
            return
        end if
        
        write(unit, iostat=iostat) content
        close(unit)
        
        if (iostat /= 0) then
            print *, "Error writing file: ", trim(filename)
        end if
    end subroutine write_file
    
    ! Helper function to add a file to ZIP
    subroutine add_file_to_zip(zip_handle, filename)
        type(c_ptr), intent(in) :: zip_handle
        character(len=*), intent(in) :: filename
        
        integer(c_int) :: error
        integer(c_int64_t) :: index
        type(c_ptr) :: source
        integer(c_signed_char), allocatable, target :: data(:)
        integer :: unit, iostat, file_size
        
        ! Read file content
        open(unit, file=filename, access='stream', form='unformatted', iostat=iostat, status='old')
        if (iostat /= 0) then
            print *, "Error opening file: ", trim(filename)
            return
        end if
        
        inquire(unit, size=file_size)
        allocate(data(file_size))
        read(unit, iostat=iostat) data
        close(unit)
        
        if (iostat /= 0) then
            print *, "Error reading file: ", trim(filename)
            deallocate(data)
            return
        end if
        
        ! Create ZIP source from file content
        source = zip_source_buffer(zip_handle, c_loc(data), int(file_size, c_size_t), 1)
        if (.not. c_associated(source)) then
            print *, "Error creating source for: ", trim(filename)
            deallocate(data)
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
        
        deallocate(data)
        print *, "Added to ZIP: ", trim(filename)
    end subroutine add_file_to_zip

    ! Helper function to add a string as a file to ZIP
    subroutine add_string_to_zip(zip_handle, filename, content)
        type(c_ptr), intent(in) :: zip_handle
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), target :: content
        
        integer(c_int) :: error
        integer(c_int64_t) :: index
        type(c_ptr) :: source
        
        ! Create ZIP source from string content
        source = zip_source_buffer(zip_handle, c_loc(content(1:1)), int(len(content), c_size_t), 0)
        if (.not. c_associated(source)) then
            print *, "Error creating source for: ", trim(filename)
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
    subroutine get_zip_entry_name(zip_handle, entry_index, name)
        type(c_ptr), intent(in) :: zip_handle
        integer(c_int64_t), intent(in) :: entry_index
        character(len=:), allocatable, intent(out) :: name
        
        type(c_ptr) :: name_ptr
        integer :: name_len, i
        character(kind=c_char), pointer :: f_ptr(:)
        
        ! Get name from ZIP
        name_ptr = zip_get_name(zip_handle, entry_index, 0)
        if (.not. c_associated(name_ptr)) then
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
            allocate(character(len=name_len) :: name)
            name = transfer(f_ptr(1:name_len), name)
        else
            name = ""
        end if
    end subroutine get_zip_entry_name

    ! Manifest creation and parsing functions
    subroutine write_manifest(gene_ids_file, expression_file, gene_to_family_file, &
                            family_ids_file, family_centroids_file, shift_vectors_file, &
                            manifest_content)
        character(len=*), intent(in) :: gene_ids_file, expression_file, gene_to_family_file, &
                                        family_ids_file, family_centroids_file, shift_vectors_file
        character(len=:), allocatable, target, intent(out) :: manifest_content
        integer :: total_length, line_count
        
        ! Calculate total length needed
        total_length = 0
        line_count = 0
        
        if (len_trim(gene_ids_file) > 0) then
            total_length = total_length + len('gene_ids=""') + len_trim(gene_ids_file)
            line_count = line_count + 1
        end if
        
        if (len_trim(expression_file) > 0) then
            total_length = total_length + len('expression=""') + len_trim(expression_file)
            line_count = line_count + 1
        end if
        
        if (len_trim(gene_to_family_file) > 0) then
            total_length = total_length + len('gene_to_family=""') + len_trim(gene_to_family_file)
            line_count = line_count + 1
        end if
        
        if (len_trim(family_ids_file) > 0) then
            total_length = total_length + len('family_ids=""') + len_trim(family_ids_file)
            line_count = line_count + 1
        end if
        
        if (len_trim(family_centroids_file) > 0) then
            total_length = total_length + len('family_centroids=""') + len_trim(family_centroids_file)
            line_count = line_count + 1
        end if
        
        if (len_trim(shift_vectors_file) > 0) then
            total_length = total_length + len('shift_vectors=""') + len_trim(shift_vectors_file)
            line_count = line_count + 1
        end if
        
        ! Add space for newline characters
        total_length = total_length + line_count
        
        ! Allocate and build the manifest content
        allocate(character(total_length) :: manifest_content)
        manifest_content = ""
        
        if (len_trim(gene_ids_file) > 0) then
            manifest_content = 'gene_ids="' // trim(gene_ids_file) // '"'
        end if
        
        if (len_trim(expression_file) > 0) then
            if (len_trim(manifest_content) > 0) then
                manifest_content = manifest_content // new_line('a') // 'expression="' // trim(expression_file) // '"'
            else
                manifest_content = 'expression="' // trim(expression_file) // '"'
            end if
        end if
        
        if (len_trim(gene_to_family_file) > 0) then
            if (len_trim(manifest_content) > 0) then
                manifest_content = manifest_content // new_line('a') // 'gene_to_family="' // trim(gene_to_family_file) // '"'
            else
                manifest_content = 'gene_to_family="' // trim(gene_to_family_file) // '"'
            end if
        end if
        
        if (len_trim(family_ids_file) > 0) then
            if (len_trim(manifest_content) > 0) then
                manifest_content = manifest_content // new_line('a') // 'family_ids="' // trim(family_ids_file) // '"'
            else
                manifest_content = 'family_ids="' // trim(family_ids_file) // '"'
            end if
        end if
        
        if (len_trim(family_centroids_file) > 0) then
            if (len_trim(manifest_content) > 0) then
                manifest_content = manifest_content // new_line('a') // 'family_centroids="' // trim(family_centroids_file) // '"'
            else
                manifest_content = 'family_centroids="' // trim(family_centroids_file) // '"'
            end if
        end if
        
        if (len_trim(shift_vectors_file) > 0) then
            if (len_trim(manifest_content) > 0) then
                manifest_content = manifest_content // new_line('a') // 'shift_vectors="' // trim(shift_vectors_file) // '"'
            else
                manifest_content = 'shift_vectors="' // trim(shift_vectors_file) // '"'
            end if
        end if
    end subroutine write_manifest

    subroutine read_manifest(manifest_content, gene_ids_file, expression_file, gene_to_family_file, &
                            family_ids_file, family_centroids_file, shift_vectors_file, status)
        character(len=*), intent(in) :: manifest_content
        character(len=:), allocatable, intent(out) :: gene_ids_file, expression_file, gene_to_family_file, &
                                                    family_ids_file, family_centroids_file, shift_vectors_file
        integer, intent(out) :: status
        
        integer :: gene_start, gene_end, expr_start, expr_end, gtf_start, gtf_end
        integer :: fid_start, fid_end, fc_start, fc_end, sv_start, sv_end
        integer :: line_break_pos
        
        status = 0
        
        ! Initialize all outputs to empty strings
        gene_ids_file = ""
        expression_file = ""
        gene_to_family_file = ""
        family_ids_file = ""
        family_centroids_file = ""
        shift_vectors_file = ""
        
        ! Parse gene_ids
        gene_start = index(manifest_content, 'gene_ids="')
        if (gene_start > 0) then
            gene_start = gene_start + len('gene_ids="')
            gene_end = index(manifest_content(gene_start:), '"')
            if (gene_end > 0) then
                gene_end = gene_start + gene_end - 2
                gene_ids_file = manifest_content(gene_start:gene_end)
            end if
        end if
        
        ! Parse expression
        expr_start = index(manifest_content, 'expression="')
        if (expr_start > 0) then
            expr_start = expr_start + len('expression="')
            expr_end = index(manifest_content(expr_start:), '"')
            if (expr_end > 0) then
                expr_end = expr_start + expr_end - 2
                expression_file = manifest_content(expr_start:expr_end)
            end if
        end if
        
        ! Parse gene_to_family
        gtf_start = index(manifest_content, 'gene_to_family="')
        if (gtf_start > 0) then
            gtf_start = gtf_start + len('gene_to_family="')
            gtf_end = index(manifest_content(gtf_start:), '"')
            if (gtf_end > 0) then
                gtf_end = gtf_start + gtf_end - 2
                gene_to_family_file = manifest_content(gtf_start:gtf_end)
            end if
        end if
        
        ! Parse family_ids
        fid_start = index(manifest_content, 'family_ids="')
        if (fid_start > 0) then
            fid_start = fid_start + len('family_ids="')
            fid_end = index(manifest_content(fid_start:), '"')
            if (fid_end > 0) then
                fid_end = fid_start + fid_end - 2
                family_ids_file = manifest_content(fid_start:fid_end)
            end if
        end if
        
        ! Parse family_centroids
        fc_start = index(manifest_content, 'family_centroids="')
        if (fc_start > 0) then
            fc_start = fc_start + len('family_centroids="')
            fc_end = index(manifest_content(fc_start:), '"')
            if (fc_end > 0) then
                fc_end = fc_start + fc_end - 2
                family_centroids_file = manifest_content(fc_start:fc_end)
            end if
        end if
        
        ! Parse shift_vectors
        sv_start = index(manifest_content, 'shift_vectors="')
        if (sv_start > 0) then
            sv_start = sv_start + len('shift_vectors="')
            sv_end = index(manifest_content(sv_start:), '"')
            if (sv_end > 0) then
                sv_end = sv_start + sv_end - 2
                shift_vectors_file = manifest_content(sv_start:sv_end)
            end if
        end if
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
        integer, intent(out) :: ierr
        
        character(len=:), allocatable :: actual_gene_ids_file, actual_expression_file, actual_gene_to_family_file, &
                                        actual_family_ids_file, actual_family_centroids_file, actual_shift_vectors_file
        logical :: gene_ids_present, expression_present, gene_to_family_present, &
                family_ids_present, family_centroids_present, shift_vectors_present
        
        ierr = 0
        
        ! Determine which arrays are present and set appropriate filenames
        gene_ids_present = present(gene_ids) .and. present(gene_ids_file)
        if (gene_ids_present) then
            actual_gene_ids_file = gene_ids_file
            call save_gene_ids(gene_ids, actual_gene_ids_file, ierr)
        else
            actual_gene_ids_file = ""
        end if
        
        expression_present = present(expression) .and. present(expression_file)
        if (expression_present) then
            actual_expression_file = expression_file
            call save_expression_vectors(expression, actual_expression_file, ierr)
        else
            actual_expression_file = ""
        end if
        
        gene_to_family_present = present(gene_to_family) .and. present(gene_to_family_file)
        if (gene_to_family_present) then
            actual_gene_to_family_file = gene_to_family_file
            call save_gene_to_family(gene_to_family, actual_gene_to_family_file, ierr)
        else
            actual_gene_to_family_file = ""
        end if
        
        family_ids_present = present(family_ids) .and. present(family_ids_file)
        if (family_ids_present) then
            actual_family_ids_file = family_ids_file
            call save_family_ids(family_ids, actual_family_ids_file, ierr)
        else
            actual_family_ids_file = ""
        end if
        
        family_centroids_present = present(family_centroids) .and. present(family_centroids_file)
        if (family_centroids_present) then
            actual_family_centroids_file = family_centroids_file
            call save_family_centroids(family_centroids, actual_family_centroids_file, ierr)
        else
            actual_family_centroids_file = ""
        end if
        
        shift_vectors_present = present(shift_vectors) .and. present(shift_vectors_file)
        if (shift_vectors_present) then
            actual_shift_vectors_file = shift_vectors_file
            call save_shift_vectors(shift_vectors, actual_shift_vectors_file, ierr)
        else
            actual_shift_vectors_file = ""
        end if
        
        ! Create the ZIP archive
        call create_zip_archive(zip_filename, actual_gene_ids_file, actual_expression_file, &
                            actual_gene_to_family_file, actual_family_ids_file, &
                            actual_family_centroids_file, actual_shift_vectors_file, ierr)
        
        ! Clean up temporary files
        if (gene_ids_present) call delete_file(actual_gene_ids_file)
        if (expression_present) call delete_file(actual_expression_file)
        if (gene_to_family_present) call delete_file(actual_gene_to_family_file)
        if (family_ids_present) call delete_file(actual_family_ids_file)
        if (family_centroids_present) call delete_file(actual_family_centroids_file)
        if (shift_vectors_present) call delete_file(actual_shift_vectors_file)
        
    contains
        subroutine delete_file(filename)
            character(len=*), intent(in) :: filename
            integer(int32) :: unit
            logical :: exists
            inquire(file=filename, exist=exists)
            if (exists) then
                open(newunit=unit, file=filename, status='old')
                close(unit, status='delete')
            end if
        end subroutine delete_file
    end subroutine save_tox_data

    subroutine read_tox_data(zip_filename, ierr, gene_ids, gene_ids_file, expression, expression_file, &
                        gene_to_family, gene_to_family_file, family_ids, family_ids_file, &
                        family_centroids, family_centroids_file, shift_vectors, shift_vectors_file)
        use array_utils, only: get_array_metadata
        use tox_data_read_write
        use iso_fortran_env, only: real64, int32
        implicit none
        
        character(len=*), intent(in) :: zip_filename
        integer, intent(out) :: ierr
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
        integer(int32) :: max_dims, ndims, dims(5), metadata_ierr, char_len
        
        ierr = 0
        max_dims = 5  ! Maximum number of dimensions supported
        
        ! Extract the ZIP archive and get file names from manifest
        call extract_zip_archive(zip_filename, extracted_gene_ids_file, extracted_expression_file, &
                                extracted_gene_to_family_file, extracted_family_ids_file, &
                                extracted_family_centroids_file, extracted_shift_vectors_file, ierr)
        if (ierr /= 0) return
        
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
            call get_array_metadata(extracted_gene_ids_file, dims, max_dims, ndims, metadata_ierr, char_len)
            if (metadata_ierr == 0 .and. ndims == 1) then
                ! Allocate array based on metadata with proper character length
                allocate(character(len=char_len) :: gene_ids(dims(1)))
                call load_gene_ids(gene_ids, extracted_gene_ids_file, ierr)
            else
                ierr = metadata_ierr
                print *, "Error getting metadata for gene_ids file"
            end if
            call delete_file(extracted_gene_ids_file)
        end if
        
        expression_requested = present(expression) .and. len_trim(extracted_expression_file) > 0
        if (expression_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_expression_file, dims, max_dims, ndims, metadata_ierr)
            if (metadata_ierr == 0 .and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(expression(dims(1), dims(2)))
                call load_expression_vectors(expression, extracted_expression_file, ierr)
            else
                ierr = metadata_ierr
                print *, "Error getting metadata for expression file"
            end if
            call delete_file(extracted_expression_file)
        end if
        
        gene_to_family_requested = present(gene_to_family) .and. len_trim(extracted_gene_to_family_file) > 0
        if (gene_to_family_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_gene_to_family_file, dims, max_dims, ndims, metadata_ierr)
            if (metadata_ierr == 0 .and. ndims == 1) then
                ! Allocate array based on metadata
                allocate(gene_to_family(dims(1)))
                call load_gene_to_family(gene_to_family, extracted_gene_to_family_file, ierr)
            else
                ierr = metadata_ierr
                print *, "Error getting metadata for gene_to_family file"
            end if
            call delete_file(extracted_gene_to_family_file)
        end if
        
        family_ids_requested = present(family_ids) .and. len_trim(extracted_family_ids_file) > 0
        if (family_ids_requested) then
            ! Get array metadata to determine size and character length
            call get_array_metadata(extracted_family_ids_file, dims, max_dims, ndims, metadata_ierr, char_len)
            if (metadata_ierr == 0 .and. ndims == 1) then
                ! Allocate array based on metadata with proper character length
                allocate(character(len=char_len) :: family_ids(dims(1)))
                call load_family_ids(family_ids, extracted_family_ids_file, ierr)
            else
                ierr = metadata_ierr
                print *, "Error getting metadata for family_ids file"
            end if
            call delete_file(extracted_family_ids_file)
        end if
        
        family_centroids_requested = present(family_centroids) .and. len_trim(extracted_family_centroids_file) > 0
        if (family_centroids_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_family_centroids_file, dims, max_dims, ndims, metadata_ierr)
            if (metadata_ierr == 0 .and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(family_centroids(dims(1), dims(2)))
                call load_family_centroids(family_centroids, extracted_family_centroids_file, ierr)
            else
                ierr = metadata_ierr
                print *, "Error getting metadata for family_centroids file"
            end if
            call delete_file(extracted_family_centroids_file)
        end if
        
        shift_vectors_requested = present(shift_vectors) .and. len_trim(extracted_shift_vectors_file) > 0
        if (shift_vectors_requested) then
            ! Get array metadata to determine size
            call get_array_metadata(extracted_shift_vectors_file, dims, max_dims, ndims, metadata_ierr)
            if (metadata_ierr == 0 .and. ndims == 2) then
                ! Allocate array based on metadata
                allocate(shift_vectors(dims(1), dims(2)))
                call load_shift_vectors(shift_vectors, extracted_shift_vectors_file, ierr)
            else
                ierr = metadata_ierr
                print *, "Error getting metadata for shift_vectors file"
            end if
            call delete_file(extracted_shift_vectors_file)
        end if
        
        ! Clean up the manifest file
        call delete_file("manifest.txt")
        
    contains
        subroutine delete_file(filename)
            character(len=*), intent(in) :: filename
            integer :: unit
            logical :: exists
            inquire(file=filename, exist=exists)
            if (exists) then
                open(newunit=unit, file=filename, status='old')
                close(unit, status='delete')
            end if
        end subroutine delete_file
    end subroutine read_tox_data

end module tox_archive