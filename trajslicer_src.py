def detect_file_type(input_file):
    """
    Detect if the input file is LAMMPS dump format or XYZ format
    
    Args:
        input_file (str): Path to input file
        
    Returns:
        str: 'lammps' or 'xyz'
    """
    with open(input_file, 'r') as f:
        first_line = f.readline().strip()
        
        # LAMMPS dump files start with "ITEM: TIMESTEP"
        if first_line.startswith("ITEM: TIMESTEP"):
            return 'lammps'
        
        # XYZ files start with a number (atom count)
        try:
            int(first_line)
            return 'xyz'
        except ValueError:
            # If neither, assume LAMMPS (default behavior)
            return 'lammps'

def convert_lammps_to_xyz(input_file, output_file, filter_type=None, sample_rate=1, 
                     atom_labels=None, index_assignments=None, start_frame=None, end_frame=None):
    """
    Convert a LAMMPS dump file to XYZ format with customizable options.
    
    Args:
        input_file (str): Path to input LAMMPS dump file
        output_file (str): Path to output XYZ file
        filter_type (list, optional): List of atom types to keep. If None, keep all atoms.
        sample_rate (int, optional): Sample rate for frames (1 = every frame, 2 = every other, etc.)
        atom_labels (dict, optional): Mapping of atom types to element labels.
                                     Example: {1: 'C', 2: 'Xe'}
        index_assignments (dict, optional): Mapping of atom indices to element labels.
                                          Example: {1: 'H', 2: 'C', 3: 'N'}
                                          Takes precedence over atom_labels when specified.
        start_frame (int, optional): Starting snapshot index (0-based). If None, start from beginning.
        end_frame (int, optional): Ending snapshot index (0-based, inclusive). If None, go to end.
    
    Note:
        Output XYZ format: element x y z atom_id
        Header includes: Timestep, Lattice, and Properties=species:S:1:pos:R:3:id:I:1
    """
    import os
    from tqdm import tqdm
    
    # Default atom labels if not provided
    if atom_labels is None:
        atom_labels = {1: 'C', 2: 'Xe'}
    
    # Count number of frames for progress bar
    total_frames = 0
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                total_frames += 1
    
    # Set default values for start_frame and end_frame
    if start_frame is None:
        start_frame = 0
    if end_frame is None:
        end_frame = total_frames - 1
    
    # Validate frame range
    if start_frame < 0:
        start_frame = 0
    if end_frame >= total_frames:
        end_frame = total_frames - 1
    if start_frame > end_frame:
        print(f"Error: start_frame ({start_frame}) cannot be greater than end_frame ({end_frame})")
        return
    
    frames_to_process = end_frame - start_frame + 1
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        frame_lines = []
        current_section = None
        current_frame = 0
        frames_written = 0
        
        # We'll collect complete frames and process them
        print(f"Converting LAMMPS dump file to XYZ format...")
        print(f"Processing snapshots {start_frame} to {end_frame} (total: {frames_to_process} frames)")
        pbar = tqdm(total=frames_to_process, unit="frames")
        
        for line in infile:
            line = line.strip()
            
            if line.startswith("ITEM: TIMESTEP"):
                # If we have collected a previous frame, process it
                if frame_lines and "ITEM: TIMESTEP" in frame_lines[0]:
                    # Check if this frame is within our desired range
                    if start_frame <= current_frame <= end_frame:
                        # Check if this frame matches our sampling rate
                        if (current_frame - start_frame) % sample_rate == 0:
                            convert_frame_to_xyz(frame_lines, outfile, filter_type, atom_labels, index_assignments)
                            frames_written += 1
                        pbar.update(1)
                    
                    current_frame += 1
                    
                    # If we've passed the end_frame, we can stop processing
                    if current_frame > end_frame:
                        break
                
                # Start a new frame
                frame_lines = [line]
                current_section = "timestep"
            else:
                frame_lines.append(line)
        
        # Process the last frame if it's within range
        if frame_lines and start_frame <= current_frame <= end_frame:
            if (current_frame - start_frame) % sample_rate == 0:
                convert_frame_to_xyz(frame_lines, outfile, filter_type, atom_labels, index_assignments)
                frames_written += 1
            pbar.update(1)
        
        pbar.close()
        print(f"Conversion complete. Output saved to {output_file}")
        print(f"Processed frames {start_frame} to {end_frame}, wrote {frames_written} frames")

def sample_xyz_file(input_file, output_file, sample_rate=1, start_frame=None, end_frame=None):
    """
    Sample an existing XYZ file with specified sampling rate and frame range.
    
    Args:
        input_file (str): Path to input XYZ file
        output_file (str): Path to output XYZ file
        sample_rate (int, optional): Sample rate for frames (1 = every frame, 2 = every other, etc.)
        start_frame (int, optional): Starting snapshot index (0-based). If None, start from beginning.
        end_frame (int, optional): Ending snapshot index (0-based, inclusive). If None, go to end.
    """
    from tqdm import tqdm
    
    # First pass: count total frames
    total_frames = 0
    with open(input_file, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            
            # First line of each frame is atom count
            try:
                atom_count = int(line.strip())
                total_frames += 1
                
                # Skip comment line
                f.readline()
                
                # Skip atom lines
                for _ in range(atom_count):
                    f.readline()
            except ValueError:
                break
    
    # Set default values for start_frame and end_frame
    if start_frame is None:
        start_frame = 0
    if end_frame is None:
        end_frame = total_frames - 1
    
    # Validate frame range
    if start_frame < 0:
        start_frame = 0
    if end_frame >= total_frames:
        end_frame = total_frames - 1
    if start_frame > end_frame:
        print(f"Error: start_frame ({start_frame}) cannot be greater than end_frame ({end_frame})")
        return
    
    frames_to_process = end_frame - start_frame + 1
    
    # Second pass: process and sample frames
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_frame = 0
        frames_written = 0
        
        print(f"Sampling XYZ file...")
        print(f"Processing snapshots {start_frame} to {end_frame} (total: {frames_to_process} frames)")
        print(f"Sample rate: every {sample_rate} frame(s)")
        pbar = tqdm(total=frames_to_process, unit="frames")
        
        while True:
            # Read atom count line
            atom_count_line = infile.readline()
            if not atom_count_line:
                break
            
            try:
                atom_count = int(atom_count_line.strip())
            except ValueError:
                break
            
            # Read comment line
            comment_line = infile.readline()
            if not comment_line:
                break
            
            # Read atom lines
            atom_lines = []
            for _ in range(atom_count):
                atom_line = infile.readline()
                if not atom_line:
                    break
                atom_lines.append(atom_line)
            
            # Check if this frame is within our desired range
            if start_frame <= current_frame <= end_frame:
                # Check if this frame matches our sampling rate
                if (current_frame - start_frame) % sample_rate == 0:
                    # Write this frame to output
                    outfile.write(atom_count_line)
                    outfile.write(comment_line)
                    for atom_line in atom_lines:
                        outfile.write(atom_line)
                    frames_written += 1
                pbar.update(1)
            
            current_frame += 1
            
            # If we've passed the end_frame, we can stop processing
            if current_frame > end_frame:
                break
        
        pbar.close()
        print(f"Sampling complete. Output saved to {output_file}")
        print(f"Processed frames {start_frame} to {end_frame}, wrote {frames_written} frames")

def chunk_xyz_file(input_file, output_base, num_chunks, sample_rate=1, start_frame=None, end_frame=None):
    """
    Split an XYZ file into multiple chunk files with specified number of chunks.
    
    Args:
        input_file (str): Path to input XYZ file
        output_base (str): Base name for output files (without extension)
        num_chunks (int): Number of chunks to create
        sample_rate (int, optional): Sample rate for frames (1 = every frame, 2 = every other, etc.)
        start_frame (int, optional): Starting snapshot index (0-based). If None, start from beginning.
        end_frame (int, optional): Ending snapshot index (0-based, inclusive). If None, go to end.
    """
    from tqdm import tqdm
    import os
    
    # First pass: count total frames and collect frame data
    total_frames = 0
    frame_data = []  # Store (atom_count_line, comment_line, atom_lines) for each frame
    
    print("Reading XYZ file...")
    with open(input_file, 'r') as f:
        while True:
            # Read atom count line
            atom_count_line = f.readline()
            if not atom_count_line:
                break
            
            try:
                atom_count = int(atom_count_line.strip())
            except ValueError:
                break
            
            # Read comment line
            comment_line = f.readline()
            if not comment_line:
                break
            
            # Read atom lines
            atom_lines = []
            for _ in range(atom_count):
                atom_line = f.readline()
                if not atom_line:
                    break
                atom_lines.append(atom_line)
            
            frame_data.append((atom_count_line, comment_line, atom_lines))
            total_frames += 1
    
    # Set default values for start_frame and end_frame
    if start_frame is None:
        start_frame = 0
    if end_frame is None:
        end_frame = total_frames - 1
    
    # Validate frame range
    if start_frame < 0:
        start_frame = 0
    if end_frame >= total_frames:
        end_frame = total_frames - 1
    if start_frame > end_frame:
        print(f"Error: start_frame ({start_frame}) cannot be greater than end_frame ({end_frame})")
        return
    
    # Apply frame range filtering and sampling
    selected_frames = []
    for i in range(start_frame, end_frame + 1):
        if (i - start_frame) % sample_rate == 0:
            selected_frames.append(frame_data[i])
    
    total_selected_frames = len(selected_frames)
    if total_selected_frames == 0:
        print("No frames selected after applying filters.")
        return
    
    # Calculate frames per chunk
    frames_per_chunk = total_selected_frames // num_chunks
    remaining_frames = total_selected_frames % num_chunks
    
    print(f"Total frames after filtering and sampling: {total_selected_frames}")
    print(f"Creating {num_chunks} chunks...")
    print(f"Base frames per chunk: {frames_per_chunk}")
    if remaining_frames > 0:
        print(f"First {remaining_frames} chunks will have {frames_per_chunk + 1} frames")
    
    # Prepare output file name template
    base_name, ext = os.path.splitext(output_base)
    if not ext:
        ext = '.xyz'
    
    # Split frames into chunks and write to separate files
    frame_idx = 0
    for chunk_num in range(num_chunks):
        # Calculate number of frames for this chunk
        chunk_size = frames_per_chunk + (1 if chunk_num < remaining_frames else 0)
        
        # Generate output filename
        output_file = f"{base_name}_chunk_{chunk_num + 1:0{len(str(num_chunks))}d}{ext}"
        
        print(f"Writing chunk {chunk_num + 1}/{num_chunks}: {chunk_size} frames -> {output_file}")
        
        with open(output_file, 'w') as outfile:
            for _ in range(chunk_size):
                if frame_idx < total_selected_frames:
                    atom_count_line, comment_line, atom_lines = selected_frames[frame_idx]
                    outfile.write(atom_count_line)
                    outfile.write(comment_line)
                    for atom_line in atom_lines:
                        outfile.write(atom_line)
                    frame_idx += 1
        
        print(f"  -> Saved {chunk_size} frames to {output_file}")
    
    print(f"\nChunking complete! Created {num_chunks} chunk files.")

def convert_frame_to_xyz(frame_lines, outfile, filter_type=None, atom_labels=None, index_assignments=None):
    """
    Convert a single LAMMPS frame to XYZ format
    
    Args:
        frame_lines (list): Lines for a single frame
        outfile (file): Output file to write to
        filter_type (list, optional): List of atom types to keep. If None, keep all atoms.
        atom_labels (dict): Mapping of atom types to element labels
        index_assignments (dict, optional): Mapping of atom indices to element labels.
                                          Takes precedence over atom_labels when specified.
    
    Note:
        Output XYZ format includes atom ID as additional column: element x y z atom_id
        Header includes Properties specification: Properties=species:S:1:pos:R:3:id:I:1
    """
    # Find indices of important sections
    atom_count_index = None
    atoms_section_index = None
    box_bounds_index = None
    timestep_value = None
    
    for i, line in enumerate(frame_lines):
        if i == 1:  # Second line is the timestep value
            timestep_value = line
        elif line.startswith("ITEM: NUMBER OF ATOMS"):
            atom_count_index = i
        elif line.startswith("ITEM: BOX BOUNDS"):
            box_bounds_index = i
        elif line.startswith("ITEM: ATOMS"):
            atoms_section_index = i
            break
    
    if atom_count_index is None or atoms_section_index is None:
        return  # Skip this frame if it doesn't have the required sections
    
    # Extract atom lines (skip the header)
    atom_header = frame_lines[atoms_section_index]
    atom_lines = frame_lines[atoms_section_index + 1:]
    
    # Get the column indices for x, y, z and type from the atom header
    header_parts = atom_header.split()
    col_indices = {part: i-2 for i, part in enumerate(header_parts[2:])}
    
    # We need positions (x, y, z), type, and optionally id
    x_idx = col_indices.get('x')
    y_idx = col_indices.get('y')
    z_idx = col_indices.get('z')
    type_idx = col_indices.get('type')
    id_idx = col_indices.get('id')  # Atom ID/index for index assignments
    
    if x_idx is None or y_idx is None or z_idx is None or type_idx is None:
        print("Error: Could not find x, y, z, or type columns in the atom header")
        return
    
    # Warn if no id column found - atom ID is important for tracking atoms
    if id_idx is None:
        print("Warning: No 'id' column found in LAMMPS dump file. Using atom type as fallback for ID column.")
        print("Available columns:", list(col_indices.keys()))
    
    # Filter atoms if filter_type is specified, otherwise keep all
    xyz_atom_lines = []
    for line in atom_lines:
        if not line or line.startswith("ITEM:"):
            continue  # Skip empty lines or new section headers
            
        parts = line.split()
        if len(parts) >= 2:
            try:
                atom_type = int(parts[2 + type_idx])
                atom_id = int(parts[2 + id_idx]) if id_idx is not None else None
                
                # Check if we should keep this atom type (for filtering)
                if filter_type is None or atom_type in filter_type:
                    # Determine element label: index assignments take precedence over type labels
                    element = None
                    if index_assignments is not None and atom_id is not None and atom_id in index_assignments:
                        element = index_assignments[atom_id]
                    elif atom_labels is not None:
                        element = atom_labels.get(atom_type, f"Type{atom_type}")
                    else:
                        element = f"Type{atom_type}"
                    
                    # XYZ format: element x y z atom_id
                    x = float(parts[2 + x_idx])
                    y = float(parts[2 + y_idx])
                    z = float(parts[2 + z_idx])
                    atom_id_to_write = atom_id if atom_id is not None else atom_type  # Fallback to type if no ID
                    xyz_atom_lines.append(f"{element} {x} {y} {z} {atom_id_to_write}")
            except (ValueError, IndexError):
                continue  # Skip lines that cause errors
    
    # Get box dimensions for XYZ comment line
    box_x = None
    box_y = None
    box_z = None
    
    if box_bounds_index is not None and box_bounds_index + 3 <= len(frame_lines):
        try:
            box_x_parts = frame_lines[box_bounds_index + 1].split()
            box_y_parts = frame_lines[box_bounds_index + 2].split()
            box_z_parts = frame_lines[box_bounds_index + 3].split()
            
            box_x = float(box_x_parts[1]) - float(box_x_parts[0])
            box_y = float(box_y_parts[1]) - float(box_y_parts[0])
            box_z = float(box_z_parts[1]) - float(box_z_parts[0])
        except (ValueError, IndexError):
            pass
    
    # Write the XYZ frame
    # First line: number of atoms
    outfile.write(f"{len(xyz_atom_lines)}\n")
    
    # Second line: comment line with box dimensions, timestep, and properties
    comment = f"Timestep={timestep_value}"
    if box_x is not None and box_y is not None and box_z is not None:
        comment += f" Lattice=\"{box_x} 0.0 0.0 0.0 {box_y} 0.0 0.0 0.0 {box_z}\""
    comment += " Properties=species:S:1:pos:R:3:id:I:1"
    outfile.write(f"{comment}\n")
    
    # Atom lines: element x y z atom_id
    for line in xyz_atom_lines:
        outfile.write(f"{line}\n")

def check_tqdm_installed():
    """Check if tqdm is installed, if not suggest installing it"""
    try:
        import tqdm
        return True
    except ImportError:
        print("The tqdm package is not installed. For a progress bar, install it using:")
        print("pip install tqdm")
        return False

if __name__ == "__main__":
    import sys
    import argparse
    
    # Check if tqdm is installed
    has_tqdm = check_tqdm_installed()
    if not has_tqdm:
        # Import a simple progress indicator as fallback
        from tqdm import tqdm
        def tqdm(iterable, **kwargs):
            return iterable
    
    # Create command-line argument parser
    parser = argparse.ArgumentParser(description='Convert LAMMPS dump file to XYZ format, sample XYZ files, or split XYZ files into chunks')
    parser.add_argument('input_file', help='Path to input file (LAMMPS dump or XYZ format)')
    parser.add_argument('output_file', help='Path to output XYZ file (or base name for chunks)')
    parser.add_argument('--filter', type=int, nargs='+', 
                        help='Atom types to keep (e.g., --filter 2 for Xe only) - LAMMPS files only')
    parser.add_argument('--sample', type=int, default=1,
                        help='Sample rate for frames (default: 1, meaning keep all frames)')
    parser.add_argument('--labels', type=str, nargs='+', 
                        help='Custom atom type labels (format: 1:C 2:Xe) - LAMMPS files only')
    parser.add_argument('--index_assignments', type=str, nargs='+',
                        help='Assign elements by atom index (format: 1:H 2:C 3:N) - LAMMPS files only')
    parser.add_argument('--start', type=int, default=None,
                        help='Starting snapshot index (0-based, default: 0)')
    parser.add_argument('--end', type=int, default=None,
                        help='Ending snapshot index (0-based, inclusive, default: last frame)')
    parser.add_argument('--chunks', type=int, default=None,
                        help='Split XYZ file into specified number of chunks - XYZ files only')
    
    args = parser.parse_args()
    
    # Detect file type
    file_type = detect_file_type(args.input_file)
    print(f"Detected file type: {file_type.upper()}")
    
    if file_type == 'xyz':
        # Handle XYZ file operations
        if args.filter is not None:
            print("Warning: --filter option is not available for XYZ files (only LAMMPS files)")
        if args.labels is not None:
            print("Warning: --labels option is not available for XYZ files (only LAMMPS files)")
        if args.index_assignments is not None:
            print("Warning: --index_assignments option is not available for XYZ files (only LAMMPS files)")
        
        if args.chunks is not None:
            # Split XYZ file into chunks
            if args.chunks <= 0:
                print("Error: Number of chunks must be greater than 0")
                sys.exit(1)
            
            chunk_xyz_file(
                args.input_file,
                args.output_file,
                args.chunks,
                sample_rate=args.sample,
                start_frame=args.start,
                end_frame=args.end
            )
        else:
            # Regular XYZ sampling
            sample_xyz_file(
                args.input_file,
                args.output_file,
                sample_rate=args.sample,
                start_frame=args.start,
                end_frame=args.end
            )
    
    else:  # LAMMPS file
        if args.chunks is not None:
            print("Warning: --chunks option is only available for XYZ files")
        
        # Process custom atom labels if provided
        atom_labels = {1: 'C', 2: 'Xe'}  # Default labels
        if args.labels:
            try:
                atom_labels = {}
                for label_pair in args.labels:
                    type_num, element = label_pair.split(':')
                    atom_labels[int(type_num)] = element
            except ValueError:
                print("Error: Labels must be in format '1:C 2:Xe'")
                sys.exit(1)
        
        # Process index assignments if provided
        index_assignments = None
        if args.index_assignments:
            try:
                index_assignments = {}
                for assignment_pair in args.index_assignments:
                    atom_index, element = assignment_pair.split(':')
                    index_assignments[int(atom_index)] = element
            except ValueError:
                print("Error: Index assignments must be in format '1:H 2:C 3:N'")
                sys.exit(1)
        
        convert_lammps_to_xyz(
            args.input_file, 
            args.output_file,
            filter_type=args.filter,
            sample_rate=args.sample,
            atom_labels=atom_labels,
            index_assignments=index_assignments,
            start_frame=args.start,
            end_frame=args.end
        )
