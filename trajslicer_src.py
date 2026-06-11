import os
import sys


def check_tqdm_installed():
    """Check if tqdm is installed, if not suggest installing it"""
    try:
        import tqdm
        return True
    except ImportError:
        print("The tqdm package is not installed. For a progress bar, install it using:")
        print("pip install tqdm")
        return False


def detect_file_type(input_file):
    """
    Detect if the input file is LAMMPS dump format or XYZ format.
    Returns 'lammps' or 'xyz'.
    """
    with open(input_file, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith("ITEM: TIMESTEP"):
            return 'lammps'
        try:
            int(first_line)
            return 'xyz'
        except ValueError:
            return 'lammps'


def convert_frame_to_xyz(frame_lines, outfile, filter_type=None, atom_labels=None, index_assignments=None):
    """
    Convert a single LAMMPS frame to XYZ format with proper triclinic box handling.

    Output XYZ format: element x y z atom_id
    Header includes: Timestep, Lattice, and Properties=species:S:1:pos:R:3:id:I:1
    """
    atom_count_index = None
    atoms_section_index = None
    box_bounds_index = None
    timestep_value = None

    for i, line in enumerate(frame_lines):
        if i == 1:
            timestep_value = line
        elif line.startswith("ITEM: NUMBER OF ATOMS"):
            atom_count_index = i
        elif line.startswith("ITEM: BOX BOUNDS"):
            box_bounds_index = i
        elif line.startswith("ITEM: ATOMS"):
            atoms_section_index = i
            break

    if atom_count_index is None or atoms_section_index is None:
        return

    atom_header = frame_lines[atoms_section_index]
    atom_lines = frame_lines[atoms_section_index + 1:]

    header_parts = atom_header.split()
    col_indices = {part: i - 2 for i, part in enumerate(header_parts[2:])}

    x_idx    = col_indices.get('x')
    y_idx    = col_indices.get('y')
    z_idx    = col_indices.get('z')
    type_idx = col_indices.get('type')
    id_idx   = col_indices.get('id')

    if x_idx is None or y_idx is None or z_idx is None or type_idx is None:
        print("Error: Could not find x, y, z, or type columns in the atom header")
        return

    if id_idx is None:
        print("Warning: No 'id' column found. Using atom type as fallback for ID column.")

    xyz_atom_lines = []
    for line in atom_lines:
        if not line or line.startswith("ITEM:"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                atom_type = int(parts[2 + type_idx])
                atom_id   = int(parts[2 + id_idx]) if id_idx is not None else None

                if filter_type is None or atom_type in filter_type:
                    element = None
                    if index_assignments is not None and atom_id is not None and atom_id in index_assignments:
                        element = index_assignments[atom_id]
                    elif atom_labels is not None:
                        element = atom_labels.get(atom_type, f"Type{atom_type}")
                    else:
                        element = f"Type{atom_type}"

                    x = float(parts[2 + x_idx])
                    y = float(parts[2 + y_idx])
                    z = float(parts[2 + z_idx])
                    atom_id_to_write = atom_id if atom_id is not None else atom_type
                    xyz_atom_lines.append(f"{element} {x} {y} {z} {atom_id_to_write}")
            except (ValueError, IndexError):
                continue

    lattice_matrix = None
    if box_bounds_index is not None and box_bounds_index + 3 <= len(frame_lines):
        try:
            box_header = frame_lines[box_bounds_index]
            is_triclinic = "xy xz yz" in box_header

            x_bounds = frame_lines[box_bounds_index + 1].split()
            y_bounds = frame_lines[box_bounds_index + 2].split()
            z_bounds = frame_lines[box_bounds_index + 3].split()

            if is_triclinic:
                xlo_bound, xhi_bound, xy = float(x_bounds[0]), float(x_bounds[1]), float(x_bounds[2])
                ylo_bound, yhi_bound, xz = float(y_bounds[0]), float(y_bounds[1]), float(y_bounds[2])
                zlo_bound, zhi_bound, yz = float(z_bounds[0]), float(z_bounds[1]), float(z_bounds[2])

                xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
                xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
                ylo = ylo_bound - min(0.0, yz)
                yhi = yhi_bound - max(0.0, yz)
                zlo = zlo_bound
                zhi = zhi_bound

                lx = xhi - xlo
                ly = yhi - ylo
                lz = zhi - zlo

                lattice_matrix = [
                    [lx,  0.0, 0.0],
                    [xy,  ly,  0.0],
                    [xz,  yz,  lz ],
                ]
            else:
                xlo, xhi = float(x_bounds[0]), float(x_bounds[1])
                ylo, yhi = float(y_bounds[0]), float(y_bounds[1])
                zlo, zhi = float(z_bounds[0]), float(z_bounds[1])

                lx = xhi - xlo
                ly = yhi - ylo
                lz = zhi - zlo

                lattice_matrix = [
                    [lx,  0.0, 0.0],
                    [0.0, ly,  0.0],
                    [0.0, 0.0, lz ],
                ]
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse box bounds: {e}")

    outfile.write(f"{len(xyz_atom_lines)}\n")

    comment = f"Timestep={timestep_value}"
    if lattice_matrix is not None:
        lattice_str = " ".join([f"{lattice_matrix[i][j]:.10f}" for i in range(3) for j in range(3)])
        comment += f" Lattice=\"{lattice_str}\""
    comment += " Properties=species:S:1:pos:R:3:id:I:1"
    outfile.write(f"{comment}\n")

    for line in xyz_atom_lines:
        outfile.write(f"{line}\n")


def convert_lammps_to_xyz(input_file, output_file, filter_type=None, sample_rate=1,
                           atom_labels=None, index_assignments=None,
                           start_frame=None, end_frame=None):
    """
    Convert a LAMMPS dump file to XYZ format.
    Streams frame by frame — memory usage is O(1 frame).
    """
    from tqdm import tqdm

    if atom_labels is None:
        atom_labels = {1: 'C', 2: 'Xe'}

    # Count frames
    total_frames = 0
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                total_frames += 1

    if start_frame is None:
        start_frame = 0
    if end_frame is None:
        end_frame = total_frames - 1

    start_frame = max(0, start_frame)
    end_frame   = min(total_frames - 1, end_frame)

    if start_frame > end_frame:
        print(f"Error: start_frame ({start_frame}) > end_frame ({end_frame})")
        return

    frames_to_process = end_frame - start_frame + 1

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        frame_lines   = []
        current_frame = 0
        frames_written = 0

        print(f"Converting LAMMPS dump to XYZ...")
        print(f"Processing snapshots {start_frame} to {end_frame} ({frames_to_process} frames)")
        pbar = tqdm(total=frames_to_process, unit="frames")

        for line in infile:
            line = line.strip()

            if line.startswith("ITEM: TIMESTEP"):
                if frame_lines and "ITEM: TIMESTEP" in frame_lines[0]:
                    if start_frame <= current_frame <= end_frame:
                        if (current_frame - start_frame) % sample_rate == 0:
                            convert_frame_to_xyz(frame_lines, outfile, filter_type,
                                                 atom_labels, index_assignments)
                            frames_written += 1
                        pbar.update(1)
                    current_frame += 1
                    if current_frame > end_frame:
                        break
                frame_lines = [line]
            else:
                frame_lines.append(line)

        # Last frame
        if frame_lines and start_frame <= current_frame <= end_frame:
            if (current_frame - start_frame) % sample_rate == 0:
                convert_frame_to_xyz(frame_lines, outfile, filter_type,
                                     atom_labels, index_assignments)
                frames_written += 1
            pbar.update(1)

        pbar.close()
        print(f"Done. Wrote {frames_written} frames to {output_file}")


def sample_xyz_file(input_file, output_file, sample_rate=1, start_frame=None, end_frame=None):
    """
    Sample an existing XYZ file with specified sampling rate and frame range.
    Streams frame by frame — memory usage is O(1 frame).
    """
    from tqdm import tqdm

    # Pass 1: count frames
    total_frames = 0
    with open(input_file, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            try:
                n_atoms = int(header.strip())
            except ValueError:
                break
            total_frames += 1
            f.readline()
            for _ in range(n_atoms):
                f.readline()

    if start_frame is None:
        start_frame = 0
    if end_frame is None:
        end_frame = total_frames - 1

    start_frame = max(0, start_frame)
    end_frame   = min(total_frames - 1, end_frame)

    if start_frame > end_frame:
        print(f"Error: start_frame ({start_frame}) > end_frame ({end_frame})")
        return

    frames_to_process = end_frame - start_frame + 1

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_frame  = 0
        frames_written = 0

        print(f"Sampling XYZ file...")
        print(f"Processing snapshots {start_frame} to {end_frame} ({frames_to_process} frames)")
        print(f"Sample rate: every {sample_rate} frame(s)")
        pbar = tqdm(total=frames_to_process, unit="frames")

        while current_frame <= end_frame:
            header = infile.readline()
            if not header:
                break
            try:
                n_atoms = int(header.strip())
            except ValueError:
                break

            comment = infile.readline()

            if start_frame <= current_frame <= end_frame:
                if (current_frame - start_frame) % sample_rate == 0:
                    outfile.write(header)
                    outfile.write(comment)
                    for _ in range(n_atoms):
                        outfile.write(infile.readline())
                    frames_written += 1
                else:
                    for _ in range(n_atoms):
                        infile.readline()
                pbar.update(1)
            else:
                for _ in range(n_atoms):
                    infile.readline()

            current_frame += 1

        pbar.close()
        print(f"Done. Wrote {frames_written} frames to {output_file}")


def chunk_xyz_file(input_file, output_base, num_chunks, sample_rate=1,
                   start_frame=None, end_frame=None, chunk_start=1):
    """
    Split an XYZ file into multiple chunk files.
    Streams frame by frame — no full file loaded into memory.
    """
    from tqdm import tqdm

    # ----------------------------------------------------------------
    # Pass 1: count frames by reading only header lines
    # ----------------------------------------------------------------
    print("Counting frames (fast pass)...")
    total_frames = 0
    with open(input_file, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            try:
                n_atoms = int(header.strip())
            except ValueError:
                break
            total_frames += 1
            f.readline()                      # skip comment
            for _ in range(n_atoms):          # skip atom lines
                f.readline()

    print(f"Total frames in file: {total_frames}")

    # ----------------------------------------------------------------
    # Resolve frame range
    # ----------------------------------------------------------------
    if start_frame is None:
        start_frame = 0
    if end_frame is None:
        end_frame = total_frames - 1

    start_frame = max(0, start_frame)
    end_frame   = min(total_frames - 1, end_frame)

    if start_frame > end_frame:
        print(f"Error: start_frame ({start_frame}) > end_frame ({end_frame})")
        return

    # Indices of frames that will actually be written
    selected_indices = [
        i for i in range(start_frame, end_frame + 1)
        if (i - start_frame) % sample_rate == 0
    ]
    total_selected = len(selected_indices)

    if total_selected == 0:
        print("No frames selected after applying filters.")
        return

    # ----------------------------------------------------------------
    # Assign selected frames to chunks
    # ----------------------------------------------------------------
    frames_per_chunk = total_selected // num_chunks
    remaining        = total_selected % num_chunks

    chunk_assignment = []
    for chunk_num in range(num_chunks):
        chunk_size = frames_per_chunk + (1 if chunk_num < remaining else 0)
        chunk_assignment.extend([chunk_num] * chunk_size)

    # O(1) lookup: global frame index -> chunk index
    frame_to_chunk = {
        selected_indices[k]: chunk_assignment[k]
        for k in range(total_selected)
    }

    print(f"Selected frames  : {total_selected}")
    print(f"Frames per chunk : {frames_per_chunk}"
          + (f" (first {remaining} chunks get one extra)" if remaining else ""))

    # ----------------------------------------------------------------
    # Open all output files upfront
    # ----------------------------------------------------------------
    base_name, ext = os.path.splitext(output_base)
    if not ext:
        ext = '.xyz'

    pad = len(str(chunk_start + num_chunks - 1))
    out_files = {}
    for chunk_num in range(num_chunks):
        fname = f"{base_name}_chunk_{chunk_num + chunk_start:0{pad}d}{ext}"
        out_files[chunk_num] = open(fname, 'w')
        print(f"  Chunk {chunk_num + chunk_start:0{pad}d} -> {fname}")

    # ----------------------------------------------------------------
    # Pass 2: stream file, write selected frames to their chunk
    # ----------------------------------------------------------------
    print("\nStreaming and writing frames...")
    frames_written = [0] * num_chunks

    with open(input_file, 'r') as f:
        pbar = tqdm(total=total_selected, unit="frames")
        current_frame = 0

        while current_frame <= end_frame:
            header = f.readline()
            if not header:
                break
            try:
                n_atoms = int(header.strip())
            except ValueError:
                break

            comment = f.readline()

            if current_frame in frame_to_chunk:
                chunk_num = frame_to_chunk[current_frame]
                out = out_files[chunk_num]
                out.write(header)
                out.write(comment)
                for _ in range(n_atoms):
                    out.write(f.readline())
                frames_written[chunk_num] += 1
                pbar.update(1)
            else:
                for _ in range(n_atoms):
                    f.readline()

            current_frame += 1

        pbar.close()

    for fh in out_files.values():
        fh.close()

    # ----------------------------------------------------------------
    # Summary
    # ----------------------------------------------------------------
    print(f"\nChunking complete! Created {num_chunks} chunk files.")
    for chunk_num in range(num_chunks):
        fname = f"{base_name}_chunk_{chunk_num + chunk_start:0{pad}d}{ext}"
        print(f"  Chunk {chunk_num + chunk_start:0{pad}d}: {frames_written[chunk_num]} frames -> {fname}")


# ====================================================================
# Entry point
# ====================================================================
if __name__ == "__main__":
    import argparse

    check_tqdm_installed()

    parser = argparse.ArgumentParser(
        description='Convert LAMMPS dump to XYZ, sample XYZ files, or split XYZ into chunks'
    )
    parser.add_argument('input_file',  help='Path to input file (LAMMPS dump or XYZ)')
    parser.add_argument('output_file', help='Path to output XYZ file (or base name for chunks)')
    parser.add_argument('--filter', type=int, nargs='+',
                        help='Atom types to keep (LAMMPS files only)')
    parser.add_argument('--sample', type=int, default=1,
                        help='Sample rate for frames (default: 1 = keep all)')
    parser.add_argument('--labels', type=str, nargs='+',
                        help='Atom type labels, format: 1:C 2:Xe (LAMMPS files only)')
    parser.add_argument('--index_assignments', type=str, nargs='+',
                        help='Element by atom index, format: 1:H 2:C (LAMMPS files only)')
    parser.add_argument('--start', type=int, default=None,
                        help='Starting snapshot index, 0-based (default: 0)')
    parser.add_argument('--end', type=int, default=None,
                        help='Ending snapshot index, 0-based inclusive (default: last frame)')
    parser.add_argument('--chunks', type=int, default=None,
                        help='Split XYZ file into N chunks (XYZ files only)')
    parser.add_argument('--chunk_start', type=int, default=1,
                        help='Starting index for chunk file numbering (default: 1)')

    args = parser.parse_args()

    file_type = detect_file_type(args.input_file)
    print(f"Detected file type: {file_type.upper()}")

    if file_type == 'xyz':
        if args.filter is not None:
            print("Warning: --filter is not available for XYZ files")
        if args.labels is not None:
            print("Warning: --labels is not available for XYZ files")
        if args.index_assignments is not None:
            print("Warning: --index_assignments is not available for XYZ files")

        if args.chunks is not None:
            if args.chunks <= 0:
                print("Error: number of chunks must be > 0")
                sys.exit(1)
            chunk_xyz_file(
                args.input_file,
                args.output_file,
                args.chunks,
                sample_rate=args.sample,
                start_frame=args.start,
                end_frame=args.end,
                chunk_start=args.chunk_start,
            )
        else:
            sample_xyz_file(
                args.input_file,
                args.output_file,
                sample_rate=args.sample,
                start_frame=args.start,
                end_frame=args.end,
            )

    else:  # LAMMPS
        if args.chunks is not None:
            print("Warning: --chunks is only available for XYZ files")

        atom_labels = {1: 'C', 2: 'Xe'}
        if args.labels:
            try:
                atom_labels = {}
                for pair in args.labels:
                    t, e = pair.split(':')
                    atom_labels[int(t)] = e
            except ValueError:
                print("Error: --labels must be in format '1:C 2:Xe'")
                sys.exit(1)

        index_assignments = None
        if args.index_assignments:
            try:
                index_assignments = {}
                for pair in args.index_assignments:
                    idx, e = pair.split(':')
                    index_assignments[int(idx)] = e
            except ValueError:
                print("Error: --index_assignments must be in format '1:H 2:C'")
                sys.exit(1)

        convert_lammps_to_xyz(
            args.input_file,
            args.output_file,
            filter_type=args.filter,
            sample_rate=args.sample,
            atom_labels=atom_labels,
            index_assignments=index_assignments,
            start_frame=args.start,
            end_frame=args.end,
        )
