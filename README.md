![Logo](./trajslicer_logo.svg)
---
📄 Author: **Ouail Zakary**
- 📧 Email: [Ouail.Zakary@oulu.fi](mailto:Ouail.Zakary@oulu.fi)
- 🔗 ORCID: [0000-0002-7793-3306](https://orcid.org/0000-0002-7793-3306)
- 🌐 Website: [cc.oulu.fi/~nmrwww/members/Ouail_Zakary.html](https://cc.oulu.fi/~nmrwww/members/Ouail_Zakary.html)  
- 📁 Personal Website: [ozakary.github.io](https://ozakary.github.io/)
---
A versatile Python tool for converting LAMMPS dump files to XYZ format and sampling molecular dynamics trajectories with precise frame control.
## What's New in v3
- **Configurable chunk numbering offset** — `chunk_xyz_file` now accepts a `--chunk_start` argument that sets the first chunk label. Instead of always starting from `001`, you can now start from any index (e.g. `--chunk_start 10` produces `010`, `011`, ...). The zero-padding width is computed from the highest label (`chunk_start + num_chunks - 1`), so output filenames stay consistently formatted regardless of the chosen offset. Default behavior (starting from `1`) is unchanged.

## What's New in v2
- **Fully streaming I/O throughout** — no function in the codebase loads an entire trajectory into memory at any point. All operations (conversion, sampling, chunking) now process one frame at a time regardless of trajectory size. This fixes `MemoryError` crashes on large trajectories (tested on multi-GB XYZ files with 15,000+ atoms per snapshot).
- **Fast frame counting** — the first pass over the file reads only the atom-count header line of each frame (1 line per frame instead of `n_atoms + 2`), making the counting step orders of magnitude faster on large trajectories.
- **Direct chunk streaming** — `chunk_xyz_file` no longer buffers all selected frames before writing. Frames are written to their destination chunk file as they are read, with all output files open simultaneously. Peak memory is one frame at a time regardless of how many chunks are requested.
- **Consistent O(1) memory guarantee** — all three main operations (`convert_lammps_to_xyz`, `sample_xyz_file`, `chunk_xyz_file`) now share the same streaming design. The old `sample_xyz_file` already partially streamed; it has been tightened to skip non-selected frames with bare `readline()` calls rather than storing atom lines.
## Features
- Process both LAMMPS dump files and existing XYZ files
- Extract every nth frame with customizable sampling rates
- Specify exact start and end snapshots (0-based indexing)
- Assign and keep only specific atom types (LAMMPS files only)
- Map LAMMPS atom types to element symbols
- Divide the MD trajectory into sub-trajectories (chunks)
- Intelligently detects input file format
- Real-time progress bars with `tqdm`
- Processes trajectories of any size with constant memory usage
## Installation
### Requirements
- Python 3.6+
- `tqdm` (optional, for progress bars)
### Setup
1. Clone this repository:
```bash
git clone https://github.com/ozakary/TrajSlicer.git
cd TrajSlicer
```
2. Install optional dependencies:
```bash
pip install tqdm
```
3. Make the script executable:
```bash
chmod +x trajslicer_src.py
```
## Usage
### Basic Syntax
```bash
python trajslicer_src.py input_file output_file [options]
```
### Command-Line Options
| Option | Description | Example |
|--------|-------------|---------|
| `--sample N` | Sample every Nth frame | `--sample 10` |
| `--start N` | Starting snapshot index (0-based) | `--start 100` |
| `--end N` | Ending snapshot index (0-based, inclusive) | `--end 999` |
| `--chunks N` | Sequentially divide the trajectory into N chunks | `--chunks 10` |
| `--chunk_start N` | Starting index for chunk file numbering (default: 1) | `--chunk_start 10` |
| `--filter TYPE [TYPE ...]` | Keep only specified atom types (LAMMPS only) | `--filter 1 2` |
| `--labels TYPE:ELEMENT [...]` | Custom element labels (LAMMPS only) | `--labels 1:C 2:Xe` |
## Examples
### LAMMPS Dump File Conversion
**Convert entire trajectory:**
```bash
python trajslicer_src.py production.dump trajectory.xyz
```
**Sample every 10th frame:**
```bash
python trajslicer_src.py production.dump sampled.xyz --sample 10
```
**Extract first 1000 snapshots:**
```bash
python trajslicer_src.py production.dump first_1000.xyz --start 0 --end 999
```
**Keep only Xenon atoms (type 2), every 5th frame:**
```bash
python trajslicer_src.py production.dump xe_only.xyz --filter 2 --sample 5
```
**Custom element labels and range selection:**
```bash
python trajslicer_src.py production.dump custom.xyz \
    --labels 1:C 2:Xe \
    --start 500 --end 1500 \
    --sample 2
```
### XYZ File Sampling
**Sample existing XYZ file (every 10th frame):**
```bash
python trajslicer_src.py large_trajectory.xyz sampled.xyz --sample 10
```
**Extract specific frame range from XYZ:**
```bash
python trajslicer_src.py trajectory.xyz subset.xyz --start 1000 --end 2000
```
**Combine range and sampling for XYZ:**
```bash
python trajslicer_src.py trajectory.xyz final.xyz \
    --start 0 --end 5000 \
    --sample 25
```
**Divide MD trajectory into chunks:**
```bash
python trajslicer_src.py trajectory.xyz final.xyz \
    --start 0 --end 5000 \
    --sample 25 \
    --chunks 10
```
**Divide MD trajectory into chunks with custom starting index:**
```bash
python trajslicer_src.py trajectory.xyz final.xyz \
    --start 0 --end 5000 \
    --sample 25 \
    --chunks 91 \
    --chunk_start 10
```
This produces `final_chunk_010.xyz` through `final_chunk_100.xyz`.
## File Format Support
### Input Formats
- **LAMMPS Dump Files**: Standard LAMMPS trajectory files with `ITEM:` headers
- **XYZ Files**: Standard XYZ molecular coordinate files
### Output Format
- **XYZ Files**: Standard XYZ format with atom counts, comments, and coordinates
- **Comment Lines**: Include timestep information and lattice parameters (when available)
### Example Output (XYZ format)
```
1000
Timestep=100000 Lattice="50.0 0.0 0.0 0.0 50.0 0.0 0.0 0.0 50.0" Properties=species:S:1:pos:R:3:id:I:1
C 10.5 20.3 15.7 4231
C 11.2 21.1 16.4 4232
Xe 25.8 30.2 25.1 9001
...
```
## Technical Details
### Frame Indexing
- All frame indices are **0-based**
- `--end` parameter is **inclusive**
- Example: `--start 0 --end 999` extracts exactly 1000 frames (indices 0–999)
### Sampling Logic
- Sampling is applied **after** frame range selection
- Formula: `(current_frame - start_frame) % sample_rate == 0`
- Example: `--start 100 --end 200 --sample 2` extracts frames 100, 102, 104, ..., 200
### Memory Usage
All operations stream the input file in two passes:
1. **Count pass** — reads only the atom-count header line of each frame to determine the total number of frames. For a trajectory with `N` frames of `M` atoms each, this reads `N` lines instead of `N × (M + 2)`.
2. **Write pass** — streams frames sequentially. Selected frames are written immediately to their output file; skipped frames are discarded with minimal I/O. Peak RAM usage is proportional to one frame, not the full trajectory.
This design makes TrajSlicer suitable for trajectories of any size, including multi-GB files that would otherwise cause `MemoryError` with in-memory approaches.
### Chunk Distribution
When `--chunks N` is used, frames are distributed as evenly as possible. If the total number of selected frames is not exactly divisible by `N`, the first `remainder` chunks each receive one extra frame. All output files are opened before the streaming pass begins to avoid repeated open/close overhead.
### Chunk Numbering
By default, chunk output files are numbered starting from `1` (e.g. `_chunk_001`, `_chunk_002`, ...). Use `--chunk_start` to offset the numbering — useful when splitting a trajectory in multiple separate runs and needing contiguous file naming across them. Zero-padding is computed automatically from the highest label so filenames sort correctly.
## Contributing
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request
## License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
## Troubleshooting
### Common Issues
**`MemoryError` on large trajectories**  
Upgrade to v2. The v1 `chunk_xyz_file` function loaded the entire trajectory into memory before writing. v2 streams all operations and has no such limitation.  
**"Error: Could not find x, y, z, or type columns"**  
Ensure your LAMMPS dump file includes position and type information. Check that the `ITEM: ATOMS` line contains the required columns.  
**"Labels must be in format '1:C 2:Xe'"**  
Use a colon to separate atom type number from element symbol. Example: `--labels 1:C 2:Xe 3:O`.  
**Output chunk files appear empty or truncated**  
Check that `--start` and `--end` fall within the actual frame range of your trajectory. The script will warn if the requested range exceeds the file contents.  
### Performance Tips  
- Use `--sample` to reduce output file size and processing time on dense trajectories
- Combine `--start` and `--end` to process only the relevant segment of a long trajectory
- Install `tqdm` for real-time progress monitoring
- For very large files (tens of GB), the count pass may take a minute — this is normal and only happens once per run
