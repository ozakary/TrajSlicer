# TrajSlicer

![Logo](./trajslicer_logo.svg)

---

📄 Author: **Ouail Zakary**  
- 📧 Email: [Ouail.Zakary@oulu.fi](mailto:Ouail.Zakary@oulu.fi)  
- 🔗 ORCID: [0000-0002-7793-3306](https://orcid.org/0000-0002-7793-3306)  
- 🌐 Website: [Personal Webpage](https://cc.oulu.fi/~nmrwww/members/Ouail_Zakary.html)  
- 📁 Portfolio: [GitHub Portfolio](https://ozakary.github.io/)

---
A versatile Python tool for converting LAMMPS dump files to XYZ format and sampling molecular dynamics trajectories with precise frame control.

## Features

- Process both LAMMPS dump files and existing XYZ files
- Extract every nth frame with customizable sampling rates
- Specify exact start and end snapshots (0-based indexing)
- Keep only specific atom types (LAMMPS files only)
- Map LAMMPS atom types to element symbols
- Intelligently detects input file format
- Real-time progress bars with `tqdm`
- Processes large trajectories without loading entire files into memory

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
| `--filter TYPE [TYPE ...]` | Keep only specified atom types (LAMMPS only) | `--filter 1 2` |
| `--labels TYPE:ELEMENT [TYPE:ELEMENT ...]` | Custom element labels (LAMMPS only) | `--labels 1:C 2:Xe` |

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
    --labels 1:Carbon 2:Xenon \
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
Timestep=100000 Lattice="50.0 0.0 0.0 0.0 50.0 0.0 0.0 0.0 50.0"
C 10.5 20.3 15.7
C 11.2 21.1 16.4
Xe 25.8 30.2 25.1
...
```

## Technical Details

### Frame Indexing
- All frame indices are **0-based**
- `--end` parameter is **inclusive**
- Example: `--start 0 --end 999` extracts exactly 1000 frames (indices 0-999)

### Sampling Logic
- Sampling is applied **after** frame range selection
- Formula: `(current_frame - start_frame) % sample_rate == 0`
- Example: `--start 100 --end 200 --sample 2` extracts frames 100, 102, 104, ..., 200

### Memory Usage
- Processes files frame-by-frame (constant memory usage)
- Suitable for trajectories of any size
- Progress tracking with minimal overhead

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Troubleshooting

### Common Issues

**"Error: Could not find x, y, z, or type columns"**
- Ensure your LAMMPS dump file includes position and type information
- Check that the `ITEM: ATOMS` line contains the required columns

**"Labels must be in format '1:C 2:Xe'"**
- Use colon to separate atom type number from element symbol
- Example: `--labels 1:C 2:Xe 3:O`

**Memory issues with large files**
- Use frame range selection to process in smaller chunks
- Consider sampling to reduce output file size

### Performance Tips
- Use `--sample` for faster processing of large trajectories
- Combine `--start` and `--end` to process specific trajectory segments
- Install `tqdm` for progress monitoring on long conversions
