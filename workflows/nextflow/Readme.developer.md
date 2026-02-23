# Understanding Nextflow Work Directory Structure

This document explains the Nextflow work directory structure and the purpose of each file, particularly the `.command.run` wrapper script.

## Work Directory Example

After running the GATK FastQC process, Nextflow creates this directory structure:

```
work/8d/4e3739526e64f6d4cea84e25ea7788/
├── .command.begin          # Timestamp marker for task start
├── .command.err            # Standard error output
├── .command.log            # Combined stdout + stderr
├── .command.out            # Standard output
├── .command.run            # Nextflow wrapper script (executes .command.sh)
├── .command.sh             # The actual process script
├── .exitcode               # Exit status (0 = success)
├── sample1_R1_fastqc.html  # Process output
├── sample1_R1_fastqc.zip   # Process output
├── sample1_R1.fastq.gz ->  # Symlink to input file
├── sample1_R2_fastqc.html  # Process output
├── sample1_R2_fastqc.zip   # Process output
├── sample1_R2.fastq.gz ->  # Symlink to input file
└── versions.yml            # Process output
```

## Key Files Explained

### `.command.sh` - The Process Script

This is the **actual command** you wrote in your Nextflow process:

```bash
#!/bin/bash -euo pipefail
# Quality Control with FastQC
# This generates HTML and ZIP files with quality metrics
fastqc \
     \
    --threads 2 \
    sample1_R1.fastq.gz sample1_R2.fastq.gz

# Create versions file
cat <<-END_VERSIONS > versions.yml
"GATK_VARIANT_CALLING:FASTQC":
    fastqc: $(fastqc --version | sed 's/FastQC v//')
END_VERSIONS
```

**Key points:**
- This is a direct translation of your process script block
- Uses `bash -euo pipefail` (from `nextflow.config`)
- Input files are referenced by their staged names
- This is what actually gets executed

### `.command.run` - The Wrapper Script

Nextflow wraps `.command.sh` with `.command.run`, which provides:

1. **Error handling and cleanup**
2. **File staging (input symlinks)**
3. **Output capture and logging**
4. **Exit code management**
5. **Signal handling (SIGTERM, etc.)**

#### Key sections of `.command.run`:

**1. Task Metadata (Lines 1-8)**
```bash
#!/bin/bash
### ---
### name: 'GATK_VARIANT_CALLING:FASTQC (sample1)'
### outputs:
### - '*.html'
### - '*.zip'
### - 'versions.yml'
### ...
```
Documents the task name and expected outputs

**2. Error Handling Setup (Lines 9-11)**
```bash
set -e
set -u
NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
```
Enables strict error checking

**3. Helper Functions (Lines 15-83)**
- `nxf_sleep()` - Sleep with fallback
- `nxf_date()` - Cross-platform timestamps
- `nxf_env()` - Environment variable logging
- `nxf_kill()` - Process tree cleanup
- File staging functions (`nxf_fs_copy`, `nxf_fs_move`, etc.)

**4. File Staging (Lines 104-111)**
```bash
nxf_stage() {
    true
    # stage input files
    rm -f sample1_R1.fastq.gz
    rm -f sample1_R2.fastq.gz
    ln -s /home/.../data/sample1_R1.fastq.gz sample1_R1.fastq.gz
    ln -s /home/.../data/sample1_R2.fastq.gz sample1_R2.fastq.gz
}
```
Creates symlinks to input files in the work directory

**5. Command Execution (Line 101)**
```bash
nxf_launch() {
    /bin/bash -euo pipefail /home/.../work/.../. command.sh
}
```
Executes the actual process script

**6. Output Capture (Lines 145-147)**
```bash
(set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
pid=$!
wait $pid || nxf_main_ret=$?
```
Captures stdout to `.command.out` and stderr to `.command.err`

**7. Exit Code Management (Lines 85-93)**
```bash
on_exit() {
    local last_err=$?
    local exit_status=${nxf_main_ret:=0}
    [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
    [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
    printf -- $exit_status > .exitcode
    set +u
    exit $exit_status
}
```
Writes the final exit code to `.exitcode` file

## Why This Matters for Validation

When comparing bash to Nextflow:

1. **The `.command.sh` is what matters** - This is your actual process logic
2. **The `.command.run` is infrastructure** - Nextflow's execution wrapper
3. **Symlinks vs copies** - Nextflow uses symlinks for efficiency (same data)
4. **Work directory isolation** - Each task runs in its own directory

## Bash vs Nextflow Execution Comparison

### Bash Approach
```bash
# Direct execution in a single directory
fastqc -o results/qc/sample1 -t 2 data/sample1_R1.fastq.gz data/sample1_R2.fastq.gz
```

### Nextflow Approach
```bash
# 1. Create isolated work directory
mkdir -p work/8d/4e3739...

# 2. Stage inputs (symlink)
ln -s /full/path/data/sample1_R1.fastq.gz work/8d/4e3739.../sample1_R1.fastq.gz

# 3. Execute in work directory
cd work/8d/4e3739...
fastqc --threads 2 sample1_R1.fastq.gz sample1_R2.fastq.gz

# 4. Publish outputs (copy to results)
cp *.html results_nextflow/qc/sample1/
cp *.zip results_nextflow/qc/sample1/
```

## Benefits of the Work Directory Pattern

1. **Isolation** - Each task has its own directory (parallel-safe)
2. **Resume** - Can restart from cached work directories (`-resume`)
3. **Debugging** - All artifacts preserved (`.command.log`, `.exitcode`)
4. **Reproducibility** - Complete execution trace
5. **Cleanup** - Can delete entire `work/` directory when done

## Debugging Tips

### View what was executed:
```bash
cat work/8d/4e3739.../. command.sh
```

### Check output logs:
```bash
cat work/8d/4e3739.../. command.log
```

### See exit code:
```bash
cat work/8d/4e3739.../. exitcode
```

### Re-run manually:
```bash
cd work/8d/4e3739...
bash .command.run
```

## Summary

- `.command.sh` = Your process script (the "what")
- `.command.run` = Nextflow's wrapper (the "how")
- Work directory = Isolated execution environment
- Symlinks = Efficient input staging
- The **scientific logic is identical** between bash and Nextflow
- The **execution infrastructure differs** (isolation, logging, caching)
