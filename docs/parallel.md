# Parallel Execution Guide

This document explains the parallel execution capabilities in HelixForge for processing large genomes efficiently on multi-core systems and HPC clusters.

## Overview

The parallel execution module (`helixforge.parallel`) provides:

1. **Flexible Genome Chunking**: Multiple strategies to partition genomes into processable units
2. **Task File Generation**: Create task files for HyperShell, GNU Parallel, or xargs
3. **Local Parallel Execution**: ProcessPoolExecutor and ThreadPoolExecutor backends
4. **Memory Monitoring**: Track and limit memory usage during processing
5. **SLURM Utilities**: Minimal helpers for users who need direct SLURM integration

## Recommended Workflow

The recommended approach for HPC parallelization uses **task files** that work with any scheduler:

```bash
# 1. Create a chunk plan
helixforge parallel plan --genome genome.fa -o chunks.json

# 2. Generate task file (include all required options for your command)
helixforge parallel tasks --chunk-plan chunks.json \
    --command 'helixforge confidence -p predictions.h5 -g predictions.gff3 --genome genome.fa --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o outputs/{chunk_id}.tsv' \
    --output tasks.txt --output-dir outputs/

# 3. Execute with HyperShell (recommended) or GNU Parallel
hs cluster tasks.txt --num-tasks 32
# OR
parallel -j 32 < tasks.txt

# 4. Aggregate results
helixforge parallel aggregate --input-dir outputs/ --pattern '*.tsv' -o combined.tsv --type merge_tsv
```

## Chunk-Aware Commands

The core HelixForge commands (`confidence`, `refine`, and `evidence`) support chunk-aware processing through these flags:

- `--region seqid:start-end` - Process only genes in this region (1-based inclusive coordinates)
- `--scaffold name` - Process only genes on this scaffold (simpler alternative to --region)
- `--chunk-id id` - Chunk identifier for logging and output organization

### Confidence Scoring by Region

```bash
# Generate tasks for confidence scoring with region filtering
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command "helixforge confidence \
        -p predictions.h5 \
        -g predictions.gff3 \
        --genome genome.fa \
        --region {seqid}:{start}-{end} \
        --chunk-id {chunk_id} \
        -o outputs/{chunk_id}_confidence.tsv" \
    --output tasks.txt \
    --output-dir outputs/
```

### Confidence Scoring by Scaffold

For genomes with many small scaffolds, scaffold-based chunking is simpler:

```bash
# Create scaffold-based chunk plan
helixforge parallel plan \
    --genome genome.fa \
    --strategy scaffold \
    --output chunks.json

# Generate tasks using --scaffold flag
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command "helixforge confidence \
        -p predictions.h5 \
        -g predictions.gff3 \
        --genome genome.fa \
        --scaffold {seqid} \
        --chunk-id {chunk_id} \
        -o outputs/{chunk_id}_confidence.tsv" \
    --output tasks.txt \
    --output-dir outputs/
```

### Gene Refinement (Main Pipeline)

```bash
# Single BAM file
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command "helixforge refine \
        -p predictions.h5 \
        -g predictions.gff3 \
        --genome genome.fa \
        --rnaseq-bam rnaseq.bam \
        --region {seqid}:{start}-{end} \
        --chunk-id {chunk_id} \
        -o outputs/{chunk_id}_refined.gff3 \
        -r outputs/{chunk_id}_report.tsv" \
    --output tasks.txt \
    --output-dir outputs/

# Multiple BAM files (multi-tissue)
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command "helixforge refine \
        -p predictions.h5 \
        -g predictions.gff3 \
        --genome genome.fa \
        --rnaseq-bam liver.bam,brain.bam,heart.bam \
        --min-tissues 2 \
        --min-reads 3 \
        --region {seqid}:{start}-{end} \
        --chunk-id {chunk_id} \
        -o outputs/{chunk_id}_refined.gff3 \
        -r outputs/{chunk_id}_report.tsv" \
    --output tasks.txt \
    --output-dir outputs/

# Using a BAM list file
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command "helixforge refine \
        -p predictions.h5 \
        -g predictions.gff3 \
        --genome genome.fa \
        --rnaseq-bam-list tissue_bams.txt \
        --min-tissues 2 \
        --region {seqid}:{start}-{end} \
        --chunk-id {chunk_id} \
        -o outputs/{chunk_id}_refined.gff3" \
    --output tasks.txt \
    --output-dir outputs/

# Using STAR SJ.out.tab junction files
helixforge parallel tasks \
    --chunk-plan chunks.json \
    --command "helixforge refine \
        -p predictions.h5 \
        -g predictions.gff3 \
        --genome genome.fa \
        --junctions-bed sample1_SJ.out.tab,sample2_SJ.out.tab,sample3_SJ.out.tab \
        --min-tissues 2 \
        --region {seqid}:{start}-{end} \
        --chunk-id {chunk_id} \
        -o outputs/{chunk_id}_refined.gff3" \
    --output tasks.txt \
    --output-dir outputs/
```

### Coordinate Convention Summary

| Context | Convention | Example |
|---------|------------|---------|
| User CLI input (--region) | 1-based inclusive | chr1:1000-2000 |
| Internal storage | 0-based half-open | (999, 2000) |
| GFF3 files | 1-based inclusive | 1000\t2000 |
| BED files | 0-based half-open | 999\t2000 |
| Task template {start}/{end} | 1-based inclusive | 1000-2000 |
| Task template {start_0}/{end_0} | 0-based half-open | 999-2000 |

The `parse_region()` function converts CLI input (1-based) to internal (0-based).
The `format_command()` function converts back when generating task commands.

### Why Task Files?

- **Scheduler-agnostic**: Works with SLURM, PBS, SGE, LSF, or local execution
- **Simple**: One command per line, easy to inspect and debug
- **Flexible**: Works with HyperShell, GNU Parallel, xargs, or custom scripts
- **HPC-friendly**: HyperShell handles distributed execution across nodes

## Chunking Strategies

HelixForge supports four chunking strategies optimized for different scenarios:

### BY_SCAFFOLD (default)

One chunk per scaffold/chromosome. Best for genomes with many scaffolds.

```python
from helixforge.parallel import GenomeChunker, ChunkStrategy

chunker = GenomeChunker(genome)
plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)
```

Use when:
- Draft assemblies with many contigs
- Small genomes where each scaffold fits in memory
- Even gene distribution across scaffolds

### BY_SIZE

Fixed base-pair windows. Best for genomes with few large chromosomes.

```python
plan = chunker.create_plan(
    ChunkStrategy.BY_SIZE,
    chunk_size=10_000_000,  # 10 Mb windows
    min_chunk_size=100_000,  # Don't create tiny chunks
)
```

Use when:
- Chromosome-level assemblies (maize, wheat)
- Need consistent memory usage per chunk
- Processing doesn't depend on gene boundaries

### BY_GENES

Fixed number of genes per chunk. Best for annotation-heavy tasks.

```python
# Requires GFF parser
chunker = GenomeChunker(genome, gff_parser)
plan = chunker.create_plan(
    ChunkStrategy.BY_GENES,
    chunk_size=100,  # 100 genes per chunk
)
```

Use when:
- Processing time scales with gene count
- Want balanced workloads across tasks
- Gene density varies significantly

### ADAPTIVE

Automatically balances chunk sizes based on genome content.

```python
plan = chunker.create_plan(
    ChunkStrategy.ADAPTIVE,
    target_chunks=50,  # Aim for ~50 chunks
    min_chunk_size=100_000,
)
```

Use when:
- Unsure which strategy is best
- Variable gene density
- Want automatic load balancing

## Using the CLI

### Create a Chunking Plan

```bash
# Scaffold-based (default)
helixforge parallel plan \
    --genome genome.fa \
    --strategy scaffold \
    --output chunks.json

# Size-based with 50 Mb chunks
helixforge parallel plan \
    --genome genome.fa \
    --strategy size \
    --chunk-size 50000000 \
    --output chunks.json

# Gene-based (requires GFF)
helixforge parallel plan \
    --genome genome.fa \
    --gff predictions.gff3 \
    --strategy genes \
    --chunk-size 200 \
    --output chunks.json
```

### Generate Task File

```bash
# Simple task file with direct commands (include all required options)
helixforge parallel tasks --chunk-plan chunks.json \
    --command 'helixforge confidence -p predictions.h5 -g predictions.gff3 --genome genome.fa --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o {output_dir}/{chunk_id}.tsv' \
    --output tasks.txt --output-dir outputs/

# With wrapper script for complex workflows
helixforge parallel tasks --chunk-plan chunks.json \
    --wrapper wrapper.sh \
    --wrapper-setup 'module load python/3.10' \
    --wrapper-setup 'conda activate helixforge' \
    --command 'helixforge refine -p predictions.h5 -g predictions.gff3 --genome genome.fa --rnaseq-bam rnaseq.bam --region ${SEQID}:${START}-${END} -o ${OUTPUT_DIR}/${CHUNK_ID}.gff3' \
    --output tasks.txt
```

Available placeholders in command templates:
- `{chunk_id}` - Unique chunk identifier
- `{seqid}` - Scaffold/chromosome name
- `{start}` - Start coordinate (1-based, for CLI --region flags)
- `{end}` - End coordinate (1-based inclusive, for CLI --region flags)
- `{start_0}` - Start coordinate (0-based, for internal use)
- `{end_0}` - End coordinate (0-based exclusive, for internal use)
- `{size}` - Chunk size in bases
- `{output_dir}` - Output directory (if --output-dir provided)

Note: `{start}` and `{end}` output 1-based coordinates suitable for CLI commands like `--region chr1:1000-2000`. The internal 0-based half-open coordinates are converted automatically.

### Execute with HyperShell

[HyperShell](https://hypershell.readthedocs.io/) is the recommended tool for HPC execution:

```bash
# Local execution with 32 parallel workers
hs cluster tasks.txt --num-tasks 32

# With timeout per task (in seconds)
hs cluster tasks.txt --num-tasks 32 --task-timeout 3600

# Submit tasks to database for distributed execution
hs submit tasks.txt
```

### Execute with GNU Parallel

```bash
# Local execution
parallel -j 32 < tasks.txt

# With progress bar
parallel -j 32 --bar < tasks.txt

# With job log
parallel -j 32 --joblog parallel.log < tasks.txt
```

### Aggregate Results

```bash
# Merge TSV files
helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern '*.tsv' \
    -o combined.tsv \
    --type merge_tsv

# Merge GFF3 files
helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern '*.gff3' \
    -o combined.gff3 \
    --type merge_gff

# Simple concatenation
helixforge parallel aggregate \
    --input-dir outputs/ \
    --pattern '*.txt' \
    -o combined.txt \
    --type concat
```

### Get Parameter Suggestions

```bash
helixforge parallel suggest \
    --genome genome.fa \
    --memory 64 \
    --workers 16
```

### Generate Example SLURM Script

For users who prefer direct SLURM submission:

```bash
# HyperShell-based script
helixforge parallel example-sbatch -o run_job.sbatch --executor hypershell

# GNU Parallel-based script
helixforge parallel example-sbatch -o run_job.sbatch --executor parallel
```

## Python API

### Task File Generation

```python
from helixforge.parallel import (
    GenomeChunker,
    ChunkStrategy,
    ChunkPlan,
    TaskGenerator,
    generate_hypershell_command,
)
from helixforge.io.fasta import GenomeAccessor

# Load genome and create chunk plan
genome = GenomeAccessor("genome.fa")
chunker = GenomeChunker(genome)
plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)

# Generate task file (include all required options in command template)
gen = TaskGenerator(plan)
task_file = gen.generate(
    command_template="helixforge confidence -p predictions.h5 -g predictions.gff3 --genome genome.fa --chunk-id {chunk_id} --region {seqid}:{start}-{end} -o outputs/{chunk_id}.tsv",
    output_path="tasks.txt",
    output_dir="outputs/",
)

print(f"Generated {task_file.n_tasks} tasks")
print(f"Preview: {task_file.preview(3)}")

# Get HyperShell command
cmd = generate_hypershell_command("tasks.txt", parallelism=32)
print(f"Execute with: {cmd}")
```

### Local Parallel Execution

```python
from helixforge.parallel import (
    GenomeChunker,
    ChunkStrategy,
    ParallelExecutor,
    ExecutorBackend,
)

# Create chunk plan
chunker = GenomeChunker(genome)
plan = chunker.create_plan(ChunkStrategy.BY_SCAFFOLD)

# Define processing function
def process_chunk(chunk):
    # Your processing logic here
    return {"chunk_id": chunk.chunk_id, "genes_processed": 10}

# Execute in parallel
executor = ParallelExecutor(
    n_workers=8,
    backend=ExecutorBackend.PROCESSES,
)
results, stats = executor.map_chunks(process_chunk, list(plan.chunks))

print(f"Processed {stats.successful}/{stats.total_tasks} chunks")
print(f"Total time: {stats.total_duration:.1f}s")
```

### With Memory Monitoring

```python
executor = ParallelExecutor(
    n_workers=4,
    backend=ExecutorBackend.PROCESSES,
    max_memory_gb=32.0,  # Will log warnings if exceeded
)

# Track progress
def on_progress(completed, total, chunk_id):
    print(f"Progress: {completed}/{total} - {chunk_id}")

executor.progress_callback = on_progress
results, stats = executor.map_chunks(process_chunk, chunks)

if stats.peak_memory_mb:
    print(f"Peak memory: {stats.peak_memory_mb:.1f} MB")
```

### Wrapper Script Generation

For complex workflows with multiple commands per chunk:

```python
from helixforge.parallel import TaskGenerator

# Generate wrapper script (include all required options)
TaskGenerator.generate_wrapper_script(
    output_path="wrapper.sh",
    chunk_plan_path="chunks.json",
    commands=[
        "helixforge refine -p predictions.h5 -g predictions.gff3 --genome genome.fa --rnaseq-bam liver.bam,brain.bam,heart.bam --min-tissues 2 --region ${SEQID}:${START}-${END} -o ${OUTPUT_DIR}/${CHUNK_ID}.gff3 -r ${OUTPUT_DIR}/${CHUNK_ID}.tsv",
    ],
    setup_commands=[
        "module load python/3.10",
        "conda activate helixforge",
    ],
    output_dir="outputs/",
)

# Generate task file that calls wrapper
gen = TaskGenerator(plan)
task_file = gen.generate_with_wrapper(
    output_path="tasks.txt",
    wrapper_script="wrapper.sh",
    chunk_plan_output="chunks.json",
)
```

### SLURM Environment Detection

For scripts that need to detect SLURM context:

```python
from helixforge.parallel import (
    is_slurm_job,
    get_slurm_task_id,
    get_slurm_resources,
    get_chunk_for_task,
)

if is_slurm_job():
    # Get allocated resources
    resources = get_slurm_resources()
    print(f"CPUs: {resources['cpus']}, Memory: {resources['memory_mb']} MB")

    # Get chunk for this array task
    task_id = get_slurm_task_id()
    if task_id is not None:
        chunk = get_chunk_for_task("chunks.json", task_id)
        print(f"Processing chunk {chunk.chunk_id}")
        print(f"Region: {chunk.seqid}:{chunk.start}-{chunk.end}")
```

## Memory Monitoring

### Track Memory Usage

```python
from helixforge.parallel import track_memory, get_memory_stats

# One-time check
stats = get_memory_stats()
print(f"Current: {stats.current_mb:.1f} MB")
print(f"Peak: {stats.peak_mb:.1f} MB")
print(f"Available: {stats.available_mb:.1f} MB")

# Context manager for tracking
with track_memory() as get_stats:
    # Do some work
    large_array = [0] * 10_000_000

    stats = get_stats()
    print(f"Peak during processing: {stats.peak_mb:.1f} MB")
```

### Background Monitoring

```python
from helixforge.parallel import MemoryMonitor

monitor = MemoryMonitor(
    warning_threshold_percent=80.0,
    critical_threshold_percent=95.0,
    sample_interval_seconds=1.0,
)

monitor.start()

# Do work...

is_ok, message = monitor.check()
if not is_ok:
    print(f"Memory issue: {message}")

final_stats = monitor.stop()
print(f"Peak memory: {final_stats.peak_mb:.1f} MB")
```

## Plan Serialization

Chunk plans can be saved and loaded:

```python
# Save plan
plan.save("chunks.json")

# Load plan
from helixforge.parallel import ChunkPlan
plan = ChunkPlan.load("chunks.json")

# Get chunk by ID
chunk = plan.get_chunk("chunk_0005")

# Get chunk by index
chunk = plan.get_chunk_by_index(5)

# Find chunk containing a region
chunk = plan.find_chunk_for_region("chr1", 1000000, 1500000)
```

## Result Aggregation

### Merge Overlapping Results

When chunks may overlap (e.g., to avoid splitting genes):

```python
from helixforge.parallel import merge_overlapping_results

# results is list of (chunk, genes) tuples
merged_genes = merge_overlapping_results(results)
```

## Best Practices

### For HPC Clusters

1. **Use HyperShell**: Provides scheduler-agnostic parallelization
2. **Generate more chunks than workers**: Enables better load balancing
3. **Use consistent chunk sizes**: BY_SIZE strategy for predictable memory usage
4. **Test locally first**: Run a few chunks to estimate memory/time

### For Local Execution

1. **Start conservative**: Begin with fewer workers and increase
2. **Monitor memory**: Use `max_memory_gb` to prevent OOM
3. **Use processes for CPU-bound**: `ExecutorBackend.PROCESSES`
4. **Use threads for I/O-bound**: `ExecutorBackend.THREADS`

### For Large Genomes (>1 Gb)

1. **Use BY_SIZE strategy**: Predictable memory usage
2. **Limit chunk size to 50-100 Mb**: Balances memory and overhead
3. **Request 16-32 GB per task**: Allows for peak usage spikes
4. **Generate more chunks than workers**: Better load balancing

## Choosing Parameters

### Optimal Worker Count

```python
from helixforge.parallel import get_optimal_workers

# Automatically determine based on CPU and memory
workers = get_optimal_workers(
    max_workers=16,
    memory_per_worker_mb=2000,  # 2 GB per worker
)
```

### Suggested Chunk Parameters

```python
from helixforge.parallel import suggest_chunk_parameters

suggestion = suggest_chunk_parameters(
    genome_size=2_300_000_000,  # 2.3 Gb (maize)
    n_genes=40_000,
    available_memory_gb=64,
    n_workers=16,
)

print(f"Strategy: {suggestion['strategy']}")
print(f"Chunk size: {suggestion['chunk_size']}")
print(f"Rationale: {suggestion['rationale']}")
```

## Troubleshooting

### Out of Memory

1. Reduce `chunk_size` or increase `target_chunks`
2. Reduce `n_workers`
3. Switch to `BY_SIZE` strategy with smaller chunks
4. Increase SLURM `--mem` allocation

### Slow Processing

1. Increase `n_workers` (if memory allows)
2. Use `BY_GENES` strategy to balance workload
3. Check for I/O bottlenecks (use local scratch on clusters)

### Task Failures

```bash
# With HyperShell, check the failed tasks
hs status

# With GNU Parallel, check the job log
cat parallel.log | grep -v "^1"  # Show non-successful tasks

# Re-run failed tasks manually or regenerate task file
```
