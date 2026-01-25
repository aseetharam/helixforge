"""SLURM job templates for HelixForge parallel execution.

Available templates:
    - submit.sh.j2: Main SBATCH submission script
    - task.sh.j2: Per-task processing script
    - aggregate.sh.j2: Result aggregation script
    - splice_job.sh.j2: Complete splice refinement job
    - confidence_job.sh.j2: Complete confidence scoring job

Usage:
    >>> from helixforge.templates import render_template
    >>> script = render_template("slurm/submit.sh.j2",
    ...     job_name="my_job",
    ...     partition="normal",
    ...     time="4:00:00",
    ...     memory_gb=16,
    ...     n_chunks=10,
    ... )
"""
