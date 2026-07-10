# `helixforge parallel example-sbatch`

Write a copy-paste SBATCH wrapper that runs a ``tasks`` file on one node.

Convenience starting point for the scheduler-agnostic path: it runs the whole
``tasks.txt`` (from ``parallel tasks``) inside a single allocation, fanning
out across ``$SLURM_CPUS_PER_TASK`` cores via HyperShell or GNU Parallel. Edit
the partition / time / resources for your cluster. For a true one-task-per-
chunk Slurm *array* driven off a plan, use ``reconcile --hpc slurm`` instead.

## Options

| Option | Description |
|---|---|
| `-o`, `--output` | Where to write the example SBATCH script. (required) |
| `--executor` | Run the task file with HyperShell or GNU Parallel. (choices: hypershell/parallel; default: `hypershell`) |

## Examples

```bash
  # HyperShell wrapper (default)
  helixforge parallel example-sbatch -o run.sbatch

  # GNU Parallel wrapper
  helixforge parallel example-sbatch -o run.sbatch --executor parallel
```
