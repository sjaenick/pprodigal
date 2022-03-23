PProdigal: Parallelized gene prediction based on Prodigal.

This is just a small wrapper around the prodigal gene prediction program
that splits input into chunks and processes them im parallel, since prodigal
does not support multithreading by itself. The wrapper supports all command
line parameters accepted by prodigal itself, with two additional parameters
that control the parallelization:

  -T TASKS, --tasks TASKS
                        number of prodigal processes to start in parallel (default: 20)
  -C CHUNKSIZE, --chunksize CHUNKSIZE
                        number of input sequences to process within a chunk (default: 2000)

Due to prodigal's self-training phase, chunks should be chosen sufficiently
large in order to avoid suboptimal results.

