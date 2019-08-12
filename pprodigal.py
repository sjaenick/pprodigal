#!/usr/bin/env python3

import argparse
import tempfile
from concurrent.futures import ThreadPoolExecutor
import subprocess
import threading

def run_prodigal(opts, workDir, currentId, chunkFile):
    cmd = ["prodigal", "-q", "-i", chunkFile.name ]

    if opts.proteins:
        cmd.append("-a")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".faa")

    if opts.closed:
        cmd.append("-c")

    if opts.nucl:
        cmd.append("-d")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".fna")

    if opts.format:
        cmd.append("-f")
        cmd.append(opts.format)

    if opts.gencode:
        cmd.append("-g")
        cmd.append(str(opts.gencode))

    if opts.mask:
        cmd.append("-m")

    if opts.nosd:
        cmd.append("-n")

    if opts.procedure:
        cmd.append("-p")
        cmd.append(opts.procedure)

    if opts.scorefile:
        cmd.append("-s")
        cmd.append(workDir.name + "/chunk" + str(currentId) + ".score")

    cmd.append("-o")
    cmd.append(workDir.name + "/chunk" + str(currentId) + ".out")

    #print(str(cmd))
    subprocess.run(cmd, shell=False, check=True)

def append_file(file, targetFile):
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            trgt.write(input.read())


def main():
    argp=argparse.ArgumentParser(description='Parallel Prodigal gene prediction')
    argp.add_argument('-a', "--proteins", type=str, help="Write protein translations to the selected file.")
    argp.add_argument('-c', "--closed", action="store_true", help="Closed ends.  Do not allow genes to run off edges.")
    argp.add_argument('-d', "--nucl", type=str, help="Write nucleotide sequences of genes to the selected file.")
    argp.add_argument('-f', "--format", type=str, help="Select output format (gbk, gff, or sco).  Default is gbk.")
    argp.add_argument('-g', "--gencode", type=int, help="Specify a translation table to use (default 11).")
    argp.add_argument('-i', "--input", type=str, help="Specify FASTA/Genbank input file (default reads from stdin).")
    argp.add_argument('-m', "--mask", action="store_true", help="Treat runs of N as masked sequence; don't build genes across them.")
    argp.add_argument('-n', "--nosd", action="store_true", help="Bypass Shine-Dalgarno trainer and force a full motif scan.")
    argp.add_argument('-o', "--output", type=str, help="Specify output file (default writes to stdout).")
    argp.add_argument('-p', "--procedure", type=str, help="Select procedure (single or meta).  Default is single.")
    argp.add_argument('-s', "--scorefile", type=str, help="Write all potential genes (with scores) to the selected file.")
    argp.add_argument('-T', "--threads", type=int, help="number of threads")
    opts = argp.parse_args()

    if opts.threads and opts.threads < 1:
        raise ValueError

    seqsPerChunk = 5000
    seqCnt = 0
    currentChunk = 1

    workDir = tempfile.TemporaryDirectory()
    executor = ThreadPoolExecutor(max_workers=opts.threads)
    currentFile = open(workDir.name + "/chunk" + str(currentChunk), 'w')

    queryFile = None
    if not opts.input:
        queryFile = "/dev/fd/0"
    else:
        queryFile = opts.input

    with open(queryFile, 'r') as fasta:
        for line in fasta:

            if line[0] == '>' and seqCnt == seqsPerChunk:
                currentFile.close()
                executor.submit(run_prodigal, opts, workDir, currentChunk, currentFile)
                currentFile = None
                seqCnt = 0
                currentChunk = currentChunk + 1

            if currentFile is None:
                currentFile = open(workDir.name + "/chunk" + str(currentChunk), 'w')

            currentFile.write(line)

            if line[0] == '>':
                seqCnt = seqCnt + 1

    if seqCnt > 0:
        currentFile.close()
        executor.submit(run_prodigal, opts, workDir, currentChunk, currentFile)

    # await completion of tasks
    executor.shutdown(wait=True)

    # collect output
    #
    proteinFile = opts.proteins
    nuclFile = opts.nucl
    outFile = opts.output
    scoreFile = opts.scorefile

    for cur in range(1, currentChunk + 1):
        if proteinFile:
            append_file(workDir.name + "/chunk" + str(cur) + ".faa", proteinFile)
        if nuclFile:
            append_file(workDir.name + "/chunk" + str(cur) + ".fna", nuclFile)
        if scoreFile:
            append_file(workDir.name + "/chunk" + str(cur) + ".score", scoreFile)

        if outFile:
            append_file(workDir.name + "/chunk" + str(cur) + ".out", outFile)
        else:
            with open(workDir.name + "/chunk" + str(cur) + ".out") as out:
                print(out.read())


if __name__== "__main__":
    main()
