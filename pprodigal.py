#!/usr/bin/env python3

#### Jianshu Zhao jianshu.zhao@gatech.edu. Modified according to a online version and fixed a lot of bugs
#### gff format output now is correct. Feel free to contact me if you find any bugs.
#### usage: python pprodigal.py -i contig.fasta -T 24 -f gff -o contig.gff -p meta


import argparse
import tempfile
from concurrent.futures import ThreadPoolExecutor
from shutil import which
import subprocess
import threading
import sys
import re


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

def append_fasta_file(file, startNum, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                if line[0] == '>':
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4)
                trgt.write(line)
    return startNum

def append_gff_file(file, startNum, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                if line[0] != '#' and "ID=" in line:
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4) + "\n"
                trgt.write(line)
    return startNum


def append_gbk_file(file, startNum, targetFile):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            for line in input:
                if line[0] == ' ' and "ID=" in line:
                    match = re.match(pattern, line)
                    if match and match.group(3) == "1":
                        startNum = startNum + 1
                    line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4)
                trgt.write(line)
    return startNum


def append_raw_file(file, targetFile):
    with open(targetFile, "a") as trgt:
        with open(file, "r") as input:
            trgt.write(input.read())


def print_gff_file(file, startNum):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(file, "r") as input:
        for line in input:
            if line[0] != '#' and "ID=" in line:
                match = re.match(pattern, line)
                if match and match.group(3) == "1":
                    startNum = startNum + 1
                line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4)
            print(line)
    return startNum


def print_gbk_file(file, startNum):
    pattern = re.compile(r"(.*ID=)(\d+)_(\d+)(.*)")
    with open(file, "r") as input:
        for line in input:
            if line[0] == ' ' and "ID=" in line:
                match = re.match(pattern, line)
                if match and match.group(3) == "1":
                    startNum = startNum + 1
                line = match.group(1) + str(startNum) + "_" + match.group(3) + match.group(4)
            print(line)
    return startNum


def print_raw_file(file):
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
    argp.add_argument('-T', "--tasks", type=int, help="number of prodigal processes to start in parallel (default: 20)")
    argp.add_argument('-C', "--chunksize", type=int, help="number of input sequences to process within a chunk (default: 2000)")
    opts = argp.parse_args()

    ### check if prodigal installed
    rc = subprocess.call(['which', 'prodigal'])
    if rc == 0:
        print ("wget installed!")
    else:
        print ("wget missing in path!")

    ### check if tasks less than 1
    tasks = 20
    if opts.tasks is not None:
        if opts.tasks < 1:
            raise ValueError
        tasks = opts.tasks

    if opts.chunksize and opts.chunksize < 1:
        raise ValueError

    seqsPerChunk = 2000
    if opts.chunksize is not None:
        if opts.chunksize < 1:
            raise ValueError
        seqsPerChunk = opts.chunksize

    seqCnt = 0
    currentChunk = 1

    workDir = tempfile.TemporaryDirectory()
    executor = ThreadPoolExecutor(max_workers=tasks)
    currentFile = open(workDir.name + "/chunk" + str(currentChunk), 'w')

    queryFile = None
    if not opts.input:
        queryFile = "/dev/fd/0"
        if sys.stdin.isatty():
            print("Cannot read sequences from STDIN.")
            exit(1)
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

    protIdStart = 0
    nuclIdStart = 0
    gffIdStart = 0
    gbkIdStart = 0
    for cur in range(1, currentChunk + 1):
        if proteinFile:
            protIdStart = append_fasta_file(workDir.name + "/chunk" + str(cur) + ".faa", protIdStart, proteinFile)
        if nuclFile:
            nuclIdStart = append_fasta_file(workDir.name + "/chunk" + str(cur) + ".fna", nuclIdStart, nuclFile)
        if scoreFile:
            append_raw_file(workDir.name + "/chunk" + str(cur) + ".score", scoreFile)

        if outFile:
            if opts.format == "gff":
                gffIdStart = append_gff_file(workDir.name + "/chunk" + str(cur) + ".out", gffIdStart, outFile)
            elif opts.format == "sco":
                append_raw_file(workDir.name + "/chunk" + str(cur) + ".out", outFile)
            else:
                gbkIdStart = append_gbk_file(workDir.name + "/chunk" + str(cur) + ".out", gbkIdStart, outFile)
        else:
            if opts.format == "gff":
                gffIdStart = print_gff_file(workDir.name + "/chunk" + str(cur) + ".out", gffIdStart)
            elif opts.format == "sco":
                print_raw_file(workDir.name + "/chunk" + str(cur) + ".out")
            else:
                gbkIdStart = print_gbk_file(workDir.name + "/chunk" + str(cur) + ".out", gbkIdStart)


if __name__== "__main__":
    main()
