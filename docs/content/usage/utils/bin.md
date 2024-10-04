---
title: bin
weight: 0
---

## Usage

```plain

Bin reads:
$ lexicmap utils bin -r report.tsv -o /tmp/outputs ~/input.fastx

Bin reads but do not make seperate output for reads uniquely binned.
$ lexicmap utils bin -r report.tsv -u -o /tmp/outputs ~/input.fastx

Distribute fasta or fastq files into separate files based on the output of a previously generated report from search.

Previously mapped sequences can go be placed into file in which the sequence is put into a file for each genome it maps too, or
into a folder containing each read it is mapped too along with a separate folder containing only reads that uniquely map to one genome.

Output files can be gzipped depending on the input format. e.g. if the input file is gzipped, output files will be gzipped too.

Output Format:

Output Directory
  ├── All
  |   ├── Genome1.sequences
  |   ├── Genome2.sequences
  |
  ├──Unique
      ├── Genome1.unique.sequences
      ├── Genome2.unique.sequences


or if unique sequences are not binned.

Output Directory
   ├── Genome1.sequences
   ├── Genome2.sequences

Usage:
  lexicmap utils bin [flags]

Flags:
  -u, --bin-unique-reads     ► Create separate reads from unique source into a separate folder.
                             (default true)
  -b, --buffer-size string   ► Size of buffer, supported unit: K, M, G. You need increase the value
                             when "bufio.Scanner: token too long" error reported (default "20M")
  -h, --help                 help for bin
  -o, --out-dir string       ► Output directory, supports and recommends a ".gz" suffix ("-" for
                             stdout). (default "binned")
  -r, --report string        ► The generated output of the lexicmap

Global Flags:
  -X, --infile-list string   ► File of input file list (one file per line). If given, they are
                             appended to files from CLI arguments.
      --log string           ► Log file.
      --quiet                ► Do not print any verbose information. But you can write them to a file
                             with --log.
  -j, --threads int          ► Number of CPU cores to use. By default, it uses all available cores.
                             (default 6)

```