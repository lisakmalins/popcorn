# Popcorn

<br /><br />
<div align="center">
  <img src="Images/popcorn_logo_transparent.png" alt="Popcorn logo" width="450px" />
</div>
<br /><br />

A kernel-based, alignment-free distance calculator for amino acid sequences.

Inspired by Stephen Smale and colleagues' 2012 paper, "Towards a mathematical foundation of immunology and amino acid chains."

## Usage
Popcorn requires the g++ compiler from the [GNU Compiler Collection](https://gcc.gnu.org/). In addition, [GNU Make](https://www.gnu.org/software/make/) is recommended for easy compiling.

First, clone the repository from GitHub. (Or obtain and unpack a zip archive of Popcorn.)
```
git clone https://github.com/lisakmalins/popcorn.git
cd popcorn
```

To compile Popcorn, use the included Makefile.
```
make
```

If you do not have GNU Make installed, you may also compile popcorn directly from the command line.
```
g++ -std=c++11 popcorn.cpp -o popcorn
g++ -std=c++11 popcorn_demo.cpp -o popcorn_demo
g++ -std=c++11 test_max_substrings.cpp -o test_max_substrings
```

The executables can be run from the command-line.

### Main program and options
The `popcorn` executable can run on any two single-fasta files. If not specified, beta defaults to 0.01 and there is no limit on k-mer length.
```
# Run standard popcorn on seq1 and seq2
./popcorn data/seq1.fa data/seq2.fa
```

Beta can be specified with option `-b` and maximum k-mer length can be specified with option `-k`.
```
# Run popcorn on seq1 and seq2 with beta = 0.1 and max k-mer length = 20
./popcorn data/seq1.fa data/seq2.fa -b 0.1 -k 20
```

### Quickstart demos
This repository also includes two executables that demonstrate kernel functions without arguments.

```
./popcorn_demo
```

`popcorn_demo` is a simple demonstration of distance calculation between three amino acid sequences.

This program shows that the distance is equal regardless of the order of the sequences. It also shows that the distance is 0 when the sequences are identical.

```
./test_max_substrings
```

`test_max_substrings` demonstrates the effect of limiting the maximum substring length.

The output is tab-delimited; each line lists the maximum substring length and the distances between all three sequence pairs.

### Error handling
While parsing fasta files, Popcorn will print an error message and exit if it detects any of the following:

- Input file missing fasta header
- Input file has multiple fasta sequences
- Unexpected whitespace
- Numbers or symbols in sequence
- Nonstandard amino acids in sequence

Examples of invalid fasta files are included in `data/invalid_fasta/`.

```
# Testing error handling
./popcorn data/invalid_fasta/missingheader.fa data/seq1.fa
./popcorn data/invalid_fasta/multifasta.fa data/seq1.fa
./popcorn data/invalid_fasta/whitespace.fa data/seq1.fa
./popcorn data/invalid_fasta/nonstandard.fa data/seq1.fa
./popcorn data/invalid_fasta/garbage.fa data/seq1.fa
```

On the other hand, Popcorn will handle fasta sequences with uppercase or lowercase amino acids without errors.

The file `data/lowercase.fa` is a copy of `data/seq1.fa` with some letters lowercase to demonstrate.
```
# Same as distance between seq1 and seq3
./popcorn data/lowercase.fa data/seq3.fa
# Distance is 0 because protein is the same
./popcorn data/lowercase.fa data/seq1.fa
```


## Collaborators
- Abigaela Boroica (UC Davis Computer Science, Class of 2021)
- Lisa Malins (UC Davis Biotechnology, Class of 2020)

## References
- [Shen WJ, Wong HS, Xiao QW, Guo X, Smale S. Towards a mathematical foundation of immunology and amino acid chains. arXiv preprint arXiv:1205.6031. 2012 May 28.](https://arxiv.org/abs/1205.6031)
