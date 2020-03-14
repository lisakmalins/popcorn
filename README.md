# Popcorn
![Popcorn logo](Images/popcorn_logo_transparent.png)
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
g++ -std=c++11 popcorn_demo.cpp -o popcorn_demo
g++ -std=c++11 test_max_substrings.cpp -o test_max_substrings
```

The executables can be run from the command-line.
```
./popcorn_demo
./test_max_substrings
```

`popcorn_demo` is a simple demonstration with three amino acid sequences. The program shows that the distance calculation is equivalent regardless of the order of the sequences. It also demonstrates that the distance is 0 when the sequences are the same.

`test_max_substrings` demonstrates the effect of limiting the maximum substring length. The output is tab-delimited; each line lists the maximum substring length and the distances between all three sequence pairs.



## Collaborators
- Abigaela Boroica (UC Davis Computer Science, Class of 2021)
- Lisa Malins (UC Davis Biotechnology, Class of 2020)

## References
- [Shen WJ, Wong HS, Xiao QW, Guo X, Smale S. Towards a mathematical foundation of immunology and amino acid chains. arXiv preprint arXiv:1205.6031. 2012 May 28.](https://arxiv.org/abs/1205.6031)
