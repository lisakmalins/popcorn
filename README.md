# Popcorn
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
g++ -std=c++11 popcorn.cpp -o run_popcorn
```

The executable will be named run_popcorn.
```
./run_popcorn
```

## Collaborators
_(In alphabetical order)_
- Abigaela Boroica (UC Davis Computer Science, Class of 2021)
- Lisa Malins (UC Davis Biotechnology, Class of 2020)

## References
- [Shen WJ, Wong HS, Xiao QW, Guo X, Smale S. Towards a mathematical foundation of immunology and amino acid chains. arXiv preprint arXiv:1205.6031. 2012 May 28.](https://arxiv.org/abs/1205.6031)
