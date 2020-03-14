#include <iostream>
#include <algorithm>
#include "kernels.h"

/* Parse a fasta line. Check if all characters are valid and cast to uppercase.
** Parameters: Amino acid sequence to check.
** Returns: Sanitized amino acid sequence.
*/
std::string valid_AA(std::string seq) {
  std::string sanitized;

  static const char amino_acids[] =
  { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };

  for (auto c : seq) {

    // First check if alphabetical
    if (!isalpha(c)) {
      throw runtime_error("ERROR: Character " + std::string(1, c) + " not a valid amino acid");
    }

    // Convert to uppercase if necessary
    if (islower(c)) {
      c = toupper(c);
    }

    // Check if in amino acid alphabet
    if (std::find(std::begin(amino_acids), std::end(amino_acids), c) < std::end(amino_acids)) {
      sanitized += c;
    } else {
      throw runtime_error("ERROR: Character " + std::string(1, c) + " not a valid amino acid");
    }
  }

  return sanitized;
}

/* Parse a valid single-fasta file.
** Throw runtime error if first line is not a header or if there is more than one header.
** Parameters: Name of fasta file.
** Returns: Sequence from fasta without headers or newlines.
*/
std::string parse_fasta(std::string filename) {
  std::string seq = "";

  // Open file
  std::ifstream input_file;
  input_file.open(filename);

  if (!input_file) {
    throw runtime_error("ERROR: Unable to open file " + filename + "\n");
  }

  // Read first line; should be fasta header
  std::string line;
  std::getline(input_file, line);

  if (line[0] != '>') {
    throw runtime_error("ERROR: Expected fasta header, instead found:\n" + line + "\n");
  }

  // Read sequence
  while (std::getline(input_file, line)) {

    // Throw runtime error if there is another fasta header
    if (line[0] == '>') {
      cout << "Unexpected input: " <<  line << endl;
      throw runtime_error("ERROR: Multifasta files not supported.\n"
      "Please choose one fasta sequence from " + filename + " and save as single fasta.\n");
    }

    // Check sequence and sanitize
    try {
      line = valid_AA(line);
    }
    catch (const runtime_error& e) {
        cout << endl << e.what();
        throw runtime_error("ERROR: Problem parsing " + filename + "\n");
    }

    // Concatenate line to sequence
    seq += line;

  }

  return seq;
}

int main(int argc, char *argv[])
{
  // Parse arguments
  std::string usage = "USAGE: \t./popcorn seq1.fa seq2.fa [ -b beta ] [ -k maximum k-mer length ]";
  double beta = 0.01;
  int max_substring_length = 0;

  try {
    // 2 fasta files are required arguments
    if (argc < 3) {
      throw runtime_error("ERROR: Popcorn requires 2 single fasta files.\n");
    }

    // Support optional arguments -b {beta} -k {max substring length}
    else if (argc == 5 || argc == 7) {

      for (int i = 3; i < argc; i += 2) {
        std::stringstream argument(argv[i+1]);

        if (std::string(argv[i]) == "-b") {
          argument >> beta;

          if (argument.fail() || beta <= 0) {
            throw runtime_error("ERROR: Beta must be a positive decimal number.\n");
          }
        }

        else if (std::string(argv[i]) == "-k") {
          argument >> max_substring_length;

          if (argument.fail() || max_substring_length <= 0) {
            throw runtime_error("ERROR: K must be a strictly positive integer.\n");
          }
        }

        else {
          if (std::string(argv[i])[0] == '-') {
            throw runtime_error("ERROR: Unrecognized option " + std::string(argv[i]) + "\n");
          }
          else {
            throw runtime_error("ERROR: Unsupported positional argument " + std::string(argv[i]) + "\n");
          }
        }
      }
    }

    // User is drunk and should go home
    else if (argc != 3) {
      throw runtime_error("ERROR: Unexpected number of arguments.\n");
    }
  }
  catch (const runtime_error& e) {
      cout << endl << e.what();
      cout << usage << endl << endl;
      return 1;
  }


  // Parse sequences
  std::string seq1;
  std::string seq2;

  try {
    seq1 = parse_fasta(argv[1]);
    seq2 = parse_fasta(argv[2]);
  }
  catch (const runtime_error& e) {
      cout << endl << e.what();
      cout << usage << endl;
      return 1;
  }


  // Echo extra options to user
  cerr << "Using beta = " << beta << endl;
  if (max_substring_length != 0) {
    cerr << "Limiting algorithm to substrings of max length " << max_substring_length << endl;
  }


  // Calculate kernels
  double calculated_distance;
  unordered_map <string, double> local_K1 = read_blosum_build_kernel(beta);

  calculated_distance = protein_distance(seq1, seq2, local_K1);
  cout << calculated_distance << endl;

  return 0;
}
