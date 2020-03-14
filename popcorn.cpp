#include <iostream>
#include "kernels.h"


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

    // Concatenate line to sequence
    seq += line;

  }

  return seq;
}

int main(int argc, char *argv[])
{
  std::string usage = "USAGE: \t./popcorn seq1.fa seq2.fa beta";
  if (argc != 3) {
    cout << "ERROR: Popcorn requires 2 single fasta files.\n" << usage << endl;
    return 1;
  }

  std::string seq1;
  std::string seq2;

  try {
    seq1 = parse_fasta(argv[1]);
    seq2 = parse_fasta(argv[2]);
  }
  catch (const runtime_error& e)
  {
      cout << endl << e.what();
      cout << usage << endl;
      return 1;
  }

  cout << "Here is seq1" << endl << seq1 << endl;
  cout << "Here is seq2" << endl << seq2 << endl;

  return 0;
}
