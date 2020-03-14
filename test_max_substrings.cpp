#include <iostream>
#include "kernels.h"

int main()
{
  double beta = 0.01;
  double calculated_distance;

  unordered_map <string, double> local_K1 = read_blosum_build_kernel(beta);

  string s1 = "EFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK";
  string s2 = "VPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEG";
  string s3 = "EFDVILKAAGANKVAVIKAVRGATGLALKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK";

  cout << "max_length\t" << "seq1_seq2\t" << "seq1_seq3\t" << "seq2_seq3" << endl;

  for(int i=10; i<=100; i+=10)
  {
    cout << i << "\t";
    calculated_distance = protein_distance(s1, s2, local_K1, i);
    cout << calculated_distance << "\t";
    calculated_distance = protein_distance(s3, s1, local_K1, i);
    cout << calculated_distance << "\t";
    calculated_distance = protein_distance(s3, s2, local_K1, i);
    cout << calculated_distance << endl;
  }
  return 0;
}
