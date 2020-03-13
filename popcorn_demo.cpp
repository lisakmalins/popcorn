#include <iostream>
#include "popcorn.cpp"

int main()
{
    double beta;
    double calculated_distance;


    /* These are the 3 sequences given with assignment */
    string s1 = "EFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK";
    string s2 = "VPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEG";
    string s3 = "EFDVILKAAGANKVAVIKAVRGATGLALKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK";

    cout << "Please type the value of beta: ";

    try
    {
        cin >> beta;
        if (cin.fail()) throw runtime_error("Input is not an a correct data type! Please try again!\n");
        if (beta < 0) throw runtime_error("Input is not a positive real number! Please try again!\n");

    }
    catch (const runtime_error& e)
    {
        cout << endl << e.what();
        return 1;
    }


    unordered_map <string, double> local_K1 = read_blosum_build_kernel(beta);

    /* Results from calculating the distance using all possible substrings */
    /* Verify results are the same regardless of order */
    calculated_distance= protein_distance(s1, s2, local_K1, 0);
    cout << "Distance between sequence 1 and sequence 2 = " << calculated_distance << endl;
    calculated_distance= protein_distance(s2, s1, local_K1, 0);
    cout << "Distance between sequence 2 and sequence 1 = " << calculated_distance << endl;
    cout << endl;

    calculated_distance = protein_distance(s1, s3, local_K1);
    cout << "Distance between sequence 1 and sequence 3 = " << calculated_distance << endl;
    calculated_distance = protein_distance(s3, s1, local_K1);
    cout << "Distance between sequence 3 and sequence 1 = " << calculated_distance << endl;
    cout << endl;

    calculated_distance = protein_distance(s2, s3, local_K1);
    cout << "Distance between sequence 2 and sequence 3 = " << calculated_distance << endl;
    calculated_distance = protein_distance(s3, s2, local_K1);
    cout << "Distance between sequence 3 and sequence 2 = " << calculated_distance << endl;
    cout << endl;


    /* Verify distance is 0 for identities */
    calculated_distance = protein_distance(s1, s1, local_K1);
    cout << "Distance between sequence 1 and sequence 1 = " << calculated_distance << endl;
    calculated_distance = protein_distance(s2, s2, local_K1);
    cout << "Distance between sequence 2 and sequence 2 = " << calculated_distance << endl;
    calculated_distance = protein_distance(s3, s3, local_K1);
    cout << "Distance between sequence 3 and sequence 3 = " << calculated_distance << endl;
    cout << endl;


    return 0;
}
