#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unordered_map>
#include <map>
#include <set>
#include <bits/stdc++.h>

using namespace std;

/*************************************************************************************
                    GLOBAL VARIABLES
*************************************************************************************/

/* this string stores the row-column codes for the blosum and k1 matrices */

std::string amino_acids = ("ARNDCQEGHILKMFPSTWYV");

/*this stores K1 matrix */
unordered_map <string, double> K1;

/*************************************************************************************
                    PROCEDURES
*************************************************************************************/

/* find a subsequence of variable length */

string subsequence_generator(string str, int subsequence, int length)
{
    string local_subsequence = "";
    for (int j = 0; j < length; j++)

        /* check if jth bit in subsequence is 1 and if it is, include it in the subsequence*/
        if (subsequence & (1 << j))

            local_subsequence += str[j];

    return local_subsequence;
}

/* take all-size generated subsequences and store them in a map */

vector<string> store_subsequences(string str){

    /* store subsequences in ascending order by length */
    map<int, set<string> > possible_subsequences;

    /* subsequences stored in sorted order for return */
    vector<string> returned_sunbsequence;

    /* maximum length of the subsequence; make it 10 now but change it later if needed */
    int max_length = str.length();

    /* limit sequence length to 10 */
    int limit = pow(2, max_length);

    /* start at i=2 because you want K2 to receive sequences of minimum 2 letters */
    for (int i = 2; i <= limit - 1; i++) {

        /* subsequence for binary pattern i */

        string subsequence = subsequence_generator(str, i, max_length);

        /* store subsequences in a map before putting them in order in a vector */
        possible_subsequences[subsequence.length()].insert(subsequence);
    }
    int i = 0;
    int j = 0;
    for (auto it : possible_subsequences) {
        // it.first is length of Subsequence
        // it.second is set<string>
        cout << "Subsequences of length = "
             << it.first << " are:" << endl;

        for (auto ii : it.second){
             // ii is iterator of type set<string>
             if(j > 0){

                returned_sunbsequence.push_back(ii);
                cout << returned_sunbsequence.at(i) << " ";
                i++;
             }

        }
        j++;


        cout << endl;
    }
    return returned_sunbsequence;
}


/* the function below reads in the BLOSUM62-2 matrix and
** raises each of is elements to the power of beta to create
** kernel K^1; kernel K1 can be accessed on this way: K1[index] where index is matrix row+column index
** i.e) index = AA or AR etc so K1[index] <=> K1["AA"] or K1["AR"]
** parameters: beta
** returns: K1
*/

unordered_map <string,double> read_blosum_build_kernel(double beta)
{
    std::ifstream file("BLOSUM62.txt");

    std::string   line = "";
    double value = 0.0;
    std::string index; //this will store the index of the matrix row-column format

    // Read one line at a time into the variable line:
    int i = 0, j = 0;
    while(getline(file, line))
    {
        for(int k = 0; k < 20; k++)
        {
            /*reset the column number for each new row */
            std::stringstream  lineStream(line);

            j = 0;

            /* read each value from a line in blosum62 */
            while(lineStream >> value)
            {

                index.push_back(amino_acids[k]);
                index.push_back(amino_acids[j]);
                /* if failing to read, throw and error */
                if(lineStream.fail())
                {

                    std::cout << "lineStream failed" << std::endl;
                }
                K1[index] = pow(value, beta); //raise each element in K to beta

                /*below is the test case to see if values get mapped correctly */
                /*if(k < 1) {
                    cout << K1[index] << " and the index is " << index << "  ";
                }*/

                index.erase();
                j++; //increment column number


            }

        }

    }

    return K1;
}

/* this function computes kernel K2
** parameters: two protein sequences (1st sequence is smaller than or equals 2nd sequence), beta
** returns: K2 of the 2 sequences
*/

double compute_K2(string seq1, string seq2, double beta) {

    double K2 = 1.0;
    unordered_map <string, double> local_k1 = read_blosum_build_kernel(beta);
    string index;

    for(int i = 0; i < seq1.length(); i++)
    {
        for(int j = 0; j < seq2.length(); j++)
        {
            index.push_back(seq1[i]);
            index.push_back(seq2[j]);
            //cout << "index = " << index << "  && k1[index] = " << local_k1[index] << "  ";
            K2 *= local_k1[index];
            index.erase();
        }

    }
return K2;
}

double compute_K3(string seq1, string seq2, double beta) {

    double K3 = 0.0;
    string index;

    unordered_map <string, double> local_k1 = read_blosum_build_kernel(beta);
    vector<string> sequence1;
    vector<string> sequence2;

    /* K3 for 1 letter sequence combination */

    for(int i = 0; i < seq1.length(); i++)
    {
        for(int j = 0; j < seq2.length(); j++)
        {
            index.push_back(seq1[i]);
            index.push_back(seq2[j]);
            //cout << "index = " << index << "  && k1[index] = " << local_k1[index] << "  ";
            K3 += local_k1[index];
            index.erase();
        }

    }
   cout << endl;
    /* K3 for 2+ letter sequence combination */

    sequence1 = store_subsequences(seq1);
    sequence2 = store_subsequences(seq2);
    for(auto seq1 = sequence1.begin(); seq1 != sequence1.end(); ++seq1){
        for(auto seq2 = sequence2.begin(); seq2 != sequence2.end(); ++seq2){

            if((*seq1).length() == (*seq2).length())
            {
               // cout << "seq 1 = " << *seq1 << " seq 2 = " << *seq2 << endl;
                K3 += compute_K2(*seq1, *seq2, beta);

            }

        }
    }
    return K3;

}




int main()
{

    unordered_map <std::string, double> blosum62;
    double beta;
    double K2;
    double K3;
    std::cout << "Please type the value of beta: ";
    std::cin >> beta;

    string s1 = "EFDVI";
    string s2 = "VPCSD";

    //K2 = compute_K2(s1,s2, beta);
    //cout << "K2 = " << K2 << endl;

    K3 = compute_K3(s1,s2,beta);
    cout << "K3 = " << K3 << endl;
     //read_blosum_build_kernel(beta);

   // Pick starting point

  // It computes all the subsequence of an string

   // store_subsequences(s2);








    return 0;
}
