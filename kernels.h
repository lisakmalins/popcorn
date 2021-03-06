#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <unordered_map>
#include <map>
#include <set>
#include <stdexcept>


using namespace std;

/*************************************************************************************
                    GLOBAL VARIABLES
*************************************************************************************/

/* this string stores the row-column codes for the blosum and k1 matrices */

string amino_acids = ("ARNDCQEGHILKMFPSTWYV");


/*************************************************************************************
                    PROCEDURES
*************************************************************************************/



/* this function will find all possible substring of different sizes
** parameters: protein sequence str, size of the sequence n, max length for substring size
** return: returns all possible substrings stored in a vector for future use
*/

vector<string> substring_generator(const string &str, int n, int substring_max_size)
{
    vector<string> possible_subsequences;
    string temp;

    /* choose the starting letter number of the subsequence*/

    for (int L = 1; L <= substring_max_size; L++)
    {
        /* choose the position in the sequence where you want to stop */
        for (int i = 0; i <= n-L; i++)
        {
            /* Find all characters between the starting letter and the end point*/
            int j = i + L - 1;

            for (int k = i; k <= j; k++)
            {

                /* use a temporary string to store all characters of substrings bigger than size 1 */
                temp.push_back(str[k]);

            }

            /* store all the characters of the substring in vector */
            possible_subsequences.push_back(temp);

            /* reset the temporary string to use in the next substring */
            temp.erase();

        }
    }


    /* return the vector containing all substrings in ascending order */
    return possible_subsequences;

}
/* the function below reads in the BLOSUM62-2 matrix and
** raises each of is elements to the power of beta to create
** kernel K^1; kernel K1 can be accessed on this way: K1[index] where index is matrix row+column index
** i.e) index = AA or AR etc so K1[index] <=> K1["AA"] or K1["AR"]
** parameters: beta
** returns: K1
*/

unordered_map <string, double> read_blosum_build_kernel(double beta)
{
    ifstream file("data/BLOSUM62.txt");
    unordered_map <string, double> K1;

    string line = "";
    double value = 0.0;
    string index; /* this will store the index of the matrix row-column format */

    /* Read one line at a time from blosum 62 into the variable line */
    int i = 0, j = 0;
    for(int k = 0; k < 20; k++)
    {
        getline(file, line);
        stringstream  lineStream(line);

        /*reset the column number for each new row */
        j = 0;

        /* read each value from a line in blosum62 */
        while(lineStream >> value)
        {

            /* attach all matrix index letters to create the index of K1 */
            index.push_back(amino_acids[k]);
            index.push_back(amino_acids[j]);

            /* if failing to read, throw an error */
            if(lineStream.fail())
            {

                cout << "lineStream failed" << endl;
            }

            /*raise each element in K to beta */
            K1.insert(std::make_pair<string,double> ((string)index,pow(value, beta)));

            index.erase();
            j++; /* increment column number */


        }

    }


    return K1;
}

/* this function computes kernel K2
** parameters: two protein sequences/subsequences, beta, kernel K1, max length for substring size
** returns: K2 of the 2 sequences
*/

double compute_K2(string &seq1, string &seq2, unordered_map <string, double> &K1)
{

    double K2 = 1.0;
    string index;

    for(int i = 0; i < seq1.length(); i++)
    {

        index.push_back(seq1[i]);
        index.push_back(seq2[i]);
        K2 *= K1.at(index);
        index.erase();


    }

    return K2;
}

/* This function computes kernel K3.
** Parameters: two protein sequences, beta, kernel K1, max length for substring size
** Substring max size = 0 (default argument) indicates no limit
*/

double compute_K3(string &seq1, string &seq2, unordered_map <string, double> &K1, int substring_max_size = 0)
{

    double K3 = 0.0;
    double val = 0.0;
    vector<string> sequence1;
    vector<string> sequence2;


    /* K3 for 1+ letter sequence combination */

    /* Choose if you want to generate substrings of unlimited or limited length */
    if(substring_max_size == 0)
    {
        sequence1 = substring_generator(seq1, seq1.length(), seq1.length());
        sequence2 = substring_generator(seq2, seq2.length(), seq2.length());
    }
    else
    {
        sequence1 = substring_generator(seq1, seq1.length(), substring_max_size);
        sequence2 = substring_generator(seq2, seq2.length(), substring_max_size);
    }


    /* compare each substring of the same length in sequence1 with one in sequence2 then compute K3 */
    for(auto s1 : sequence1)
    {
        for(auto s2 : sequence2)
        {
            /* only add K2 of sequences of the same length */
            if(s1.length() == s2.length())
            {
                K3 += compute_K2(s1, s2, K1);

            }
        }
    }
    return K3;

}

/* this function computes correlation kernel normalized from kernel K3
** parameters: two protein sequences, beta, kernel K1, max length for substring size
** returns: correlation K3
*/
double correlation_kernel_K3(string &seq1, string &seq2, unordered_map <string, double> &K1, int substring_max_size = 0)
{

    double correlation_K3 = 0.0;

    double K3_f_g = compute_K3(seq1, seq2, K1, substring_max_size);
    double K3_f_f = compute_K3(seq1, seq1, K1, substring_max_size);
    double K3_g_g = compute_K3(seq2, seq2, K1, substring_max_size);

    correlation_K3 = K3_f_g / sqrt(K3_f_f * K3_g_g);

    return correlation_K3;

}

/* this function computes the final metric distance between 2 protein sequences
** parameters: two protein sequences (1st sequence is smaller than or equals 2nd sequence), beta
** returns: the metric distance between 2 sequences of protein
*/
double protein_distance(string &seq1, string &seq2, unordered_map <string, double> &K1, int substring_max_size = 0)
{
    double distance = 0.0;
    double correlation_K3 = correlation_kernel_K3(seq1, seq2, K1, substring_max_size);

    distance = sqrt(2 * (1 - correlation_K3));

    return distance;
}
