#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
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


/*************************************************************************************
                    PROCEDURES
*************************************************************************************/



/* this function will find all possible substring of different sizes
** parameters: protein sequence str, size of the sequence n
** return: returns all possible substrings stored in a vector for future use
*/

vector<string> substring_generator(const string &str, int n)
{
    vector<string> possible_subsequences;
    string temp;

    /* choose the starting letter numbers of the subsequence; only create lengths of 10 max */

    for (int L = 1; L <= 10; L++)
    {
        /* choose the position in the sequence where you want to stop */
        for (int i = 0; i <= n-L; i++)
        {
            /* Find all characters between the starting letter and the end point*/
            int j = i + L - 1;

            for (int k = i; k <= j; k++)  {

               //cout << str[k];
               /* use a temporary string to store all characters of substrings bigger than size 1 */
                temp.push_back(str[k]);

            }

            /* store all the characters of the substring in vector */
             possible_subsequences.push_back(temp);
            //cout << endl;

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
    std::ifstream file("BLOSUM62.txt");
    unordered_map <string, double> K1;

    std::string   line = "";
    double value = 0.0;
    std::string index; /* this will store the index of the matrix row-column format */

    /* Read one line at a time from blosum 62 into the variable line */
    int i = 0, j = 0;
    while(getline(file, line))
    {
        for(int k = 0; k < 20; k++)
        {

            std::stringstream  lineStream(line);

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

                    std::cout << "lineStream failed" << std::endl;
                }
                //cout << "value at " << index << " = " << value << " ";

                /*raise each element in K to beta */
                K1.insert(make_pair<string,double> ((string)index ,pow(value, beta)));

               //cout << "K1[" << index << "] = " << K1[index] << " | value =  " << value << " ";

                /*below is the test case to see if values get mapped correctly */
               /* if(k < 2) {
                    cout << K1[index] << " and the index is " << index << "  ";
                } */

                index.erase();
                j++; /* increment column number */


            }

        }

    }

    /* the following writes all values of K1 in a file called K1 - for the purpose of validating if
     * elements of blosum were raised correctly to the power of beta */

    ofstream os;

    os.open("K1.txt");
    for(int k = 0; k < 20; k++){
        for(int j = 0; j < 20; j++){
             index.push_back(amino_acids[k]);
             index.push_back(amino_acids[j]);
             os << "K1[" << index << "] = " << K1.at(index) << endl;
             index.erase();
        }

    }
    os.close();
return K1;
}

/* this function computes kernel K2
** parameters: two protein sequences (1st sequence is smaller than or equals 2nd sequence), beta
** returns: K2 of the 2 sequences
*/

double compute_K2(string &seq1, string &seq2, double beta) {

    double K2 = 1.0;
   unordered_map <string, double> local_K1 = read_blosum_build_kernel(beta);
    string index;

    for(int i = 0; i < seq1.length(); i++)
    {

            index.push_back(seq1[i]);
            index.push_back(seq2[i]);
            K2 *= local_K1.at(index);
            index.erase();


    }

return K2;
}

/* This function computer kernel K3.
** parameters: two protein sequences (the first one shorter than or equals the 2nd one), beta
** returns: kernel K3
*/

double compute_K3(string &seq1, string &seq2, double beta) {

    double K3 = 0.0;
    vector<string> sequence1;
    vector<string> sequence2;


    /* K3 for 1+ letter sequence combination */
     sequence1 = substring_generator(seq1, seq1.length());
     sequence2 = substring_generator(seq2, seq2.length());

    /* compare each substring of the same length in sequence1 with one in sequence2 then compute K3 */
    for(auto s1 : sequence1){
        for(auto s2 : sequence2){
            if(s1.length() == s2.length())
            {
                cout << "add to K3: K2(" << s1 << "," << s2 << ")" << endl;
                K3 += compute_K2(s1, s2, beta);

            }
        }
    }
    cout << endl;
    return K3;

}

/* this function computes correlation kernel normalized from kernel K3
** parameters: two protein sequences (1st sequence is smaller than or equals 2nd sequence), beta
** returns: correlation K3
*/
 double correlation_kernel_K3(string &seq1, string &seq2, double beta){

    double correlation_K3 = 0.0;
    double K3_f_g = compute_K3(seq1, seq2, beta);
    cout << "K3(f,g) = " << K3_f_g << endl;
    double K3_f_f = compute_K3(seq1, seq1, beta);
    cout << "K3(f,f) = " << K3_f_f << endl;
    double K3_g_g = compute_K3(seq2, seq2, beta);
    cout << "K3(g,g) = " << K3_g_g << endl;

    correlation_K3 = K3_f_g / sqrt(K3_f_f * K3_g_g);

    return correlation_K3;

 }

 /* this function computes the final metric distance between 2 protein sequences
** parameters: two protein sequences (1st sequence is smaller than or equals 2nd sequence), beta
** returns: the metric distance between 2 sequences of protein
*/
 double protein_distance(string &seq1, string &seq2, double beta) {
    double distance = 0.0;
    double correlation_K3 = correlation_kernel_K3(seq1, seq2, beta);

    distance = sqrt(2 * (1 - correlation_K3));

    return distance;
 }




int main()
{

    double beta;
    double calculated_distance;
    string index;

    std::cout << "Please type the value of beta: ";
    std::cin >> beta;

    /* these are the first 2 sequences given in project 4 */
   // string s1 = "EFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK";
    //string s2 = "VPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEG";

    /* these are the test sequences */
    string s1 = "AR";
    string s2 = "ARN";


    calculated_distance = protein_distance(s1, s2, beta);
    cout << endl << "final distance = " << calculated_distance << endl;


    return 0;
}
