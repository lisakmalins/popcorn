#include "/usr/include/c++/7/iostream"
#include "/usr/include/c++/7/string"
#include "/usr/include/c++/7/fstream"
#include "/usr/include/c++/7/sstream"
#include "/usr/include/c++/7/cmath"


//the function below reads in the BLOSUM62-2 matrix and
//raises each of is elements to the power of beta to create
// kernel K^1 

 void read_blosum_build_kernel(double K[][20], double beta)
{
    std::ifstream file("BLOSUM62.txt");

    std::string   line = "";
    // Read one line at a time into the variable line:
    int i = 0, j = 0;
    while(getline(file, line))
    {
         //cout << line << endl;
         /*reset the column number for each new row */
        std::stringstream  lineStream(line);

         j = 0;

         double value = 0.0; 

        /* read each value from a line in blosum62 */
        while(lineStream >> value)
        {
            //cout << value << endl;

            /* if failing to read, throw and error */
            if(lineStream.fail()){

                std::cout << "lineStream failed" << std::endl;
            }
            K[i][j] = pow(value, beta); //raise each element in K to beta
            j++; //increment column number


        }

        i++; //increment row number

    }
}




int main()
{

    double blosum62[20][20];
    double beta;
    std::cout << "Please type the value of beta: ";
    std::cin >> beta;

    read_blosum_build_kernel(blosum62, beta);



   for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            std::cout << blosum62[i][j] << std::endl;
        }
    }




    return 0;
}
