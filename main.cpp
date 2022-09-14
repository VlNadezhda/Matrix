#include <iostream>
#include "Matrix.h"
#include <vector>


int main(int argc, char *argv[]){

    Matrix<int> M{{3,4,5},
                  {3,5,6}};
    Matrix<int> N{{3,3},
                  {3,3},
                 {3,3}};

    std::cout<<M*N<<"\n";
//    {{"dhg","gjhg","hjy"},{"dytj","fjyt","gjyt"},{"dj","gjkg","khjg"}};
//    std::vector<std::vector<int>> vec {{3,5,6},{3.24,56.66}};
//    M(1,1) = 3;
    return 0;
}
