#include <iostream>
#include <vector>
#include "neighbourfind2.hpp"

using namespace std;

int main()
{
    int Nxj = 100;
    int Nyj = 100;
    double h = 1;

    int Nj = Nxj*Nyj;
    vector<double> xj(Nj);
    vector<double> yj(Nj);

    int k = 0;
    for (int i = 0; i < Nxj; i++)
    {
        for (int j = 0; j < Nyj; j++)
        {
            xj[k] = i*h;
            yj[k] = j*h;
            k = k + 1;
        }
    }
}