#include <iostream>
#include <vector>
#include "neighbourfind2.hpp"

using namespace std;

int main()
{
    int Nxj = 100;
    int Nyj = 100;
    double h = 1.0;

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

    int Nxi = 20;
    int Nyi = 20;
    h = 2.33;

    int Ni = Nxi*Nyi;

    vector<double> xi(Ni);
    vector<double> yi(Ni);

    k = 0;
    for (int i = 0; i < Nxi; i++)
    {
        for (int j = 0; j < Nyi; j++)
        {
            xi[k] = i*h;
            yi[k] = j*h;
            k = k + 1;
        }
    }

    vector<vector<int>> neighbour;
    vector<vector<double>> range;
    vector<int> isboundary;

    neighbourfind_limited_boundary(neighbour,range,isboundary,xi,yi,xj,yj,4.601*1.0,1.0,20,20);
}