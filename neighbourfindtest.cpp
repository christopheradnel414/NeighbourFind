#include <iostream>
#include <vector>
#include "neighbourfind2.hpp"

using namespace std;

void printcoordinates  (vector<double> x,
                        vector<double> y,
                        vector<vector<int>> neighbour);

int main()
{
    int Nxj = 10;
    int Nyj = 10;
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

    int Nxi = 10;
    int Nyi = 10;
    h = 1.0;

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

    //neighbourfind_limited_boundary(neighbour,range,isboundary,xi,yi,xj,yj,4.601*1.0,1.0,20,20);
    neighbourfind_boundary(neighbour,range,isboundary,xi,yi,xj,yj,2.101*1.0,1.0,20);

    printcoordinates(xi,yi,neighbour);

}

void printcoordinates  (vector<double> x,
                        vector<double> y,
                        vector<vector<int>> neighbour)
{
    cout << "No." << "\t" << "x" << "\t" << "y"<< "\t" << "Nb_num" << "\t" <<"Neighbours" <<"\n";
    for (int i = 0; i < x.size(); i++)
    {
        cout << i << "\t" << x[i] << "\t" << y[i]<< "\t" << neighbour[i].size() << "\t";
        for (int j = 0; j < neighbour[i].size(); j++)
        {
            cout << neighbour[i][j] <<",";
        }
        cout << endl;
    }
}