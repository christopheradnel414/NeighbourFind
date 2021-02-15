#ifndef NEIGHBOURFIND_HPP
#define NEIGHBOURFIND_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <chrono>

using namespace std;

void zip   (vector<int> sorted,
            vector<double> sorter,
            vector<pair<int,double>> &zipped)
{
    for (int i = 0; i < sorted.size(); i++)
    {
        zipped.push_back(make_pair(sorted[i],sorter[i]));
    }
}

void unzip (vector<int> &sorted,
            vector<double> &sorter,
            vector<pair<int,double>> zipped)
{
    for (int i = 0; i < sorted.size(); i++)
        {
            sorted[i] = zipped[i].first;
            sorter[i] = zipped[i].second;
        }
}

void create_grid   (vector<vector<int>>& cellneighbouroutput,
                    vector<int>& cellneighbour_numoutput,
                    vector<double>& xstartoutput,
                    vector<double>& xendoutput,
                    vector<double>& ystartoutput,
                    vector<double>& yendoutput,
                    vector<double> x,
                    vector<double> y,
                    double rc,
                    double h)
{
    int size = x.size();

    double x_min = *min_element(x.begin(),x.end()) - 0.0000001;
    double x_max = *max_element(x.begin(),x.end()) + 0.0000001;
    double y_min = *min_element(y.begin(),y.end()) - 0.0000001;
    double y_max = *max_element(y.begin(),y.end()) + 0.0000001;
    
    double dY = 6.0*h; // Can be adjusted for efficiency
    double dX = 6.0*h; // Can be adjusted for efficiency

    // Checking (dX & dY must be larger than rc)
    if (dX < rc || dY < rc)
    {
        cout << "WARNING: dY or dX is set lower than rc" << endl;
    }

    int NcellX = ceil((x_max-x_min)/dX);
    int NcellY = ceil((y_max-y_min)/dY);

    int Ncell = NcellX*NcellY;

    vector<vector<int>> cellneighbour(Ncell, vector<int> (9));
    vector<int> cellneighbour_num(Ncell);
    vector<double> xstart(Ncell);
    vector<double> xend(Ncell);
    vector<double> ystart(Ncell);
    vector<double> yend(Ncell);

    // Dividing Domain
    int counter1 = 0;
    int counter2 = 0;
    for (int i = 0; i < NcellX; i++)
    {
        for (int j = 0; j < NcellY; j++)
        {
            // Cell Neighbour Search
            counter2 = 0;
            if ((i != 0) && (j != 0))
            {
                cellneighbour[counter1][counter2] = (counter1-NcellY-1);
                counter2 = counter2+1;
            }
            if ((i != 0) && (j != NcellY-1))
            {
                cellneighbour[counter1][counter2] = (counter1-NcellY+1);
                counter2 = counter2+1;
            }
            if ((i != NcellX-1) && (j != 0))
            {
                cellneighbour[counter1][counter2] = (counter1+NcellY-1);
                counter2 = counter2+1;
            }
            if ((i != NcellX-1) && (j != NcellY-1))
            {
                cellneighbour[counter1][counter2] = (counter1+NcellY+1);
                counter2 = counter2+1;
            }
            if (i != 0)
            {
                cellneighbour[counter1][counter2] = (counter1-NcellY);
                counter2 = counter2+1;
            }
            if (i != NcellX-1)
            {
                cellneighbour[counter1][counter2] = (counter1+NcellY);
                counter2 = counter2+1;
            }
            if (j != 0)
            {
                cellneighbour[counter1][counter2] = (counter1-1);
                counter2 = counter2+1;
            }
            if (j != NcellY-1)
            {
                cellneighbour[counter1][counter2] = (counter1+1);
                counter2 = counter2+1;
            }

            cellneighbour[counter1][counter2] = counter1; // Self neighbour
            counter2 = counter2+1;

            cellneighbour_num[counter1] = counter2;

            // Cell Starting and Ending Coordinates
            xstart[counter1] = x_min + i*dX;
            xend[counter1] = x_min + (i+1)*dX;
            ystart[counter1] = y_min + j*dY;
            yend[counter1] = y_min + (j+1)*dY;

            counter1 = counter1 + 1;
        }
    }

    cellneighbouroutput = cellneighbour;
    cellneighbour_numoutput = cellneighbour_num;
    xstartoutput = xstart;
    xendoutput = xend;
    ystartoutput = ystart;
    yendoutput = yend;
}

void assign_grid   (vector<int>& cellidoutput,
                    vector<vector<int>>& cell_particleoutput,
                    vector<int>& cell_particle_numoutput,
                    vector<double> x,
                    vector<double> y,
                    vector<double> xstart,
                    vector<double> xend,
                    vector<double> ystart,
                    vector<double> yend,
                    int Ncell,
                    int part_num)
{
    vector<int> cellid(part_num,-2);
    vector<vector<int>> cell_particle(Ncell);
    vector<int> cell_particle_num(Ncell);
    //#pragma omp parallel for
    for (int k = 0; k < part_num; k++)
    {
        for (int i = 0; i < Ncell; i++)
        {
            if ((x[k] >= xstart[i]) && (x[k] <= xend[i]) && (y[k] >= ystart[i]) && (y[k] <= yend[i]))
            {
                //#pragma omp critical
                cell_particle[i].push_back(k);
                cell_particle_num[i] = cell_particle_num[i] + 1;
                cellid[k] = i;
                break;
            }
        }
    }
    cellidoutput = cellid;
    cell_particleoutput = cell_particle;
    cell_particle_numoutput = cell_particle_num;
}

void find_neighbours_grid  (vector<vector<int>>& neighbourdata,
                            vector<vector<double>>& rangedata,
                            vector<double> x_i,
                            vector<double> y_i,
                            vector<double> x_j,
                            vector<double> y_j,
                            vector<int> cellid_i,
                            vector<int> cellid_j)
{

}























#endif