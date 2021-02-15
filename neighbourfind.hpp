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

void neighbourfind (vector<vector<int>>& neighbour,
                    vector<double> x,
                    vector<double> y,
                    double rc,
                    double h)

{   
    cout << "Neighbour Finding ...";
    auto start = chrono::steady_clock::now();

    int part_num = x.size();

    double xmin = *min_element(x.begin(),x.end()) - 0.0000001;
    double xmax = *max_element(x.begin(),x.end()) + 0.0000001;
    double ymin = *min_element(y.begin(),y.end()) - 0.0000001;
    double ymax = *max_element(y.begin(),y.end()) + 0.0000001;

    //double rc_per_h = rc/h;

    //double dY = 6.0/(rc_per_h)*rc;
    //double dX = 6.0/(rc_per_h)*rc;

    double dY = 6.0*h;
    double dX = 6.0*h;

    //double maxPartInCell = 4.0*dY*dX/(h*h);

    int NcellX = ceil((xmax-xmin)/dX);
    int NcellY = ceil((ymax-ymin)/dY);

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
            xstart[counter1] = xmin + i*dX;
            xend[counter1] = xmin + (i+1)*dX;
            ystart[counter1] = ymin + j*dY;
            yend[counter1] = ymin + (j+1)*dY;

            counter1 = counter1 + 1;

        }
    }

    // Assigning Particles into Cells and Assigning Cell Number to Each Particles
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

    double dx;
    double dy;
    double dr;
    int nbi;
    vector<vector<int>> neighbourdata(part_num);
    for (int k = 0; k < part_num; k++)
    {
        for (int i = 0; i < cellneighbour_num[cellid[k]]; i++)
        {
            for (int j = 0; j < cell_particle_num[cellneighbour[cellid[k]][i]]; j++)
            {
                nbi = cell_particle[cellneighbour[cellid[k]][i]][j];
                dx = x[nbi] - x[k];
                dy = y[nbi] - y[k];
                dr = sqrt(dx*dx + dy*dy);

                if ((dr < rc) && k != nbi)
                {
                    neighbourdata[k].push_back(nbi);
                }
            }
        }
    }

    neighbour = neighbourdata;
    
    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}

void neighbourfindlimited (vector<vector<int>>& neighbour,
                    vector<double> x,
                    vector<double> y,
                    double rc,
                    double h,
                    int neighbourlimit)

{   
    cout << "Neighbour Finding ...";
    auto start = chrono::steady_clock::now();

    int part_num = x.size();

    double xmin = *min_element(x.begin(),x.end()) - 0.0000001;
    double xmax = *max_element(x.begin(),x.end()) + 0.0000001;
    double ymin = *min_element(y.begin(),y.end()) - 0.0000001;
    double ymax = *max_element(y.begin(),y.end()) + 0.0000001;

    //double rc_per_h = rc/h;

    //double dY = 6.0/(rc_per_h)*rc;
    //double dX = 6.0/(rc_per_h)*rc;

    double dY = 6.0*h;
    double dX = 6.0*h;

    //double maxPartInCell = 4.0*dY*dX/(h*h);

    int NcellX = ceil((xmax-xmin)/dX);
    int NcellY = ceil((ymax-ymin)/dY);

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
            xstart[counter1] = xmin + i*dX;
            xend[counter1] = xmin + (i+1)*dX;
            ystart[counter1] = ymin + j*dY;
            yend[counter1] = ymin + (j+1)*dY;

            counter1 = counter1 + 1;

        }
    }

    // Assigning Particles into Cells and Assigning Cell Number to Each Particles
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

    double dx;
    double dy;
    double dr;
    int nbi;

    vector<vector<int>> neighbourdata(part_num);


    vector<int> temp_neighbourdata;
    vector<double> temp_drdata;
    vector<pair<int,double>> zipped;

    for (int k = 0; k < part_num; k++)
    {

        temp_neighbourdata.clear();
        temp_drdata.clear();

        for (int i = 0; i < cellneighbour_num[cellid[k]]; i++)
        {
            
            for (int j = 0; j < cell_particle_num[cellneighbour[cellid[k]][i]]; j++)
            {
                nbi = cell_particle[cellneighbour[cellid[k]][i]][j];
                dx = x[nbi] - x[k];
                dy = y[nbi] - y[k];
                dr = sqrt(dx*dx + dy*dy);

                if ((dr < rc) && k != nbi)
                {
                    temp_neighbourdata.push_back(nbi);
                    temp_drdata.push_back(dr);
                }
            }
        }

        if (temp_neighbourdata.size() < neighbourlimit)
        {
            cout<<"Neighbour Search Error: Not enough neighbour"<<endl;
            exit(EXIT_FAILURE);
        }

        zip(temp_neighbourdata,temp_drdata,zipped);

        sort(begin(zipped), end(zipped), 
        [&](const auto& a, const auto& b)
        {
            return a.second < b.second;
        });

        unzip(temp_neighbourdata,temp_drdata,zipped);
        zipped.clear();

        for (int i = 0; i < neighbourlimit; i++)
        {
            neighbourdata[k].push_back(temp_neighbourdata[i]);
        }
    }

    neighbour = neighbourdata;
    
    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}

void neighbourfindlimitedboundary  (vector<vector<int>>& neighbour,
                                    vector<int>& isboundary,
                                    vector<double> x,
                                    vector<double> y,
                                    double rc,
                                    double h,
                                    int neighbourlimit,
                                    int boundlimit)

{   
    int part_num = x.size();

    double xmin = *min_element(x.begin(),x.end()) - 0.0000001;
    double xmax = *max_element(x.begin(),x.end()) + 0.0000001;
    double ymin = *min_element(y.begin(),y.end()) - 0.0000001;
    double ymax = *max_element(y.begin(),y.end()) + 0.0000001;

    //double rc_per_h = rc/h;

    //double dY = 6.0/(rc_per_h)*rc;
    //double dX = 6.0/(rc_per_h)*rc;

    double dY = 6.0*h;
    double dX = 6.0*h;

    //double maxPartInCell = 4.0*dY*dX/(h*h);

    int NcellX = ceil((xmax-xmin)/dX);
    int NcellY = ceil((ymax-ymin)/dY);

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
            xstart[counter1] = xmin + i*dX;
            xend[counter1] = xmin + (i+1)*dX;
            ystart[counter1] = ymin + j*dY;
            yend[counter1] = ymin + (j+1)*dY;

            counter1 = counter1 + 1;

        }
    }

    // Assigning Particles into Cells and Assigning Cell Number to Each Particles
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

    double dx;
    double dy;
    double dr;
    int nbi;

    vector<vector<int>> neighbourdata(part_num);
    vector<int> isboundarytemp(part_num);


    vector<int> temp_neighbourdata;
    vector<double> temp_drdata;
    vector<pair<int,double>> zipped;

    for (int k = 0; k < part_num; k++)
    {

        temp_neighbourdata.clear();
        temp_drdata.clear();

        for (int i = 0; i < cellneighbour_num[cellid[k]]; i++)
        {
            
            for (int j = 0; j < cell_particle_num[cellneighbour[cellid[k]][i]]; j++)
            {
                nbi = cell_particle[cellneighbour[cellid[k]][i]][j];
                dx = x[nbi] - x[k];
                dy = y[nbi] - y[k];
                dr = sqrt(dx*dx + dy*dy);

                if ((dr < rc) && k != nbi)
                {
                    temp_neighbourdata.push_back(nbi);
                    temp_drdata.push_back(dr);
                }
            }
        }

        if (temp_neighbourdata.size() < neighbourlimit)
        {
            cout<<"Neighbour Search Error: Not enough neighbour"<<endl;
            exit(EXIT_FAILURE);
        }

        if (temp_neighbourdata.size() < boundlimit)
        {
            isboundarytemp[k] = 1;
        }
        else
        {
            isboundarytemp[k] = 0;
        }
        

        zip(temp_neighbourdata,temp_drdata,zipped);

        sort(begin(zipped), end(zipped), 
        [&](const auto& a, const auto& b)
        {
            return a.second < b.second;
        });

        unzip(temp_neighbourdata,temp_drdata,zipped);
        zipped.clear();

        for (int i = 0; i < neighbourlimit; i++)
        {
            neighbourdata[k].push_back(temp_neighbourdata[i]);
        }
    }

    neighbour = neighbourdata;
    isboundary = isboundarytemp;
}

void neighbourfindboundary  (vector<vector<int>>& neighbour,
                                    vector<int>& isboundary,
                                    vector<double> x,
                                    vector<double> y,
                                    double rc,
                                    double h,
                                    int boundlimit)

{   
    cout << "Neighbour Finding ...";
    auto start = chrono::steady_clock::now();

    int part_num = x.size();

    double xmin = *min_element(x.begin(),x.end()) - 0.0000001;
    double xmax = *max_element(x.begin(),x.end()) + 0.0000001;
    double ymin = *min_element(y.begin(),y.end()) - 0.0000001;
    double ymax = *max_element(y.begin(),y.end()) + 0.0000001;

    //double rc_per_h = rc/h;

    //double dY = 6.0/(rc_per_h)*rc;
    //double dX = 6.0/(rc_per_h)*rc;

    double dY = 6.0*h;
    double dX = 6.0*h;

    //double maxPartInCell = 4.0*dY*dX/(h*h);

    int NcellX = ceil((xmax-xmin)/dX);
    int NcellY = ceil((ymax-ymin)/dY);

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
            xstart[counter1] = xmin + i*dX;
            xend[counter1] = xmin + (i+1)*dX;
            ystart[counter1] = ymin + j*dY;
            yend[counter1] = ymin + (j+1)*dY;

            counter1 = counter1 + 1;

        }
    }

    // Assigning Particles into Cells and Assigning Cell Number to Each Particles
    vector<int> cellid(part_num);
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

    double dx;
    double dy;
    double dr;
    int nbi;

    vector<vector<int>> neighbourdata(part_num);
    vector<int> isboundarytemp(part_num);


    vector<int> temp_neighbourdata;
    vector<double> temp_drdata;
    vector<pair<int,double>> zipped;

    for (int k = 0; k < part_num; k++)
    {

        temp_neighbourdata.clear();
        temp_drdata.clear();

        for (int i = 0; i < cellneighbour_num[cellid[k]]; i++)
        {
            
            for (int j = 0; j < cell_particle_num[cellneighbour[cellid[k]][i]]; j++)
            {
                nbi = cell_particle[cellneighbour[cellid[k]][i]][j];
                dx = x[nbi] - x[k];
                dy = y[nbi] - y[k];
                dr = sqrt(dx*dx + dy*dy);

                if ((dr < rc) && k != nbi)
                {
                    temp_neighbourdata.push_back(nbi);
                    temp_drdata.push_back(dr);
                }
            }
        }

        if (temp_neighbourdata.size() < boundlimit)
        {
            isboundarytemp[k] = 1;
        }
        else
        {
            isboundarytemp[k] = 0;
        }

        for (int i = 0; i < temp_neighbourdata.size(); i++)
        {
            neighbourdata[k].push_back(temp_neighbourdata[i]);
        }
    }

    neighbour = neighbourdata;
    isboundary = isboundarytemp;
    
    auto end = chrono::steady_clock::now();
    cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}


void neighbourfindwithcollision(vector<vector<int>>& neighbour,
                                vector<vector<int>>& neighbour_colls,
                                vector<double> x,
                                vector<double> y,
                                double rc,
                                double r_collision,
                                double h)

{   
    //cout << "Neighbour Finding ...";
    //auto start = chrono::steady_clock::now();

    int part_num = x.size();

    double xmin = *min_element(x.begin(),x.end()) - 0.0000001;
    double xmax = *max_element(x.begin(),x.end()) + 0.0000001;
    double ymin = *min_element(y.begin(),y.end()) - 0.0000001;
    double ymax = *max_element(y.begin(),y.end()) + 0.0000001;

    //double rc_per_h = rc/h;

    //double dY = 6.0/(rc_per_h)*rc;
    //double dX = 6.0/(rc_per_h)*rc;

    double dY = 6.0*h;
    double dX = 6.0*h;

    //double maxPartInCell = 4.0*dY*dX/(h*h);

    int NcellX = ceil((xmax-xmin)/dX);
    int NcellY = ceil((ymax-ymin)/dY);

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
            xstart[counter1] = xmin + i*dX;
            xend[counter1] = xmin + (i+1)*dX;
            ystart[counter1] = ymin + j*dY;
            yend[counter1] = ymin + (j+1)*dY;

            counter1 = counter1 + 1;

        }
    }

    // Assigning Particles into Cells and Assigning Cell Number to Each Particles
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

    double dx;
    double dy;
    double dr;
    int nbi;
    vector<vector<int>> neighbourdata(part_num);
    vector<vector<int>> collisiondata(part_num);
    for (int k = 0; k < part_num; k++)
    {
        for (int i = 0; i < cellneighbour_num[cellid[k]]; i++)
        {
            for (int j = 0; j < cell_particle_num[cellneighbour[cellid[k]][i]]; j++)
            {
                nbi = cell_particle[cellneighbour[cellid[k]][i]][j];
                dx = x[nbi] - x[k];
                dy = y[nbi] - y[k];
                dr = sqrt(dx*dx + dy*dy);

                if ((dr < rc) && k != nbi)
                {
                    neighbourdata[k].push_back(nbi);

                    if (dr < r_collision)
                    {
                        collisiondata[k].push_back(nbi);
                    }
                }
            }
        }
    }

    neighbour = neighbourdata;
    neighbour_colls = collisiondata;
    
    //auto end = chrono::steady_clock::now();
    //cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}

void neighbourfindold (vector<vector<int>>& neighbour,
                    vector<double> x,
                    vector<double> y,
                    double rc,
                    double h)
{
    auto start = chrono::steady_clock::now();
    double dX;
    double dY;
    double dR;

    int particlenumber = x.size();
    vector<vector<int>> neighbourdata(particlenumber);
    
    for (int i = 0; i < particlenumber-1; i++)
    {
        for (int j = i+1; j < particlenumber; j++)
        {
            dX = x[i] - x[j];
            dY = y[i] - y[j];
            dR = sqrt(pow(dX,2) + pow(dY,2));
            if (dR < rc)
            {
                neighbourdata[i].push_back(j);
                neighbourdata[j].push_back(i);
            }

        }
    }
    neighbour = neighbourdata;
    auto end = chrono::steady_clock::now();
    //cout << "Finished (" <<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" ms)"<< endl;
}
#endif