#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<gsl/gsl_math.h>
#include<gsl/gsl_complex_math.h>

//compile:  g++ Sk.cpp -lgsl -lgslcblas

using namespace std;

int N;  //number of particles,
//THIS IS USEFUL for scaling purpose, has to be made consistent with the actual number of particles in the system, although maynot be exactly same...
int cell_counter = 0; //the acutally number of cells within a square box

const double Kmax = 80; //K-space cutoff
const int maxN = 75;  //counting threshold; greater than Kmax
const int mesh = 500;  //k-space mesh


double GetMin(double a, double b)
{
    if(a<b) return a;
    else return b;
}


//NOTES:  It is computationally more efficient to compute everything with complex variables, requiring only O(N) operations instead of O(N^2).
int main()
{
    gsl_complex Sksum = gsl_complex_rect(0, 0);
    
    double dkx, dky, kmag;
    
    double* x;
    double* y; //for the coordinates...
    fstream inputFile;
    
    /*
     inputFile.open("coords", ios::in);
     if(!inputFile)
     {
     cout<<"Can not open coords! Abort!"<<endl;
     exit(1);
     }
     */
    
    double MAXX, MAXY;
    
    cin >> MAXX; cin >> MAXY;
    cin >> N;
    
    double L = GetMin(MAXX, MAXY);  //box length, the shortest edge length if a rectangular box is used...
    
    double deltaK = Kmax/mesh*(2.0*M_PI/L); //spacing for k-space mesh
    
    
    x = new double[N];
    y = new double[N];
    
    double tempx, tempy;
    cell_counter = 0;
    
    for(int i=0; i<N; i++)
    {
        cin >> tempx;
        cin >> tempy;
        
        if(tempx<L && tempy<L)
        {
            x[cell_counter] = tempx;
            y[cell_counter] = tempy;
            
            cell_counter++;
        }
    }
    
    N = cell_counter;
    
    double rho = N/(L*L);  //density
    
    //inputFile.close();  //input coordinates
    double count[mesh];
    double S[mesh];
    double k[mesh];
    double bin[mesh]; //for averaging over wavevectors
    double Sk;
    
    for(int i=0; i<mesh; i++)
    {
        S[i] = 0;
        k[i] = (i+0.5)*deltaK;
        bin[i] = i*deltaK;
        count[i] = 0;
    }
    
    for(int ki=0; ki<=2*maxN; ki++)
    {
        for(int kj=0; kj<=2*maxN; kj++)
        {
            dkx = ki-maxN;
            dky = kj-maxN;
            
            if(dkx==0 && dky==0)
            {
                continue;
            }
            if(pow(dkx, 2.0)+pow(dky, 2.0)> pow(Kmax, 2.0))
            {
                continue;
            }
            
            dkx*=2.0*M_PI/L;
            dky*=2.0*M_PI/L;
            kmag = sqrt(pow(dkx, 2.0)+pow(dky, 2.0));
            
            Sksum = gsl_complex_rect(0,0);
            
            for(int i=0; i<N; i++)
            {
                Sksum = gsl_complex_add(Sksum, gsl_complex_exp(gsl_complex_rect(0, dkx*x[i]+dky*y[i])));
            }  //sum over all collective coordinates
            
            Sk = gsl_complex_abs2(Sksum)/N;  //square-modulus over number of particles
            
            for(int i=0; i<mesh; i++)
            {
                if(bin[i]<=kmag && kmag<bin[i]+deltaK)
                {
                    count[i]+=1.0;
                    S[i]+=Sk;
                    break;
                }
            }  //bin Sk and keep track of binned wavevectors
        }
    }
    
    fstream outputFile;
    outputFile.open("Sk_bin.txt", ios::out);
    
    for(int i=0; i<mesh; i++)
    {
        if(count[i]!=0)
        {
            outputFile << k[i]/(2*M_PI*sqrt(rho)) << " " << S[i]/count[i] << endl;
        }
    } //output radially-averaged structure factor
    
    
    return 0;
}
