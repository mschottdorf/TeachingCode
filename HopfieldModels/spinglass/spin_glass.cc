
//////////////////////////////////////////////////////////////////
//
// This counts the metastable states of the SK spin glass model.
//
// The number of metastable states is well described by an exponential
// N = a*exp(b*L)
// with the number of spins "L", and the coefficients
// a = 0.199(1), b = of 1.04(1)
//
// For the theory, see F Tanakat and S F Edwards: "Analytic theory of the
// ground state properties of a  spin glass: I. Ising spin glass" J. Phys.
// F: Met. Phys. 10 2769 1980
// They found that a=0.19923 and b close to 1.
/////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
using namespace std;

#include "rnd.h"
#include "cmath"


int N=25;							//Initialisation parameters
double samples=100;

double bm(void){						//Box-Mueller-Transformation for mu=0; sigma=1
    return sqrt(-2.*log(rnd()))*cos(2.*3.141592653587*rnd());	//to get an gaussean distribution
}

int* increase(int* array, int n){				// Increase function, to initialise a new spin-configuration
	for(int i=0;i<n;i++){					// principaly bool +1:   0000 ->0001->0010->0011->0100 etc
		if(array[i]==1)					// for 1->1 and 0->-1
			array[i]=-1;
		else{
			array[i]=1;
			i=n;
		}
	}
return array;
}

int metacheck(int* spins, double** couple, int n){		//Is *spin metastable? true-> return 1, false ->return 0.

int var=1;
double sum;

for(int i=0;i<n&&var;i++){
	sum=0;
	for(int j=0;j<n;j++)
		sum+=couple[i][j]*spins[j];
	if(spins[i]*sum<=0)
		var=0;
}

return var;
}

//----------------------------------------------------------------------------//
//----------------------------MAIN--------------------------------------------//
//----------------------------------------------------------------------------//


int main(){
    
    InitRandom(1);
    ofstream os1("./results_6.dat");
    os1<<"Number of Spins: "<<N<<endl;

    int *spin=new int[N];					//Array for the Spins
    double** couple = new double*[N];			//Array for Coupling-Parameters
    for(int i=0;i<N;i++)
        couple[i] = new double[N];
    double limit=pow(2,N-1);				//The Limit is N-1 from the Symmetry of the Problem.
                                //Later the result is multiplied by 2 to obtain the correct result.
    int counter=0;

    for(int n=0;n<samples;n++){

        for(int j=0;j<N;j++){
            for(int i=j+1;i<N;i++){
                double h=bm();			//Symmetric initialisation of the J-Matrix
                couple[i][j]=h;
                couple[j][i]=h;
            }
            couple[j][j]=0;				//J_{ii}=0 because of no interaction of a Spin with itself.
            spin[j]=-1;				//And all Spins are up
        }


        for(int j=0;j<limit;j++){			//Checking half the Number of all Spins
            counter+=metacheck(spin,couple,N);
            increase(spin,N);
        }
        
        cout<<"Sample: "<<n+1<<" of "<<samples<<endl;        //Statusausgabe
        os1<<n+1<<"	"<<counter*2<<endl;		//And Multiplying the Number of metastable States with 2.
        counter=0;
    }

return 0;
} 



