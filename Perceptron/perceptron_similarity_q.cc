#include <standard>
#include <fstream>
#include "rnd.h"


const int N=500;
double alphamax=3;
int Pmax=alphamax*N;


const int samples=100;


int compare(const void *a,const void *b){
 return (*(double*)a-*(double*)b);
}


double median(double* arr, int size)
{
int middle = size/2;
double a;
if (size%2==0) 
	a=(arr[middle-1]+arr[middle])/2.;
else
	a=(double)arr[middle];

return a;
}




int main(){

InitRandom();
ofstream os;
os.open("corr_0.2.dat");

os.precision(14);



double ** w = new double*[N];				//Neuron viele Gewichtungsparameter
for(int j=0;j<N;j++)
	w[j]=new double[samples];

double * O = new double[Pmax];				//Patterns viele Output-Ergebnisse

double ** x= new double*[N];				//Patterns viele x-Vektoren mit N Elementen
for(int j=0;j<N;j++)
	x[j]=new double[Pmax];

double * datm= new double[samples];			//Table to store information for median


int recmax=5000;					//Max Recursions starting

double alpha=0.2;

	int P=alpha*N;
	for(int i=0;i<samples;i++)
		datm[i]=0;

	for(int i=0;i<P;i++){
		O[i]=rndsign();
		for(int j=0;j<N;j++){
			x[j][i]=rndsign();
		}
	}

	for(int c=0;c<samples;c++){

		for(int i=0;i<P;i++){					//Init
			for(int j=0;j<N;j++){
				w[j][c]=2*rnd()-1;
			}
		}
		bool flag=true;
		int escape=0;
		double foundratio=0;					//Number of Patterns that are found
		double value=0;
	
		while(flag  &&  escape<recmax){				//Perceptron Learning Algorithm
			foundratio=0;
			flag=false;
			for(int i=0;i<P;i++){				//Run all patterns
				value=0;
				for(int j=0;j<N;j++)
					value+=w[j][c]*x[j][i];
				if(O[i]*value<0){
					for(int j=0;j<N;j++){
						w[j][c]+=x[j][i]*O[i];
						flag=true;
					}
				}else					//Pattern identified
					foundratio++;
			}
			if(foundratio==P)
				datm[c]=escape;
			if(escape==recmax-1)
				datm[c]=recmax;
			escape++;
		}
		cout<<c<<endl;
	}

	qsort(datm,samples,sizeof(double),compare);				//Alg. for rec. depth
	if(alpha>0.5)
		recmax=10*median(datm,samples)*(alpha-0.1);
	else
		recmax=1000;

	for(int i=0;i<samples;i++){
		for(int j=0;j<samples;j++){
			if(j==i)
				j++;
			double z=0;
			double n1=0;
			double n2=0;
			for(int k=0;k<N;k++){
				z+=w[k][i]*w[k][j];
				n1+=w[k][i]*w[k][i];
				n2+=w[k][j]*w[k][j];
			}
			os<<z/sqrt(n1*n2)<<endl;
		}
	}


os.close();

return 0;
}

