#include <standard>
#include <fstream>
#include "rnd.h"


const int N=5000;
double samples=5;
double agreement=0.8;
double alpha=0.07;



int main(){

InitRandom();
ofstream os;
os.open("hopfield_5000_0.138_0.140_critcheck.dat");
os.precision(14);


double ** w= new double*[N];
for(int j=0;j<N;j++)
	w[j]=new double[N];

double ** xmu = new double*[N];
for(int j=0;j<N;j++)					//P<N
	xmu[j]=new double[N];


double * u = new double[N];
double * s = new double[N];

for(double alpha=0.138;alpha<=0.140;alpha+=0.001){			//Alpha is the ration of patterns to neurons.
cout<<"Alpha: "<<alpha<<endl;
	int P=N*alpha;	//Patterns = neurons * alpha
	double r=1./N;
	double x=0;

	for(int run=0;run<samples;run++){
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++)
					w[i][j]=0;

		for(int mu=0;mu<P;mu++){				//Calculation W_ij from N random patterns
			for(int j=0;j<N;j++)
				xmu[mu][j]=rndsign();

			for(int i=0;i<N;i++)
				for(int j=i+1;j<N;j++)
					w[j][i]=w[i][j]+=xmu[mu][i]*xmu[mu][j]*r;
		}
		
		
		cout<<"Init done"<<endl;				//From here: Can all patterns be reconstructed?
		for(int mu=0;mu<P;mu++){
			bool change=true;
			double m=0;
			cout<<"Mu ist "<<mu<<endl;
			for(int i=0;i<N;i++){				//Take pattern mu
				u[i]=0;
				s[i]=xmu[mu][i];
				}
		
			for(int counter=0;counter<N*(1-agreement);){	//And add noise by random flipping bits
				int i=rndint(N);
				if(u[i]==0){
					s[i]=-s[i];
					u[i]=1;
					counter++;
				}
			}

			while(change){					//In here are the Hebbian updates. The loop runs as long as a static point is reached.
				change=false;
				for(int i=0;i<N;i++){
					double h=0;
					for(int j=0;j<N;j++)
						h+=w[i][j]*s[j];	//local field
		
					if(h>0){
						if(s[i]==-1){
							s[i]=1;
							change=true;
						}
					}else{
						if(s[i]==1){
							s[i]=-1;
							change=true;
						}
					}
				}
			}
		
			for(int i=0;i<N;i++)				//Compare the minimum with the original pattern
				m+=s[i]*xmu[mu][i];
			if(m*r>0.9)
				x++;
		}
	cout<<run<<endl;
	}
	os<<alpha<<"\t"<<(x/samples)/P<<endl;
	cout<<"Save: "<<alpha<<"\t"<<x/samples/P<<endl;		//Save number of patterns and the ratio of patterns identified with >90%
}
os.close();


return 0;
}

