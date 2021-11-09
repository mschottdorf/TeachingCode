///////////////////////////////////////////////////////////////////////////
//
// This calculates the average fraction of correctly classified patterns for a
// perceptron of length N unter load alpha. The result is stored as a textfile
// --->  "performance.dat"
//
// in the format
// alpha | classification accuracy
//  0.1     1
//  ...     ...
//  2.5     0
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
using namespace std;
//Comfortable and fast random number generator
#include "rnd.h"


const int N = 500;         // Size of the perceptron
double alphamax = 3;       // Max load to test
double alpha_step = 0.2;   // In steps of 0.1.
const int samples = 100;     // Number of samples per load

int recmax = 2000;         // Max Recursions for perceptron learning algorithm
int Pmax=alphamax*N;       // Max patterns

int main(){

    //Randomseed
    InitRandom();
    
    //parameter for output
    ofstream osmedian;
    ofstream osconver;
    osconver.open("performance.dat");
    osmedian.precision(14);
    osconver.precision(14);

    //Perceptron data arrays and initialization
    double *w = new double[N];				    //Weights
    double *O = new double[Pmax];				//Output patterns
    double ** x= new double*[N];				//Input patterns
    for(int j=0;j<N;j++){
        x[j]=new double[Pmax];
    }
    
    
    double *separability_per_sample = new double[samples];

    //Calculate performance for different load values
    for(double alpha = alpha_step;alpha<alphamax;alpha+=alpha_step){
        int P=alpha*N;
        
        //Average over samples
        double found_patterns = 0; //counter for correctly identified patterns across samples
        for(int s=0;s<samples;s++){
            
            //Initialization for the weights and P patterns to be learned
            for(int i=0;i<P;i++){
                O[i]=rndsign();
                for(int j=0;j<N;j++){
                    x[j][i]=rndsign();
                    w[j]=rndint(40)-20;
                }
            }
            
            //Train the perceptron with all pairs of inputs|outputs
            int steps_performed=0;
            bool training_needed = true;
            while(training_needed  &&  steps_performed < recmax){
                training_needed = false;
                steps_performed++;
                for(int i=0;i<P;i++){				            // Run all patterns
                    double value=0;
                    for(int j=0;j<N;j++){
                        value += w[j]*x[j][i];                  // calculate output of perceptron for pairing x|O
                    }
                    if(O[i]*value<0){                           // if it is not correct, do a learning step
                        for(int j=0;j<N;j++){
                            w[j] += x[j][i]*O[i];               // Perceptron learning algorithm
                            training_needed = true;
                        }
                    }
                }
                
            }
            
            // Test performance
            double number_of_found_patters=0;
            for(int i=0;i<P;i++){
                double value=0;
                for(int j=0;j<N;j++){
                    value += w[j]*x[j][i];
                }
                if(O[i]*value>0){
                    number_of_found_patters++;
                }
            }
            if(number_of_found_patters == P){
                separability_per_sample[s] = 1;
            }else{
                separability_per_sample[s] = 0;
            }
        }
        
        //Calculate the average performance across samples
        double av = 0;
        for(int i=0;i<samples;i++){
            av += separability_per_sample[i];
        }
        
        //And save it
        cout<<alpha<<"\t"<<av/samples<<endl;
        osconver<<alpha<<"\t"<<av/samples<<endl;     //Save result for this load value
    }

    osconver.close();
    return 0;
}
