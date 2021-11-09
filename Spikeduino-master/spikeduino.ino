//
// This small program is a event based simulation of a neural network 
// on a random graph, coded for an arduino.
//

const double N_PI = 3.1415;
const int NUMBER_OF_PINRONS = 57;   //Number of Pins is number of Neurons
int K = 5;
double J0 = 1;
double tauM = 0.01;


boolean calculate_next_spike = false;
double I0 = 0.5*tauM*sqrt(K);
double w = 2.*sqrt( I0 * sqrt(K));
double C = -1.* J0 / (sqrt(K)*sqrt(I0*sqrt(K)));
int A[NUMBER_OF_PINRONS][NUMBER_OF_PINRONS];
int phi[NUMBER_OF_PINRONS];
long blinktime = 0;

int ledPIN = 0;
int ledPIN_old = 0;
double phi_j = -N_PI;
int j = 0;

void setup() {
    randomSeed(42);  //To make the random matrix and the phases reproducible.
    Serial.begin(9600);
    for(int ledPin=0; ledPin<NUMBER_OF_PINRONS; ledPin++){
        pinMode(ledPin, OUTPUT); 
        digitalWrite(ledPin, LOW);
        phi[ledPIN] = random(-(int) N_PI*10000,(int)N_PI*10000)/10000.;   //The phase is initialized as a random number between -pi and pi.
        for(int target=0; target<NUMBER_OF_PINRONS; target++){  //The connectivity matrix is 0|1
            if(random(0,100) < (100.*K/NUMBER_OF_PINRONS) ){
                A[target][ledPin] = 1;
            }else{
                A[target][ledPin] = 0;
            }
        }
    }
}


void loop(){
    if(calculate_next_spike){
        //This finds the next spiking neuron "j" with phase "phi_j"
        phi_j = -N_PI;
        j = 0;
        for(int i=0; i<NUMBER_OF_PINRONS; i++){
            if(phi_j < phi[i]){
                phi_j = phi[i];
                j = i;
            }
        }
        //After the spike, this forwards all the phases given the phase velocity.
        //the postsynapic cells are advanced even more corresponding to the phase response curve.
        double dt = ( N_PI - phi_j )/w;
        for(int i=0; i<NUMBER_OF_PINRONS; i++){
            phi[i] += w*dt;
            if(A[i][j] == 1){    //Only the postsynaptic cells
                phi[i] = 2.*atan( tan( phi[i]/2.) + C );
            }
        }
        //Finally spiking cell is reset
        phi[j] = -N_PI; 

        calculate_next_spike = false;
        blinktime = millis() + dt*tauM;  //this is the next spiketime
    }

 
    if(blinktime < millis()) {           //Once next spiketime has been reached...
        blinktime = millis(); 
        ledPIN = j; //This neuron "fired" and the output is switched on.
        digitalWrite(ledPIN, HIGH);
        digitalWrite(ledPIN_old, LOW); //the last neuron is switched off.
        
        /* Health Check; firing cells' indices are sent via serial port */
        Serial.print(ledPIN);
        Serial.print("\n");
        
        ledPIN_old = ledPIN;
        calculate_next_spike = true;
    }

}
