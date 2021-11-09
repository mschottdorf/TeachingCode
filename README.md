# TeachingCode
 This repository contains various toy examples for influential comp neuro models and ideas. 

# NetworkModels

The file `network_models.ipynb` contains an implementation of the Ring model, and an event based simulation of a spiking RNN to plot a Poincare section.

# KeplerLSTM
This contains a jupyter notebook to fit an LSTM to Tycho Brahe's Mars data. The original data is somewhat unyielding, so instead the data is obtained through JPL's [HORIZONS system](https://ssd.jpl.nasa.gov/horizons.cgi?s_time=1).

# The XY model

An implementation of the XY model in C. Read summary [on wiki](https://en.wikipedia.org/wiki/.Classical_XY_model). Output with OpenCV.

# Spinglasses

Counting the number of ground states of a Spinglas.

# HopfieldModel

Implements, and numerically calculates the capacity of a Hopfield network. Evaluates to a=0.138.

# Perceptron

Implements, and numerically calculates the capacity of the perceptron. Finds a=2.

# Spikeduino

An event based neural network implementation for an arduino. The C Program contains the code to be run on the chip. The Python Program allows to read the serial data stream, and plot it as an animation. This code is:

1. an implementation of a (small) balanced networks. I.e. synapses scale like 1/sqrt(K) where K is the indegree. The network is a random graph that is initialized with the build-in random number generator. The typical activity pattern is asynchronous and irregular.
1. The Code generates Voltage-Pulses on the digital pins of the controller. Each pin corresponding to one neuron.
1. The serial interface sends a continuous stream of numbers, corresponding to the index of the spiking cell. 



