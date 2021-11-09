#
#  This program allows to read the Serial data stream from the Arduino, and 
#  plots the spiking neurons on a grid.
#

import serial
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as animation

ser = serial.Serial('/dev/ttyACM2', 9600)

fig = pl.figure()
rr = np.zeros((7,7))
img = pl.imshow(rr, cmap='gist_gray_r', vmin=0, vmax=1)

def animate(i):
    by = ser.readline()
    if len(by)>2:
        neuron_index = float(by)
        if neuron_index<49:
            print neuron_index
            rr = np.zeros((7,7))
            rr = rr.reshape((1,49))
            rr[0,int(neuron_index)]=1
            rr = rr.reshape(7,7)
            img.set_data(rr)
    else:
        img.set_data(np.zeros((7,7)))
    return img

def init():
    img.set_data(np.zeros((7,7)))

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=200000, interval=10)

pl.show()
