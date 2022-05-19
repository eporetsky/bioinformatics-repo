import os
import serial
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import argrelmax
from scipy.signal import savgol_filter

#arduino_port = "/dev/cu.usbmodem1420" #serial port of Arduino
arduino_port = "/dev/cu.usbserial-1420" #serial port of Arduino
baud = 9600 #arduino uno runs at 9600 baud
print("Connected to Arduino port:" + arduino_port)

experiment_name = "220511" #name of the CSV file generated
try:
	os.mkdir(experiment_name)
	os.mkdir(experiment_name+"/plots")
except:
	None

output = open(experiment_name+"/data.raw.csv", "a")
print("Writing raw data to: "+experiment_name+"/data.raw.csv")
output_integrated = open(experiment_name+"/data.integrated.csv", "a")
print("Writing Integrated data to: "+experiment_name+"/data.integrated.csv")

ser = serial.Serial(arduino_port, baud)
dsa = 0
ln_count = 1
data_ls = []
time_ls, time = [], 1

while True:
	dat = ser.readline().decode('ascii').split("\r")[0]

	if dat == "start":
		print("Starting: sample "+str(ln_count))
		data_ls = ["sample_"+str(ln_count)]
		time_ls = []
		time =  1
		continue
		
	# Write to file 
	if dat == "stop":
		print("Complete")
		output = open(experiment_name+"/data.raw.csv", "a")
		output.write(",".join(data_ls) + "\n") #write data with a newline
		output.close()
		
		y_raw = [int(val) for val in  data_ls[1:]]
		y = savgol_filter(y_raw, 51, 2) # window size 51, polynomial order 3
		x = range(len(y))
		maxi = argrelmax(y, order=200)[0]
		maxi_val = [y[n] for n in maxi]
		solvent_peak = maxi[maxi_val.index(max(maxi_val))]

		# The solvent peak is approximately this distance from the etyhlene peak
		strt, end = solvent_peak+1100 ,solvent_peak+1600
		trapz = np.trapz(y[strt:end]-y[strt])
		
		output_integrated = open(experiment_name+"/data.integrated.csv", "a")
		output_integrated.write(",".join([data_ls[0], str(trapz)]) + "\n") #write data with a newline
		output_integrated.close()

		# Draw and save the plot
		ax = sns.lineplot(y=y_raw, x=x)
		ax = sns.lineplot(y=y, x=x)
		plt.axvline(strt, color='blue')
		plt.axvline(end, color='blue')
		plt.savefig(experiment_name+"/plots/sample_"+str(ln_count)+".png")
		plt.close()

		ln_count += 1
		continue
	
	data_ls.append(dat)
	time_ls.append(time)
	time += 1