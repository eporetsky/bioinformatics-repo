import serial
import matplotlib.pyplot as plt


#arduino_port = "/dev/cu.usbmodem1420" #serial port of Arduino
arduino_port = "/dev/cu.usbserial-1420" #serial port of Arduino
baud = 9600 #arduino uno runs at 9600 baud
fileName="analog-data.csv" #name of the CSV file generated

ser = serial.Serial(arduino_port, baud)
print("Connected to Arduino port:" + arduino_port)
output = open(fileName, "a")
print("Created file")
dsa = 0
ln_count = 1
data_ls = []
time_ls, time = [], 1

while True:
	dat = ser.readline().decode('ascii').split("\r")[0]

	if dat == "start":
		print("Starting: sample "+str(ln_count))
		data_ls = ["sample_"+str(ln_count)]
		time_ls, time = [], 1

		continue
		
	# Write to file 
	if dat == "stop":
		print("Complete")
		output = open(fileName, "a")
		output.write(",".join(data_ls) + "\n") #write data with a newline
		plt.plot(time_ls, data_ls[1:])
		plt.savefig("sample_"+str(ln_count)+".png")
		plt.close()

		ln_count += 1
		continue
	
	data_ls.append(dat)
	time_ls.append(time)
	time += 1