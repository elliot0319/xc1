import serial
import time

ser = serial.Serial('/dev/ttyACM0', 115200)

while(True):
	if(ser.in_waiting > 0):
		myData = ser.readline()
		print(myData)
