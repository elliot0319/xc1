import serial

ser = serial.Serial('/dev/ttyACM0',9600)

def serial_check():
	if(ser.in_waiting > 0):
		line = ser.readline()
		print("Arduino reports,", line)
