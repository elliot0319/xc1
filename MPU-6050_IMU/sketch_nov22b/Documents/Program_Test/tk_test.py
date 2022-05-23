from tkinter import *
import sys
import activater

ver = str(sys.version)

root = Tk()
root.title("ISAC")
root.geometry("320x320")

# This is where I put the doge pic
#intro_photo = PhotoImage(file="/home/pi/Documents/Program_Test/download.png")
#label1 = Label(root, image=Intro_photo)
#label1.pack()

label2 = Label(root, text=ver)
label2.pack()

def btncmd1():
	activater.activate()
btn1 = Button(root, text="PTT", command=btncmd1)
btn1.pack()

def btncmd2():
	root.quit()
btn2 = Button(root, text="Exit", command=btncmd2)
btn2.pack()

root.mainloop()
