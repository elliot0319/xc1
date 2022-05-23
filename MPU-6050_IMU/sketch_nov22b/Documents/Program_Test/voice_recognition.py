import speech_recognition as sr

def take_command():
	try:
		with sr.Microphone() as source:
			sr.Recognizer().adjust_for_ambient_noise(source)
			print("Adjusted to ambient noise")
			print("Please Speak!")
			audio = sr.Recognizer().listen(source)
			print("Audio recorded")
		return("[Sphinx] You said" + sr.Recognizer().recognize_sphinx(audio))
	except:
		return("Error")
		print("Sphinx did not take the input.")
