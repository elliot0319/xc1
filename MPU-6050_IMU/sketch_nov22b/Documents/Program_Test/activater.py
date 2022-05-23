import voice_recognition as vr
import voice_output as vo
#import qualify as qua

def activate():
	command = vr.take_command()
	vo.talk(str(command))
	#qua.check(command)

