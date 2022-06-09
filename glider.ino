/*Arduino 101 Glider Stabilizer
*By: Kemal Ficici and Tarik Agcayazi
*
*Current code only sabilizes glider back to 0,0,0 degrees X,Y,Z.
*TO DO - TEST OUT PID WITH SERVOS ATTACHED TO GLIDER FLAPS
*TO DO - Fix PID mixing(line 130,131)
*TO DO - Maneuvers(Initial Dive, Turns, Spiral)
*TO DO - Navigate to GPS Coordinate, Dive till set altitude, release mechanism
*
*https://playground.arduino.cc/Code/PIDLibrary
*https://www.arduino.cc/en/Tutorial/Genuino101CurieIMUOrientationVisualiser
*/

#include <CurieIMU.h>
#include <MadgwickAHRS.h>
#include <PID_v1.h>
#include <Servo.h>


///////////////////////////////IMU Stuff///////////////////////////////////
Madgwick filter;
unsigned long microsPerReading, microsPrevious;
float accelScale, gyroScale;
float roll, pitch, heading;
///////////////////////////////IMU Stuff///////////////////////////////////



/////////////////////////////////PID Stuff////////////////////////////////////////
double Kp=1, Ki=1, Kd=0.07;  //PID Tuning Parameters (NOTE: idk if I tuned it properly)
//     Roll PID stuff
double rollSetpoint, rollInput, rollOutput;
PID rollPID(&rollInput, &rollOutput, &rollSetpoint, Kp, Ki, Kd, DIRECT);


//   PitchPID Stuff
double pitchSetpoint, pitchInput, pitchOutput;
PID pitchPID(&pitchInput, &pitchOutput, &pitchSetpoint, Kp, Ki, Kd, DIRECT);


//Specify the links and initial tuning parameters
/////////////////////////////////PID Stuff////////////////////////////////////////


Servo leftServo, rightServo; //Servos for flaps
int leftOutput, rightOutput; //Output for flaps

void setup(){
    Serial.begin(9600);
    ///////////////PID Setup/////////////
      rollInput = roll; //Sets the PID inputs to the gyro's roll and pitch axes
      pitchInput = pitch;
    
      rollSetpoint = 0; 
      pitchSetpoint = 0;
    
      rollPID.SetOutputLimits(-90, 90); //limits the PID loop to output numbers in between -90 and 90. I used values from -90 to 90 in order to better visualize the glider position(when flaps are flat, the PID output is 0)
      pitchPID.SetOutputLimits(-90,90);
    
      rollPID.SetMode(AUTOMATIC);  //Turns PID on
      pitchPID.SetMode(AUTOMATIC);
    ///////////////PID Setup/////////////
    
    
    ///////////////IMU Setup/////////////
      // start the IMU and filter
      CurieIMU.begin();
      CurieIMU.setGyroRate(25);
      CurieIMU.setAccelerometerRate(25);
      filter.begin(25);
    
      // Set the accelerometer range to 2G
      CurieIMU.setAccelerometerRange(2);
      // Set the gyroscope range to 250 degrees/second
      CurieIMU.setGyroRange(250);
    
      // initialize variables to pace updates to correct rate
      microsPerReading = 1000000 / 25;
      microsPrevious = micros();
    ///////////////IMU Setup/////////////
    
      leftServo.attach(9);
      rightServo.attach(10); 
}//setup()


void loop() {
  rollInput = roll; //sets PID inputs to Roll and Pitch
  pitchInput = pitch;
  rollPID.SetTunings(Kp, Ki, Kd); //PID Tunings
  pitchPID.SetTunings(Kp, Ki, Kd);

  int aix, aiy, aiz;  //IMU Stuff
  int gix, giy, giz;
  float ax, ay, az;
  float gx, gy, gz;

  unsigned long microsNow;

  // check if it's time to read data and update the filter
  microsNow = micros();
  if (microsNow - microsPrevious >= microsPerReading) {

    // read raw data from CurieIMU
    CurieIMU.readMotionSensor(aix, aiy, aiz, gix, giy, giz);

    // convert from raw data to gravity and degrees/second units
    ax = convertRawAcceleration(aix);
    ay = convertRawAcceleration(aiy);
    az = convertRawAcceleration(aiz);
    gx = convertRawGyro(gix);
    gy = convertRawGyro(giy);
    gz = convertRawGyro(giz);

    // update the filter, which computes orientation
    filter.updateIMU(gx, gy, gz, ax, ay, az);

    // print the heading, pitch and roll
    roll = filter.getRoll();
    pitch = filter.getPitch();
    heading = filter.getYaw();



    //Calculates PID
    rollPID.Compute();
    pitchPID.Compute();

    //Combines the outputs and adds 90 to make it work with the servos
    leftOutput = ((rollOutput + pitchOutput) / 2) + 90;
    rightOutput = ((rollOutput - pitchOutput) / 2) + 90;

    
    /*  ^^^^^^^^^^NOTE TO TARIK: PLEASE FIX MIXING PIDS ^^^^^^^^^^^^^^^
     *  One issue with the code is that the way I mix the Roll and Pitch PID loops
     * is not the best way to mix PID loops. For instance, if the pitchPID was 90, and 
     * rollPID was 0, both flaps would not trim up, rather, rightServo would be at -45,
     * while leftServo trims at 135.  
     * 
     * ^^^^^^^^^^^NOTE TO TARIK: TEST OUT CODE WITH SERVOS ATTACHED TO GLIDER FLAPS^^^^^^^^^^^^
     */

     
    //writes to the servos
    leftServo.write(leftOutput);
    rightServo.write(rightOutput);
    Serial.println(pitch);
    // increment previous time, so we keep proper pace
    microsPrevious = microsPrevious + microsPerReading;
  } //if
    
} //loop()

float convertRawAcceleration(int aRaw) {
  // since we are using 2G range
  // -2g maps to a raw value of -32768
  // +2g maps to a raw value of 32767
  
  float a = (aRaw * 2.0) / 32768.0;
  return a;
}

float convertRawGyro(int gRaw) {
  // since we are using 250 degrees/seconds range
  // -250 maps to a raw value of -32768
  // +250 maps to a raw value of 32767
  
  float g = (gRaw * 250.0) / 32768.0;
  return g;
}
  
