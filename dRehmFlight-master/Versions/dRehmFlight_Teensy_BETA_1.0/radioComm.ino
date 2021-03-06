//Arduino/Teensy Flight Controller - dRehmFlight
//Author: Nicholas Rehm
//Project Start: 1/6/2020
//Version: Beta 1.0

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unsigned long rising_edge_start_1, rising_edge_start_2, rising_edge_start_3, rising_edge_start_4, rising_edge_start_5, rising_edge_start_6; 
unsigned long channel_1_raw, channel_2_raw, channel_3_raw, channel_4_raw, channel_5_raw, channel_6_raw;
int ppm_counter = 0;
unsigned long time_ms = 0;


void readPPM_setup(int pin) {
  // DESCRIPTION: Initialize software interrupts on radio channel pins
  // Declare interrupt pins
  pinMode(pin, INPUT_PULLUP);
  delay(20);
  //Attach interrupt and point to corresponding functions
  attachInterrupt(digitalPinToInterrupt(pin), getPPM, CHANGE);
}


void readPWM_setup(int ch1, int ch2, int ch3, int ch4, int ch5, int ch6) {
  //DESCRIPTION: Initialize software interrupts on radio channel pins
  //Declare interrupt pins 
  pinMode(ch1, INPUT_PULLUP);
  pinMode(ch2, INPUT_PULLUP);
  pinMode(ch3, INPUT_PULLUP);
  pinMode(ch4, INPUT_PULLUP);
  pinMode(ch5, INPUT_PULLUP);
  pinMode(ch6, INPUT_PULLUP);
  delay(20);
  //Attach interrupt and point to corresponding functions
  attachInterrupt(digitalPinToInterrupt(ch1), getCh1, CHANGE);
  attachInterrupt(digitalPinToInterrupt(ch2), getCh2, CHANGE);
  attachInterrupt(digitalPinToInterrupt(ch3), getCh3, CHANGE);
  attachInterrupt(digitalPinToInterrupt(ch4), getCh4, CHANGE);
  attachInterrupt(digitalPinToInterrupt(ch5), getCh5, CHANGE);
  attachInterrupt(digitalPinToInterrupt(ch6), getCh6, CHANGE);
  delay(20);
}

unsigned long getRadioPWM(int ch_num) {
  //DESCRIPTION: Get current radio commands from interrupt routines 
  unsigned long returnPWM = 0;
  
  if (ch_num == 1) {
    returnPWM = channel_1_raw;
  }
  else if (ch_num == 2) {
    returnPWM = channel_2_raw;
  }
  else if (ch_num == 3) {
    returnPWM = channel_3_raw;
  }
  else if (ch_num == 4) {
    returnPWM = channel_4_raw;
  }
  else if (ch_num == 5) {
    returnPWM = channel_5_raw;
  }
  else if (ch_num == 6) {
    returnPWM = channel_6_raw;
  }
  
  return returnPWM;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//INTERRUPT SERVICE ROUTINES 

void getPPM() {
  unsigned long dt_ppm;
  int trig = digitalRead(PPM_Pin);
  
  if (trig==1) { //only care about rising edge
    dt_ppm = micros() - time_ms;
    time_ms = micros();

    if (dt_ppm > 5000) { //ms
      ppm_counter = 0;
    }
  
    if (ppm_counter == 1) {
      channel_2_raw = dt_ppm;
    }
  
    if (ppm_counter == 2) {
      channel_3_raw = dt_ppm;
    }
  
    if (ppm_counter == 3) {
      channel_1_raw = dt_ppm;
    }
  
    if (ppm_counter == 4) {
      channel_4_raw = dt_ppm;
    }
  
    if (ppm_counter == 5) {
      channel_6_raw = dt_ppm;
    }
  
    if (ppm_counter == 6) {
      channel_5_raw = dt_ppm;
    }
    
    ppm_counter = ppm_counter + 1;
  }
}

void getCh1() {
  int trigger = digitalRead(ch1Pin);
  if(trigger == 1) {
    rising_edge_start_1 = micros();
  }
  else if(trigger == 0) {
    channel_1_raw = micros() - rising_edge_start_1;
  }
}

void getCh2() {
  int trigger = digitalRead(ch2Pin);
  if(trigger == 1) {
    rising_edge_start_2 = micros();
  }
  else if(trigger == 0) {
    channel_2_raw = micros() - rising_edge_start_2;
  }
}

void getCh3() {
  int trigger = digitalRead(ch3Pin);
  if(trigger == 1) {
    rising_edge_start_3 = micros();
  }
  else if(trigger == 0) {
    channel_3_raw = micros() - rising_edge_start_3;
  }
}

void getCh4() {
  int trigger = digitalRead(ch4Pin);
  if(trigger == 1) {
    rising_edge_start_4 = micros();
  }
  else if(trigger == 0) {
    channel_4_raw = micros() - rising_edge_start_4;
  }
}

void getCh5() {
  int trigger = digitalRead(ch5Pin);
  if(trigger == 1) {
    rising_edge_start_5 = micros();
  }
  else if(trigger == 0) {
    channel_5_raw = micros() - rising_edge_start_5;
  }
}

void getCh6() {
  int trigger = digitalRead(ch6Pin);
  if(trigger == 1) {
    rising_edge_start_6 = micros();
  }
  else if(trigger == 0) {
    channel_6_raw = micros() - rising_edge_start_6;
  }
}
