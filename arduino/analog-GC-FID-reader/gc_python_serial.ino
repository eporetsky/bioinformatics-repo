#include <Wire.h>
#include <Adafruit_ADS1015.h>
int var = 0;
int integrated = 0;

const int buttonPin = 2;     // the number of the pushbutton pin
const int ledPin =  8;      // the number of the LED pin
// variables will change:
int buttonCurrentState = 0;         // variable for reading the pushbutton status
int buttonLastState = 0;
int buttonTimer = 0;

Adafruit_ADS1115 ads;  /* Use this for the 16-bit version */
// Adafruit_ADS1015 ads;     /* Use thi for the 12-bit version */

void setup(void) 
{
  Serial.begin(9600);
  // initialize the LED pin as an output:
  pinMode(ledPin, OUTPUT);
  // initialize the pushbutton pin as an input:
  pinMode(buttonPin, INPUT);
  
  ads.setGain(GAIN_SIXTEEN);    // 16x gain  +/- 0.256V  1 bit = 0.125mV  0.0078125mV
  ads.begin();
}

void loop(void) 
{

  buttonCurrentState = digitalRead(buttonPin);

  // check if the pushbutton is pressed. If it is, the buttonState is HIGH:
  if (buttonCurrentState == LOW & buttonLastState == HIGH) {
    // turn LED on:
    digitalWrite(ledPin, HIGH);
    Serial.println("start");
    buttonTimer = 6000;
  }
  
  int16_t adc0; //, adc1, adc2, adc3;
  if (buttonTimer == 1){
    Serial.println("stop");
  }
  
  if (buttonTimer > 0){
    adc0 = ads.readADC_SingleEnded(0);
    Serial.println(adc0);
    buttonTimer = buttonTimer - 1;
  }
  

  digitalWrite(ledPin, LOW);
  buttonLastState = buttonCurrentState;
}
