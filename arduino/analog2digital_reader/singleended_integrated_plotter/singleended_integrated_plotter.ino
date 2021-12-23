#include <Wire.h>
#include <Adafruit_ADS1015.h>
int var = 0;
int integrated = 0;
int line_num = 0;

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
  Serial.begin(115200);
  //Serial.begin(9600);
  
  // initialize the LED pin as an output:
  pinMode(ledPin, OUTPUT);
  // initialize the pushbutton pin as an input:
  pinMode(buttonPin, INPUT);
  
  // Serial.println("Getting single-ended readings from AIN0..3");
  // Serial.println("ADC Range: +/- 6.144V (1 bit = 3mV/ADS1015, 0.1875mV/ADS1115)");
  
  // The ADC input range (or gain) can be changed via the following
  // functions, but be careful never to exceed VDD +0.3V max, or to
  // exceed the upper and lower limits if you adjust the input range!
  // Setting these values incorrectly may destroy your ADC!
  //                                                                ADS1015  ADS1115
  //                                                                -------  -------
  // ads.setGain(GAIN_TWOTHIRDS);  // 2/3x gain +/- 6.144V  1 bit = 3mV      0.1875mV (default)
  // ads.setGain(GAIN_ONE);        // 1x gain   +/- 4.096V  1 bit = 2mV      0.125mV
  // ads.setGain(GAIN_TWO);        // 2x gain   +/- 2.048V  1 bit = 1mV      0.0625mV
  // ads.setGain(GAIN_FOUR);       // 4x gain   +/- 1.024V  1 bit = 0.5mV    0.03125mV
  // ads.setGain(GAIN_EIGHT);      // 8x gain   +/- 0.512V  1 bit = 0.25mV   0.015625mV
  ads.setGain(GAIN_SIXTEEN);    // 16x gain  +/- 0.256V  1 bit = 0.125mV  0.0078125mV
  
  ads.begin();
}

void loop(void) 
{

  buttonCurrentState = digitalRead(buttonPin);
  digitalWrite(ledPin, HIGH);
  // check if the pushbutton is pressed. If it is, the buttonState is HIGH:
  if (buttonCurrentState == LOW & buttonLastState == HIGH) {
    // turn LED on:
    
    Serial.println("");
    Serial.print(line_num);
    Serial.print(",");
    line_num++;
    buttonTimer = 5400;
    //buttonTimer = 600;
  }
  
  int16_t adc0; //, adc1, adc2, adc3;
  
  if (buttonTimer > 0){
    digitalWrite(ledPin, LOW);
    adc0 = ads.readADC_SingleEnded(0);
    
    //if(var == 10){
    if(var == 120){
        //if(integrated<0){integrated=1000;}
        //Serial.println(integrated); // Use of you want to plot, and comment next two lines
        Serial.print(integrated);
        Serial.print(",");
        var = 0;
        integrated = 0;
        //delay(20);
     }
    var++;
    integrated = integrated + adc0;    
    buttonTimer = buttonTimer - 1;
  }
 
 buttonLastState = buttonCurrentState;
}
