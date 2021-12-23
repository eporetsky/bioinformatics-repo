#include <Wire.h>
#include <Adafruit_ADS1015.h>

Adafruit_ADS1115 ads;  /* Use this for the 16-bit version */
const int BUTTON_PIN = 2; // the number of the pushbutton pin

int var = 0;
//int last_state= HIGH;
int current_state=LOW; 

// Adafruit_ADS1015 ads;     /* Use thi for the 12-bit version */

void setup(void) 
{
  pinMode(BUTTON_PIN, INPUT);
  Serial.begin(115200);
  //Serial.begin(9600);
  //Serial.println("Hello!");
  
  //Serial.println("Getting single-ended readings from AIN0..3");
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

void loop(){
  //Serial.println(current_state);
  current_state = digitalRead(BUTTON_PIN);
  int16_t adc0; //, adc1, adc2, adc3;
  //Serial.println(current_state);
  if (current_state == HIGH){    
      var = 0;
      while(var < 100){                      //A while loop to flash the LED2 on and off
        adc0 = ads.readADC_SingleEnded(0);
        Serial.println(adc0);
        //Serial.println("The state changed from LOW to HIGH");
        var++;
        delay(20);
     }
  //last_state = current_state;

      //delay(120000);                        //A two minute delay before the button can be pressed again
  }
}
  
  // adc_integrated = 0;
  //adc1 = ads.readADC_SingleEnded(1);
  //adc2 = ads.readADC_SingleEnded(2);
  //adc3 = ads.readADC_SingleEnded(3);
  //Serial.print("AIN0: "); 
  
  //Serial.println(analogRead(A0));
  //Serial.print("AIN1: "); Serial.println(adc1);
  //Serial.print("AIN2: "); Serial.println(adc2);
  //Serial.print("AIN3: "); Serial.println(adc3);
  //Serial.println(" ");
//}
