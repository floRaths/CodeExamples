#include <FastLED.h>
#include <Ramp.h>
#include <OneButton.h>

#include <Adafruit_DotStar.h>
Adafruit_DotStar strip(DOTSTAR_NUM, PIN_DOTSTAR_DATA, PIN_DOTSTAR_CLK, DOTSTAR_BRG);

// Pin definitions
#define LED_PIN 5
#define BTN_PIN 7

// Button Initialization
OneButton btn = OneButton(BTN_PIN, true, true);

// ################## matrix ###################
// define physical layout of LED matrix
const uint8_t kMatrixWidth  = 10;
const uint8_t kMatrixHeight = 13;

#define NUM_LEDS kMatrixWidth * kMatrixHeight
CRGB leds[NUM_LEDS];

boolean coil          = true;  // is it coiled or zig-zag
boolean flip          = false; // invert it
boolean ser_col       = true;  // serpentine layout?
boolean prototyping   = true;  // when prototyping flip row order
boolean reporting     = false; // report function output to serial monitor
boolean dataSmoothing = true;  // smooth noise data to reduce flicker

// ################## config ###################
uint8_t hurry = 6; // speed multiplier

// three brightness values to choose via button
uint8_t Bri1 = 255;
uint8_t Bri2 = 86;
uint8_t Bri3 = 86;

// prameters for initial palette selection
uint8_t base_hue1 = 30;  // first hue
uint8_t base_hue2 = 50; // second hue
uint8_t base_hue3 = base_hue2; // second hue
uint8_t range = 10;       // fluctuation
uint8_t sat_range = 55;
uint8_t bri_range = 155;

#include "declarations/declarations.h"
#include "auxFns/auxFns.h"

// #############################################
// ################## SETUP ####################
void setup() {
  
  // turn off onboard LED
  strip.begin(); 
  strip.setBrightness(0);
  strip.show(); // Turn all LEDs off ASAP

  delay(1000); // startup safety delay
  Serial.begin(115200);
  randomSeed(analogRead(0));

  FastLED.addLeds < WS2812B, LED_PIN, GRB > (leds, NUM_LEDS);
  FastLED.setBrightness(0);
  FastLED.setCorrection(TypicalLEDStrip);
  FastLED.setTemperature(Tungsten40W);

  btn.attachClick(brightnessAreaButton);
  btn.attachLongPressStart(paletteButton);
  btn.setPressMs(250);

  lumRampX.go(500, 0, LINEAR);
  lumRampY.go(500, 0, LINEAR);
  colRampX.go(500, 0, LINEAR);
  colRampY.go(500, 0, LINEAR);

  // random xy values for the noise field to ensure different starting points
  for (int i = 0; i < 4; i++)
  {
    noiRampMin[i] = 2000;
    noiRampMax[i] = 12000;
    xyVals[i]     = random(10000);
  }

  changeScales(10000);
  
  buildPalette(range, true, false, sat_range, bri_range);
  for (uint8_t i = 0; i < 4; i++)
  { runCol[i] = newCol[i];
    oldCol[i] = newCol[i];
  }

  bri_speed = 4500; lo_speed = 5000; up_speed = 5000;
  TargetBri = Bri1;
  brightnessAreaButton();

  Serial.println("Hello Lamp");
}



// #############################################
// ################## MAIN #####################
void loop() {

  buttonSwitches();
  rollColors();

  makeNoise();

  EVERY_N_MILLISECONDS(200)
  {
    if (indexDrift == true)
    {
        paletteIndex ++;
    }
  }

  EVERY_N_SECONDS(55)
  {
      changeScales(35000);
  }

  EVERY_N_SECONDS(34)
  {
    if (palRamp2.isFinished() == 1 && palette_changed == false)
    {
      if (rainbow)
      {
        newHues(30);
        range = random(5, 15);
      }

      buildPalette(range, false, triple, sat_range, bri_range);
      triggerRoll(30000);
    }
  }

  moveRange(lowerRamp.update(), upperRamp.update(), 8);
  fadeToBlackBy(leds, NUM_LEDS, 64);
  FastLED.setBrightness(briRamp.update());
  FastLED.show();
  btn.tick();
}
