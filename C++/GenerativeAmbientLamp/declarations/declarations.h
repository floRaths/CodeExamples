// ################## inits ###################
// booleans for conditional execution
boolean palette_changed = false;
boolean new_colors = false;
boolean pressed = false;
boolean rainbow = false;
boolean rolling = false;
boolean triple = false;
boolean indexDrift = false;

uint8_t paletteIndex, CurrentBri, TargetBri, stored_bri;
int speed1, speed2;
uint32_t xyVals[4];
int noiRampMin[4];
int noiRampMax[4];

CRGB runCol[4];  
CRGB newCol[4];  
CRGB newCol2[4]; // 
CRGB oldCol[4];  // 
CRGB tmpCol[4];  // 

// parameter for moving the lit area
uint16_t lower = 0;            // lower end of lights
uint16_t upper = NUM_LEDS; // upper end of lights
uint16_t up_speed, lo_speed, bri_speed;

// initializing smoothing functions
rampInt briRamp, palRamp1, palRamp2;             // smooth palette blending 1
ramp lowerRamp, upperRamp;                       // smooth area blending 1
rampLong lumRampX, lumRampY, colRampX, colRampY; // smooth luminance scale blending on X

// arrays holding the raw, intermediate and final noise data for colors and luminosity
uint8_t lumNoise[kMatrixHeight][kMatrixHeight];
uint8_t lumData[kMatrixHeight][kMatrixHeight];
uint8_t scaled_lumData[kMatrixHeight][kMatrixHeight];

uint8_t colNoise[kMatrixHeight][kMatrixHeight];
uint8_t colData[kMatrixHeight][kMatrixHeight];
//uint8_t scaled_colData[kMatrixHeight][kMatrixHeight];

int maxLumValue = 0; // Assume the first element is the maximum

uint8_t list[] = {0, 1, 2, 3};

uint8_t switchBrightness = 0;
uint8_t switchArea = 2;
uint8_t switchIndex = 4;