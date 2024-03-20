void makeNoise()
{
    memset(lumNoise, 0, NUM_LEDS);
    fill_raw_2dnoise16into8(
        (uint8_t *)lumNoise,
        kMatrixHeight,     // width
        kMatrixHeight,     // height
        1,                 // octaves
        xyVals[0],         // x
        lumRampX.update(), // scalex
        xyVals[1],         // y
        lumRampY.update(), // scaley
        millis() * hurry   // timeVal
    );

    memset(colNoise, 0, NUM_LEDS);
    fill_raw_2dnoise16into8(
        (uint8_t *)colNoise,
        kMatrixHeight,       // width
        kMatrixHeight,       // height
        1,                   // octaves
        xyVals[2],           // x
        colRampX.update(),   // scalex
        xyVals[3],           // y
        colRampY.update(),   // scalex
        millis() * 4 // timeVal
    );

    CRGBPalette16 runPal = CRGBPalette16(runCol[0], runCol[1], runCol[2], runCol[3]);

    for (int x = 0; x < kMatrixWidth; x++)
    {
        for (int y = 0; y < kMatrixHeight; y++)
        {
            // Find the greatest brightness value in the frame and store it
            if (lumNoise[x][y] > maxLumValue)
            { maxLumValue = lumNoise[x][y]; }

            // upscale the data realtive to the max value
            scaled_lumData[x][y] = qadd8(lumNoise[x][y], scale8(lumNoise[x][y], maxLumValue));

            if (dataSmoothing)
            {
                uint8_t  old_lumData = lumData[x][y];
                uint16_t new_lumData = (uint16_t)old_lumData * (256 - 128) + (uint16_t)scaled_lumData[x][y] * 128;
                lumData[x][y] = new_lumData / 256; // Convert back to uint8_t

                uint8_t old_colData = colData[x][y];
                uint16_t new_colData = (uint16_t)old_colData * (256 - 128) + (uint16_t)colNoise[x][y] * 128;
                colData[x][y] = new_colData / 256; // Convert back to uint8_t
            }
            else
            {
                lumData[x][y] = scaled_lumData[x][y];
                colData[x][y] = colNoise[x][y];
            }

            leds[mtx(x, y)] = ColorFromPalette(runPal,
                                               // noiseCols[(y * kMatrixWidth) + x], // when used with 1D colors
                                               colData[x][y] + paletteIndex,
                                               brighten8_lin(lumData[x][y])
                                               );
        }
    }
}