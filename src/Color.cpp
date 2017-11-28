/*

Excited States software: KGS
Contributors: See CONTRIBUTORS.txt
Contact: kgs-contact@simtk.org

Copyright (C) 2009-2017 Stanford University

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

This entire text, including the above copyright notice and this permission notice
shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

*/



#include <iostream>
#include <sstream>
#include <limits.h>
#include <math.h>
#include <iomanip>


#include "Color.h"

using namespace std;

// initialize static variables
unsigned int Color::numberOfBins = 0;
bool Color::useSpectrum = false;

std::map<ColorNumber, SiteID> Color::mapFromColorNumberToMeanSiteID;
std::map<ColorNumber, SiteID> Color::mapFromColorNumberOtRMSDfromMeanSiteID;

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string Color::next_rcd_color_as_name( bool reset ){

  static int new_color = -1;

  if( reset ){
    new_color = -1;
    return("");
  }

  new_color++;

  return( getPymolColor(new_color) );
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
unsigned int Color::convertFractionToInt(float fractionValue, unsigned int maxValue) {

  unsigned int intValue = (unsigned int)floor(fractionValue*maxValue);

  return intValue;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string Color::getHexString(unsigned int n) {

  stringstream hexString;

  RGBintegerComponents rgbIntegerComponents = getRGBintegerComponents(n);

  hexString << setfill('0') << setbase(16)
	    << setw(2) << rgbIntegerComponents.red
	    << setw(2) << rgbIntegerComponents.green
	    << setw(2) << rgbIntegerComponents.blue;

  return hexString.str();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string Color::getPymolColor(unsigned int n) {

  stringstream pymolColor;

  pymolColor << "0x" << getHexString(n);

  return pymolColor.str();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void Color::setNumberOfBins(unsigned int value) {
  numberOfBins = value;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void Color::setUseSpectrum(bool value) {
  useSpectrum = value;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string Color::getJmolColor(unsigned int n) {
	stringstream jmolColor;

	jmolColor << "[x" << getHexString(n) << "]";


	return jmolColor.str();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string Color::getSVGColor(unsigned int n) {

  stringstream svgColor;

  svgColor << "#" << getHexString(n);

  return svgColor.str();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string Color::getPSColor( unsigned int n){

  stringstream psHSBColor;
  HSBfractionComponents hsbFractionComponents = getHSBfractionComponents( n );

  psHSBColor << showpoint << fixed << setprecision(4);

  psHSBColor << hsbFractionComponents.hue << " "
	     << hsbFractionComponents.saturation << " "
	     << hsbFractionComponents.brightness << " "
	     << "sethsbcolor" << endl;

  return( psHSBColor.str() );
}

////////////////////////////////////////////////////////////////////////////////
RGBfractionComponents Color::getRGBfractionComponents(ColorNumber colorNumber) {

  HSBfractionComponents hsbFractionComponents = getHSBfractionComponents(colorNumber);

  RGBfractionComponents rgbFractionComponents = convertFromHSBtoRGB(hsbFractionComponents);

  return rgbFractionComponents;
}

////////////////////////////////////////////////////////////////////////////////
RGBintegerComponents Color::getRGBintegerComponents(ColorNumber colorNumber) {

  RGBfractionComponents rgbFractionComponents = getRGBfractionComponents(colorNumber);

  RGBintegerComponents rgbIntegerComponents;

  rgbIntegerComponents.red   = convertFractionToInt(rgbFractionComponents.red,   255);
  rgbIntegerComponents.green = convertFractionToInt(rgbFractionComponents.green, 255);
  rgbIntegerComponents.blue  = convertFractionToInt(rgbFractionComponents.blue,  255);

  return rgbIntegerComponents;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is the definitive source for colors (all others are generated based on this)
////////////////////////////////////////////////////////////////////////////////
HSBfractionComponents Color::getHSBfractionComponents(ColorNumber colorNumber) {
  HSBfractionComponents hsbFractionComponents;

  if (useSpectrum) {
    if (colorNumber == UINT_MAX) {
      hsbFractionComponents.hue = 0.;
      hsbFractionComponents.brightness = 0.5;
      hsbFractionComponents.saturation = 0.;
    } else {

      if (numberOfBins == 0) {
        numberOfBins = 1;
      }

      SiteID meanSiteID = mapFromColorNumberToMeanSiteID[colorNumber];
      float fractionalMeanSiteID = (float)meanSiteID / (float)numberOfBins;

      SiteID rmsdFromMeanSiteID = mapFromColorNumberOtRMSDfromMeanSiteID[colorNumber];

      float fractionalRMSDfromMeanSiteID = (float)rmsdFromMeanSiteID / ((float)numberOfBins * 0.288675);

      float deltaFractionalMeanSiteID = (fractionalMeanSiteID - 0.5);
      float hue = deltaFractionalMeanSiteID + 2.0/3.0;

      hsbFractionComponents.hue = hue - floor(hue); // ensure that hue is in the interval [0,1)

      float saturation = pow(fractionalRMSDfromMeanSiteID, 0.5f);
      hsbFractionComponents.saturation = (saturation-floor(saturation));
      hsbFractionComponents.saturation = 1.0 - 0.5 * hsbFractionComponents.saturation;

      hsbFractionComponents.saturation = 1.0;

      float luminosity = sqrt((.6-.4*fractionalRMSDfromMeanSiteID) * (.5+.5*cos(10.*fractionalMeanSiteID)));
      luminosity = .55 - .4 * pow(fractionalRMSDfromMeanSiteID, 0.5f);

      //luminosity = .7 - .4 * (1.0 - pow(fractionalRMSDfromMeanSiteID, 0.5f));

      //luminosity = .2 + .4 * pow(fractionalRMSDfromMeanSiteID, 0.5f);

      //luminosity = .7 - .4 * pow(fractionalRMSDfromMeanSiteID, 0.5f);
      //hsbFractionComponents.saturation = .2+.8*(luminosity - .3)/.4;
      hsbFractionComponents.brightness = (luminosity - floor(luminosity));
    }

  } else {
    if (colorNumber == UINT_MAX) {
      hsbFractionComponents.hue = 0.;
      hsbFractionComponents.brightness = 0.5;
      hsbFractionComponents.saturation = 0.;

    } else {
      float x = (float)colorNumber;

      // simple irrational model to ensure that all hues are different
      float hue = 2./3. + 95.*x*sqrt(3./11.);

      hsbFractionComponents.hue = hue - floor(hue); // ensure that hue is in the interval [0,1)

      float saturation = 3593. * x * sqrt(7./23.);
      hsbFractionComponents.saturation = 1.-.7*(saturation-floor(saturation));

      float brightness = 2741.*x*sqrt(11./13.);
      hsbFractionComponents.brightness = .35+.4*(brightness - floor(brightness));
    }
  }

  return hsbFractionComponents;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Based on http://www.w3.org/TR/css3-color/
//
// NOTE / FIXME - this is actually HSL, not HSB; TODO - change all names accordingly
////////////////////////////////////////////////////////////////////////////////
RGBfractionComponents Color::convertFromHSBtoRGB(HSBfractionComponents hsbFractionComponents) {

  RGBfractionComponents rgbFractionComponents;

  float m2 = 0.0;
  if (hsbFractionComponents.brightness < 0.5) {
    m2 = hsbFractionComponents.brightness * (hsbFractionComponents.saturation +1.);

  } else {
    m2 = hsbFractionComponents.brightness +hsbFractionComponents.saturation - hsbFractionComponents.brightness * hsbFractionComponents.saturation;
  }

  float m1 = hsbFractionComponents.brightness * 2. - m2;

  rgbFractionComponents.red   = HUEtoRGB(m1, m2, hsbFractionComponents.hue + 1. / 3.);
  rgbFractionComponents.green = HUEtoRGB(m1, m2, hsbFractionComponents.hue);
  rgbFractionComponents.blue  = HUEtoRGB(m1, m2, hsbFractionComponents.hue - 1. / 3.);

  return rgbFractionComponents;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
float Color::HUEtoRGB(float m1, float m2, float hue) {

  float h = hue;

  if (h < 0.) {
    h += 1.;
  }

  if (h > 1.) {
    h -= 1.;
  }

  if (h * 6. < 1.) {
    return m1 + (m2 - m1) * h * 6.;
  }

  if (h * 2. < 1.) {
    return m2;
  }

  if (h * 3. < 2.) {
    return m1 + (m2 - m1) * (2. / 3. - h) * 6.;
  }

  return m1;
}

