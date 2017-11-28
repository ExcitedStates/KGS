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


#ifndef COLOR_H_
#define COLOR_H_
#include <string>
#include <map>
#include "SiteID.h"


static const unsigned int MIN_CLUSTER_SIZE = 1;

struct RGBfractionComponents{float red; float green; float blue;};
struct RGBintegerComponents{int red; int green; int blue;};

struct HSBfractionComponents {
  float hue;
  float saturation;
  float brightness;
};

typedef unsigned int ColorNumber;

class Color {
 public:
	static std::string next_rcd_color_as_name( bool reset=false );

	static std::string getHexString(unsigned int n);

	static std::string getPymolColor(unsigned int n);
	static std::string getJmolColor(unsigned int n);
	static std::string getSVGColor(unsigned int n);
	static std::string getPSColor(unsigned int n);

	static RGBfractionComponents getRGBfractionComponents(ColorNumber colorNumber);
	static RGBintegerComponents getRGBintegerComponents(ColorNumber colorNumber);

	static HSBfractionComponents getHSBfractionComponents(ColorNumber colorNumber);

	static unsigned int convertFractionToInt(float fractionValue, unsigned int maxValue);

	static RGBfractionComponents convertFromHSBtoRGB(HSBfractionComponents hsbFractionComponents);

	static float HUEtoRGB(float m1, float m2, float hue);

	static void setUseSpectrum(bool);
	static void setNumberOfBins(unsigned int);

	static std::map<ColorNumber, SiteID> mapFromColorNumberToMeanSiteID;
	static std::map<ColorNumber, SiteID> mapFromColorNumberOtRMSDfromMeanSiteID;

 private:
	static bool useSpectrum;
	static unsigned int numberOfBins;
};

#endif /* COLOR_H_ */
