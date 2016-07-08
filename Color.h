/*
 * Color.h
 *
 *  Created on: 27.10.2014
 *      Author: StDoBudd
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
