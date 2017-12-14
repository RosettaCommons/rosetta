// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/color_util.cc
/// @author Sam DeLuca
/// @details RGB <-> HSV conversion utilities based off of http://www.cs.rit.edu/~ncs/color/t_convert.html

#include <numeric/xyzVector.hh>
#include <numeric/color_util.hh>

namespace numeric {

numeric::xyzVector<platform::Real> rgb_to_hsv(platform::Real r, platform::Real g, platform::Real b)
{
	numeric::xyzVector<platform::Real> rgb_values(r,g,b);

	return rgb_to_hsv(rgb_values);
}

numeric::xyzVector<platform::Real> rgb_to_hsv(numeric::xyzVector<platform::Real> rgb_triplet)
{
	//x = Hue
	//y = Saturation
	//z = Value
	numeric::xyzVector<platform::Real> hsv_values;
	platform::Real rgb_min = rgb_triplet.minimum_value();
	platform::Real rgb_max = rgb_triplet.maximum_value();

	platform::Real red = rgb_triplet.x();
	platform::Real green = rgb_triplet.y();
	platform::Real blue = rgb_triplet.z();
	//assert( red < 0 || red > 1 || green < 0 || green > 1 || blue < 0 || blue > 1 );

	//Set the Value
	hsv_values.z(rgb_max);

	platform::Real delta = rgb_max - rgb_min;

	//Set the saturation
	if ( rgb_max != 0 ) {
		hsv_values.y(delta/rgb_max);
	} else {
		//if rgb_max == 0, the object is black.  set saturation to 0
		hsv_values.y(0.0);
	}


	platform::Real hue = 0.0;
	//if the object is black, set hue to zero
	//if the object is grey, delta=0 and hue is undefined, so set it to 0
	if ( rgb_max == 0.0 || delta == 0.0 ) {
		hue = 0.0;
	} else if ( red == rgb_max ) {
		hue = (green - blue) / delta;
	} else if ( green == rgb_max ) {
		hue = 2.0+(blue - red) / delta;
	} else {
		hue = 4.0 + (red - green) / delta;
	}

	//convert hue to 0-360
	hue *= 60;
	if ( hue < 0 ) {
		hue += 360;
	}

	//set the hue
	hsv_values.x(hue);

	return hsv_values;
}

numeric::xyzVector<platform::Real> hsv_to_rgb(platform::Real h, platform::Real s, platform::Real v)
{
	numeric::xyzVector<platform::Real> hsv_values(h,s,v);

	return hsv_to_rgb(hsv_values);
}

numeric::xyzVector<platform::Real> hsv_to_rgb(numeric::xyzVector<platform::Real> hsv_triplet)
{
	numeric::xyzVector<platform::Real> rgb_values;

	platform::Real hue = hsv_triplet.x();
	platform::Real saturation = hsv_triplet.y();
	platform::Real value = hsv_triplet.z();

	//assert( hue < 0 || hue > 360 || saturation < 0 || saturation > 1 || value < 0 || value > 1 );

	//special case for grey colors
	if ( saturation == 0.0 ) {
		rgb_values.x(value);
		rgb_values.y(value);
		rgb_values.z(value);
	}

	// Look at http://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV for an explanation of whats going on here
	hue /= 60;
	auto i = static_cast<platform::Size>(std::floor(hue));
	platform::Real f = hue - i;
	platform::Real p = value * ( 1 - saturation );
	platform::Real q = value * ( 1 - saturation * f );
	platform::Real t = value * ( 1 - saturation * ( 1 - f ) );

	switch(i)
			{
			case 0 :
				rgb_values.x(value);
				rgb_values.y(t);
				rgb_values.z(p);
				break;
			case 1 :
				rgb_values.x(q);
				rgb_values.y(value);
				rgb_values.z(p);
				break;
			case 2 :
				rgb_values.x(p);
				rgb_values.y(value);
				rgb_values.z(t);
				break;
			case 3 :
				rgb_values.x(p);
				rgb_values.y(q);
				rgb_values.z(value);
				break;
			case 4 :
				rgb_values.x(t);
				rgb_values.y(p);
				rgb_values.z(value);
				break;
			default :
				rgb_values.x(value);
				rgb_values.y(p);
				rgb_values.z(q);
				break;
			}

	return rgb_values;

}


}

