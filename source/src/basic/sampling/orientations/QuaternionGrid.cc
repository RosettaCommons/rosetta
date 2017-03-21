// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/sampling/orientations/QuaternionGrid.cc
/// @brief  provide an evenly spaced set of rotations
/// @author Will Sheffler <willsheffler@gmail.com>

// much of this comes from http://charles.karney.info/orientation/

// unit headers
#include <basic/sampling/orientations/QuaternionGrid.hh>

// basic headers
#include <basic/database/open.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// C++ headers
#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <limits>

namespace basic {
namespace sampling {
namespace orientations {

using namespace std;

QuaternionGridManager::QuaternionGridManager(){
	fill_metadata();
}

bool QuatDBMetatada_cmpsize  (QuatDBMetadata i,QuatDBMetadata j) { return i.N > j.N; }
bool QuatDBMetatada_cmpradius(QuatDBMetadata i,QuatDBMetadata j) { return i.radius > j.radius; }
bool QuatDBMetatada_cmpcover (QuatDBMetadata i,QuatDBMetadata j) { return i.cover  > j.cover ; }

void QuaternionGridManager::fill_metadata() {
	//name                    N          radius          cover          delta
	by_radius_.push_back(QuatDBMetadata("c48u1",                  24,        62.80 ,       1.57514,         0.70000  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48n9",                 216,        36.47 ,       2.89689,         0.26091  ) );//        7.00
	by_radius_.push_back(QuatDBMetadata("c48u27",                648,        20.83 ,       1.64091,         0.33582  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u83",               1992,        16.29 ,       2.42065,         0.25970  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u157",              3768,        14.49 ,       3.22614,         0.20710  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u181",              4344,        12.29 ,       2.27013,         0.19415  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48n309",              7416,         9.72 ,       1.91567,         0.15167  ) );//        1.86
	by_radius_.push_back(QuatDBMetadata("c48u519",             12456,         9.05 ,       2.60257,         0.13807  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48n527",             12648,         8.17 ,       1.94334,         0.12599  ) );//        1.86
	by_radius_.push_back(QuatDBMetadata("c48u815",             19560,         7.40 ,       2.23719,         0.11607  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u1153",            27672,         6.60 ,       2.23735,         0.10330  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u1201",            28824,         6.48 ,       2.20918,         0.09999  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u1641",            39384,         5.75 ,       2.10646,         0.08993  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u2219",            53256,         5.27 ,       2.20117,         0.08249  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u2867",            68808,         5.24 ,       2.79649,         0.07531  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u2947",            70728,         4.71 ,       2.07843,         0.07359  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u3733",            89592,         4.37 ,       2.11197,         0.06836  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u4701",           112824,         4.22 ,       2.39041,         0.06372  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u4749",           113976,         4.00 ,       2.05300,         0.06248  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u5879",           141096,         3.74 ,       2.07325,         0.05837  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u7111",           170664,         3.53 ,       2.11481,         0.05514  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u8649",           207576,         3.26 ,       2.02898,         0.05094  ) );//        0.00
	by_radius_.push_back(QuatDBMetadata("c48u10305",          247320,         3.102,       2.08130,         0.048456 ) );//        2/41.2973
	by_radius_.push_back(QuatDBMetadata("c48u12083",          289992,         3.096,       2.42678,         0.046023 ) );//
	by_radius_.push_back(QuatDBMetadata("c48u12251",          294024,         2.903,       2.02950,         0.045354 ) );//        2s/18.2657
	by_radius_.push_back(QuatDBMetadata("c48u14251",          342024,         2.767,       2.04269,         0.043215 ) );//        2/46.2973
	by_radius_.push_back(QuatDBMetadata("c48u16533",          396792,         2.655,       2.09385,         0.041421 ) );//        2/48.2973
	by_radius_.push_back(QuatDBMetadata("c48u19181",          460344,         2.497,       2.02149,         0.039000 ) );//        2s/21.2450
	by_radius_.push_back(QuatDBMetadata("c48u21863",          524712,         2.403,       2.05419,         0.037534 ) );//        2/53.2973
	by_radius_.push_back(QuatDBMetadata("c48u25039",          600936,         2.282,       2.01458,         0.035641 ) );//        2s/23.2450
	by_radius_.push_back(QuatDBMetadata("c48u28329",          679896,         2.197,       2.03407,         0.034313 ) );//        2/58.2973
	by_radius_.push_back(QuatDBMetadata("c48u31793",          763032,         2.162,       2.17361,         0.033137 ) );//
	by_radius_.push_back(QuatDBMetadata("c48u32081",          769944,         2.116,       2.05852,         0.032786 ) );//
	by_radius_.push_back(QuatDBMetadata("c48u35851",          860424,         2.024,       2.01113,         0.031601 ) );//        2/63.2973
	by_radius_.push_back(QuatDBMetadata("c48u40003",          960072,         1.962,       2.04420,         0.030633 ) );//        2/65.2973
	by_radius_.push_back(QuatDBMetadata("c48u44709",         1073016,         1.877,       2.00081,         0.029307 ) );//        2s/28.2657
	by_radius_.push_back(QuatDBMetadata("c48u49397",         1185528,         1.822,       2.02304,         0.028453 ) );//        2/70.2973
	by_radius_.push_back(QuatDBMetadata("c48u54799",         1315176,         1.753,       1.99776,         0.027370 ) );//        2s/30.2657
	by_radius_.push_back(QuatDBMetadata("c48u60279",         1446696,         1.701,       2.00892,         0.026563 ) );//        2/75.2973
	by_radius_.push_back(QuatDBMetadata("c48u65985",         1583640,         1.657,       2.03291,         0.025876 ) );//        2/77.2973
	by_radius_.push_back(QuatDBMetadata("c48u72521",         1740504,         1.596,       1.99529,         0.024918 ) );//        2s/33.2450
	by_radius_.push_back(QuatDBMetadata("c48u79099",         1898376,         1.557,       2.01914,         0.024303 ) );//        2/82.2973
	by_radius_.push_back(QuatDBMetadata("c48u86451",         2074824,         1.505,       1.99648,         0.023504 ) );//        2s/35.2450
	by_radius_.push_back(QuatDBMetadata("c48u93701",         2248824,         1.467,       2.00411,         0.022911 ) );//        2/87.2973
	by_radius_.push_back(QuatDBMetadata("c48u101477",        2435448,         1.447,       2.07920,         0.022389 ) );//
	by_radius_.push_back(QuatDBMetadata("c48u101917",        2446008,         1.444,       2.07768,         0.022222 ) );//
	by_radius_.push_back(QuatDBMetadata("c48u110143",        2643432,         1.388,       1.99316,         0.021669 ) );//        2/92.2973
	by_radius_.push_back(QuatDBMetadata("c48u118647",        2847528,         1.358,       2.01352,         0.021210 ) );//        2/94.2973
	by_radius_.push_back(QuatDBMetadata("c48u128249",        3077976,         1.318,       1.98655,         0.020574 ) );//        2s/40.2657
	by_radius_.push_back(QuatDBMetadata("c48u137809",        3307416,         1.290,       2.00301,         0.020142 ) );//        2/99.2973
	by_radius_.push_back(QuatDBMetadata("c48u148395",        3561480,         1.255,       1.98744,         0.019600 ) );//        2s/42.2657
	by_radius_.push_back(QuatDBMetadata("c48u158763",        3810312,         1.228,       1.99130,         0.019176 ) );//        2/104.2973
	by_radius_.push_back(QuatDBMetadata("c48u169757",        4074168,         1.205,       2.01122,         0.018815 ) );//        2/106.2973
	by_radius_.push_back(QuatDBMetadata("c48u181909",        4365816,         1.173,       1.98631,         0.018310 ) );//        2s/45.2450
	by_radius_.push_back(QuatDBMetadata("c48u193767",        4650408,         1.151,       2.00013,         0.017970 ) );//        2/111.2973
	by_radius_.push_back(QuatDBMetadata("c48u207023",        4968552,         1.123,       1.98553,         0.017535 ) );//        2s/47.2450
	by_radius_.push_back(QuatDBMetadata("c48u220121",        5282904,         1.102,       1.99143,         0.017197 ) );//        2/116.2973
	by_radius_.push_back(QuatDBMetadata("c48u233569",        5605656,         1.083,       2.00765,         0.016906 ) );//        2/118.2973
	by_radius_.push_back(QuatDBMetadata("c48u248571",        5965704,         1.056,       1.98203,         0.016488 ) );//        2/121.2973
	by_radius_.push_back(QuatDBMetadata("c48u263339",        6320136,         1.039,       1.99944,         0.016221 ) );//        2/123.2973
	by_radius_.push_back(QuatDBMetadata("c48u279565",        6709560,         1.015,       1.98032,         0.015850 ) );//        2s/52.2657
	by_radius_.push_back(QuatDBMetadata("c48u295333",        7087992,         0.999,       1.99038,         0.015589 ) );//        2/128.2973
	by_radius_.push_back(QuatDBMetadata("c48u312831",        7507944,         0.978,       1.97997,         0.015266 ) );//        2s/54.2657
	by_radius_.push_back(QuatDBMetadata("c48u330023",        7920552,         0.961,       1.98309,         0.015004 ) );//        2/133.2973
	by_radius_.push_back(QuatDBMetadata("c48u347617",        8342808,         0.947,       1.99747,         0.014782 ) );//        2/135.2973
	by_radius_.push_back(QuatDBMetadata("c48u367113",        8810712,         0.927,       1.97956,         0.014472 ) );//        2s/57.2450
	by_radius_.push_back(QuatDBMetadata("c48u386211",        9269064,         0.913,       1.99027,         0.014255 ) );//        2/140.2973
	by_radius_.push_back(QuatDBMetadata("c48u407099",        9770376,         0.896,       1.98011,         0.013983 ) );//        2s/59.2450
	by_radius_.push_back(QuatDBMetadata("c48u427333",       10255992,         0.882,       1.98284,         0.013765 ) );//        2/145.2973
	by_radius_.push_back(QuatDBMetadata("c48u448437",       10762488,         0.870,       1.99711,         0.013578 ) );//        2/147.2973
	by_radius_.push_back(QuatDBMetadata("c48u471503",       11316072,         0.852,       1.97662,         0.013307 ) );//        2/150.2973
	by_radius_.push_back(QuatDBMetadata("c48u493799",       11851176,         0.841,       1.98949,         0.013132 ) );//        2/152.2973
	by_radius_.push_back(QuatDBMetadata("c48u518377",       12441048,         0.826,       1.97564,         0.012891 ) );//        2s/64.2657
	by_radius_.push_back(QuatDBMetadata("c48u542361",       13016664,         0.814,       1.98354,         0.012715 ) );//        2/157.2973
	by_radius_.push_back(QuatDBMetadata("c48u566819",       13603656,         0.814,       2.06674,         0.012551 ) );//
	by_radius_.push_back(QuatDBMetadata("c48u568499",       13643976,         0.807,       2.02337,         0.012499 ) );//
	by_radius_.push_back(QuatDBMetadata("c48u593755",       14250120,         0.789,       1.97681,         0.012323 ) );//        2/162.2973
	by_radius_.push_back(QuatDBMetadata("c48u619981",       14879544,         0.780,       1.98967,         0.012173 ) );//        2/164.2973
	by_radius_.push_back(QuatDBMetadata("c48u648549",       15565176,         0.766,       1.97599,         0.011964 ) );//        2s/69.2450
	by_radius_.push_back(QuatDBMetadata("c48u676103",       16226472,         0.757,       1.98293,         0.011813 ) );//        2/169.2973
	by_radius_.push_back(QuatDBMetadata("c48u706351",       16952424,         0.747,       1.98930,         0.011627 ) );//        2s/71.2450
	by_radius_.push_back(QuatDBMetadata("c48u735777",       17658648,         0.735,       1.97798,         0.011475 ) );//        2/174.2973
	by_radius_.push_back(QuatDBMetadata("c48u765729",       18377496,         0.727,       1.98881,         0.011344 ) );//        2/176.2973
	by_radius_.push_back(QuatDBMetadata("c48u798587",       19166088,         0.715,       1.97265,         0.011155 ) );//        2/179.2973
	by_radius_.push_back(QuatDBMetadata("c48u830491",       19931784,         0.707,       1.98390,         0.011032 ) );//        2/181.2973
	by_radius_.push_back(QuatDBMetadata("c48u865149",       20763576,         0.696,       1.97317,         0.010863 ) );//        2s/76.2657
	by_radius_.push_back(QuatDBMetadata("c48u898517",       21564408,         0.688,       1.97768,         0.010735 ) );//        2/186.2973
	by_radius_.push_back(QuatDBMetadata("c48u932999",       22391976,         0.683,       2.01538,         0.010620 ) );//
	by_radius_.push_back(QuatDBMetadata("c48u970447",       23290728,         0.670,       1.97320,         0.010455 ) );//        2/191.2973
	by_radius_.push_back(QuatDBMetadata("c48u1006449",      24154776,         0.663,       1.98364,         0.010347 ) );//        2/193.2973
	by_radius_.push_back(QuatDBMetadata("c48u1045817",      25099608,         0.653,       1.97289,         0.010197 ) );//        2s/81.2450
	by_radius_.push_back(QuatDBMetadata("c48u1083955",      26014920,         0.646,       1.97878,         0.010086 ) );//        2/198.2973
	by_size_ = by_radius_;
	by_cover_ = by_radius_;
	std::sort(by_size_  .begin(),by_size_  .end(),QuatDBMetatada_cmpsize  );
	std::sort(by_radius_.begin(),by_radius_.end(),QuatDBMetatada_cmpradius);
	std::sort(by_cover_ .begin(),by_cover_ .end(),QuatDBMetatada_cmpcover );
	// for(long i = 1; i <= (long)by_radius_.size(); ++i) cout << by_radius_[i].name << " " << by_radius_[i].N << " " << by_radius_[i].radius << " " << by_radius_[i].cover << endl;
	// cout << endl;
	// for(long i = 1; i <= (long)by_size_  .size(); ++i) cout << by_size_  [i].name << " " << by_size_  [i].N << " " << by_size_  [i].radius << " " << by_size_  [i].cover << endl;
	// cout << endl;
	// for(long i = 1; i <= (long)by_cover_ .size(); ++i) cout << by_cover_ [i].name << " " << by_cover_ [i].N << " " << by_cover_ [i].radius << " " << by_cover_ [i].cover << endl;
	// cout << endl;
}

QuaternionGridCOP QuaternionGridManager::request_by_name(std::string const & name){
	if ( grids_.find(name)==grids_.end() ) {
		utility::io::izstream infile;
		basic::database::open(infile,"sampling/orientations/orientgridall/data/"+name+".grid");
		grids_[name] = utility::pointer::shared_ptr<const class basic::sampling::orientations::QuaternionGrid>( utility::pointer::shared_ptr<const class basic::sampling::orientations::QuaternionGrid>( new QuaternionGrid(name,infile) ) );
		infile.close();
	}
	return grids_[name];
}

QuaternionGridCOP QuaternionGridManager::request_by_size(long target_size){
	for ( numeric::Size i = 1; i <= by_size_.size(); ++i ) {
		if ( target_size >= by_size_[i].N ) {
			return request_by_name(by_size_[i].name);
		}
	}
	return request_by_name(by_size_.back().name);
}
QuaternionGridCOP QuaternionGridManager::request_by_radius(numeric::Real target_radius){
	for ( numeric::Size i = 1; i <= by_radius_.size(); ++i ) {
		if ( target_radius >= by_radius_[i].radius ) {
			return request_by_name(by_radius_[i].name);
		}
	}
	return request_by_name(by_radius_.back().name);
}

// The rotational symmetries of the cube.  (Not normalized, since
// QuatSet.Add does this.)
static numeric::Real CubeSyms[24][4] = {
{1, 0, 0, 0},
// 180 deg rotations about 3 axes
{0, 1, 0, 0},
{0, 0, 1, 0},
{0, 0, 0, 1},
// +/- 120 degree rotations about 4 leading diagonals
{1, 1, 1, 1},
{1, 1, 1,-1},
{1, 1,-1, 1},
{1, 1,-1,-1},
{1,-1, 1, 1},
{1,-1, 1,-1},
{1,-1,-1, 1},
{1,-1,-1,-1},
// +/- 90 degree rotations about 3 axes
{1, 1, 0, 0},
{1,-1, 0, 0},
{1, 0, 1, 0},
{1, 0,-1, 0},
{1, 0, 0, 1},
{1, 0, 0,-1},
// 180 degree rotations about 6 face diagonals
{0, 1, 1, 0},
{0, 1,-1, 0},
{0, 1, 0, 1},
{0, 1, 0,-1},
{0, 0, 1, 1},
{0, 0, 1,-1},
};

// Convert from index to position.  The sinh scaling tries to compensate
// for the bunching up that occurs when [1 x y z] is projected onto the
// unit sphere.
numeric::Real pind(numeric::Real ind, numeric::Real delta, numeric::Real sigma) {
	return (sigma == 0) ? ind * delta : sinh(sigma * ind * delta) / sigma;
}

void Quaternion::Print(ostream& s) const {
	s << fixed << setprecision(9) << setw(12) << w << " ";
	s << setw(12) << x << " ";
	s << setw(12) << y << " ";
	s << setw(12) << z;
}

void QuatSet::Print(std::ostream& s, bool euler, size_t prec ) const {
	for ( size_t i = 0; i < Number(); ++i ) {
		if ( euler ) {
			m_v[i].PrintEuler(s);
		} else {
			m_v[i].Print(s);
		}
		s << " " << fixed << setprecision(prec) << setw(prec + 2) << m_w[i] << endl;
	}
}

void Quaternion::Normalize() {
	numeric::Real t = w*w + x*x + y*y + z*z;
	debug_assert(t > 0);
	t = 1.0/std::sqrt(t);
	w *= t;    x *= t;    y *= t;    z *= t;
	return;
}

numeric::xyzVector<numeric::Real> Quaternion::euler() const {
	// Print out orientation as a set of Euler angles, following the
	// convention given in
	//
	//    http://www.mhl.soton.ac.uk/research/help/Euler/index.html
	//
	// Rotation by Euler angles [a,b,c] is defined as rotation by c about
	// z axis, followed by rotation by b about y axis. followed by
	// rotation by a about z axis (again).
	//
	// Convert to rotation matrix (assume quaternion is already
	// normalized)
	numeric::Real
		// m00 = 1 - 2*y*y - 2*z*z,
		m01 =     2*x*y - 2*z*w,
		m02 =     2*x*z + 2*y*w,
		// m10 =     2*x*y + 2*z*w,
		m11 = 1 - 2*x*x - 2*z*z,
		m12 =     2*y*z - 2*x*w,
		m20 =     2*x*z - 2*y*w,
		m21 =     2*y*z + 2*x*w,
		m22 = 1 - 2*x*x - 2*y*y;
	// Taken from Ken Shoemake, "Euler Angle Conversion", Graphics Gems
	// IV, Academic 1994.
	//
	//    http://vered.rose.utoronto.ca/people/david_dir/GEMS/GEMS.html
	numeric::Real sy = sqrt(m02*m02 + m12*m12);
	//  numeric::Real sy = sqrt(m10*m10 + m20*m20);
	numeric::Real a,  b, c;
	b = atan2(sy, m22);
	if ( sy > 16 * std::numeric_limits<numeric::Real>::epsilon() ) {
		a = atan2(m12, m02);
		c = atan2(m21, -m20);
	} else {
		a = atan2(-m01, m11);
		c = 0;
	}
#if !defined(NDEBUG)
	// Sanity check.  Convert from Euler angles back to a quaternion, q
	Quaternion q = Quaternion(cos(a/2), 0, 0, sin(a/2)). // a about z
		Times(Quaternion(cos(b/2), 0, sin(b/2), 0). // b about y
		Times(Quaternion(cos(c/2), 0, 0, sin(c/2)))); // c about z
	// and check that q is parallel to *this.
	numeric::Real t = std::abs(q.w * w + q.x * x + q.y * y + q.z * z);
	debug_assert(t > 1 - 16 * numeric_limits<numeric::Real>::epsilon());
#endif
	return numeric::xyzVector<numeric::Real>(a,b,c);
}

numeric::xyzMatrix<numeric::Real> Quaternion::rotation_matrix() const {
	return numeric::xyzMatrix<numeric::Real>::cols(
		1.0 - 2.0*y*y - 2.0*z*z,       2.0*x*y - 2.0*z*w,       2.0*x*z + 2.0*y*w,
		2.0*x*y + 2.0*z*w, 1.0 - 2.0*x*x - 2.0*z*z,       2.0*y*z - 2.0*x*w,
		2.0*x*z - 2.0*y*w,       2.0*y*z + 2.0*x*w, 1.0 - 2.0*x*x - 2.0*y*y
	);
}


void Quaternion::PrintEuler(ostream& s) const {
	numeric::xyzVector<numeric::Real> e = euler();
	s << fixed << setprecision(9) << setw(12) << e.x() << " "
		<< setw(12) << e.y() << " " << setw(12) << e.z();
}

QuaternionGrid::QuaternionGrid(std::string const & name, std::istream & in): name_(name) {
	debug_assert(in.good());
	string line;
	while ( in.peek() == '#' ) {
		getline(in, line);
		// cout << line << endl;
	}
	debug_assert(in.good());
	getline(in, line);
	debug_assert(line == "format grid");
	in >> delta >> sigma >> ntot >> ncell >> nent >> maxrad_ >> coverage;
	// cout << ntot << " " << fixed
	// << setprecision(3) << maxrad_ << " "
	// << setprecision(5) << coverage << endl;
	size_t ncell1 = 0;
	for ( size_t n = 0; n < nent; ++n ) {
		long k, l, m;
		size_t mult;
		numeric::Real r, w;
		debug_assert(in.good());
		in >> k >> l >> m >> w >> r >> mult;
		Permute p(Triple(k, l, m));
		debug_assert(mult == p.Number());
		for ( size_t i = 0; i < mult; ++i ) {
			Triple t = p.Member(i);
			s.Add(Quaternion(1.0,
				pind(0.5 * t.a, delta, sigma),
				pind(0.5 * t.b, delta, sigma),
				pind(0.5 * t.c, delta, sigma)),
				w);
		}
		ncell1 += mult;
	}
	debug_assert(in.good());
	debug_assert(ncell1 == ncell);
	{
		size_t nc = s.Number();
		debug_assert(nc == ncell);
		for ( size_t n = 1; n < 24; ++n ) {
			Quaternion q(CubeSyms[n][0], CubeSyms[n][1], CubeSyms[n][2], CubeSyms[n][3] );
			for ( size_t i = 0; i < nc; ++i ) {
				s.Add(q.Times(s.Orientation(i)), s.Weight(i));
			}
		}
		debug_assert(s.Number() == ntot);
	}
}

std::ostream & operator<<(std::ostream & out, QuaternionGrid const & q){
	out << "QuaternionGrid " << q.name_ << ", nsamp: " << q.ntot << ", covering radius: " << q.maxrad_ << " degrees";
	return out;
}


void QuaternionGrid::print() const {

	debug_assert(s.Number() == ntot);
	// s.Print(cout, false, 7);
	for ( long i = 1; i <= num_samples(); ++i ) {
		quaternion(i).Print(cout);

		// cout << " " << quaternion(i).euler();
		// numeric::xyzMatrix<numeric::Real> rot = quaternion(i).rotation_matrix();
		// numeric::HomogeneousTransform<numeric::Real> ht(rot,numeric::xyzVector<numeric::Real>(0,0,0));
		// cout << endl << ht.euler_angles_rad() << endl << endl;

		cout << " " << weight(i) << endl;
	}

	cout << "DONE!!!!!" << endl;
}


}
}
}
