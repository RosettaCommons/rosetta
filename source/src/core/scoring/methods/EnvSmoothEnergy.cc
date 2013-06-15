// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/EnvSmoothEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/EnvSmoothEnergy.hh>
#include <core/scoring/methods/EnvSmoothEnergyCreator.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Atom.hh>
// AUTO-REMOVED #include <core/scoring/EnvPairPotential.hh>
//#include <core/scoring/CachedData.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>



// Utility headers



// C++


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the EnvSmoothEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
EnvSmoothEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new EnvSmoothEnergy;
}

ScoreTypes
EnvSmoothEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( envsmooth );
	return sts;
}


//const core::Real start_sig =  9.9999;
//const core::Real end_sig   = 10.001;

Distance const start_sig = 9.8;
Distance const end_sig   = 10.2;

DistanceSquared const start_sig2 = start_sig*start_sig;
DistanceSquared const end_sig2   = end_sig*end_sig;


// This is temporary - this data will soon be read from a file in minirosetta_database
///
///core::Real envdata[20][40] = { { 0.215,  0.211,  0.199,  0.184,  0.168,  0.155,  0.152,  0.152,  0.156,  0.154,  0.161,  0.187,  0.222,  0.264,  0.290,  0.306,  0.290,  0.254,  0.191,  0.117,  0.037, -0.032, -0.098, -0.162, -0.241, -0.318, -0.393, -0.467, -0.547, -0.634, -0.721, -0.807, -0.894, -0.989, -1.089, -1.179, -1.241, -1.273, -1.282, -1.284 },
///{1.947,  1.943,  1.936,  1.933,  1.932,  1.937,  1.934,  1.916,  1.842,  1.699,  1.509,  1.313,  1.124,  0.940,  0.740,  0.561,  0.369,  0.196,  0.008, -0.147, -0.277, -0.371, -0.448, -0.510, -0.567, -0.606, -0.644, -0.669, -0.700, -0.721, -0.734, -0.744, -0.758, -0.810, -0.867, -0.933, -0.976, -1.027, -1.064, -1.088 },
///{-0.935, -0.923, -0.892, -0.847, -0.798, -0.752, -0.713, -0.668, -0.620, -0.564, -0.514, -0.466, -0.416, -0.351, -0.272, -0.172, -0.068,  0.043,  0.147,  0.264,  0.391,  0.526,  0.643,  0.741,  0.816,  0.863,  0.884,  0.886,  0.888,  0.900,  0.925,  0.951,  0.972,  0.983,  0.991,  0.997,  1.006,  1.013,  1.016,  1.016 },
///{-0.286, -0.309, -0.346, -0.398, -0.442, -0.506, -0.562, -0.624, -0.652, -0.651, -0.620, -0.582, -0.534, -0.469, -0.377, -0.262, -0.129,  0.021,  0.190,  0.372,  0.561,  0.736,  0.907,  1.060,  1.204,  1.310,  1.389,  1.440,  1.473,  1.492,  1.509,  1.540,  1.573,  1.602,  1.615,  1.621,  1.619,  1.618,  1.616,  1.617 },
///{ 0.760,  0.762,  0.772,  0.798,  0.863,  0.966,  1.085,  1.168,  1.192,  1.167,  1.104,  1.016,  0.896,  0.744,  0.566,  0.361,  0.166, -0.008, -0.149, -0.271, -0.371, -0.441, -0.476, -0.484, -0.470, -0.442, -0.400, -0.349, -0.277, -0.158,  0.009,  0.199,  0.353,  0.450,  0.497,  0.519,  0.533,  0.541,  0.548,  0.551 },
///{-0.041, -0.050, -0.056, -0.055, -0.037, -0.028, -0.029, -0.051, -0.076, -0.096, -0.104, -0.106, -0.098, -0.078, -0.040,  0.007,  0.049,  0.080,  0.099,  0.120,  0.134,  0.150,  0.153,  0.150,  0.126,  0.072, -0.006, -0.104, -0.214, -0.341, -0.464, -0.579, -0.678, -0.788, -0.888, -1.000, -1.133, -1.311, -1.489, -1.579 },
///{ 0.020,  0.018,  0.016,  0.013,  0.018,  0.044,  0.088,  0.140,  0.172,  0.181,  0.166,  0.137,  0.092,  0.033, -0.036, -0.104, -0.151, -0.179, -0.177, -0.159, -0.116, -0.062, -0.007,  0.040,  0.087,  0.128,  0.172,  0.209,  0.242,  0.263,  0.272,  0.277,  0.279,  0.282,  0.280,  0.275,  0.270,  0.271,  0.278,  0.284 },
///{ 1.167,  1.170,  1.178,  1.193,  1.230,  1.270,  1.315,  1.326,  1.309,  1.247,  1.151,  1.038,  0.899,  0.746,  0.577,  0.407,  0.244,  0.091, -0.041, -0.159, -0.262, -0.362, -0.449, -0.517, -0.551, -0.557, -0.543, -0.514, -0.468, -0.401, -0.331, -0.268, -0.190, -0.082,  0.046,  0.144,  0.200,  0.224,  0.238,  0.244 },
///{-0.026, -0.049, -0.092, -0.150, -0.214, -0.289, -0.359, -0.416, -0.456, -0.488, -0.510, -0.519, -0.508, -0.481, -0.438, -0.377, -0.282, -0.154,  0.018,  0.226,  0.454,  0.683,  0.908,  1.122,  1.328,  1.497,  1.638,  1.730,  1.800,  1.870,  1.969,  2.081,  2.168,  2.219,  2.244,  2.260,  2.266,  2.262,  2.254,  2.249 },
///{ 1.128,  1.132,  1.141,  1.156,  1.178,  1.212,  1.251,  1.270,  1.244,  1.167,  1.060,  0.936,  0.791,  0.624,  0.449,  0.281,  0.131, -0.007, -0.130, -0.242, -0.333, -0.406, -0.452, -0.476, -0.471, -0.444, -0.396, -0.323, -0.229, -0.120, -0.010,  0.108,  0.222,  0.374,  0.552,  0.783,  0.986,  1.139,  1.208,  1.232 },
///{ 0.756,  0.758,  0.766,  0.784,  0.807,  0.841,  0.874,  0.914,  0.910,  0.867,  0.766,  0.668,  0.557,  0.452,  0.333,  0.230,  0.133,  0.053, -0.028, -0.097, -0.169, -0.238, -0.317, -0.381, -0.423, -0.438, -0.443, -0.434, -0.423, -0.388, -0.387, -0.381, -0.412, -0.418, -0.438, -0.445, -0.454, -0.461, -0.466, -0.471 },
///{-0.510, -0.515, -0.528, -0.545, -0.569, -0.579, -0.568, -0.513, -0.457, -0.397, -0.361, -0.322, -0.290, -0.250, -0.206, -0.152, -0.083, -0.006,  0.081,  0.162,  0.246,  0.327,  0.404,  0.462,  0.502,  0.529,  0.548,  0.561,  0.567,  0.576,  0.582,  0.594,  0.606,  0.621,  0.634,  0.645,  0.653,  0.662,  0.670,  0.675 },
///{-1.240, -1.229, -1.194, -1.123, -1.022, -0.910, -0.810, -0.720, -0.620, -0.510, -0.394, -0.296, -0.204, -0.126, -0.054, -0.002,  0.042,  0.083,  0.137,  0.191,  0.247,  0.302,  0.359,  0.408,  0.430,  0.431,  0.409,  0.383,  0.354,  0.342,  0.336,  0.335,  0.338,  0.348,  0.360,  0.365,  0.370,  0.375,  0.380,  0.383 },
///{-0.072, -0.067, -0.061, -0.055, -0.050, -0.066, -0.103, -0.174, -0.242, -0.309, -0.354, -0.386, -0.390, -0.372, -0.337, -0.284, -0.220, -0.130, -0.028,  0.099,  0.227,  0.361,  0.488,  0.602,  0.704,  0.786,  0.852,  0.898,  0.947,  1.006,  1.083,  1.158,  1.240,  1.335,  1.446,  1.548,  1.633,  1.702,  1.766,  1.796 },
///{ 0.190,  0.188,  0.183,  0.177,  0.177,  0.173,  0.154,  0.111,  0.047, -0.032, -0.119, -0.200, -0.271, -0.324, -0.358, -0.363, -0.346, -0.297, -0.213, -0.097,  0.039,  0.189,  0.352,  0.519,  0.681,  0.831,  0.982,  1.116,  1.252,  1.380,  1.509,  1.606,  1.668,  1.699,  1.713,  1.719,  1.720,  1.720,  1.720,  1.720 },
///{-0.321, -0.323, -0.327, -0.338, -0.351, -0.365, -0.371, -0.367, -0.350, -0.321, -0.283, -0.237, -0.186, -0.129, -0.072, -0.019,  0.030,  0.076,  0.121,  0.163,  0.202,  0.233,  0.252,  0.258,  0.249,  0.224,  0.179,  0.118,  0.056,  0.000, -0.057, -0.112, -0.159, -0.184, -0.196, -0.203, -0.211, -0.216, -0.220, -0.221 },
///{ 0.073,  0.073,  0.073,  0.073,  0.073,  0.073,  0.066,  0.039, -0.004, -0.058, -0.102, -0.136, -0.149, -0.149, -0.139, -0.121, -0.087, -0.039,  0.020,  0.073,  0.124,  0.164,  0.197,  0.209,  0.199,  0.163,  0.118,  0.066,  0.015, -0.038, -0.083, -0.120, -0.147, -0.167, -0.182, -0.190, -0.197, -0.202, -0.206, -0.208 },
///{ 1.137,  1.138,  1.138,  1.142,  1.143,  1.148,  1.145,  1.133,  1.094,  1.027,  0.942,  0.838,  0.729,  0.608,  0.489,  0.367,  0.253,  0.143,  0.036, -0.078, -0.183, -0.277, -0.353, -0.422, -0.482, -0.532, -0.569, -0.589, -0.592, -0.588, -0.577, -0.573, -0.567, -0.577, -0.585, -0.599, -0.602, -0.606, -0.605, -0.605 },
///{ 1.430,  1.430,  1.430,  1.413,  1.374,  1.311,  1.240,  1.182,  1.134,  1.088,  1.022,  0.916,  0.773,  0.597,  0.407,  0.228,  0.051, -0.097, -0.239, -0.340, -0.410, -0.431, -0.428, -0.408, -0.381, -0.336, -0.263, -0.153, -0.014,  0.128,  0.247,  0.357,  0.474,  0.607,  0.729,  0.823,  0.877,  0.902,  0.912,  0.917 },
///{ 1.037,  1.040,  1.046,  1.051,  1.051,  1.051,  1.050,  1.047,  1.021,  0.973,  0.896,  0.804,  0.684,  0.534,  0.362,  0.182,  0.018, -0.133, -0.254, -0.343, -0.397, -0.422, -0.424, -0.396, -0.343, -0.263, -0.178, -0.068,  0.049,  0.200,  0.367,  0.550,  0.725,  0.874,  0.991,  1.067,  1.115,  1.142,  1.158,  1.165 } };


core::Real const envdata_coarse[20][40] = {
{ 0.86, 3.30, 3.30, 2.79, 1.75, 0.99, 0.91, 0.60, 0.46, 0.33, 0.19, 0.25, 0.20, 0.22, 0.18, 0.18, 0.12, 0.08,-0.03,-0.07,-0.18,-0.23,-0.26,-0.25,-0.29,-0.29,-0.28,-0.28,-0.29,-0.32,-0.35,-0.45,-0.59,-0.59,-0.74,-1.05,-0.71,-0.79,-0.71,-1.00},
{-0.28, 2.17, 2.17, 2.42, 2.20, 2.35, 1.99, 2.21, 2.12, 1.68, 1.48, 1.39, 1.26, 0.96, 0.84, 0.60, 0.40, 0.17,-0.00,-0.22,-0.37,-0.46,-0.53,-0.57,-0.61,-0.69,-0.69,-0.73,-0.72,-0.84,-0.68,-0.74,-0.67,-0.88,-0.93,-0.41,-1.15,-1.28,-1.28,-1.28},
{-0.50, 2.55, 2.55, 1.09,-0.12,-0.35,-0.42,-0.54,-0.48,-0.54,-0.52,-0.51,-0.42,-0.41,-0.31,-0.22,-0.11,-0.03, 0.12, 0.28, 0.39, 0.59, 0.65, 0.75, 0.93, 1.08, 1.20, 1.15, 1.18, 1.16, 1.07, 1.23, 1.26, 1.41, 1.85, 1.88, 1.02, 1.02, 1.02, 1.02},
{-0.17,-0.03,-0.03,-0.79,-0.83,-0.72,-0.87,-0.82,-0.81,-0.69,-0.65,-0.54,-0.45,-0.35,-0.25,-0.12, 0.07, 0.25, 0.40, 0.54, 0.77, 0.90, 1.10, 1.17, 1.35, 1.43, 1.50, 1.66, 1.72, 1.71, 1.68, 1.80, 1.86, 1.79, 1.53, 2.29, 2.06, 2.06, 1.01, 0.39},
{ 0.83, 0.33, 0.01, 0.05, 0.34, 0.49, 0.77, 0.88, 1.07, 1.10, 1.07, 0.99, 0.89, 0.76, 0.58, 0.45, 0.28, 0.19, 0.10,-0.03,-0.16,-0.27,-0.32,-0.44,-0.53,-0.58,-0.65,-0.74,-0.80,-0.78,-0.74,-0.68,-0.58,-0.32,-0.41,-0.25,-0.39, 0.73, 0.51, 0.51},
{ 0.36, 2.30, 2.30, 3.06, 2.20, 1.23, 0.75, 0.41, 0.19, 0.09, 0.00,-0.11,-0.11,-0.15,-0.16,-0.12,-0.11,-0.11,-0.12,-0.05,-0.07,-0.01, 0.02, 0.04, 0.12, 0.11, 0.09, 0.05, 0.03,-0.09,-0.13,-0.26,-0.33,-0.55,-0.62,-0.80,-0.83,-1.30,-0.94,-1.68},
{-0.42,-0.22,-0.38,-0.72,-0.47,-0.38,-0.18,-0.17,-0.07,-0.03,-0.04,-0.06,-0.15,-0.10,-0.15,-0.19,-0.20,-0.12,-0.01,-0.09, 0.02, 0.07, 0.20, 0.27, 0.34, 0.33, 0.37, 0.36, 0.31, 0.34, 0.23, 0.38, 0.38, 0.37, 0.17, 0.03, 0.60, 0.18, 0.18, 0.18},
{-0.63, 2.00, 2.00, 1.03, 0.94, 1.32, 1.35, 1.50, 1.30, 1.27, 1.06, 1.05, 1.02, 0.87, 0.74, 0.60, 0.48, 0.39, 0.18, 0.09,-0.09,-0.19,-0.28,-0.42,-0.55,-0.65,-0.78,-0.86,-0.97,-1.01,-1.07,-1.05,-1.03,-0.97,-1.02,-0.90,-0.84,-0.72,-0.76,-0.47},
{ 0.09,-0.92,-1.60,-1.63,-1.51,-1.35,-1.24,-1.07,-0.87,-0.70,-0.59,-0.50,-0.41,-0.29,-0.17, 0.01, 0.21, 0.44, 0.63, 0.87, 1.18, 1.31, 1.51, 1.75, 1.89, 2.04, 2.20, 2.14, 2.20, 2.22, 2.24, 1.96, 2.04, 2.29, 2.42, 2.14, 2.61, 2.61, 2.61, 2.61},
{ 0.38, 3.38, 3.38, 2.13, 1.48, 1.56, 1.44, 1.52, 1.28, 1.13, 1.10, 1.05, 0.85, 0.69, 0.50, 0.30, 0.15,-0.04,-0.12,-0.26,-0.30,-0.40,-0.51,-0.57,-0.61,-0.62,-0.60,-0.54,-0.43,-0.38,-0.26,-0.18,-0.09, 0.06, 0.41, 0.60, 0.01, 0.84, 0.62, 0.62},
{ 0.04, 0.24, 0.69, 0.58, 0.65, 0.52, 0.75, 0.75, 0.73, 0.80, 0.84, 0.80, 0.88, 0.61, 0.52, 0.40, 0.30, 0.18, 0.12, 0.03,-0.11,-0.14,-0.22,-0.34,-0.42,-0.52,-0.62,-0.68,-0.80,-0.85,-0.96,-0.90,-0.91,-0.89,-0.70,-0.79,-0.33,-0.46,-0.28,-0.28},
{ 0.48, 1.95, 1.95, 1.05,-0.21,-0.34,-0.35,-0.33,-0.33,-0.40,-0.34,-0.30,-0.31,-0.27,-0.21,-0.17,-0.14,-0.12, 0.00, 0.08, 0.16, 0.25, 0.38, 0.53, 0.60, 0.61, 0.75, 0.79, 0.78, 0.71, 0.71, 0.72, 0.85, 0.96, 1.09, 0.65, 0.93, 0.93, 0.57, 0.57},
{ 0.96, 2.02, 2.02, 1.58, 0.24,-0.18,-0.37,-0.46,-0.36,-0.33,-0.26,-0.20,-0.18,-0.13,-0.10,-0.11,-0.07,-0.04,-0.00,-0.00, 0.08, 0.14, 0.23, 0.27, 0.42, 0.45, 0.43, 0.48, 0.55, 0.63, 0.62, 0.55, 0.56, 0.73, 0.44, 0.54, 0.78, 0.87, 0.64, 0.64},
{-0.19, 0.46, 0.46,-0.20,-0.30,-0.25,-0.49,-0.45,-0.51,-0.50,-0.41,-0.36,-0.36,-0.33,-0.27,-0.21,-0.18,-0.09, 0.04, 0.16, 0.32, 0.51, 0.58, 0.73, 0.89, 0.86, 0.98, 1.06, 0.97, 1.06, 1.02, 0.86, 0.80, 0.65, 0.93, 1.40, 2.15, 0.63, 0.41, 0.41},
{-0.67,-2.65,-2.20,-1.66,-1.17,-0.94,-0.63,-0.49,-0.44,-0.42,-0.34,-0.37,-0.41,-0.39,-0.31,-0.25,-0.19,-0.01, 0.07, 0.25, 0.40, 0.59, 0.76, 0.90, 1.16, 1.34, 1.54, 1.70, 1.84, 1.75, 1.82, 2.14, 2.21, 2.24, 2.09, 2.41, 1.77, 1.77, 1.77, 1.77},
{-0.06, 1.39, 2.54, 2.80, 1.26, 0.45, 0.39, 0.11,-0.08,-0.12,-0.19,-0.22,-0.21,-0.17,-0.17,-0.15,-0.14,-0.12,-0.11,-0.03, 0.02, 0.04, 0.11, 0.21, 0.17, 0.26, 0.31, 0.30, 0.29, 0.33, 0.21, 0.28, 0.18, 0.06,-0.07, 0.31,-0.02,-0.40,-0.92, 0.25},
{-0.47, 2.89, 2.89, 2.86, 1.39, 0.90, 0.82, 0.57, 0.28, 0.11, 0.01,-0.14,-0.19,-0.19,-0.24,-0.27,-0.28,-0.25,-0.18,-0.13,-0.08,-0.03, 0.09, 0.10, 0.17, 0.22, 0.20, 0.21, 0.29, 0.24, 0.25, 0.23, 0.18, 0.10, 0.43, 0.36,-0.07,-0.05,-0.28, 0.20},
{-0.09, 3.86, 3.86, 4.12, 2.10, 1.67, 1.69, 1.62, 1.44, 1.20, 1.05, 0.85, 0.73, 0.54, 0.40, 0.22, 0.05,-0.06,-0.17,-0.32,-0.38,-0.47,-0.52,-0.52,-0.54,-0.50,-0.47,-0.42,-0.36,-0.28,-0.25,-0.26,-0.17,-0.33,-0.18, 0.11, 0.01,-0.06,-0.00,-0.62},
{-0.00,-0.00, 1.55, 0.42, 0.76, 0.85, 0.96, 0.97, 1.05, 1.05, 0.95, 0.87, 0.65, 0.47, 0.19, 0.05,-0.09,-0.19,-0.32,-0.41,-0.47,-0.41,-0.41,-0.46,-0.34,-0.36,-0.29,-0.22,-0.25,-0.10, 0.10, 0.28, 0.39, 0.66, 0.34, 1.86, 1.86, 1.86, 1.86, 1.86},
{ 1.40,-0.20, 0.01,-0.29, 0.24, 0.38, 0.55, 0.63, 0.70, 0.65, 0.58, 0.44, 0.30, 0.12,-0.01,-0.15,-0.18,-0.31,-0.33,-0.31,-0.28,-0.28,-0.26,-0.16,-0.14,-0.11,-0.05, 0.02,-0.02, 0.11, 0.16, 0.12, 0.17, 0.38, 0.45, 0.12, 0.74, 0.74, 0.38, 0.38} };


core::Real const envdata[20][40] = {
{ 3.46, 3.33, 3.05, 2.54, 1.89, 1.30, 0.90, 0.65, 0.48, 0.35, 0.27, 0.23, 0.21, 0.21, 0.18, 0.16, 0.11, 0.06,-0.01,-0.09,-0.16,-0.21,-0.25,-0.26,-0.28,-0.28,-0.28,-0.29,-0.30,-0.33,-0.39,-0.46,-0.55,-0.63,-0.70,-0.74,-0.76,-0.81,-0.89,-1.00 },
{ 2.60, 2.54, 2.47, 2.39, 2.29, 2.22, 2.18, 2.11, 1.97, 1.76, 1.55, 1.37, 1.20, 1.01, 0.81, 0.60, 0.40, 0.19,-0.01,-0.19,-0.33,-0.44,-0.51,-0.57,-0.62,-0.66,-0.69,-0.73,-0.74,-0.75,-0.73,-0.74,-0.76,-0.84,-0.93,-1.04,-1.14,-1.22,-1.27,-1.28 },
{-0.50,-0.47,-0.45,-0.43,-0.43,-0.43,-0.45,-0.48,-0.50,-0.52,-0.51,-0.48,-0.44,-0.38,-0.30,-0.22,-0.11, 0.00, 0.13, 0.27, 0.41, 0.54, 0.66, 0.79, 0.91, 1.02, 1.09, 1.13, 1.16, 1.18, 1.20, 1.24, 1.29, 1.35, 1.41, 1.46, 1.49, 1.50, 1.47, 1.43 },
{-0.81,-0.82,-0.82,-0.80,-0.80,-0.80,-0.81,-0.80,-0.77,-0.71,-0.63,-0.54,-0.45,-0.35,-0.23,-0.09, 0.07, 0.23, 0.40, 0.57, 0.74, 0.91, 1.06, 1.19, 1.32, 1.42, 1.53, 1.62, 1.68, 1.71, 1.74, 1.78, 1.82, 1.85, 1.88, 1.89, 1.88, 1.87, 1.86, 1.88 },
{ 0.37, 0.36, 0.36, 0.38, 0.45, 0.56, 0.72, 0.88, 1.00, 1.05, 1.04, 0.97, 0.87, 0.74, 0.59, 0.45, 0.31, 0.19, 0.08,-0.03,-0.14,-0.25,-0.34,-0.43,-0.51,-0.59,-0.66,-0.72,-0.76,-0.76,-0.72,-0.64,-0.54,-0.43,-0.37,-0.32,-0.30,-0.26,-0.25,-0.24 },
{ 2.37, 2.32, 2.26, 2.10, 1.78, 1.32, 0.86, 0.49, 0.26, 0.11, 0.00,-0.07,-0.11,-0.14,-0.14,-0.13,-0.12,-0.11,-0.10,-0.07,-0.05,-0.02, 0.02, 0.06, 0.09, 0.09, 0.08, 0.05, 0.00,-0.08,-0.15,-0.26,-0.37,-0.51,-0.64,-0.76,-0.84,-0.92,-0.99,-1.06 },
{-0.42,-0.35,-0.32,-0.32,-0.31,-0.28,-0.21,-0.15,-0.09,-0.06,-0.06,-0.08,-0.11,-0.13,-0.15,-0.17,-0.16,-0.14,-0.10,-0.05, 0.01, 0.09, 0.18, 0.25, 0.31, 0.34, 0.35, 0.35, 0.34, 0.34, 0.35, 0.36, 0.37, 0.36, 0.34, 0.33, 0.30, 0.26, 0.22, 0.18 },
{ 1.39, 1.37, 1.35, 1.33, 1.32, 1.32, 1.32, 1.31, 1.28, 1.23, 1.15, 1.07, 0.97, 0.86, 0.74, 0.61, 0.48, 0.35, 0.21, 0.07,-0.06,-0.18,-0.30,-0.42,-0.54,-0.65,-0.76,-0.86,-0.94,-1.00,-1.04,-1.04,-1.03,-1.00,-0.95,-0.90,-0.85,-0.79,-0.73,-0.64 },
{-1.55,-1.60,-1.61,-1.57,-1.49,-1.36,-1.22,-1.05,-0.89,-0.73,-0.61,-0.50,-0.40,-0.28,-0.14, 0.03, 0.22, 0.43, 0.66, 0.89, 1.12, 1.33, 1.52, 1.71, 1.87, 2.00, 2.08, 2.14, 2.18, 2.21, 2.22, 2.23, 2.27, 2.32, 2.41, 2.49, 2.56, 2.60, 2.61, 2.61 },
{ 2.77, 2.62, 2.40, 2.11, 1.83, 1.62, 1.48, 1.38, 1.29, 1.19, 1.10, 0.99, 0.85, 0.68, 0.50, 0.32, 0.15, 0.00,-0.12,-0.23,-0.32,-0.41,-0.49,-0.55,-0.59,-0.60,-0.57,-0.52,-0.44,-0.36,-0.27,-0.17,-0.05, 0.08, 0.23, 0.36, 0.47, 0.55, 0.60, 0.62 },
{ 0.99, 0.98, 0.96, 0.95, 0.94, 0.92, 0.91, 0.89, 0.87, 0.85, 0.82, 0.77, 0.70, 0.61, 0.51, 0.40, 0.30, 0.20, 0.11, 0.02,-0.07,-0.15,-0.24,-0.33,-0.42,-0.52,-0.61,-0.69,-0.78,-0.85,-0.90,-0.91,-0.90,-0.87,-0.81,-0.72,-0.60,-0.49,-0.40,-0.33 },
{-0.17,-0.21,-0.22,-0.24,-0.26,-0.30,-0.33,-0.34,-0.35,-0.35,-0.34,-0.32,-0.29,-0.26,-0.22,-0.18,-0.13,-0.08,-0.01, 0.08, 0.17, 0.27, 0.38, 0.49, 0.58, 0.65, 0.71, 0.75, 0.76, 0.74, 0.74, 0.75, 0.79, 0.81, 0.84, 0.86, 0.89, 0.89, 0.86, 0.83 },
{ 0.06, 0.04, 0.02,-0.02,-0.09,-0.20,-0.31,-0.37,-0.37,-0.32,-0.26,-0.22,-0.17,-0.14,-0.11,-0.09,-0.07,-0.04,-0.01, 0.03, 0.08, 0.15, 0.22, 0.30, 0.37, 0.42, 0.46, 0.50, 0.55, 0.59, 0.61, 0.61, 0.61, 0.61, 0.62, 0.63, 0.63, 0.64, 0.64, 0.64 },
{-0.19,-0.14,-0.15,-0.19,-0.27,-0.33,-0.41,-0.46,-0.48,-0.46,-0.42,-0.38,-0.35,-0.31,-0.27,-0.22,-0.15,-0.07, 0.04, 0.18, 0.32, 0.47, 0.60, 0.72, 0.82, 0.90, 0.96, 1.00, 1.02, 1.01, 0.96, 0.88, 0.79, 0.71, 0.66, 0.63, 0.60, 0.55, 0.48, 0.41 },
{-2.73,-2.48,-2.12,-1.70,-1.28,-0.95,-0.71,-0.55,-0.46,-0.41,-0.38,-0.38,-0.38,-0.36,-0.31,-0.24,-0.15,-0.03, 0.10, 0.25, 0.41, 0.58, 0.76, 0.94, 1.14, 1.34, 1.52, 1.66, 1.75, 1.82, 1.92, 2.05, 2.14, 2.21, 2.17, 2.11, 1.95, 1.84, 1.77, 1.77 },
{ 1.44, 1.28, 1.08, 0.86, 0.66, 0.48, 0.32, 0.14,-0.01,-0.11,-0.17,-0.19,-0.20,-0.18,-0.17,-0.15,-0.14,-0.12,-0.08,-0.04, 0.01, 0.06, 0.11, 0.17, 0.21, 0.25, 0.28, 0.30, 0.29, 0.29, 0.26, 0.22, 0.16, 0.09, 0.03,-0.01,-0.03,-0.05,-0.04,-0.02 },
{ 1.36, 1.32, 1.26, 1.18, 1.06, 0.92, 0.75, 0.55, 0.34, 0.15, 0.01,-0.10,-0.16,-0.20,-0.23,-0.25,-0.26,-0.23,-0.18,-0.13,-0.07,-0.01, 0.06, 0.11, 0.16, 0.19, 0.21, 0.23, 0.25, 0.25, 0.24, 0.21, 0.17, 0.13, 0.08, 0.03,-0.03,-0.07,-0.11,-0.14 },
{ 2.92, 2.80, 2.64, 2.39, 2.11, 1.85, 1.69, 1.55, 1.41, 1.23, 1.05, 0.87, 0.71, 0.55, 0.39, 0.23, 0.08,-0.06,-0.18,-0.29,-0.38,-0.45,-0.50,-0.52,-0.52,-0.50,-0.46,-0.41,-0.36,-0.30,-0.26,-0.25,-0.24,-0.22,-0.16,-0.10,-0.04,-0.03,-0.03,-0.05 },
{ 0.67, 0.69, 0.71, 0.73, 0.78, 0.85, 0.93, 0.98, 1.01, 1.00, 0.93, 0.81, 0.64, 0.45, 0.24, 0.07,-0.08,-0.19,-0.30,-0.38,-0.42,-0.43,-0.42,-0.41,-0.37,-0.34,-0.29,-0.24,-0.18,-0.06, 0.08, 0.22, 0.30, 0.34, 0.36, 0.37, 0.39, 0.41, 0.42, 0.43 },
{ 0.01, 0.01, 0.04, 0.11, 0.23, 0.38, 0.51, 0.60, 0.64, 0.62, 0.55, 0.43, 0.29, 0.14, 0.00,-0.11,-0.20,-0.27,-0.30,-0.30,-0.29,-0.27,-0.23,-0.19,-0.14,-0.09,-0.05,-0.01, 0.03, 0.08, 0.12, 0.17, 0.24, 0.32, 0.39, 0.42, 0.41, 0.40, 0.39, 0.38 } };

////////////////////////////////////////////////////////////////////////////


/// c-tor
EnvSmoothEnergy::EnvSmoothEnergy() :
	parent( new EnvSmoothEnergyCreator )
{}


/// clone
EnergyMethodOP
EnvSmoothEnergy::clone() const
{
	return new EnvSmoothEnergy();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

inline Real sqr ( Real x ) {
	return x*x;
}


/// @details stores dScore/dNumNeighbors so that when neighbor atoms on adjacent
/// residues move, their influence on the score of the surrounding residues is
/// rapidly computed.
void
EnvSmoothEnergy::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const &
) const
{
	pose.update_residue_neighbors();
	Size nres( pose.total_residue() );

	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
		dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

 	residue_N_.clear();
 	residue_E_.clear();
 	residue_dEdN_.clear();

	// iterate over all the residues in the protein and count their neighbours
	// and save values of E, N, and dEdN
	for ( Size i = 1; i <= nres; ++i ) {
		if( symm_info && !symm_info->bb_is_independent( i ) ) {
			residue_E_.push_back(0);
			residue_N_.push_back(0);
			residue_dEdN_.push_back(0);
			continue;
		}

		// get the appropriate residue from the pose.
		conformation::Residue const & rsd( pose.residue(i) );
		// currently this is only for protein residues
		if( !rsd.is_protein() || rsd.aa() == chemical::aa_unk ) {
			residue_E_.push_back(0);
			residue_N_.push_back(0);
			residue_dEdN_.push_back(0);
			continue; //return;
		}

		Size const atomindex_i = rsd.atom_index( representative_atom_name( rsd.aa() ));

		core::conformation::Atom const & atom_i = rsd.atom(atomindex_i);

		const Energies & energies( pose.energies() );
		const TwelveANeighborGraph & graph ( energies.twelveA_neighbor_graph() );

		Real countN    =  0.0;

		// iterate across neighbors within 12 angstroms
		for ( graph::Graph::EdgeListConstIter
				ir  = graph.get_node(i)->const_edge_list_begin(),
				ire = graph.get_node(i)->const_edge_list_end();
				ir != ire; ++ir ) {
			Size const j( (*ir)->get_other_ind( i ) );
			conformation::Residue const & rsd_j( pose.residue(j) );
			Size atomindex_j( rsd_j.type().nbr_atom() );

			core::conformation::Atom const & atom_j = rsd_j.atom(atomindex_j);

			Real sqdist = atom_i.xyz().distance_squared(atom_j.xyz());
			countN += sigmoidish_neighbor( sqdist );
		}

		Real score = 0;
		Real dscoredN = 0;

		calc_energy( countN, rsd.aa(), score, dscoredN );

		residue_N_.push_back( countN );
		residue_E_.push_back( score );
		residue_dEdN_.push_back( dscoredN );

		//std::cout << "ENV:  " << i << "  " << score << std::endl;
	}

	// symmetrize
	if( symm_info ) {
		for ( Size i = 1; i <= nres; ++i ) {
			if( !symm_info->bb_is_independent( i ) ) {
				Size master_i = symm_info->bb_follows( i );
				residue_N_[i] = residue_N_[master_i];
				residue_E_[i] = residue_E_[master_i];
				residue_dEdN_[i] = residue_dEdN_[master_i];
			}
		}
	}
}

void
EnvSmoothEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const
{
	pose.update_residue_neighbors();
}

/// @details counts the number of nbr atoms within a given radius of the for the input
/// residue.  Because the representative atom on the input residue may be in a different
/// location than the representative atom on the same residue when scoring_begin() is called,
/// these neighbor counts cannot be reused; therefore, scoring_begin does not keep
/// neighbor counts.
void
EnvSmoothEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ){
			return;
	}

	// currently this is only for protein residues
	if ( ! rsd.is_protein() ) return;
	if( rsd.aa() == chemical::aa_unk ) return;

	TwelveANeighborGraph const & graph ( pose.energies().twelveA_neighbor_graph() );
	Size const atomindex_i = rsd.atom_index( representative_atom_name( rsd.aa() ));

	core::conformation::Atom const & atom_i = rsd.atom(atomindex_i);

	Real countN    =  0.0;
	// iterate across neighbors within 12 angstroms
	for ( graph::Graph::EdgeListConstIter
			ir  = graph.get_node( rsd.seqpos() )->const_edge_list_begin(),
			ire = graph.get_node( rsd.seqpos() )->const_edge_list_end();
			ir != ire; ++ir ) {
		Size const j( (*ir)->get_other_ind( rsd.seqpos() ) );
		conformation::Residue const & rsd_j( pose.residue(j) );

		// if virtual residue, don't score
		if (rsd_j.aa() == core::chemical::aa_vrt) continue;

		Size atomindex_j( rsd_j.nbr_atom() );

		core::conformation::Atom const & atom_j = rsd_j.atom(atomindex_j);

		Real sqdist = atom_i.xyz().distance_squared( atom_j.xyz() );
		countN += sigmoidish_neighbor( sqdist );
	}

	Real score = 0, dscoredN = 0;

	calc_energy( countN, rsd.aa(), score, dscoredN );

	emap[ envsmooth ] += score;
}


/// @details Special cases handled for when an atom is both the representative
/// atom for an amino acid, and its nbr_atom.
void
EnvSmoothEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & ,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( pose.residue( atom_id.rsd() ).has_variant_type( core::chemical::REPLONLY )){
		return;
	}

	conformation::Residue const & rsd = pose.residue( atom_id.rsd() );

	if (! rsd.is_protein() ) return;
	if ( rsd.aa() == chemical::aa_unk ) return;

	Size const i = rsd.seqpos();
	Size const i_nbr_atom = rsd.type().nbr_atom();
	Size const i_rep_atom = rsd.atom_index( representative_atom_name( rsd.aa() ));

	// forces act only on the nbr atom (CB or CA) or the representative atom
	if( i_nbr_atom != (Size) atom_id.atomno() && i_rep_atom != (Size) atom_id.atomno() ) return;

	core::conformation::Atom const & atom_i = rsd.atom( atom_id.atomno() );

	TwelveANeighborGraph const & graph ( pose.energies().twelveA_neighbor_graph() );

	// its possible both of these are true
	bool const input_atom_is_nbr( i_nbr_atom == Size (atom_id.atomno()) );
	bool const input_atom_is_rep( i_rep_atom == Size ( atom_id.atomno() ));

	Vector f1(0.0), f2(0.0);

	for ( graph::Graph::EdgeListConstIter
			ir  = graph.get_node(i)->const_edge_list_begin(),
			ire = graph.get_node(i)->const_edge_list_end();
			ir != ire; ++ir ) {
		Size const j( (*ir)->get_other_ind( i ) );
		conformation::Residue const & rsd_j( pose.residue(j) );

		// if virtual residue, don't score
		if (rsd_j.aa() == core::chemical::aa_vrt ) continue;

		if ( input_atom_is_nbr && input_atom_is_rep && (rsd_j.is_protein() && rsd_j.aa()<=core::chemical::num_canonical_aas) ) {
			Size const resj_rep_atom = rsd_j.atom_index( representative_atom_name( rsd_j.aa() ));
			Size const resj_nbr_atom = rsd_j.nbr_atom();
			if ( resj_rep_atom == resj_nbr_atom ) {
				/// two birds, one stone
				increment_f1_f2_for_atom_pair(
					atom_i, rsd_j.atom( resj_rep_atom ),
					weights[ envsmooth ] * ( residue_dEdN_[ j ] + residue_dEdN_[ i ] ),
					F1, F2 );
			} else {
				increment_f1_f2_for_atom_pair(
					atom_i, rsd_j.atom( resj_rep_atom ),
					weights[ envsmooth ] * ( residue_dEdN_[ j ] ),
					F1, F2 );

				increment_f1_f2_for_atom_pair(
					atom_i, rsd_j.atom( resj_nbr_atom ),
					weights[ envsmooth ] * ( residue_dEdN_[ i ] ),
					F1, F2 );
			}
		} else if ( input_atom_is_nbr && (rsd_j.is_protein() && rsd_j.aa()<=core::chemical::num_canonical_aas) ) {
			Size const resj_rep_atom = rsd_j.atom_index( representative_atom_name( rsd_j.aa() ));

			increment_f1_f2_for_atom_pair(
				atom_i, rsd_j.atom( resj_rep_atom ),
				weights[ envsmooth ] * ( residue_dEdN_[ j ] ),
				F1, F2 );

		} else {
			Size const resj_nbr_atom = rsd_j.nbr_atom();
			increment_f1_f2_for_atom_pair(
				atom_i, rsd_j.atom( resj_nbr_atom ),
				weights[ envsmooth ] * ( residue_dEdN_[ i ] ),
				F1, F2 );

		}

	}

}

/// @details returns const & to static data members to avoid expense
/// of string allocation and destruction.  Do not call this function
/// on non-canonical aas
std::string const &
EnvSmoothEnergy::representative_atom_name( chemical::AA const aa ) const
{
	// assert( aa >= 1 && aa <= chemical::num_canonical_aas );

	static std::string const cbeta_string(  "CB"  );
	static std::string const sgamma_string( "SG"  );
	static std::string const cgamma_string( "CG"  );
	static std::string const cdelta_string( "CD"  );
	static std::string const czeta_string(  "CZ"  );
	static std::string const calpha_string( "CA"  );
	static std::string const ceps_1_string( "CE1" );
	static std::string const cdel_1_string( "CD1" );
	static std::string const ceps_2_string( "CE2" );
	static std::string const sdelta_string( "SD"  );

	switch ( aa ) {
		case ( chemical::aa_ala ) : return cbeta_string;  break;
		case ( chemical::aa_cys ) : return sgamma_string; break;
		case ( chemical::aa_asp ) : return cgamma_string; break;
		case ( chemical::aa_glu ) : return cdelta_string; break;
		case ( chemical::aa_phe ) : return czeta_string;  break;
		case ( chemical::aa_gly ) : return calpha_string; break;
		case ( chemical::aa_his ) : return ceps_1_string; break;
		case ( chemical::aa_ile ) : return cdel_1_string; break;
		case ( chemical::aa_lys ) : return cdelta_string; break;
		case ( chemical::aa_leu ) : return cgamma_string; break;
		case ( chemical::aa_met ) : return sdelta_string; break;
		case ( chemical::aa_asn ) : return cgamma_string; break;
		case ( chemical::aa_pro ) : return cgamma_string; break;
		case ( chemical::aa_gln ) : return cdelta_string; break;
		case ( chemical::aa_arg ) : return czeta_string;  break;
		case ( chemical::aa_ser ) : return cbeta_string;  break;
		case ( chemical::aa_thr ) : return cbeta_string;  break;
		case ( chemical::aa_val ) : return cbeta_string;  break;
		case ( chemical::aa_trp ) : return ceps_2_string; break;
		case ( chemical::aa_tyr ) : return czeta_string;  break;
		default :
			utility_exit_with_message( "ERROR: Failed to find amino acid " + chemical::name_from_aa( aa ) + " in EnvSmooth::representative_atom_name" );
		break;
	}

	// unreachable
	return calpha_string;
}

/// @brief EnvSmoothEnergy distance cutoff
Distance
EnvSmoothEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}

/// @brief EnvSmoothEnergy
void
EnvSmoothEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ twelve_A_neighbor_graph ] = true;
}

void
EnvSmoothEnergy::calc_energy(
	Real const neighbor_count,
	chemical::AA const aa,
	Real & score,
	Real & dscore_dneighbor_count
) const
{
	Size low_bin = static_cast< Size > ( floor(neighbor_count));
	Size high_bin = static_cast< Size > ( ceil(neighbor_count));
	Real inter = neighbor_count - low_bin;

	int const aa_as_int = static_cast< int > (aa);

	if (high_bin < 40 ){
		score =    envdata[ aa_as_int - 1 ][ low_bin ]  * (1.0-inter) +
			envdata[ aa_as_int - 1 ][ high_bin ] * (inter);
		dscore_dneighbor_count = envdata[ aa_as_int - 1 ][ high_bin ] -
			envdata[ aa_as_int - 1 ][ low_bin ];
	}else{
		score = envdata[ aa_as_int - 1 ][ 39 ];
		dscore_dneighbor_count = 0;
	}
	score *= 2.019; // this factor is from rosetta++ and fuck knows where it came from originally :-)
	dscore_dneighbor_count *= 2.019;
}

Real
EnvSmoothEnergy::sigmoidish_neighbor( DistanceSquared const sqdist ) const
{
	if( sqdist > end_sig2 ) {
		return 0.0;
	} else if( sqdist < start_sig2 ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0  - sqr( (dist - start_sig) / (end_sig - start_sig) ) );
	}
}


void
EnvSmoothEnergy::increment_f1_f2_for_atom_pair(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real weighted_dScore_dN,
	Vector & F1,
	Vector & F2
) const
{
	DistanceSquared dist2 = atom1.xyz().distance_squared(atom2.xyz());
	Distance dist( 0.0 ); // only used if start_sig2 <= dist2 <= end_sig2

	Real dNdd = 0;
	if ( dist2 > end_sig2 ) {
		dNdd = 0;
	} else if ( dist2 < start_sig2 ) {
		dNdd = 0.0;
	} else {
		dist = sqrt( dist2 );
		Real x = (dist - start_sig)/ (end_sig - start_sig);
		dNdd = 4*x*(-1 + x*x) / (end_sig - start_sig);
	}

	Real dscoredd = ( weighted_dScore_dN ) * dNdd;
	if ( dscoredd != 0 ) {

		Vector const f1( cross( atom1.xyz(), atom2.xyz() ));
		Vector const f2( atom1.xyz() - atom2.xyz() );

		dscoredd /= dist;
		F1 += dscoredd * f1;
		F2 += dscoredd * f2;
	}


}

core::Size
EnvSmoothEnergy::version() const
{
	return 1; // Initial versioning
}
}
}
}

