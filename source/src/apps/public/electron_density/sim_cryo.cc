// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/public/electron_density/sim_cryo.cc
/// @brief An app that will somewhat poorly simmulate cryoEM data
/// @author Danny Farrell


// basic
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/testing.OptionKeys.gen.hh>
#include <basic/options/keys/cryst.OptionKeys.gen.hh>

// core
#include <core/chemical/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/util.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

// devel
#include <devel/init.hh>

//numeric
#include <numeric/random/random.hh>
#include <numeric/UniformRotationSampler.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>


// objexx
#include <ObjexxFCL/format.hh>

// protocols
#include <protocols/electron_density/util.hh>
#include <core/types.hh>

// apps
#include <protocols/electron_density/StarFile.hh>
#include <protocols/electron_density/PerlinNoise.hh>

// std c++ headers
#include <iostream>
#include <iosfwd>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <memory>
#include <time.h>

// utility
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/vector1.hh>


#include <basic/Tracer.hh>
static basic::Tracer TR("sim_cryo");

using namespace basic::options;
using namespace basic::options::OptionKeys;

////////////////////////////////////////////// Options                   //////////////////////////////////////////////
OPT_KEY( Integer, NR )
OPT_KEY( Real, noise_range )
OPT_KEY( Real, gauss_multiplier )
OPT_KEY( Real, gauss_offset )
OPT_KEY( Real, resolution )
OPT_KEY( Real, pixel_size )
OPT_KEY( Real, perlin_noise_multiplier )
OPT_KEY( Real, perlin_noise_offset )
OPT_KEY( Real, perlin_noise_frequency )
OPT_KEY( Real, perlin_noise_octaves )
OPT_KEY( Real, particle_wt_offset )
OPT_KEY( Real, atom_gauss_random_multiplier )
OPT_KEY( String, dupe_chains )
OPT_KEY( Real, dupe_chains_times )
OPT_KEY( Real, box_size_multiplier )

////////////////////////////////////////////// setting default namespaces //////////////////////////////////////////////

class SimulateCryoMover {
public:
	SimulateCryoMover() {
		num_rotations_to_project_from_ = 100;
		noise_range_ = 0;
		gauss_multiplier_ = 0;
		gauss_offset_ = 0;
		pixel_size_ = 3.0;
		dist_to_detector_ = 1000.0;
		resolution_ = 1.0;
		perlin_noise_multiplier_ = 0;
		perlin_noise_offset_ = 0;
		box_size_multiplier_ = 1.5;
		integration_test_ = false;
		utility::vector1< core::Real > bfactors_;
		utility::vector1< core::Real > scatter_factors_;
	}

	void
	apply(core::pose::Pose pose);

	void
	recenter_pose_at_origin( core::pose::Pose & );

	void
	set_extent( core::pose::Pose const &, core::Real );

	void
	pose_to_xyz( core::pose::Pose const &,
		utility::vector1< core::Real > & all_x,
		utility::vector1< core::Real > & all_y,
		utility::vector1< core::Real > & all_z);

	void
	downsample_xyz( core::Real,
		utility::vector1< core::Real > &,
		utility::vector1< core::Real > &,
		utility::vector1< core::Real > &);

	void
	truncate_xyz(
		utility::vector1< core::Real > &,
		utility::vector1< core::Real > &,
		utility::vector1< core::Real > &,
		utility::vector1< core::Size > &,
		utility::vector1< core::Size > &,
		utility::vector1< core::Size > &);

	void
	recenter_to_img_coordinates(
		utility::vector1< core::Real > &,
		utility::vector1< core::Real > &,
		utility::vector1< core::Real > &);

	void
	add_flat_noise( ObjexxFCL::FArray3D< float > &, core::Real );

	void
	add_gaussian_noise( ObjexxFCL::FArray3D< float > &, core::Real, core::Real );

	void
	add_perlin_noise( ObjexxFCL::FArray3D< float > &, core::Real, core::Real );

	ObjexxFCL::FArray3D< float >
	three_projections_from_xyz(
		utility::vector1< core::Size > &,
		utility::vector1< core::Size > &,
		utility::vector1< core::Size > &);

	ObjexxFCL::FArray3D< float >
	three_fine_projections_from_xyz(
		utility::vector1< core::Real > &,
		utility::vector1< core::Real > &,
		utility::vector1< core::Real > &,
		core::Real, core::Real );

	void
	append_to_density(
		ObjexxFCL::FArray3D< float > const &,
		ObjexxFCL::FArray3D< float > &,
		int, int);

	void
	writeStarFile( std::string,
		utility::vector1< protocols::electron_density::StarFile > &);

	void
	normalize_projections( ObjexxFCL::FArray3D< float > & );

	void
	offset_projections( ObjexxFCL::FArray3D< float > &, core::Real );

	void
	set_bfactors( core::pose::Pose pose );

	void
	set_image_stack( core::scoring::electron_density::ElectronDensity const & stack_in ) {
		image_stack_ = stack_in;
	}

	core::scoring::electron_density::ElectronDensity
	get_image_stack() {
		return image_stack_;
	}

	void set_NR( core::Size NR ) { num_rotations_to_project_from_ = NR; }
	core::Size get_NR() { return num_rotations_to_project_from_; }

	void set_noise_range( core::Real noise_range ) { noise_range_ = noise_range; }
	core::Real get_noise_range() { return noise_range_; }

	void set_gauss_multiplier( core::Real gauss_multiplier ) { gauss_multiplier_ = gauss_multiplier; }
	core::Real get_gauss_multiplier() { return gauss_multiplier_; }

	void set_gauss_offset( core::Real gauss_offset ) { gauss_offset_ = gauss_offset; }
	core::Real get_gauss_offset() { return gauss_offset_; }

	void set_pixel_size( core::Real pixel_size ) { pixel_size_ = pixel_size; }
	core::Real get_pixel_size() { return pixel_size_; }

	void set_resolution( core::Real resolution ) { resolution_ = resolution; }
	core::Real get_resolution() { return resolution_; }

	void set_box_size_multiplier( core::Real box_size_multiplier ) { box_size_multiplier_ = box_size_multiplier; }
	core::Real get_box_size_multiplier() { return box_size_multiplier_; }
	void set_perlin_noise_multiplier( core::Real perlin_noise_multiplier ) { perlin_noise_multiplier_ = perlin_noise_multiplier; }
	core::Real get_perlin_noise_multiplier() { return perlin_noise_multiplier_; }

	void set_perlin_noise_offset( core::Real perlin_noise_offset ) { perlin_noise_offset_ = perlin_noise_offset; }
	core::Real get_perlin_noise_offset() { return perlin_noise_offset_; }

	void set_integration_test( bool const integration_test ) { integration_test_ = integration_test; }
	core::Real get_integration_test() const { return integration_test_; }

private:
	core::scoring::electron_density::ElectronDensity image_stack_;
	numeric::xyzVector< core::Size > extents_;
	core::Size num_rotations_to_project_from_;
	core::Real gauss_offset_, gauss_multiplier_, noise_range_, pixel_size_, dist_to_detector_, resolution_, box_size_multiplier_, perlin_noise_multiplier_, perlin_noise_offset_;
	bool integration_test_;
	utility::vector1< core::Real > bfactors_;
	utility::vector1< core::Real > scatter_factors_;
};




/* Recenter pose at 0, 0, 0
*
* @param pose The pose to move to 0, 0, 0
*/
void
SimulateCryoMover::recenter_pose_at_origin( core::pose::Pose & pose ) {
	numeric::xyzVector< core::Real > COM = core::pose::get_center_of_mass(pose);
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector< core::Real > > positions;
	numeric::xyzVector< core::Real > xyz_rot;

	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		core::conformation::Residue const & rsd( pose.residue(res_num) );
		for ( core::Size atm_num = 1; atm_num <= rsd.natoms(); ++atm_num ) {
			core::conformation::Atom const & atom( rsd.atom(atm_num) );
			numeric::xyzVector< core::Real > atom_xyz = atom.xyz();
			ids.push_back( core::id::AtomID( atm_num, res_num ) );
			atom_xyz = atom_xyz - COM;
			positions.push_back( atom_xyz );
		}
	}
	pose.batch_set_xyz( ids, positions );
}

/* Set the pixel size of the images based on pose measurements
*
* @param pose that you set the extent against
*/
void
SimulateCryoMover::set_extent( core::pose::Pose const & pose, core::Real pixel_size ) {
	core::Real xmin = 0, xmax = 0;
	core::Real ymin = 0, ymax = 0;
	core::Real zmin = 0, zmax = 0;
	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		core::conformation::Residue const & rsd( pose.residue(res_num) );
		for ( core::Size atm_num = 1; atm_num <= rsd.natoms(); ++atm_num ) {
			core::conformation::Atom const & atm( rsd.atom(atm_num) );
			if ( atm.xyz()[0] < xmin ) xmin = atm.xyz()[0];
			if ( atm.xyz()[0] > xmax ) xmax = atm.xyz()[0];
			if ( atm.xyz()[1] < ymin ) ymin = atm.xyz()[1];
			if ( atm.xyz()[1] < ymin ) ymin = atm.xyz()[1];
			if ( atm.xyz()[2] > zmax ) zmax = atm.xyz()[2];
			if ( atm.xyz()[2] > zmax ) zmax = atm.xyz()[2];
		}
	}
	// TODO multiplier?
	core::Size box_size = std::round(std::max({xmax - xmin, ymax - ymin, zmax - zmin}) * box_size_multiplier_ / pixel_size);
	box_size = std::nearbyint(box_size *0.5f)*2.0f; // Must be an even number
	// normalize to square
	extents_[0] = box_size;
	extents_[1] = box_size;
	extents_[2] = box_size;
}

/* Convert pose to 2 xyz vectors
*
* @param pose: to convert
* @param all_x: all x values (reference)
* @param all_y: all y values (reference)
* @param all_z: all z values (reference)
*/
void
SimulateCryoMover::pose_to_xyz( core::pose::Pose const & pose,
	utility::vector1< core::Real > & all_x,
	utility::vector1< core::Real > & all_y,
	utility::vector1< core::Real > & all_z) {
	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
		core::conformation::Residue const & rsd( pose.residue(res_num) );
		char chain = pose.pdb_info()->chain(res_num);
		std::size_t search(std::string::npos);
		if ( option[dupe_chains].user() ) {
			search = option[dupe_chains]().find(chain);
		}
		for ( core::Size atm_num = 1; atm_num <= rsd.natoms(); ++atm_num ) {
			core::conformation::Atom const & atm( rsd.atom(atm_num) );

			all_x.push_back(atm.xyz()[0]);
			all_y.push_back(atm.xyz()[1]);
			all_z.push_back(atm.xyz()[2]);
			if ( search != std::string::npos ) {
				for ( core::Size x = 1; x <= option[dupe_chains_times](); ++x ) {
					all_x.push_back(atm.xyz()[0]);
					all_y.push_back(atm.xyz()[1]);
					all_z.push_back(atm.xyz()[2]);
				}
			}
		}
	}
}


/* Downsample xyz coordinates by a divisor
*
* @param divisor: what to divide coordinates by
* @param all_x: all x coordinate values (reference)
* @param all_y: all y coordinate values (reference)
* @param all_z: all z coordinate values (reference)
*/
void
SimulateCryoMover::downsample_xyz( core::Real divisor,
	utility::vector1< core::Real > & all_x,
	utility::vector1< core::Real > & all_y,
	utility::vector1< core::Real > & all_z) {
	for ( core::Size i = 1; i <= all_x.size(); ++i ) {
		all_x[i] /= divisor;
		all_y[i] /= divisor;
		all_z[i] /= divisor;
	}
}


/* Truncate xyz coordinates to round numbers
*
* @param all_x: all x coordinate values (reference)
* @param all_y: all y coordinate values (reference)
* @param all_z: all z coordinate values (reference)
*
*/
void
SimulateCryoMover::truncate_xyz(
	utility::vector1< core::Real > & all_x,
	utility::vector1< core::Real > & all_y,
	utility::vector1< core::Real > & all_z,
	utility::vector1< core::Size > & all_x_rnd,
	utility::vector1< core::Size > & all_y_rnd,
	utility::vector1< core::Size > & all_z_rnd) {
	for ( core::Size i = 1; i <= all_x.size(); ++i ) {
		all_x_rnd.push_back(std::round(all_x[i]));
		all_y_rnd.push_back(std::round(all_y[i]));
		all_z_rnd.push_back(std::round(all_z[i]));
	}
}


/* Recenter xyz coordinates to the center of a box (no negatives)
*
* @param all_x: all x coordinate values (reference)
* @param all_y: all y coordinate values (reference)
* @param all_z: all z coordinate values (reference)
*/
void
SimulateCryoMover::recenter_to_img_coordinates(
	utility::vector1< core::Real > & all_x,
	utility::vector1< core::Real > & all_y,
	utility::vector1< core::Real > & all_z) {
	core::Real center = extents_[0]*0.5;

	for ( core::Size i = 1; i <= all_x.size(); ++i ) {
		all_x[i] += center;
		all_y[i] += center;
		all_z[i] += center;
	}
}


/* Randomly add a small amount of noise to each pixel
*
* @param three_projections the three projections of our pose (reference)
*/
void
SimulateCryoMover::add_flat_noise( ObjexxFCL::FArray3D< float > & three_projections, core::Real noise_range) {
	numeric::random::RandomGenerator & rg(numeric::random::rg());
	if ( integration_test_ ) rg.set_seed(1);
	else { rg.set_seed(time(NULL)); }

	for ( int z = 1; z <= three_projections.u3(); ++z ) {
		for ( int y = 1; y <= three_projections.u2(); ++y ) {
			for ( int x = 1; x <= three_projections.u1(); ++x ) {
				core::Real noise_to_add = rg.uniform() * noise_range;
				three_projections(x, y, z) += noise_to_add;
			}
		}
	}

}


/* Randomly add gaussian noise to various pixels
*
* @param three_projections the three projections of our pose (reference)
* @param noise_range the sigma to simulate
*/
void
SimulateCryoMover::add_gaussian_noise( ObjexxFCL::FArray3D< float > & three_projections, core::Real gauss_multiplier, core::Real gauss_offset ) {
	for ( int z = 1; z <= three_projections.u3(); ++z ) {
		for ( int y = 1; y <= three_projections.u2(); ++y ) {
			for ( int x = 1; x <= three_projections.u1(); ++x ) {
				core::Real noise_to_add = numeric::random::gaussian() * gauss_multiplier + gauss_offset;
				three_projections(x, y, z) += noise_to_add;
			}
		}
	}
}

///* Add perlin noise to all pixels
// *
// * @param three_projections the three projections of our pose (reference)
// * @param noise_range the sigma to simulate
// */
void
SimulateCryoMover::add_perlin_noise( ObjexxFCL::FArray3D< float > & three_projections, core::Real perlin_noise_multiplier, core::Real perlin_noise_offset ) {
	if ( integration_test_ ) srand(1);
	else { srand(time(NULL)); }

	// must be between 0.1 and 64
	double const frequency = option[ perlin_noise_frequency ]();

	// must be between 1 and 16
	int const octaves = option[ perlin_noise_octaves ]();

	double const fx = three_projections.u1() / frequency;
	double const fy = three_projections.u2() / frequency;
	for ( int z = 1; z <= three_projections.u3(); ++z ) {
		protocols::electron_density::PerlinNoise perlin(rand());
		//core::Real z_const = 0.5; // could set this seperately for any proj or just leave at 0.5
		for ( int y = 1; y <= three_projections.u2(); ++y ) {
			for ( int x = 1; x <= three_projections.u1(); ++x ) {

				double noise = perlin.octaveNoise0_1(x/fx, y/fy, octaves );

				three_projections(x, y, z) += (noise * perlin_noise_multiplier + perlin_noise_offset);
				//three_projections(x, y, z) +=  (noise * perlin_noise_multiplier + perlin_noise_offset) *(1- three_projections(x, y, z));
				//three_projections(x, y, z) = ((noise * perlin_noise_multiplier + perlin_noise_offset) *(1- three_projections(x, y, z)) + three_projections(x,y,z))/2;
			}
		}
	}
}

/* Generate 3 fine images (xy, xz, yz) based on our xyz vectors
*
* @param all_x: all x coordinate values (reference)
* @param all_y: all y coordinate values (reference)
* @param all_z: all z coordinate values (reference)
* @param dist_to_detector: distance from mean to detector
* @param resolution: wanted resolution
*
* @return three_projections: three projections in FArray format (density)
*/
ObjexxFCL::FArray3D< float >
SimulateCryoMover::three_fine_projections_from_xyz(
	utility::vector1< core::Real > & all_x,
	utility::vector1< core::Real > & all_y,
	utility::vector1< core::Real > & all_z,
	core::Real resolution, core::Real pixel_size ) {

	numeric::random::RandomGenerator & rg(numeric::random::rg());
	if ( integration_test_ ) rg.set_seed(1);
	else { rg.set_seed(time(NULL)); }

	ObjexxFCL::FArray3D< float > three_projections;
	three_projections.dimension( extents_[0], extents_[0], 3, 0.0 );
	core::Real sigma = resolution / M_PI; // TODO hardcoded
	int rx2 = std::ceil( 4 / pixel_size );

	for ( core::Size i = 1; i <= all_x.size(); ++i ) {
		for ( int j = -rx2; j <= rx2; j = j+pixel_size ) {
			for ( int k = -rx2; k <= rx2; k = k+pixel_size ) {

				core::Real const distance = std::sqrt(std::pow(j,2) + std::pow(k,2));
				core::Real const gauss_f = 1 / (2*M_PI*sigma) * exp(-pow( distance/sigma, 2) );

				if ( option[ atom_gauss_random_multiplier ].user() ) {
					core::Real mp = 0.0;
					if ( basic::options::option[ basic::options::OptionKeys::cryst::crystal_refine ]() ) {
						mp = bfactors_[i] / 120.0;
						mp = mp / scatter_factors_[i];
					}
					else {
						mp = option[ atom_gauss_random_multiplier ]();
					}
					// randomly pick a rotation around 360
					// 'randomly' pick a distance around
					core::Real r_ang = rg.uniform() * M_PI * 2;
					core::Real dis = rg.uniform() * mp;
					core::Real ichange = cos(r_ang)*dis;
					core::Real jchange = sin(r_ang)*dis;
					if ( three_projections.contains(std::floor(all_x[i]+ichange+j), std::floor(all_y[i]+jchange+k), 1)  ) {
						three_projections(std::floor(all_x[i]+ichange+j), std::floor(all_y[i]+jchange+k), 1) += gauss_f;
					} else TR.Warning << "Attempted to access out of range FArray3d" << std::endl;
					r_ang = rg.uniform()*M_PI*2;
					dis = rg.uniform()*mp;
					ichange = cos(r_ang)*dis;
					jchange = sin(r_ang)*dis;
					if ( three_projections.contains(std::floor(all_x[i]+ichange+j), std::floor(all_z[i]+jchange+k), 2)  ) {
						three_projections(std::floor(all_x[i]+ichange+j), std::floor(all_z[i]+jchange+k), 2) += gauss_f;
					} else TR.Warning << "Attempted to access out of range FArray3d" << std::endl;
					r_ang = rg.uniform()*M_PI*2;
					dis = rg.uniform()*mp;
					ichange = cos(r_ang)*dis;
					jchange = sin(r_ang)*dis;
					if ( three_projections.contains(std::floor(all_y[i]+ichange+j), std::floor(all_z[i]+jchange+k), 3)  ) {
						three_projections(std::floor(all_y[i]+ichange+j), std::floor(all_z[i]+jchange+k), 3) += gauss_f;
					} else TR.Warning << "Attempted to access out of range FArray3d" << std::endl;

				} else {
					if ( three_projections.contains(std::floor(all_x[i]+j), std::floor(all_y[i]+k), 1)  ) {
						three_projections(std::floor(all_x[i]+j), std::floor(all_y[i]+k), 1) += gauss_f;
					} else TR.Warning << "Attempted to access out of range FArray3d" << std::endl;
					if ( three_projections.contains(std::floor(all_x[i]+j), std::floor(all_z[i]+k), 2) ) {
						three_projections(std::floor(all_x[i]+j), std::floor(all_z[i]+k), 2) += gauss_f;
					} else TR.Warning << "Attempted to access out of range FArray3d" << std::endl;
					if ( three_projections.contains(std::floor(all_y[i]+j), std::floor(all_z[i]+k), 3) ) {
						three_projections(std::floor(all_y[i]+j), std::floor(all_z[i]+k), 3) += gauss_f;
					} else TR.Warning << "Attempted to access out of range FArray3d" << std::endl;
				}
			}
		}
	}
	return three_projections;
}


/* Take 3 images (xy, xz, yz) based on our xyz vectors
*
* @param all_x: all x coordinate values (reference)
* @param all_y: all y coordinate values (reference)
* @param all_z: all z coordinate values (reference)
*
* @return three_projections: three projections in FArray format (density)
*/
ObjexxFCL::FArray3D< float >
SimulateCryoMover::three_projections_from_xyz(
	utility::vector1< core::Size > & all_x,
	utility::vector1< core::Size > & all_y,
	utility::vector1< core::Size > & all_z) {
	ObjexxFCL::FArray3D< float > three_projections;
	// TODO: are we sure that we want three?
	three_projections.dimension( extents_[0], extents_[0], 3, 0.0 );
	for ( core::Size i = 1; i <= all_x.size(); ++i ) {
		++three_projections(all_x[i], all_y[i], 1);
		++three_projections(all_x[i], all_z[i], 2);
		++three_projections(all_y[i], all_z[i], 3);
	}
	return three_projections;
}

/* Find max value in each image, and divide by the max (norm between 0 - 1)
*
* @param three_projections projections to normalize
*/
void
SimulateCryoMover::normalize_projections( ObjexxFCL::FArray3D< float > & three_projections ) {
	for ( int z = 1; z <= three_projections.u3(); ++z ) {
		float max = 0, min = 1000;
		for ( int y = 1; y <= three_projections.u2(); ++y ) {
			for ( int x = 1; x <= three_projections.u1(); ++x ) {
				if ( three_projections(x, y, z) > max ) max = three_projections(x, y, z);
				if ( three_projections(x, y, z) < min ) min = three_projections(x, y, z);
			}
		}
		for ( int y = 1; y <= three_projections.u2(); ++y ) {
			for ( int x = 1; x <= three_projections.u1(); ++x ) {
				three_projections(x, y, z) /= max;
			}
		}
	}
}

void
SimulateCryoMover::offset_projections( ObjexxFCL::FArray3D< float > & three_projections, core::Real particle_offset  ) {
	for ( int z = 1; z <= three_projections.u3(); ++z ) {
		for ( int y = 1; y <= three_projections.u2(); ++y ) {
			for ( int x = 1; x <= three_projections.u1(); ++x ) {
				three_projections(x, y, z) += particle_offset;
			}
		}
	}
}


/* Append images to density data based on given index (x3)
*
* @param images_in: FArray3D with images stored in first 3 z stacks (const)
* @param density: Final total stack of images (reference)
* @param index: offset to place images in density
*/
void
SimulateCryoMover::append_to_density(
	ObjexxFCL::FArray3D< float > const & three_projections,
	ObjexxFCL::FArray3D< float > & density,
	int offset, int proj_num) {
	for ( int y = 1; y <= three_projections.u2(); ++y ) {
		for ( int x = 1; x <= three_projections.u1(); ++x ) {
			density(x, y, offset ) = three_projections( x, y, proj_num );
		}
	}
}

void
SimulateCryoMover::set_bfactors( core::pose::Pose pose ) {
	utility::vector1< core::Real > bfactors;
	utility::vector1< core::Real > scatter_factors;

	std::map< std::string, core::Real > scatter_map;
	scatter_map["C"] = 6.0;
	scatter_map["N"] = 5.28737;
	scatter_map["O"] = 4.74213;
	scatter_map["NA"] = 11.42607;
	scatter_map["MG"] = 12.45197;
	scatter_map["P"] = 13.12395;
	scatter_map["S"] = 12.34197;
	scatter_map["K"] = 21.48425;
	scatter_map["CA"] = 23.70586;
	scatter_map["FE"] = 17.13431;
	scatter_map["CO"] = 15.70905;
	scatter_map["NI"] = 15.70905;
	scatter_map["ZN"] = 15.70905;

	for ( core::Size res_num = 1; res_num <= pose.size(); ++res_num ) {
                core::conformation::Residue const & rsd( pose.residue(res_num) );
                char chain = pose.pdb_info()->chain(res_num);
                std::size_t search(std::string::npos);
                if ( option[dupe_chains].user() ) {
                        search = option[dupe_chains]().find(chain);
                }
                for ( core::Size atm_num = 1; atm_num <= rsd.natoms(); ++atm_num ) {
                        core::conformation::Atom const & atm( rsd.atom(atm_num) );

			core::Real bfactor = pose.pdb_info()->temperature( res_num, atm_num );
			bfactors.push_back(pose.pdb_info()->temperature( res_num, atm_num ));
			
			if ( scatter_map.count( rsd.atom_type(atm_num).element() ) ) {
				core::Real scatter_factor = scatter_map[ rsd.atom_type(atm_num).element() ];
				scatter_factors.push_back( scatter_factor / 6.0 );

				core ::Real mp = bfactor / 120.0;
                                mp = mp / ( scatter_factor / 6.0 );
				TR << "Max pertubation at " << rsd.name() << " " << res_num << " atom " << rsd.atom_type(atm_num).name() << " is: " << mp << std::endl;
			}
			else {
				scatter_factors.push_back( 1.0 );
			}

                        if ( search != std::string::npos ) {
                                for ( core::Size x = 1; x <= option[dupe_chains_times](); ++x ) {
					bfactors.push_back(pose.pdb_info()->temperature( res_num, atm_num ));

					if ( scatter_map.count( rsd.atom_type(atm_num).element() ) ) {
						core::Real scatter_factor = scatter_map[ rsd.atom_type(atm_num).element() ];
						scatter_factors.push_back( scatter_factor / 6.0 );
					}
					else {
						scatter_factors.push_back( 1.0 );
					}
                                }
                        }
                }
        }

	bfactors_ = bfactors;
	scatter_factors_ = scatter_factors;
}

/* Generate images in mrcs format and save them
*
* @param pose that we will Rotate and take pictures of
*/
void
SimulateCryoMover::apply( core::pose::Pose pose ) {
	std::string out_base_name = "test";
	if ( option[ OptionKeys::out::prefix ].user() ) out_base_name = option[ out::prefix ]();
	SimulateCryoMover::recenter_pose_at_origin(pose);
	SimulateCryoMover::set_extent(pose, pixel_size_);
	ObjexxFCL::FArray3D< float > density;
	density.dimension( extents_[0], extents_[0], num_rotations_to_project_from_ * 3, 0.0 );
	numeric::UniformRotationSampler urs( 3 ); // 3 degree rotation step
	numeric::xyzVector< core::Real > no_trans_vector = { 0, 0, 0 };
	numeric::xyzMatrix< core::Real > rot_matrix;
	utility::vector1< protocols::electron_density::StarFile > starfile_information;
	core::Size urs_size = urs.nrots();

	TR << "Before B factors" << std::endl;
	if ( basic::options::option[ basic::options::OptionKeys::cryst::crystal_refine ]() ) {
		TR << "Setting B factors" << std::endl;
		set_bfactors( pose );
	}

	protocols::electron_density::StarFile s;
	for ( core::Size i = 1; i <= num_rotations_to_project_from_; ++i ) {
		if ( i % 100 == 0 ) TR << "At: " << i << std::endl;

		core::Size select_pos = numeric::random::random_range( 1, urs_size );
		urs.get(select_pos, rot_matrix);
		pose.apply_transform_Rx_plus_v(rot_matrix, no_trans_vector);
		utility::vector1< core::Real > all_x;
		utility::vector1< core::Real > all_y;
		utility::vector1< core::Real > all_z;

		SimulateCryoMover::pose_to_xyz( pose, all_x, all_y, all_z );
		SimulateCryoMover::recenter_to_img_coordinates( all_x, all_y, all_z );
		SimulateCryoMover::downsample_xyz( pixel_size_, all_x, all_y, all_z );
		ObjexxFCL::FArray3D< float > three_fine_projections = SimulateCryoMover::three_fine_projections_from_xyz( all_x, all_y, all_z, resolution_, pixel_size_ );

		SimulateCryoMover::normalize_projections( three_fine_projections );
		SimulateCryoMover::offset_projections( three_fine_projections, option[ particle_wt_offset ]()  );

		if ( noise_range_ != 0 ) SimulateCryoMover::add_flat_noise( three_fine_projections, noise_range_ );
		if ( gauss_multiplier_ != 0 ) SimulateCryoMover::add_gaussian_noise( three_fine_projections, gauss_multiplier_, gauss_offset_ );
		if ( perlin_noise_multiplier_ != 0 ) SimulateCryoMover::add_perlin_noise( three_fine_projections, perlin_noise_multiplier_, perlin_noise_offset_ );

		// append to density one projection at a time to coordinate with StarFile
		// TODO: just copy add the same starfile 3 times, since these will all have similar options
		for ( core::Size j = 1; j <= 3; ++j ) {
			append_to_density( three_fine_projections, density, (i-1)*3 + j, j);
			s.append_one();
			s.set(s.size(),"rlnImageName", (std::string)(utility::to_string((i-1)*3 + j) + "@" + out_base_name + ".mrcs"));
			s.set(s.size(),"rlnDetectorPixelSize", (float)5);
			s.set(s.size(),"rlnMagnification", (float)5);
			s.set(s.size(),"rlnPixelSize", (float)pixel_size_);
			s.set(s.size(),"rlnVoltage", (float)200);
			s.set(s.size(),"rlnAmplitudeContrast", (float)0.1);
			s.set(s.size(),"rlnDefocusU", (float)0.0);
			s.set(s.size(),"rlnDefocusV", (float)0.0);
			s.set(s.size(),"rlnDefocusAngle", (float)0.0);
			s.set(s.size(),"rlnSphericalAberration", (float)0.0);
		}
		// pose.dump_pdb(utility::to_string(i) + ".pdb");
	}
	image_stack_ = core::scoring::electron_density::ElectronDensity(density);
	image_stack_.writeMRC(out_base_name + ".mrcs");
	s.writeStarFile(out_base_name + ".star");
}

void
set_options( SimulateCryoMover & GIFP ) {
	GIFP.set_NR( option[ NR ]() );
	GIFP.set_noise_range( option[ noise_range ]() );
	GIFP.set_gauss_multiplier( option[ gauss_multiplier ]() );
	GIFP.set_gauss_offset( option[ gauss_offset ]() );
	GIFP.set_resolution( option[ resolution ]() );
	GIFP.set_pixel_size( option[ pixel_size ]() );
	GIFP.set_perlin_noise_multiplier( option[ perlin_noise_multiplier ]() );
	GIFP.set_perlin_noise_offset( option[ perlin_noise_offset ]() );
	GIFP.set_box_size_multiplier( option[ box_size_multiplier ]() );
	GIFP.set_integration_test( option[ testing::INTEGRATION_TEST ]() );
}

int main(int argc, char* argv[]) {
	try {
		NEW_OPT(NR, "Number of rotations to project from", 10 );
		NEW_OPT(noise_range, "Multiplier of noise.", 0 );
		NEW_OPT(gauss_multiplier, "Multiplier of gaussian noise.", 0 );
		NEW_OPT(gauss_offset, "Add this to gaussian noise.", 0 );
		NEW_OPT(resolution, "Map resolution expected.", 6 );
		NEW_OPT(pixel_size, "What to set the pixel size to.", 3 );
		NEW_OPT(perlin_noise_multiplier, "Multiplier of perlin noise.", 0 );
		NEW_OPT(perlin_noise_offset, "What to offset (add) to perlin noise.", 0 );
		NEW_OPT(perlin_noise_frequency, "What perlin noise frequency to use.", 0 );
		NEW_OPT(perlin_noise_octaves, "What perlin noise octaves to use.", 0 );
		NEW_OPT(particle_wt_offset, "offset particle by this much", 0 );
		NEW_OPT(atom_gauss_random_multiplier, "Multiplier of normal gaussian to randomize atoms", 0 );
		NEW_OPT(dupe_chains, "List of chains to duplicate", "" );
		NEW_OPT(dupe_chains_times, "Times to add dupe chains.", 0);
		NEW_OPT(box_size_multiplier, "Size of the box of the image.", 1.5);

		devel::init( argc, argv ); // initilize options system

		core::pose::PoseOP poseOP = core::import_pose::pose_from_file( option[ in::file::s ]()[1] );
		SimulateCryoMover GIFP;
		set_options(GIFP);
		GIFP.apply(*poseOP);
	} catch ( std::runtime_error const & e ) {
		std::cout << "Caught exception " << e.what() << std::endl;
		return -1;
	}
	return 0;
}
