// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/pilot/jkleman/mp_lipid_acc.cc
/// @brief  compute lipid accessibility from a membrane protein structure
/// @author  JKLeman (julia.koehler.leman@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/util.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <protocols/membrane/util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/membrane/AddMembraneMover.hh>
// I have no idea why this is needed??? but it doesn't compile without it
#include <protocols/membrane/TranslationRotationMover.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/scoring/sasa/util.hh>

// Utility Headers
#include <core/types.hh>
#include <core/pose/util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <numeric/numeric.functions.hh>
#include <utility/numbers.hh>

// C++ headers
#include <iostream>
#include <cstdlib>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.jkleman.mp_lipid_acc" );

using namespace core;
using namespace core::pose;
using namespace numeric;
using namespace protocols::moves;
using namespace core::conformation::membrane;
using namespace protocols::simple_moves;
using namespace protocols::membrane;
using namespace basic::options;

/// @brief Load Membrane Mover: Load and Initialize a Membrane Protein
/// using RosettaMP
class MPLipidAccessibilityMover : public Mover {

public:

	/// @brief Default Constructor
	MPLipidAccessibilityMover() : Mover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MPLipidAccessibilityMover"; }

	/////////////////////////////////////////

	/// @brief apply
	void apply( Pose & pose ) {

		using namespace basic::options;

		// set defaults
		core::Real user_radius_angstrom( 10 );
		core::Real user_rel_sasa( 0.2 );
		bool user( false );
		core::Real angle_cutoff( 65.0 );
		core::Real relative_radius( 0.65 );

		// init from cmd
		if ( option[ OptionKeys::mp::lipid_acc::radius_cutoff ].user() ) {
			user_radius_angstrom = option[ OptionKeys::mp::lipid_acc::radius_cutoff ]();
			user = true;
		}
		if ( option[ OptionKeys::mp::lipid_acc::rel_sasa_cutoff ].user() ) {
			user_rel_sasa = option[ OptionKeys::mp::lipid_acc::rel_sasa_cutoff ]();
			user = true;
		}
		if ( option[ OptionKeys::mp::lipid_acc::angle_cutoff ].user() ) {
			angle_cutoff = option[ OptionKeys::mp::lipid_acc::angle_cutoff ]();
			user = true;
		}

		////////// SWITCH TO FULL-ATOM, REQUIRED FOR VECTOR DEFINITIONS ////////
		using namespace protocols::simple_moves;
		SwitchResidueTypeSetMoverOP full_atom( new SwitchResidueTypeSetMover( "fa_standard" ) );
		full_atom->apply( pose );

		// add membrane
		using namespace protocols::membrane;
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );

		//////////////////////// CHECKS ////////////////////////////////////////

		// get minimum and maximum z-coordinate
		// check whether structure is transformed into membrane, get CA position
		core::Real min_z( 10000 );
		core::Real max_z( -10000 );

		for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

			// skip for membrane residue
			if ( pose.residue( i ).name3() == "MEM" ) {
				continue;
			}

			if ( pose.residue( i ).xyz( "CA" ).z() < min_z ) {
				min_z = pose.residue( i ).xyz( "CA" ).z();
			} else if ( pose.residue( i ).xyz( "CA" ).z() > max_z ) {
				max_z = pose.residue( i ).xyz( "CA" ).z();
			}
		}

		// crude check: if no CA atoms with positive and negative z-coordinates,
		// then probably wrong
		if ( !( min_z < 0 && max_z > 0 ) ) {
			utility_exit_with_message( "Your protein is not transformed into membrane coordinates! Cannot compute lipid accessibility. Quitting." );
		}

		////////////////// SET B-FACTOR TO ZERO FOR ALL ATOMS //////////////////

		for ( core::Size r = 1; r <= nres_protein( pose ); ++r ) {
			for ( core::Size a = 1; a <= pose.residue( r ).natoms(); ++a ) {
				pose.pdb_info()->bfactor( r, a, 0 );
			}
		}

		//////////////////////// COMPUTE SASA //////////////////////////////////

		// compute relative SASA for the protein,
		// be aware: the SASA vector goes over all residues
		utility::vector1< core::Real > rel_res_sasa = core::scoring::sasa::rel_per_res_sc_sasa( pose );

		//////////////////////// COMPUTE SLICES ////////////////////////////////

		// define variables, the outer vector goes through the slices
		// inner vector goes through residues in each slice
		utility::vector1< core::Real > slice_z_cntr;
		utility::vector1< utility::vector1< core::Size > > resi;
		utility::vector1< utility::vector1< core::Vector > > ca_coord, cb_coord;
		utility::vector1< core::Vector > slice_com;

		// go from negative to positive z-direction, define slices
		core::Real slice_width = 5.0;
		core::Real iter = - pose.conformation().membrane_info()->membrane_thickness();

		// go through protein in the membrane and get slices and arrays for them
		while ( iter <= pose.conformation().membrane_info()->membrane_thickness() ) {

			slice_z_cntr.push_back( iter );

			// temp utility vectors for each slice
			utility::vector1< core::Size > slice_res;
			utility::vector1< core::Vector > slice_ca, slice_cb;

			// go through protein residues
			for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

				// skip membrane residue
				if ( pose.residue( i ).name3() == "MEM" ) {
					continue;
				}

				// get CA and CB coordinates
				core::Vector ca = pose.residue( i ).xyz( "CA" );
				core::Vector cb;

				// use 2HA instead of CB for glycine
				if ( pose.residue( i ).name3() == "GLY" ) {
					cb = pose.residue( i ).xyz( "2HA" );
				} else {
					cb = pose.residue( i ).xyz( "CB" );
				}

				// if CA-coordinate is within slice, add resi, CA and CB to arrays
				if ( ( ca.z() >= iter && ca.z() < iter + slice_width ) or
						( cb.z() >= iter && cb.z() < iter + slice_width ) ) {

					slice_res.push_back( i );
					slice_ca.push_back( ca );
					slice_cb.push_back( cb );

				}
			}

			// push back into slices
			resi.push_back( slice_res );
			ca_coord.push_back( slice_ca );
			cb_coord.push_back( slice_cb );

			// go to next slice
			iter += slice_width;

		} // get arrays of slices

		// go through slices and compute COMs
		for ( core::Size s = 1; s < slice_z_cntr.size(); ++s ) {

			core::Vector com( 0, 0, 0 );
			core::Real radius = 0;

			// go through residues, compute COM per 5Å slice
			for ( core::Size r = 1; r <= resi[ s ].size(); ++r ) {
				com += ca_coord[ s ][ r ];
			}
			com /= resi[ s ].size();

			// compute maximum radius for that slice
			for ( core::Size r = 1; r <= resi[ s ].size(); ++r ) {

				// compute distance CA-COM
				core::Vector ca_com = ca_coord[ s ][ r ] - com;

				// get maximum radius
				if ( ca_com.length() > radius ) {
					radius = ca_com.length();
				}
			}

			// go through residues, get angle and set B-factor
			for ( core::Size r = 1; r <= resi[ s ].size(); ++r ) {

				// compute distance CA-COM
				core::Vector ca_com = ca_coord[ s ][ r ] - com;

				// compute angle between COM-CA-CB
				// if angle < angle_cutoff, then sidechain face COM
				// if angle >= angle_cutoff, then sidechain faces away from COM
				core::Real angle = numeric::angle_degrees( com, ca_coord[ s ][ r ], cb_coord[ s ][ r ] );

				// go through atoms in residue and set B-factor
				for ( core::Size a = 1; a <= pose.residue( resi[ s ][ r ] ).natoms(); ++a ) {

					// IF X THEN LIPID ACCESSIBLE (B-FACTOR = 50)
					// if in the membrane
					//     if ( ca_coord[ s ][ r ].z() >= - pose.conformation().membrane_info()->membrane_thickness() && ca_coord[ s ][ r ].z() <= pose.conformation().membrane_info()->membrane_thickness() ) {

					// debug
					if ( a == 2 ) {
						TR << "residue " << resi[ s ][ r ] << " angle " << angle << " accessibility " << rel_res_sasa[ resi[ s ][ r ] ] << " ca_com " << ca_com.length() << " radius " << radius*relative_radius << std::endl;
					}

					// user-defined, i.e. uses SASA and radius in Å instead of relative radius
					if ( user == true ) {

						// facing outwards and larger than certain radius with user-defined
						// radius and accessibility
						if ( angle > angle_cutoff && ca_com.length() >= user_radius_angstrom && rel_res_sasa[ resi[ s ][ r ] ] >= user_rel_sasa ) {
							pose.pdb_info()->bfactor( resi[ s ][ r ], a, 50.0 );
						} else if ( angle > angle_cutoff && ca_com.length() >= user_radius_angstrom && rel_res_sasa[ resi[ s ][ r ] ] == 0 && pose.residue( resi[ s ][ r ] ).name3() == "GLY" ) {
							// same, just for GLY
							pose.pdb_info()->bfactor( resi[ s ][ r ], a, 50.0 );
						}
					} else {
						// doesn't use SASA (too finicky) and uses relative radius, default app

						// facing outwards and larger than certain radius
						if ( angle > angle_cutoff && ca_com.length() >= radius * relative_radius ) {
							pose.pdb_info()->bfactor( resi[ s ][ r ], a, 50.0 );
						}

						//  && rel_res_sasa[ resi[ s ][ r ] ] > 0.2
						// for single helices and 2 helices, mostly lipid exposed
						if ( pose.conformation().membrane_info()->spanning_topology()->nspans() <= 2 ) {
							pose.pdb_info()->bfactor( resi[ s ][ r ], a, 50.0 );
						}
					} // user or not
					//     } // in membrane
				} // atoms
			} // residues
		} // slices
	} // apply

};

/////////////////////////////////////////

typedef utility::pointer::shared_ptr< MPLipidAccessibilityMover > MPLipidAccessibilityMoverOP;

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		using namespace protocols::moves;
		using namespace protocols::membrane;

		protocols::jd2::register_options();

		// Create and kick off a new load membrane mover
		MPLipidAccessibilityMoverOP lipid_acc( new MPLipidAccessibilityMover() );
		protocols::jd2::JobDistributor::get_instance()->go( lipid_acc );

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

