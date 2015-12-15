// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/MgMonteCarlo.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/MgMonteCarlo.hh>
#include <protocols/magnesium/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/magnesium/util.hh>
#include <numeric/random/random.hh>
#include <numeric/UniformRotationSampler.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.magnesium.MgMonteCarlo" );

///////////////////////////////////////////////////////////////
// QUICK HACK! Does not obey detailed balance! Although
// we could certainly make that happen.
///////////////////////////////////////////////////////////////

using namespace core;
using utility::vector1;
typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace magnesium {

//Constructor
MgMonteCarlo::MgMonteCarlo():
	cycles_( 10000 ),
	temperature_( 1.0 ),
	add_delete_frequency_( 0.1 ),
	output_pdb_( false )
{}

//Destructor
MgMonteCarlo::~MgMonteCarlo()
{}

void
MgMonteCarlo::apply( pose::Pose & pose ) {

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::magnesium;
	using namespace core::conformation;
	using namespace protocols::moves;
	using namespace protocols::magnesium;
	using namespace protocols::rigid;
	using numeric::random::rg;

	// set up fold tree
	setup_mg_water_fold_tree( pose );

	static numeric::UniformRotationSampler urs( 5.0 );

	ScoreFunctionOP scorefxn = get_mg_scorefxn();
	MonteCarlo monte_carlo( pose, *scorefxn, temperature_  );
	Real const basic_add_delete_frequency( add_delete_frequency_ );
	bool const add_at_random_in_shell( false );
	Distance const box_radius( 3.1 );

	for ( Size n = 1; n <= cycles_; n++ ) {
		if ( n % 10 == 0 ) {
			TR << "Doing cycle " << n << " out of " << cycles_ << std::endl;
			if ( output_pdb_ ) pose.dump_pdb( "cycle_" + ObjexxFCL::lead_zero_string_of( n, 6 ) + ".pdb" );
		}
		// pick a residue to move
		std::string tag;
		Real const add_delete_probability = ( pose.num_jump() == 0 ) ? 1.0 : basic_add_delete_frequency; // must add if no waters.
		if ( numeric::random::rg().uniform() <= add_delete_probability ) {
			//   pose::Pose pose_working = pose;
			vector1< Size > mg_res = get_mg_res( pose );
			Size const q( numeric::random::rg().random_range( 1, mg_res.size() ) );
			Size const i( mg_res[ q ] );
			vector1< core::id::AtomID > ligands = get_mg_ligands( pose, i );
			// check full_shell
			bool const full_shell = ( ligands.size() >= 6 );
			Real add_probability = full_shell ? 0.0 : 0.8;
			if ( numeric::random::rg().uniform() <= add_probability || pose.total_residue() == 1 ) {
				// add
				if ( add_at_random_in_shell ) {
					Vector xyz_water;
					while ( xyz_water.length() < 1.9 || xyz_water.length() > 2.2 ) {
						xyz_water = Vector( box_radius * rg().uniform(), box_radius * rg().uniform(), box_radius * rg().uniform() );
					}
					ResidueOP rsd_hoh = conformation::ResidueFactory::create_residue( pose.residue( 1 ).residue_type_set()->name_map("HOH") );
					Matrix R;
					urs.get( rg().random_range(1, urs.nrots() ), R );
					for ( Size j = 1; j <= rsd_hoh->natoms(); j++ ) {
						rsd_hoh->set_xyz( j, /*R * */ rsd_hoh->xyz( j ) + xyz_water );
					}
					pose.append_residue_by_jump( *rsd_hoh, i );
				} else {
					vector1< bool > already_coordinated( 6, false );
					for ( Size m = 1; m <= ligands.size(); m++ ) {
						already_coordinated[ get_closest_orbital_axis( pose.residue( i ), pose.xyz( ligands[ m ]) ) ] = true;
					}
					vector1< Size > pick_orbitals;
					for ( Size m = 1; m <= 6; m++ ) {
						if ( !already_coordinated[ m ] ) pick_orbitals.push_back( m );
					}
					Size const m = rg().random_element( pick_orbitals );
					//     Distance hoh_distance = 2.1 + 0.2 * rg().gaussian();
					instantiate_water_at_octahedral_vertex( pose, i, m );
				}
				tag = "add_water";
			} else {
				// delete
				Size const water_res = numeric::random::rg().random_element( get_water_res( pose ) );
				vector1< Size > slice_res;
				for ( Size n = 1; n <= pose.total_residue(); n++ ) {
					if ( n != water_res ) slice_res.push_back( n );
				}
				pdbslice( pose, slice_res );
				tag = "delete_water";
			}
			update_numbers_in_pdb_info( pose, true /* reset_waters*/ );
			setup_mg_water_fold_tree( pose );
			//   pose = pose;
		} else { // resample existing residue
			Size const njump = pose.fold_tree().num_jump();
			Size const j( numeric::random::rg().random_range( 1, njump ) );
			Real const rot_mag( 10.0 ), trans_mag( 0.2 );
			RigidBodyPerturbMover perturb_mover( j, rot_mag, trans_mag );
			perturb_mover.apply( pose );
			tag = "move";
		}

		monte_carlo.boltzmann( pose, tag );
	}

	monte_carlo.show_counters();
	monte_carlo.recover_low( pose );
}



///////////////////////////////////////////
void
MgMonteCarlo::setup_mg_water_fold_tree( pose::Pose & pose ) const {
	using namespace core::kinematics;
	FoldTree f( pose.total_residue() );
	vector1< Size > mg_res = get_mg_res( pose );
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( pose.residue( n ).name3() == "HOH" ) {
			f.new_jump( n, 1, n - 1 /* cutpoint */ );
		}
	}
	pose.fold_tree( f );
}

} //magnesium
} //protocols

