// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/polar_hydrogens/PolarHydrogenPacker.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/modeler/polar_hydrogens/PolarHydrogenPacker.hh>
#include <protocols/stepwise/modeler/polar_hydrogens/util.hh>
#include <core/pose/variant_util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.modeler.polar_hydrogens.PolarHydrogenPacker" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Specialized packer for Polar Hydrogens, including sampling of bond lengths, angles, and torsions.
//
// - Motivated by cases like the 5 H-bonds in the UUCG RNA tetraloop which simply do not get scored properly with
//   ideal bond lengths & angles. 'Non-planar' hydrogens appear ubiquitous in RNA, and are critical to model since
//   almost all hydrogen bond donors are such hydrogens.
//
// - Could be part of normal core::pack framework, but that would be really slow -- and subject to unnecessary
//    rotamer explosion. A single residue like guanosine has 4 polar hydrogens that are, for the most part, independent
//    of each other -- current standard packer would have to instantiation residues with all combinations.
//
// - For scoring, currently just evaluates potential for forming H-bonds.
//
// - Allows for virtualizing/instantiation of 2'-OH; could be generalized to any proton chi torsion, actually.
//
// - TODO: Not guaranteed to give best scoring configuration (currently is greedy, traversing through hydrogens in
//    sequence order), but we could solve this exactly through dynamic programming. See, e.g.,
//
//    Leaver-Fay, et al. (2008), "Faster placement of hydrogens in protein structures by dynamic programming", JEA.
//    http://www.cs.amherst.edu/~ccm/cs34/papers/a2_5-leaver-fay-1.pdf
//
//  -- rhiju, 2014
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace core;
using namespace core::scoring;
using namespace core::scoring::hbonds;
using utility::vector1;
using namespace basic::options;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace polar_hydrogens {

//Constructor
PolarHydrogenPacker::PolarHydrogenPacker():
	Mover(),
	disallow_vary_geometry_proton_chi_( option[ OptionKeys::stepwise::polar_hydrogens::disallow_vary_geometry_proton_chi ]() ) // for debugging
{
	init();
}

//Destructor
PolarHydrogenPacker::~PolarHydrogenPacker() = default;

///////////////////////////////////////////////////////////////////////////////////
void
PolarHydrogenPacker::init(){
	hbond_options_ = HBondOptionsOP( new HBondOptions() );
	hbond_options_->use_hb_env_dep( false );
	hbond_database_ = HBondDatabase::get_database( hbond_options_->params_database_tag() );
}

void
PolarHydrogenPacker::get_best_hxyz(
	conformation::Residue const & residue,
	utility::vector1< Vector > const & ideal_hydrogen_xyz_positions,
	Size const abase,
	Size const abase2,
	Vector const & donor_xyz,
	Real & best_score,
	Vector & best_hydrogen_xyz
) {
	using namespace core::kinematics;

	for ( Vector const & ideal_hydrogen_xyz : ideal_hydrogen_xyz_positions ) {
		if ( possible_hbond_acceptors_.size() == 0 ) continue;

		// might be better to create a 'pseudo' score function penalizing deviation from ideal.
		check_hbond_score( ideal_hydrogen_xyz, donor_xyz, best_score, best_hydrogen_xyz );

		Stub const ideal_hydrogen_stub( residue.xyz( abase ) /*center*/, ideal_hydrogen_xyz, residue.xyz( abase ), residue.xyz( abase2 ) );

		// generate 'rotamers' around ideal H location.
		Vector const ideal_local_xyz = ideal_hydrogen_stub.global2local( ideal_hydrogen_xyz );
		Distance const ideal_H_length = ideal_local_xyz.x();
		int const ndist( 5 );
		Distance const length_deviation( 0.05 );
		Size const ntheta( 2 );
		Real const theta_deviation( 10.0 ); // in degrees
		Size const nphi( 6 ); // deviation will be 360.0 / nphi
		for ( int d = -(ndist - 1)/2; d <= ( ndist - 1 )/2; d++ ) {
			Real const bond_length = ideal_H_length + length_deviation * Real( d );

			for ( Size q = 1; q <= ntheta; q++ ) { // theta (angle)
				Real const theta = q * numeric::conversions::radians( theta_deviation );

				for ( Size f = 0; f < nphi; f++ ) { // phi (azimuthal)
					Real const phi = f * numeric::conversions::radians( 360.0 / Real(nphi) );
					Real const x = bond_length * cos( theta );
					Real const y = bond_length * sin( theta )  * sin( phi );
					Real const z = bond_length * sin( theta )  * cos( phi );
					Vector const H_xyz = ideal_hydrogen_stub.local2global( Vector( x,y,z) );
					//     TR << residue.name1() << residue.seqpos() << " " << residue.atom_name( j ) << " current deviation " << ( ideal_hydrogen_xyz - H_xyz ).length() << std::endl;

					check_hbond_score( H_xyz, donor_xyz, best_score, best_hydrogen_xyz );
				}
			}
		}
	} // ideal H_xyz
}

utility::vector1< Vector >
PolarHydrogenPacker::get_ideal_hxyz_positions( conformation::Residue const & residue, conformation::Residue const & ideal_res, Size const j, kinematics::Stub const & current_input_stub, kinematics::Stub const & ideal_input_stub ) {
	utility::vector1< Vector > ideal_hydrogen_xyz_positions;

	Size proton_chi_no = check_if_proton_chi_atom( residue.type(), j );
	if ( proton_chi_no == 0 ) {
		Vector const ideal_hydrogen_xyz = current_input_stub.local2global( ideal_input_stub.global2local( ideal_res.xyz( j ) ) );
		ideal_hydrogen_xyz_positions.push_back( ideal_hydrogen_xyz );
	} else {
		if ( disallow_vary_geometry_proton_chi_ ) return ideal_hydrogen_xyz_positions;
		// i.e., 10, 20  to model -40 -50 -60 -70 -80 etc.
		utility::vector1< Real > samples = residue.type().proton_chi_samples( proton_chi_no );
		utility::vector1< Real > const & extra_samples = residue.type().proton_chi_extra_samples( proton_chi_no );
		for ( Real const sample : samples ) {
			for ( Real const extra_sample : extra_samples ) {
				samples.push_back( sample - extra_sample );
				samples.push_back( sample + extra_sample );
			}
		}
		for ( Real const sample : samples ) {
			//Vector ideal_hydrogen_local = ideal_input_stub.global2local( ideal_res.xyz( j ) );
			Vector const ideal_hydrogen_local = rotation_matrix( Vector( 1.0, 0.0, 0.0 ), sample ) * ideal_input_stub.global2local( ideal_res.xyz( j ) ); //rotation about x
			//Vector const ideal_hydrogen_xyz = current_input_stub.local2global( ideal_hydrogen_local );
			ideal_hydrogen_xyz_positions.push_back( current_input_stub.local2global( ideal_hydrogen_local ) );
		}
	}

	return ideal_hydrogen_xyz_positions;
}

void
PolarHydrogenPacker::virtualize_poor_scoring_o2prime_hydrogens( pose::Pose & pose, Size const i, Size const j, Real const best_score ) {
	if ( !pose.residue_type( i ).is_RNA() ) return;
	if ( pose.residue_type( i ).RNA_info().ho2prime_index() != j  ) return;
	//     TR << "2'OH score for residue " << i << " ==> " << best_score << std::endl;
	if ( best_score > -0.1 /* cutoff */ ) {
		add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_O2PRIME_HYDROGEN, i  );
	} else {
		remove_variant_type_from_pose_residue( pose, chemical::VIRTUAL_O2PRIME_HYDROGEN, i );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
void
PolarHydrogenPacker::apply(
#ifdef GL_GRAPHICS
    core::pose::Pose & pose_to_visualize
#else
	core::pose::Pose & pose
#endif
) {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::kinematics;
#ifdef GL_GRAPHICS
	core::pose::Pose pose = pose_to_visualize; // make a local copy -- otherwise crashing graphics builds
#endif
	core::Size const nres( pose.total_residue() );

	for ( core::Size i = 1; i <= nres; i++ )  {

		ResidueOP concrete_res_op = remove_variant_type_from_residue( pose.residue( i ), chemical::VIRTUAL_O2PRIME_HYDROGEN, pose );
		Residue const & residue = *concrete_res_op;

		ResidueOP ideal_res_op = conformation::ResidueFactory::create_residue( residue.type() );
		Residue const & ideal_res = *ideal_res_op;

		for ( core::Size j = 1; j <= residue.natoms(); j++ ) {

			if ( !residue.Hpos_polar().has_value( j ) )  continue; // just polar hydrogens

			// move this to a different function!
			// define a coordinate system
			Size const j1 = residue.atom_base( j );
			Size const j2 = residue.atom_base( j1 );
			Size const j3 = residue.atom_base( j2 );
			runtime_assert( !residue.atom_is_hydrogen( j1 ) );
			runtime_assert( !residue.atom_is_hydrogen( j2 ) );
			runtime_assert( !residue.atom_is_hydrogen( j3 ) );

			Stub const current_input_stub( residue.xyz( j1 ), residue.xyz( j1 ), residue.xyz( j2 ), residue.xyz( j3 ) );
			Vector const current_hydrogen_xyz = residue.xyz( j );
			Vector const donor_xyz            = residue.xyz( j1 );

			Real best_score( 0.0 );
			Vector best_hydrogen_xyz( 0.0 );

			// find potential hydrogen bond acceptors.
			get_possible_hbond_acceptors( pose, i, j1 /* donor atm*/ );

			check_hbond_score( current_hydrogen_xyz, donor_xyz, best_score, best_hydrogen_xyz );

			Stub const ideal_input_stub( ideal_res.xyz( j1 ), ideal_res.xyz( j1 ), ideal_res.xyz( j2 ), ideal_res.xyz( j3 ) );

			utility::vector1< Vector > ideal_hydrogen_xyz_positions = get_ideal_hxyz_positions( residue, ideal_res, j, current_input_stub, ideal_input_stub );
			if ( disallow_vary_geometry_proton_chi_ && ideal_hydrogen_xyz_positions.empty() ) continue;

			get_best_hxyz( residue, ideal_hydrogen_xyz_positions, j1, j2, donor_xyz, best_score, best_hydrogen_xyz );
			//    TR << residue.name1() << residue.seqpos() << " " << residue.atom_name( j ) << " current deviation " << ( best_hydrogen_xyz - current_hydrogen_xyz).length() << " NUM ACCEPTORS " << possible_hbond_acceptors_.size() << " " << best_score << std::endl;
			// Have to supply name -- residue type may have changed...
			pose.set_xyz( id::NamedAtomID( residue.atom_name( j ), i ), best_hydrogen_xyz );

			// move this to a different function!
			if ( allow_virtual_o2prime_hydrogens_ ) {
				virtualize_poor_scoring_o2prime_hydrogens( pose, i, j, best_score );
			}
		} // j
	} // i

#ifdef GL_GRAPHICS
	pose_to_visualize = pose;
#endif
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
PolarHydrogenPacker::check_hbond_score( Vector const & H_xyz,
	Vector const & D_xyz,
	Real & best_score,
	Vector & best_hydrogen_xyz ){

	Real total_score( 0.0 );
	for ( Size n = 1; n <= possible_hbond_acceptors_.size(); n++ ) {
		Vector const & A_xyz  = possible_hbond_acceptors_[n][1];
		Vector const & B_xyz  = possible_hbond_acceptors_[n][2];
		Vector const & B2_xyz = possible_hbond_acceptors_[n][3];
		Real energy( 0.0 );
		hb_energy_deriv( *hbond_database_, *hbond_options_,
			hb_eval_tuples_[n], D_xyz, H_xyz, A_xyz, B_xyz, B2_xyz,
			energy );
		if ( energy < 0.0 ) total_score += energy;
	}
	if ( total_score < best_score || ( best_score == 0.0 && best_hydrogen_xyz == Vector( 0.0 ) ) ) {
		best_score = total_score;
		best_hydrogen_xyz = H_xyz;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
PolarHydrogenPacker::get_possible_hbond_acceptors( pose::Pose const & pose, Size const moving_res, Size const atomno ) {

	using namespace core::conformation;
	Distance contact_distance_cutoff_( 3.5 );
	possible_hbond_acceptors_.clear();
	hb_eval_tuples_.clear();

	Vector const moving_xyz = pose.residue( moving_res ).xyz( atomno );

	core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( i == moving_res ) continue;

		// this is silly, trying to prevent i, i+1 backbone H-bonds in RNA.
		if ( i > 1 &&
				pdb_info->number( i )-1 == pdb_info->number( moving_res ) &&
				pdb_info->chain( i ) == pdb_info->chain( moving_res ) ) continue;

		Residue const & rsd = pose.residue( i );
		if ( rsd.is_virtual_residue() ) continue;
		for ( Size j = 1; j <= rsd.nheavyatoms(); j++ ) {
			if ( pose.residue_type( i ).is_virtual( j ) ) continue;
			if ( !pose.residue_type( i ).heavyatom_is_an_acceptor( j ) ) continue;

			Distance dist = ( rsd.xyz( j ) - moving_xyz ).length();
			if ( dist >= contact_distance_cutoff_ ) continue;

			// check if there's already a H-bond "back" from this acceptor -- happens in, e.g.,
			//  2'-OH/2'-OH pairs for RNA. Without this check, get a 'Mexican standoff':
			//
			//     O-H ... H-O
			//
			if ( pose.residue_type( i ).heavyatom_has_polar_hydrogens( j ) ) {
				bool H_pointing_back( false );
				for ( Size jj = pose.residue_type( i ).attached_H_begin( j );
						jj <= pose.residue_type( i ).attached_H_end( j ); jj++ ) {
					if ( angle_of( moving_xyz, rsd.xyz( j ), rsd.xyz( jj ) ) < 1.0 /* radians, so ~60 degrees */ ) {
						//TR << "WATCH OUT! " << pose.pdb_info()->number(i) << " " << rsd.atom_name( j )  << "already has an H that points back to " << pose.pdb_info()->number(moving_res) << " " << pose.residue( moving_res ).atom_name( atomno ) << std::endl;
						H_pointing_back = true; break;
					}
				}
				if ( H_pointing_back ) continue;
			}

			vector1< Vector > acceptor_xyz_info;
			acceptor_xyz_info.push_back( rsd.xyz( j ) );
			acceptor_xyz_info.push_back( rsd.xyz( rsd.atom_base(j) ) );
			acceptor_xyz_info.push_back( rsd.xyz( rsd.abase2( j ) ) );
			possible_hbond_acceptors_.push_back( acceptor_xyz_info );
			hb_eval_tuples_.push_back( HBEvalTuple( atomno,
				pose.residue( moving_res ),
				j, pose.residue( i) ) );
		}
	}
}


} //polar_hydrogens
} //modeler
} //stepwise
} //protocols
