// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file metrics
/// @brief protocols that are specific to docking low resolution
/// @details
/// @author Brian Weitzner

// Unit Headers

// Package Headers

// Project Headers

// Utility Headers

// Numeric Headers and ObjexxFCL Headers

// C++ headers

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <protocols/membrane/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/metrics.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

//Auto Headers
#include <protocols/simple_filters/DdgFilter.hh>

#if (defined WIN32) && (!defined WIN_PYROSETTA)
#undef interface
#endif

static basic::Tracer TR( "protocols.docking.DockingProtocol.metrics" );

using namespace core;

namespace protocols {
namespace docking {

/// @brief Compute interface score between the input complex and complex with
///   the two partners at 500A from each other
core::Real calc_intf_score( const core::pose::Pose & pose, const core::scoring::ScoreFunctionOP scorefxn, int const jump ) {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::rigid;
	using namespace protocols::membrane;

	// initialize pose
	Pose unbound_pose = pose;

	// score pose, this is bound score
	core::Real bound_score = ( *scorefxn )( unbound_pose );

	// initialize new pose and move axis
	core::Vector axis = membrane_axis( unbound_pose, jump );

	// move partners apart
	RigidBodyTransMoverOP mover( new RigidBodyTransMover( axis, jump ) );
	mover->step_size( 500 );
	mover->apply( unbound_pose );

	// score pose, this is unbound score
	core::Real unbound_score = ( *scorefxn )( unbound_pose );

	// compute interface score
	core::Real intf_score = bound_score - unbound_score;

	return intf_score;
}


core::Real
calc_interaction_energy( const core::pose::Pose & pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps ) {
	using namespace scoring;
	using namespace basic::options;

	Real interaction_energy( 0 );

	core::scoring::ScoreFunctionOP docking_scorefxn;
	core::pose::Pose complex_pose = pose;

	docking_scorefxn = dock_scorefxn->clone() ;
	docking_scorefxn->set_weight( core::scoring::atom_pair_constraint, 0.0 );
	// docking_scorefxn->set_weight( core::scoring::dslf_ss_dst, 0.0 );
	/*
	if ( pose.is_fullatom() ){
	docking_scorefxn = new core::scoring::ScoreFunction( *docking_score_high_ ) ;
	} else {
	docking_scorefxn = new core::scoring::ScoreFunction( *docking_score_low_ ) ;
	}
	*/

	// calculate energy of complexed pose
	Real const bound_energy = (*docking_scorefxn)( complex_pose );

	// calculate energy of separated pose over each movable jump
	// ddG is the "right" way to do this, to properly penalize strained rotamers
	// but aroop reports that I_sc yields better results for antibodies
	for ( Size rb_jump : movable_jumps ) {
		/*
		Real const threshold = 100000; // dummy threshold
		Size const repeats = 3;
		protocols::simple_filters::DdgFilter ddg = protocols::simple_filters::DdgFilter( threshold, docking_scorefxn, rb_jump, repeats );
		interaction_energy += ddg.compute( pose );
		*/
		core::pose::Pose unbound_pose = complex_pose;
		Real trans_magnitude = 1000;

		// create new translation mover
		rigid::RigidBodyTransMoverOP translate_away;

		// for membrane proteins
		if ( option[ OptionKeys::mp::setup::spanfiles ].user() ) {

			using namespace protocols::membrane;

			// get membrane axis
			core::Vector trans_axis( membrane_axis( unbound_pose, rb_jump ) );
			TR << "trans_axis: " << trans_axis.to_string() << std::endl;

			// create new translation mover
			translate_away = rigid::RigidBodyTransMoverOP( new rigid::RigidBodyTransMover(trans_axis, rb_jump) );
		} else {
			// for non-membrane proteins
			translate_away = rigid::RigidBodyTransMoverOP( new rigid::RigidBodyTransMover( unbound_pose, rb_jump ) );
		}

		// set translation magnitude and translate away
		translate_away->step_size( trans_magnitude );
		translate_away->apply( unbound_pose );

		// calculate unbound energy
		Real const unbound_energy = ( *docking_scorefxn )( unbound_pose );

		// add score from this jump to total interaction energy
		interaction_energy += ( bound_energy - unbound_energy );

		TR << "unbound pose: " << std::endl;
		docking_scorefxn->show( TR, unbound_pose );
		TR << "bound pose: " << std::endl;
		docking_scorefxn->show( TR, complex_pose );
		TR << "unbound energy: " << unbound_energy << std::endl;
		TR << "bound energy: " << bound_energy << std::endl;
		TR << "interaction energy: " << interaction_energy << std::endl;
		TR << "rb_jump: " << rb_jump << std::endl;
	}

	return interaction_energy;
}

core::Real
calc_Lrmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, DockJumps const movable_jumps ){
	using namespace scoring;
	Real Lrmsd( 0 );

	for ( Size rb_jump : movable_jumps ) {
		ObjexxFCL::FArray1D_bool temp_part( pose.size(), false );
		ObjexxFCL::FArray1D_bool superpos_partner( pose.size(), false );
		/// this gets the wrong partner, therefore it is stored in a temporary
		/// array and then the opposite is put in the actualy array that is used
		/// for superpositioning.  there is probably a better way to do this
		/// need to check TODO
		pose.fold_tree().partition_by_jump( rb_jump, temp_part );
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( temp_part( i ) ) superpos_partner( i )=false;
			else superpos_partner( i ) = true;
		}

		//Lrmsd += rmsd_no_super_subset( native_pose, pose, superpos_partner, is_protein_CA );
		using namespace core::scoring;
		Lrmsd += core::scoring::rmsd_no_super_subset( native_pose, pose, superpos_partner, is_protein_backbone );
	}
	return Lrmsd;
}

core::Real
calc_P1rmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, DockJumps const movable_jumps ){
	using namespace core::scoring;
	Real P1rmsd( 0 );

	for ( Size rb_jump : movable_jumps ) {
		ObjexxFCL::FArray1D_bool temp_part( pose.size(), false );
		ObjexxFCL::FArray1D_bool superpos_partner( pose.size(), false );

		/// this gets the wrong partner, therefore it is stored in a temporary
		/// array and then the opposite is put in the actualy array that is used
		/// for superpositioning.  there is probably a better way to do this
		/// need to check TODO
		pose.fold_tree().partition_by_jump( rb_jump, temp_part );
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( temp_part( i ) ) superpos_partner( i )=true;
			else superpos_partner( i ) = false;
		}

		P1rmsd += core::scoring::rmsd_with_super_subset( native_pose, pose, superpos_partner, is_protein_backbone );
	}
	return P1rmsd;
}

core::Real
calc_P2rmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, DockJumps const movable_jumps ){
	using namespace core::scoring;
	Real P2rmsd( 0 );

	for ( Size rb_jump : movable_jumps ) {
		ObjexxFCL::FArray1D_bool temp_part( pose.size(), false );
		ObjexxFCL::FArray1D_bool superpos_partner( pose.size(), false );

		/// this gets the wrong partner, therefore it is stored in a temporary
		/// array and then the opposite is put in the actualy array that is used
		/// for superpositioning.  there is probably a better way to do this
		/// need to check TODO
		pose.fold_tree().partition_by_jump( rb_jump, temp_part );
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( temp_part( i ) ) superpos_partner( i )=false;
			else superpos_partner( i ) = true;
		}

		P2rmsd += core::scoring::rmsd_with_super_subset( native_pose, pose, superpos_partner, is_protein_backbone );
	}
	return P2rmsd;
}

core::Real
calc_Irmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps ) {
	using namespace scoring;
	Real Irmsd( 0 );

	for ( core::Size rb_jump : movable_jumps ) {
		if ( !pose.is_fullatom() || !native_pose.is_fullatom() ) {
			TR << "Irmsd calc called with non-fullatom pose!!!" << std::endl;
			return 0.0;
		}

		core::pose::Pose native_docking_pose = native_pose;
		using namespace kinematics;
		FoldTree ft( pose.fold_tree() );
		native_docking_pose.fold_tree(ft);

		// score to set up interface object
		// scoring only happened here to update the residue neighbors
		core::scoring::ScoreFunctionOP scorefxn = dock_scorefxn->clone() ;
		( *scorefxn )( native_docking_pose );
		//  dock_scorefxn->show( TR );

		//  native_docking_pose.update_residue_neighbors();

		protocols::scoring::Interface interface( rb_jump );
		interface.distance( 8.0 );
		interface.calculate( native_docking_pose );
		ObjexxFCL::FArray1D_bool is_interface ( pose.size(), false );

		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( interface.is_interface( i ) ) is_interface( i ) = true;
		}

		//Irmsd += rmsd_with_super_subset(native_docking_pose, pose, is_interface, is_heavyatom);
		using namespace core::scoring;
		Irmsd += core::scoring::rmsd_with_super_subset( native_docking_pose, pose, is_interface, is_protein_backbone );
	}
	return Irmsd;
}

core::Real
calc_CA_Irmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps ) {
	using namespace scoring;
	Real Irmsd( 0 );

	for ( core::Size rb_jump : movable_jumps ) {
		//   if (!pose.is_fullatom() || !native_pose.is_fullatom()){
		//    TR << "Irmsd calc called with non-fullatom pose!!!"<<std::endl;
		//    return 0.0;
		//   }

		core::pose::Pose native_docking_pose = native_pose;
		using namespace kinematics;
		FoldTree ft( pose.fold_tree() );
		native_docking_pose.fold_tree( ft );

		//score to set up interface object
		core::scoring::ScoreFunctionOP scorefxn = dock_scorefxn->clone();
		( *scorefxn )( native_docking_pose );

		protocols::scoring::Interface interface( rb_jump );
		interface.distance( 8.0 );
		interface.calculate( native_docking_pose );
		ObjexxFCL::FArray1D_bool is_interface ( pose.size(), false );

		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( interface.is_interface( i ) ) is_interface( i ) = true;
		}

		//Irmsd += rmsd_with_super_subset(native_docking_pose, pose, is_interface, is_heavyatom);
		using namespace core::scoring;
		Irmsd += core::scoring::rmsd_with_super_subset( native_docking_pose, pose, is_interface, is_protein_CA );
	}
	return Irmsd;
}

core::Real
calc_Fnat( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps ) {
	using namespace scoring;
	using namespace conformation;
	Real Fnat( 0 );

	for ( core::Size rb_jump : movable_jumps ) {
		if ( !pose.is_fullatom() || !native_pose.is_fullatom() ) {
			TR << "Fnat calc called with non-fullatom pose!!!" << std::endl;
			return 0.0;
		}

		core::pose::Pose native_docking_pose = native_pose;
		using namespace kinematics;
		FoldTree ft( pose.fold_tree() );
		native_docking_pose.fold_tree( ft );
		Real cutoff = 5.0;

		//score to set up interface object
		core::scoring::ScoreFunctionOP scorefxn = dock_scorefxn->clone();
		(*scorefxn)( native_docking_pose );

		ObjexxFCL::FArray1D_bool temp_part ( pose.size(), false );
		pose.fold_tree().partition_by_jump( rb_jump, temp_part );

		utility::vector1< Size > partner1;
		utility::vector1< Size > partner2;

		protocols::scoring::Interface interface( rb_jump );
		interface.distance( 8.0 );
		interface.calculate( native_docking_pose );

		//generate list of interface residues for partner 1 and partner 2
		for ( Size i = 1; i <= pose.size(); i++ ) {
			if ( interface.is_interface( i ) ) {
				if ( !temp_part( i ) ) partner1.push_back( i );
				if ( temp_part( i ) ) partner2.push_back( i );
			}
		}

		//create contact pair list
		ObjexxFCL::FArray2D_bool contact_list( partner1.size(), partner2.size(), false );

		//identify native contacts across the interface
		//this will probably be changed to use PoseMetrics once I figure out how to use them - Sid
		for ( Size i = 1; i <= partner1.size(); i++ ) {
			ResidueOP rsd1( new Residue( native_docking_pose.residue( partner1[ i ] ) ) );
			for ( Size j = 1; j <= partner2.size(); j++ ) {
				ResidueOP rsd2( new Residue( native_docking_pose.residue( partner2[ j ] ) ) );
				contact_list( i, j ) = calc_res_contact( rsd1, rsd2, cutoff );
			}
		}

		Real native_ncontact = 0;
		Real decoy_ncontact = 0;
		//  utility::io::ozstream out( "native_contact.txt");
		//identify which native contacts are recovered in the decoy
		for ( Size i = 1; i <= partner1.size(); i++ ) {
			for ( Size j = 1; j <= partner2.size(); j++ ) {
				if ( contact_list( i, j ) ) {
					//  out << partner1[i] << "  " << partner2[j] << "\n";
					native_ncontact++;
					ResidueOP rsd1( new Residue( pose.residue( partner1[ i ] ) ) );
					ResidueOP rsd2( new Residue( pose.residue( partner2[ j ] ) ) );
					if ( calc_res_contact( rsd1, rsd2, cutoff ) ) decoy_ncontact++;
				}
			}
		}

		Fnat += decoy_ncontact / native_ncontact;
	}
	return Fnat;
}

core::Real
calc_Fnat( const core::pose::Pose & pose, std::string const& list_file, DockJumps const movable_jumps ) {
	using namespace scoring;
	using namespace conformation;
	Real Fnat( 0 );
	Real cutoff = 5.0;

	for ( auto it = movable_jumps.begin(); it != movable_jumps.end(); ++it ) {
		if ( !pose.is_fullatom() ) {
			TR << "Fnat calc called with non-fullatom pose!!!" << std::endl;
			return 0.0;
		}
		//  ObjexxFCL::FArray2D_bool contact_list(partner1.size(), partner2.size(), false);
		utility::vector1< Size > partner1;
		utility::vector1< Size > partner2;
		utility::io::izstream in( list_file );
		if ( !list_file.size() ) {
			TR.Error << "no input file" << std::endl;
			return 0.0;
		}
		if ( !in.good() ) {
			TR.Error << "cannot open file " << list_file << std::endl;
			return 0.0;
		}

		std::string line;
		while ( getline( in, line ) ) {
			std::istringstream line_stream( line );
			Size res_n1;
			Size res_n2;
			line_stream >> res_n1 >> res_n2;
			partner1.push_back( res_n1 );
			partner2.push_back( res_n2 );
			//   contact_list( res_n1, res_n2 ) = true;
		}
		//create contact pair list
		if ( partner1.size() != partner2.size() ) {
			TR.Error << "two columns in the list should be of same length!" << std::endl;
			return 0.0;
		}
		Real native_ncontact = partner1.size();
		Real decoy_ncontact = 0;

		//identify which native contacts are recovered in the decoy
		for ( Size i = 1; i <= partner1.size(); i++ ) {
			TR.Debug << partner1[ i ] << "  " << partner2[ i ] << std::endl;
			ResidueOP rsd1( new Residue( pose.residue( partner1[ i ] ) ) );
			ResidueOP rsd2( new Residue( pose.residue( partner2[ i ] ) ) );
			if ( calc_res_contact( rsd1, rsd2, cutoff ) ) decoy_ncontact++;
		}

		Fnat += decoy_ncontact/native_ncontact;
	}
	return Fnat;
}

core::Real
calc_Fnonnat( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps ){
	using namespace scoring;
	using namespace conformation;
	Real Fnonnat( 0 );

	for ( core::Size rb_jump : movable_jumps ) {
		if ( !pose.is_fullatom() || !native_pose.is_fullatom() ) {
			TR << "Fnat calc called with non-fullatom pose!!!" << std::endl;
			return 0.0;
		}

		core::pose::Pose native_docking_pose = native_pose;
		using namespace kinematics;
		FoldTree ft( pose.fold_tree() );
		native_docking_pose.fold_tree(ft);
		Real cutoff = 5.0;

		//score to set up interface object
		core::scoring::ScoreFunctionOP scorefxn = dock_scorefxn->clone();
		(*scorefxn)( native_docking_pose );

		ObjexxFCL::FArray1D_bool temp_part ( pose.size(), false );
		pose.fold_tree().partition_by_jump( rb_jump, temp_part );

		utility::vector1< Size > partner1;
		utility::vector1< Size > partner2;

		protocols::scoring::Interface interface( rb_jump );
		interface.distance( 15.0 ); // 8.0 as in Fnat is too limited. It can actually exclude some true native contacts. any special reason?
		interface.calculate( native_docking_pose );

		//generate list of interface residues for partner 1 and partner 2
		Size cutpoint = 0;
		for ( Size i = 1; i < pose.size(); i++ ) {
			if ( !temp_part( i ) ) {
				cutpoint = i;
				break;
			}
		}
		TR.Debug << "start residue id of second binding partner " << cutpoint << std::endl;

		for ( Size i = 1; i <= pose.size(); i++ ) {
			if ( interface.is_interface( i ) ) {
				if ( !temp_part( i ) ) partner1.push_back( i );
				if ( temp_part( i ) ) partner2.push_back( i );
			}
		}

		//create contact pair list
		std::list< std::pair< core::Size, core::Size > > contact_list;
		//identify native contacts across the interface
		//this will probably be changed to use PoseMetrics once I figure out how to use them - Sid
		for ( Size i = 1; i <= partner1.size(); i++ ) {
			ResidueOP rsd1( new Residue( native_docking_pose.residue( partner1[ i ] ) ) );
			for ( Size j = 1; j <= partner2.size(); j++ ) {
				ResidueOP rsd2( new Residue( native_docking_pose.residue( partner2[ j ] ) ) );
				if ( calc_res_contact( rsd1, rsd2, cutoff) ) {
					std::pair< core::Size, core::Size > elem( partner1[ i ], partner2[ j ] );
					contact_list.push_back( elem );
				}
			}
		}

		Real native_ncontact = contact_list.size();
		Real decoy_n_noncontact = 0;

		//identify which native contacts are recovered in the decoy

		for ( Size i = 1; i < cutpoint-1; i++ ) {
			for ( Size j = cutpoint; j <= pose.size(); j++ ) {
				ResidueOP rsd1( new Residue( pose.residue( i ) ) );
				ResidueOP rsd2( new Residue( pose.residue( j ) ) );
				if ( calc_res_contact( rsd2, rsd1, cutoff ) ) {
					std::pair< core::Size, core::Size > elem( j, i );
					//     for ( std::list< std::pair< core::Size, core::Size> >::iterator it = contact_list.begin(); it != contact_list.end(); ++it ) {
					auto it = std::find( contact_list.begin(), contact_list.end(), elem );
					if ( it != contact_list.end() ) { // found in the contact list
						TR.Debug << "true positive contact pair " << j << " " << i << std::endl;
					} else {
						decoy_n_noncontact++;
						TR.Debug << "false positie contacts pairs " << j << " " << i << std::endl;
					}
				}
			}
		}

		Fnonnat += decoy_n_noncontact / native_ncontact;
	}
	return Fnonnat;
}


core::Real
calc_Fnonnat( const core::pose::Pose & pose, std::string const& list_file, DockJumps const movable_jumps ) {
	using namespace scoring;
	using namespace conformation;
	Real Fnonnat( 0 );
	Real cutoff = 5.0;

	for ( core::Size rb_jump : movable_jumps ) {
		if ( !pose.is_fullatom() ) {
			TR << "Fnonnat calc called with non-fullatom pose!!!" << std::endl;
			return 0.0;
		}

		ObjexxFCL::FArray1D_bool temp_part ( pose.size(), false );
		pose.fold_tree().partition_by_jump( rb_jump, temp_part );
		Size cutpoint = 0;
		for ( Size i = 1; i < pose.size(); i++ ) {
			if ( !temp_part( i ) ) {
				cutpoint = i;
				break;
			}
		}
		TR.Debug << "start residue id of second binding partner " << cutpoint << std::endl;
		//   utility::vector1< Size > partner1;
		//   utility::vector1< Size > partner2;

		std::list< std::pair< core::Size, core::Size > > contact_list;

		utility::io::izstream in( list_file );
		if ( !list_file.size() ) {
			TR.Error << "no input file" << std::endl;
			return 0.0;
		}
		if ( !in.good() ) {
			TR.Error << "cannot open file " << list_file << std::endl;
			return 0.0;
		}

		std::string line;
		while ( getline( in, line ) ) {
			std::istringstream line_stream( line );
			Size res_n1;
			Size res_n2;
			line_stream >> res_n1 >> res_n2;
			TR.Debug << "res_n1 " << res_n1 << " res_n2 " << res_n2 << std::endl;
			std::pair< core::Size, core::Size > elem( res_n1, res_n2 );
			contact_list.push_back( elem );
			//   contact_list( res_n1, res_n2 ) = true;
		}
		//create contact pair list
		//   if ( partner1.size() != partner2.size() ) {
		//    TR.Error << "two columns in the list should be of same length!" << std::endl;
		//    return 0.0;
		//   }

		//  ObjexxFCL::FArray2D_bool contact_list( pose.size()-cutpoint+1, cutpoint-1, false);
		//   std::list< std::pair< core::Size, core::Size> > contact_list;
		//   for ( Size i=1; i<=partner1.size(); i++ ) {
		//    contact_list( partner1[i], partner2[i] ) = true;
		//   }
		Real native_ncontact = contact_list.size();
		Real decoy_n_noncontact = 0;

		//identify which native contacts are recovered in the decoy

		for ( Size i = 1; i < cutpoint - 1; i++ ) {
			for ( Size j = cutpoint; j <= pose.size(); j++ ) {
				ResidueOP rsd1( new Residue( pose.residue( i ) ) );
				ResidueOP rsd2( new Residue( pose.residue( j ) ) );
				TR.Debug << "distance between residue pair " << i << " " << j << std::endl;
				if ( calc_res_contact( rsd2, rsd1, cutoff ) ) {
					std::pair< core::Size, core::Size > elem( j, i );
					//     for ( std::list< std::pair< core::Size, core::Size> >::iterator it = contact_list.begin(); it != contact_list.end(); ++it ) {
					auto it = std::find( contact_list.begin(), contact_list.end(), elem );
					if ( it != contact_list.end() ) { // found in the contact list
						TR.Debug << "true positive contact pair " << j << " " << i << std::endl;
					} else {
						decoy_n_noncontact++;
						TR.Debug << "false positie contacts pairs " << j << " " << i << std::endl;
					}
				}
			}
		}

		Fnonnat += decoy_n_noncontact / native_ncontact;
	}
	return Fnonnat;
}

bool calc_res_contact(
	conformation::ResidueOP rsd1,
	conformation::ResidueOP rsd2,
	Real dist_cutoff
)
{
	Real dist_cutoff2 = dist_cutoff * dist_cutoff;
	double dist2 = 9999.0;

	for ( Size m = 1; m <= rsd1->nheavyatoms(); m++ ) {
		for ( Size n = 1; n <= rsd2->nheavyatoms(); n++ ) {
			/* double */dist2 = rsd1->xyz( m ).distance_squared( rsd2->xyz( n ) );  //Is there a reason this is a double?
			if ( dist2 <= dist_cutoff2 ) { TR.Debug << "return true " << dist2 << std::endl; return true; }
		}
	}
	TR.Debug << "return false " << dist2 << std::endl;
	return false;
}

}//docking
}//protocols
