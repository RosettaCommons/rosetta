// Package Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

//Project Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <devel/init.hh>
#include <protocols/moves/PyMolMover.hh>


// Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

// iostream
#include <iostream>


/// parse the annotated sequence.
static void parse_sequence(
	std::string const & sequence_in,
	utility::vector1< std::string > & fullname_list) {
	
		fullname_list.clear();
		if ( sequence_in.empty() ) return;
		
		
		std::string fullname;
		
		core::Size last_index = 0;
	
		bool in_bracket = false; // currently whether scanning fullname in bracket or not.

		for ( core::Size seqpos = 0; seqpos < sequence_in.length(); ++seqpos ) {
		// inside the bracket will be the base name of this residue;
			char aa = sequence_in[ seqpos ];
		
			if ( aa == '[' ) { // bracket starts, turn on flag and reset fullname string
				in_bracket = true;
				fullname = "";
				continue;
			} else if ( aa == ']' ) { // bracket ends, save fullname and map its index
				in_bracket = false;
				fullname_list.push_back( fullname );
				last_index = fullname_list.size();
				continue;
			}

			if ( in_bracket ) { // in bracket, get fullname one char at a time
				fullname += aa;
				continue;
			} else { 
				std::cout << "Please check your peptoid sequence input format. This is going to crash in 3...2...1" << std::endl;
				continue;
			}	
		} // finish reading in the whole sequence.
}

/// @details Given a peptoid sequence where each three digit code represents an amino
/// acid, and a ResidueTypeSet, return the residue types that match the
/// sequence.
static core::chemical::ResidueTypeCOPs residue_types_from_sequence(
	std::string const & sequence_in,
	core::chemical::ResidueTypeSet const & residue_set ){
	
		core::chemical::ResidueTypeCOPs requested_types;

		if ( sequence_in.empty() ) return requested_types;

		utility::vector1< std::string > fullname_list; // a vector of non-standard full names
		parse_sequence( sequence_in, fullname_list );

		// setup the pose by appending the appropriate residues
		for ( core::Size seqpos = 1; seqpos <= fullname_list.size(); ++seqpos ) {
		
			if ( seqpos == 1){
			
				fullname_list[ seqpos ] += ":NtermPeptoidFull";
				
			} else if ( seqpos == fullname_list.size() ){
			
				fullname_list[ seqpos ] += ":CtermPeptoidFull";
			}
			requested_types.push_back( residue_set.name_map( fullname_list[ seqpos ] ).get_self_ptr() );
	}

	return requested_types;
}

int main(int argc, char ** argv){
	devel::init( argc, argv ); // Initializing Rosetta

	// Creating a blank pose
	core::pose::PoseOP peptoid_pose( new core::pose::Pose() );
	
	// Peptoid sequence to make pose from
	std::string const & sequence_in = "[450][450][450][450][450][450]";
	
	utility::vector1< std::string > fullname_list; // a vector of non-standard full names
	
	// Connecting to PyMOL listener	
	protocols::moves::AddPyMolObserver( *peptoid_pose, true, 0 );
	
	core::chemical::ResidueTypeSetCOP residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );
	
	core::chemical::ResidueTypeCOPs requested_types = residue_types_from_sequence( sequence_in, *residue_set );
	
	for ( core::Size i = 1, ie = requested_types.size(); i <= ie; ++i ) {
	
		core::chemical::ResidueType const & rsd_type = *requested_types[ i ];
	
		core::conformation::ResidueOP new_rsd( NULL );

		new_rsd = core::conformation::ResidueFactory::create_residue( rsd_type );
	
		if ( i == 1 ) {

			peptoid_pose->append_residue_by_jump( *new_rsd, 1, "", "", true );	
			
		} else {
			
			peptoid_pose->append_residue_by_bond(*new_rsd, true );		
		}
	}
	
	// Starting out with an extended structure
	for ( core::Size i=1; i<=(*peptoid_pose).n_residue() ; ++i){
		
		( *peptoid_pose ).set_phi( i, 180 );
		
		( *peptoid_pose ).set_psi( i, 180 );
		
		( *peptoid_pose ).set_omega( i, 180 );
	}
	
	std::cout << "The sequence of the pose is: " << ( *peptoid_pose ).annotated_sequence() << "." << std::endl;
		 
	// Scoring using the MM score function
	core::scoring::ScoreFunctionOP sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "mm_std" );
	core::Real score = ( *sfxn )( *peptoid_pose );
	std::cout << "The starting score is " << score << "." << std::endl;
	
	
	// Rotamer packing
	core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *peptoid_pose );
	repack_task->restrict_to_repacking();
	core::pack::pack_rotamers( *peptoid_pose, *sfxn, repack_task );
	std::cout << "The score after packing is " << score << "." << std::endl;
	
	
	return 0;
}

