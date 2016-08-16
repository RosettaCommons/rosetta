// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#if 0
void steal_constant_length_frag_set_from_pose ( pose::Pose const& pose, ConstantLengthFragSet& fragset ) {
  Size nbb ( 3 ); // three backbone torsions for Protein
  Size len = fragset.max_frag_length();
  kinematics::MoveMap dummy_movemap; //is ignored right now
  for ( Size pos = 1; pos <= pose.total_residue() - len + 1; ++pos ) {
    FragDataOP frag_raw = new FragData;
    for ( Size i = 1; i<= len; i++ ) {
      frag_raw->add_residue( new BBTorsionSRFD( nbb, pose.secstruct(pos+i-1), oneletter_code_from_aa( pose.residue( pos+i-1 ).aa() ) ) );
    };
    FrameList flist;
    FrameOP frame;
    bool frame_existed;
    if ( ( frame_existed = fragset.region_simple( (core::Size) pos, (core::Size) pos+len-1, flist )) ) {
      frame = flist.front();
    } else {
      frame = new Frame( pos, len );
    }

    frag_raw->steal( pose, *frame );
    frame->add_fragment ( frag_raw );
    if ( !frame_existed ) {
      fragset.add( frame );
    }
  };
}


void
core::pose::make_pose_from_sequence_(
	std::string sequence,
	chemical::ResidueTypeSet const& residue_set,
	pose::Pose& pose
) {
	using namespace chemical;
	// clear all of the old data in the pose
	pose.clear();

	// setup the pose by appending the appropriate residues residues
	for ( Size seqpos = 1; seqpos <= sequence.length(); ++seqpos ) {
		char aa = sequence[seqpos-1]; // string indexing is zero-based!
		AA my_aa = aa_from_oneletter_code( aa );
		ResidueTypeCOPs const & rsd_type_list( residue_set.aa_map( my_aa ) );
		Size best_index = 1;
		ResidueType const & rsd_type( *(rsd_type_list[ best_index ]) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
		if ( seqpos == 1 ) {
			pose.append_residue_by_jump( *new_rsd, 1 );
		} else {
			pose.append_residue_by_bond( *new_rsd, true );
		}
	} // for seqpos
	// pose.conformation().insert_chain_ending( pose.total_residue() - 1 );		// probably not necessary
} // make_pose_match_sequence_
#endif

#if 0
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details the function allows a pose to use a different residue_type_set to represent all its residues,
///such as from fullatom residues to centroid residues, or vice versa. During the switch, corresponding atoms
///will be copied. Redundant atoms will be removed (in case from fullatom to centroid) and missing atoms will be
///built by ideal geometry (in the case from centroid to fullatom).
void
core::util::switch_to_residue_type_set(
				 pose::Pose & pose,
				 std::string const & type_set_name
)
{
	using namespace core::chemical;
	using namespace core::conformation;

	// retrieve proper residue_type_set
	ResidueTypeSetCAP target_residue_type_set( ChemicalManager::get_instance()->residue_type_set( type_set_name ) );
	// loop each position and find new type that matches from the new type set
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		// in future we may have a conformation using mixed type set, so check this by residue
		std::string const & current_type_set_name ( rsd.type().residue_type_set().database_directory() );
		if ( current_type_set_name.find( type_set_name ) != std::string::npos ) {
			std::cerr << "core::util::switch_to_residue_type_set: residue " << i << "already in " << type_set_name
								<< " residue_type_set" << '\n';
			continue;
		}
		// get all residue types with same AA
		ResidueTypeCOPs const & rsd_types( target_residue_type_set->aa_map( rsd.aa() ) );
		ResidueOP new_rsd( 0 );
		// now look for a rsdtype with same variants
		for ( Size j=1; j<= rsd_types.size(); ++j ) {
			ResidueType const & new_rsd_type( *rsd_types[j] );
			if ( rsd.type().variants_match( new_rsd_type ) ) {
	new_rsd = ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
	break;
			}
		}
		if ( ! new_rsd ) {
			std::cerr << "can not find a residue type that matches the residue " << rsd.name()
		<< "at position " << i << '\n';
			utility_exit_with_message( "core::util::switch_to_residue_type_set fails\n" );
		}
		// switch to corresponding residue type in the new set.
		pose.replace_residue( i, *new_rsd, false );
	}
}
#endif
