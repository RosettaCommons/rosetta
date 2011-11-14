// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/annotated_sequence.cc
/// @brief  utility functions for making poses from sequences
/// @author P. Douglas Renfrew
/// @author Sam Deluca

// Unit Headers
#include <core/pose/annotated_sequence.hh>

// Package Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

namespace core {
namespace pose {

static basic::Tracer tr("core.pose");

using namespace core;
using namespace core::conformation;

////////////////////////////////////////////////////////////////////////////////
/// @details Given a protein sequence where each character represents an amino
/// acid, and a ResidueTypeSet, return the residue types that match the
/// sequence. NOTE: support making residue types from a fully annotated sequence
/// now, that is, for each residue variant or ligand which cannot be deduced
/// from one letter code directly, a [] is added directly following the one
/// letter code containig the residue's fullname, for example
/// K[lys_p:NtermProteinFull]ADFGCH[HIS_D]QNVE[glu_p:CtermProteinFull]Z[ZN].
/// This allows a pose to be constructed with full features from a silent output
/// file, such as with distiguished HIS tautomers, various chain termini and
/// cutpoint variants etc. Currently not working with disulfide variant CYD, but
/// this is on to-do list.
chemical::ResidueTypeCAPs residue_types_from_sequence(
	std::string const & sequence_in,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini /* true */
)
{
	chemical::ResidueTypeCAPs requested_types;

	using namespace core::chemical;

	if ( !sequence_in.size() ) return requested_types;

	// deal with the sequence read in; any non-standard protein AA name including lig should be put within a bracket[]
	// following the one-letter AA character. X for aa_vrt and Z for aa_unk
	std::string fullname;
	utility::vector1< std::string > fullname_list; // a vector of non-standard full names
	std::vector< Size > oneletter_to_fullname_index; // for each one-letter sequence, zero means no fullname given

	// we start with the first character in sequence and that should be a standard AA.
	std::string one_letter_sequence = sequence_in.substr( 0,1 );
	Size last_index = 0; // zero means this one-letter name does not have a fullname sepcified in bracket.
	bool in_bracket = false; // currently whether scanning fullname in bracket or not.

	for ( Size seqpos = 1; seqpos < sequence_in.length(); ++seqpos ) {
		// inside the bracket will be the base name of this residue;
		char aa = sequence_in[ seqpos ];

		// note that a full-name aa will also have its one-letter code present e.g. C[CYS]
		// hence the seqpos-count is not messed up
		if ( aa == '[' ) { // bracket starts, turn on flag and reset fullname string
			in_bracket = true;
			fullname = "";
			continue;
		} else if ( sequence_in[ seqpos ] == ']' ) { // bracket ends, save fullname and map its index
			in_bracket = false;
			fullname_list.push_back( fullname );
			last_index = fullname_list.size();
			continue;
		}

		if ( in_bracket ) { // in bracket, get fullname one char at a time
			fullname += aa;
			continue;
		} else { // outside bracket, save regular one-letter sequence.
			one_letter_sequence += aa;
			oneletter_to_fullname_index.push_back( last_index );
			last_index = 0;
		}
	} // finish reading in the whole sequence.

	oneletter_to_fullname_index.push_back( last_index );
	tr.Debug << "one_letter: " << one_letter_sequence << std::endl;
	tr.Debug << "seq_in: " << sequence_in << std::endl;

	// setup the pose by appending the appropriate residues residues
	for ( Size seqpos = 1; seqpos <= one_letter_sequence.length(); ++seqpos ) {
		char aa = one_letter_sequence[ seqpos-1 ]; // string indexing is zero-based!
		chemical::AA my_aa = chemical::aa_from_oneletter_code( aa );

		bool is_lower_terminus(false), is_upper_terminus(false);

		// is there an annotated fullname defined for this one-letter code?
		Size index = oneletter_to_fullname_index[ seqpos-1 ];
		if ( index ) { // fullname defined and get it directly from name_map
			// The next call requires reference -> CAP because ResidueTypeSet's
			// methods are not yet consistent in handing out ref vs CAP.
			requested_types.push_back( &residue_set.name_map( fullname_list[ index ] ) );
			is_lower_terminus = ( *requested_types.back() ).has_variant_type( chemical::LOWER_TERMINUS );
			is_upper_terminus = ( *requested_types.back() ).has_variant_type( chemical::UPPER_TERMINUS );
		} else {
			// use aa_map to find list of possible ResidueTypes
			chemical::ResidueTypeCAPs const & rsd_type_list( residue_set.aa_map( my_aa ) );
			// for non-annotated sequence, assume single chain for now
			is_lower_terminus = auto_termini && ( seqpos == 1 );
			is_upper_terminus = auto_termini && ( seqpos == one_letter_sequence.length() );
			bool const is_terminus( is_lower_terminus || is_upper_terminus ); // redundant, but for convenience

			Size best_index = 0;
			// iterate over rsd_types, pick one.
			for ( Size j = 1; j <= rsd_type_list.size(); ++j ) {
				chemical::ResidueType const & rsd_type( *(rsd_type_list[ j ]) );

				bool const is_polymer( rsd_type.is_polymer() );
				// pick a ResidueType
				Size nvariants = rsd_type.variant_types().size();
				if ( is_polymer && ( is_terminus && ( nvariants == 0 ) ) ) continue;
				if ( is_polymer && ( is_lower_terminus != rsd_type.has_variant_type( chemical::LOWER_TERMINUS ) ||
														 is_upper_terminus != rsd_type.has_variant_type( chemical::UPPER_TERMINUS ) ) ) continue;

				best_index = j;
				break;
			}
			if ( !best_index ) utility_exit_with_message( " can't find residue type at pos " + ObjexxFCL::string_of(seqpos) +
				"in sequence "+ sequence_in);
			// add the ResidueTypeCAP
			requested_types.push_back( rsd_type_list[ best_index ] );
		}

		tr.Trace << "residue_types_from_sequence():  seqpos: " << seqpos << " aa " << aa << " " << my_aa << std::endl;

	} // for seqpos

	return requested_types;
}


////////////////////////////////////////////////////////////////////////////////
/// @details Given a Pose, a protein sequence where each character represents an
/// amino acid, and a ResidueTypeSet, give the Pose a conformation of covalently
/// linked residues that match the sequence. NOTE: support making pose from a
/// fully annotated sequence now, that is, for each residue variant or ligand
/// which cannot be deduced from one letter code directly, a [] is added
/// directly following the one letter code containig the residue's fullname, e.g.
/// K[lys_p:NtermProteinFull]ADFGCH[HIS_D]QNVE[glu_p:CtermProteinFull]Z[ZN].
/// This allows a pose to be constructed with full features from a silent output
/// file, such as with distiguished HIS tautomers, various chain termini and
/// cutpoint variants etc. Currently not working with disulfide variant CYD, but
/// this is on to-do list.
void make_pose_from_sequence(
	pose::Pose & pose,
	std::string const & sequence_in,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini /* true */
)
{
	typedef core::Size Size;

	// grab residue types
	chemical::ResidueTypeCAPs requested_types = core::pose::residue_types_from_sequence( sequence_in, residue_set, auto_termini );
	assert( core::pose::annotated_to_oneletter_sequence( sequence_in ).length() == requested_types.size() );

	// clear the pose
	pose.clear();

	// make the pose
	bool jump_to_next = false;
	for ( Size i = 1, ie = requested_types.size(); i <= ie; ++i ) {
		// grab the new residue
		chemical::ResidueType const & rsd_type = *requested_types[ i ];
		core::conformation::ResidueOP new_rsd( NULL );
		new_rsd = conformation::ResidueFactory::create_residue( rsd_type );

		// yab 20090219: The following error check was in the original
		// code prior to the split into core::pose::residue_types_from_sequence()
		// and this function, but it doesn't appear to be triggerable
		// because ResidueFactory always returns a residue.  I leave it
		// in for now, but consider taking it out.
		if ( !new_rsd ) {
			std::cerr << "cannot create a residue that matches the residue type "
				<< rsd_type.name1() << " " << rsd_type.name() << " at position " << i << '\n';
			utility_exit_with_message( "make_pose_from_sequence fails\n" );
		}

		tr.Trace << "make_pose_from_sequence():  seqpos: " << i << " " << new_rsd->aa() << std::endl;

		// do the actual append
		if ( rsd_type.has_variant_type( chemical::LOWER_TERMINUS ) ||
				 rsd_type.has_variant_type( chemical::N_ACETYLATION ) ||
			new_rsd->aa() == chemical::aa_unk || new_rsd->aa() == chemical::aa_vrt ||
				jump_to_next ) {
			if ( new_rsd->aa() == chemical::aa_unk  || new_rsd->aa() == chemical::aa_vrt ) {
				//fpd tr.Warning << "found unknown aminoacid or X in sequence at position " << i <<  std::endl;
				//fpd if ( i< ie ) {
				//fpd 	utility_exit_with_message( "found unknown aminoacid or X in sequence\n this leads to a seg-fault if we keep going...\n");
				//fpd }

				// if you don't think so ... make the code more stable and remove this
				// but only if this sequence doesn't seg-fault: KPAFGTNQEDYASYIXNGIIK" );

				///fpd ^^^ the problem is that the residue following the X should be connected by a jump as well.
				///     it should be of LOWER_TERMINUS variant type, but if not, we'll recover & spit out a warning for now.
				///     same thing for ligands???
				jump_to_next = true;
			} else if ( jump_to_next ) {
				jump_to_next = false;
				if ( !rsd_type.has_variant_type( chemical::LOWER_TERMINUS ) )
					tr.Warning << "Residue following X, Z, or an upper terminus is _not_ a lower terminus type!  Continuing ..." << std::endl;
			}
			pose.append_residue_by_jump( *new_rsd, 1, "", "", true ); // each time this happens, a new chain should be started
		} else {
			pose.append_residue_by_bond( *new_rsd, true );

			//fpd If res i is an upper terminus but (i+1) is not a lower terminus, the code exits on a failed assertion
			//fpd Don't let this happen; always jump in these cases
			if (rsd_type.has_variant_type( chemical::UPPER_TERMINUS )) jump_to_next = true;
		}
	}

	tr.Debug << "sequence in pose: " << pose.sequence() << std::endl;
	tr.Debug << "annotated seq: " << pose.annotated_sequence() << std::endl;

} // core::pose::make_pose_from_sequence

////////////////////////////////////////////////////////////////////////////////
/// overloaded version of previous mak_pose_from_sequence, does the same
/// function, but reads in a string of the residue type set instead of a
/// ResidueTypeSet object.  Made for PyRosetta.
/// olange: DONT DUPLICATE CODE sid!  --- I removed the duplication by calling the original "core::pose::make_pose_from_sequence"
void make_pose_from_sequence(
	pose::Pose & pose,
	std::string const & sequence_in,
	std::string const & type_set_name,
	//chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini /* true */
) {
	chemical::ResidueTypeSetCAP residue_set( chemical::ChemicalManager::get_instance()->residue_type_set( type_set_name ) );
	core::pose::make_pose_from_sequence( pose, sequence_in, *residue_set, auto_termini );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string annotated_to_oneletter_sequence(
	std::string const & annotated_seq
) {
	bool add( true );
	std::string oneletter_seq;
	for ( Size i = 0; i < annotated_seq.length(); ++i ) {
		if ( annotated_seq.at(i) == '[' ) add = false;
		if ( add ) oneletter_seq += annotated_seq.at(i);
		if ( annotated_seq.at(i) == ']' ) add = true;
	}

	return oneletter_seq;
}


} // namespace core
} // namespace pose
