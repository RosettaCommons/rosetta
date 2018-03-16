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
/// @author Labonte (carbohydrate versions)

// Unit Headers
#include <core/pose/annotated_sequence.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>

// Project Headers
#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/io/carbohydrates/pose_io.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <list>
#include <sstream>


namespace core {
namespace pose {

static basic::Tracer TR( "core.pose" );

using namespace core;
using namespace core::conformation;

////////////////////////////////////////////////////////////////////////////////
/// parse the annotated sequence.
void parse_sequence(
	std::string const & sequence_in,
	utility::vector1< std::string > & fullname_list,
	std::vector< Size > & oneletter_to_fullname_index,
	std::string & one_letter_sequence
) {
	fullname_list.clear();
	oneletter_to_fullname_index.clear();
	if ( sequence_in.empty() ) return;

	// deal with the sequence read in; any non-standard protein AA name including lig should be put within a bracket[]
	// following the one-letter AA character. X for aa_vrt and Z for aa_unk
	one_letter_sequence = sequence_in.substr( 0,1 );
	std::string fullname;

	// we start with the first character in sequence and that should be a standard AA.
	Size last_index = 0; // zero means this one-letter name does not have a fullname specified in bracket.
	bool in_bracket = false; // currently whether scanning fullname in bracket or not.

	bool old_cyd_found = false;
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
			// CYD -> CYS - this requires no support from other parts of the code
			if ( fullname.substr(0,3) == "CYD" ) {
				fullname[2] = 'S';
				old_cyd_found = true;
			}

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
	// Warn if an old CYD was found
	if ( old_cyd_found ) {
		TR.Warning << "Found old-style CYD nomenclature!  ";
		TR.Warning << "Renaming as CYS for Conformation::detect_disulfides to resolve" << std::endl;
	}

	oneletter_to_fullname_index.push_back( last_index );
}

////////////////////////////////////////////////////////////////////////////////
/// Get the length of the annotated sequence
Size get_sequence_len( std::string const & sequence_in ) {
	utility::vector1< std::string > fullname_list;
	std::vector< Size > oneletter_to_fullname_index;
	std::string one_letter_sequence;
	parse_sequence( sequence_in, fullname_list, oneletter_to_fullname_index, one_letter_sequence );
	return one_letter_sequence.size();
}


////////////////////////////////////////////////////////////////////////////////
/// @details Given a protein sequence where each character represents an amino
/// acid, and a ResidueTypeSet, return the residue types that match the
/// sequence. NOTE: support making residue types from a fully annotated sequence
/// now, that is, for each residue variant or ligand which cannot be deduced
/// from one letter code directly, a [] is added directly following the one
/// letter code containing the residue's fullname, for example
/// K[lys:NtermProteinFull]ADFGCH[HIS_D]QNVE[glu:CtermProteinFull]Z[ZN].
/// This allows a pose to be constructed with full features from a silent output
/// file, such as with distinguished HIS tautomers, various chain termini and
/// cutpoint variants etc.
chemical::ResidueTypeCOPs residue_types_from_sequence(
	std::string const & sequence_in,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini /* true */
) {
	chemical::ResidueTypeCOPs requested_types;

	using namespace core::chemical;

	if ( sequence_in.empty() ) return requested_types;

	utility::vector1< std::string > fullname_list; // a vector of non-standard full names
	std::vector< Size > oneletter_to_fullname_index; // for each one-letter sequence, zero means no fullname given
	std::string one_letter_sequence;
	parse_sequence( sequence_in, fullname_list, oneletter_to_fullname_index, one_letter_sequence );

	TR.Debug << "one_letter: " << one_letter_sequence << std::endl;
	TR.Debug << "seq_in: " << sequence_in << std::endl;

	// setup the pose by appending the appropriate residues
	for ( Size seqpos = 1; seqpos <= one_letter_sequence.length(); ++seqpos ) {
		char aa = one_letter_sequence[ seqpos-1 ]; // string indexing is zero-based!

		if ( aa == '-' ) {
			requested_types.push_back( nullptr );
			continue;
		}

		if ( aa == '/' ) continue;  //fpd: force a chainbreak

		chemical::AA my_aa = chemical::aa_from_oneletter_code( aa );

		bool is_lower_terminus(false), is_upper_terminus(false);

		// is there an annotated fullname defined for this one-letter code?
		Size index = oneletter_to_fullname_index[ seqpos-1 ];
		if ( index ) { // fullname defined and get it directly from name_map
			ResidueTypeCOP rsd_type;
			rsd_type = residue_set.name_mapOP( fullname_list[ index ] );

			if ( rsd_type == nullptr ) {
				utility_exit_with_message( " can't find residue type for " + fullname_list[ index ] + " at pos " + ObjexxFCL::string_of(seqpos) +
					" in sequence "+ sequence_in);
			}

			requested_types.push_back( rsd_type );

		} else {
			// for non-annotated sequence, assume single chain for now
			is_lower_terminus = auto_termini && ( seqpos == 1 || one_letter_sequence[ seqpos-2 ]=='/');
			is_upper_terminus = auto_termini && ( seqpos == one_letter_sequence.length() || one_letter_sequence[ seqpos ]=='/' );

			ResidueTypeCOP rsd_type = get_rsd_type_from_aa( residue_set, my_aa, is_lower_terminus, is_upper_terminus );

			if ( rsd_type == nullptr ) {
				utility_exit_with_message( " can't find residue type for " + std::string(1,aa) + " at pos " + ObjexxFCL::string_of(seqpos) +
					" in sequence "+ sequence_in);
			}

			// add the ResidueTypeCOP
			requested_types.push_back( rsd_type );
		}

		TR.Trace << "residue_types_from_sequence():  seqpos: " << seqpos << " aa " << aa << " " << my_aa << std::endl;

	} // for seqpos

	return requested_types;
}

/// @brief Internal implementation detail of residue_types_from_saccharide_sequence()
/// Returns the reorganized (but not reversed) set of residuetypes
chemical::ResidueTypeCOPs
residue_types_from_saccharide_sequence_recursive( std::string const & sequence, chemical::ResidueTypeSet const & residue_set );

// Return a list of carbohydrate ResidueTypes corresponding to an annotated, linear, IUPAC polysaccharide sequence.
/// @param[in] <sequence>: an annotated IUPAC polysaccharide sequence,
/// e.g., "alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp"
/// @param[in] <residue_set>: the desired residue set
/// @return    a 1-indexed vector of ResidueType owning pointers, in Rosetta order.
/// That is, the first residue will be the reducing end (the last residue in the sequence),
/// and all branches will be consecutive, with earlier branches earlier in set of residue types.
/// @details   Format for <sequence>:\n
/// Prefixes apply to the residue to which they are attached, below indicated by residue n.\n
/// Residues are listed from N to 1, where N is the total number of residues in the saccharide.\n
/// The sequence is parsed by reading to the next hyphen, so hyphens are crucial.\n
/// Linkage indication: "(a->x)-" specifies the linkage of residue n, where a is the anomeric carbon number of residue
/// (n+1) and x is the oxygen number of residue n.  The first residue listed in the annotated sequence (residue N)
/// need not have the linkage prefix.  A ->4) ResidueType will automatically be assigned by default if not specified.\n
/// Anomer indication: The strings "alpha-" or "beta-" are supplied next, which determines the stereochemistry of the
/// anomeric carbon of the residue to which it is prefixed.  An alpha ResidueType will automatically be assigned by
/// default.\n
/// Stereochemical indication: "L-" or "D-" specifies whether residue n is an L- or D-sugar.  The default is "D-".\n
/// 3-Letter code: A three letter code (in sentence case) MUST be supplied next.  This specifies the "base sugar name",
/// e.g., Glc is for glucose.  (A list of all recognized 3-letter codes for sugars can be found in the database.)\n
/// 1-Letter suffix: If no suffix follows, residue n will be linear.  If a letter is present, it indicates the ring
/// size, where "f" is furanose, "p" is pyranose, and "s" is septanose.\n
/// Branches are indicated using nested brackets and are best explained by example:\n
/// beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc is:\n
/// beta-D-Galp-(1->4)-D-GlcpNAc\n
///                       |\n
///     alpha-L-Fucp-(1->3)
/// @note make_pose_from_saccharide_sequence() will generate a pose with a proper lower terminus.
/// glycosylate_pose() will append the fragment by bond.
chemical::ResidueTypeCOPs
residue_types_from_saccharide_sequence( std::string const & sequence, chemical::ResidueTypeSet const & residue_set )
{
	core::chemical::ResidueTypeCOPs residue_types =
		residue_types_from_saccharide_sequence_recursive( sequence, residue_set );

	std::reverse( residue_types.begin(), residue_types.end() );

	return residue_types;
}

chemical::ResidueTypeCOPs
residue_types_from_saccharide_sequence_recursive(
	std::string const & sequence,
	chemical::ResidueTypeSet const & residue_set )
{
	using namespace std;
	using namespace utility;
	using namespace chemical;
	using namespace chemical::carbohydrates;

	ResidueTypeCOPs main_residue_types;
	ResidueTypeCOPs branch_residue_types;

	if ( ! sequence.size() ) {
		return main_residue_types;
	}

	// Add delimiter to end of sequence.
	string const sequence_with_hyphen( sequence + '-' );

	// Loop through sequence one character at a time, form affixes and 3-letter codes, and assign ResidueTypes.
	uint const sequence_end( sequence_with_hyphen.length() );
	string morpheme( "" );
	stringstream residue_type_name( stringstream::out );
	utility::vector1< uint > branch_points;
	bool linkage_assigned( false );
	bool anomer_assigned( false );
	bool L_or_D_assigned( false );
	for ( uint chr_num( 0 ); chr_num < sequence_end; ++chr_num ) {
		char character = sequence_with_hyphen[ chr_num ];

		if ( character == '[' ) {
			// Branch: Stop what we are doing, recursively call self with whatever's between the brackets.
			string branch_sequence( "" );
			Size open_brackets( 1 );
			while ( open_brackets ) {
				++chr_num;
				character = sequence_with_hyphen[ chr_num ];
				if ( character == '[' ) {
					++open_brackets;
				} else if ( character == ']' ) {
					--open_brackets;
					if ( ! open_brackets ) {  // We've reached the closing bracket of this branch.
						// Peek backwards to grab the connectivity of this branch.
						// The previous character must be a right parenthesis.
						if ( sequence_with_hyphen[ chr_num - 1 ] != ')' ) {
							utility_exit_with_message( "Saccharide sequence input error: "
								"a branch must specify its connectivity to its parent chain." );
						}
						branch_points.push_back( atoi( &sequence_with_hyphen[ chr_num - 2 ] ) );
						// Set the proper VariantType below.
						break;
					}
				}
				branch_sequence += character;
			}
			// Now recurse.
			branch_residue_types.append( residue_types_from_saccharide_sequence_recursive(
				branch_sequence, residue_set ) );
		} else if ( character != '-' ) {  // '-' is the morpheme delimiter
			morpheme += character;
		} else {  // Hyphen: The morpheme is complete; interpret it....
			// Linkage indication, first half (ignored)
			if ( morpheme[ 0 ] == '(' ) {
				morpheme = "";  // The "(_" information is not needed; continue on to the next morpheme.

			} else if ( morpheme == "" ) {
				// We just came out of a branch, the user typed two hyphens in a row,
				// or we are at the beginning of an arrow; ignore.
				;

				// Linkage indication, second half
			} else if ( morpheme[ 0 ] == '>' ) {
				if ( linkage_assigned ) {
					utility_exit_with_message( "Saccharide sequence input error: "
						"the linkage notation cannot be specified twice!" );
				} else if ( anomer_assigned || L_or_D_assigned ) {
					utility_exit_with_message( "Saccharide sequence input error: "
						"the linkage notation must precede other prefixes." );
				} else {
					residue_type_name << '-' << morpheme << '-';
					linkage_assigned = true;
					morpheme = "";
				}

				// Anomer indication
			} else if ( morpheme == "alpha" || morpheme == "beta" || morpheme == "a" || morpheme == "b" ) {
				if ( anomer_assigned ) {
					utility_exit_with_message( "Saccharide sequence input error: "
						"the anomer notation cannot be specified twice!" );
				} else if ( L_or_D_assigned ) {
					utility_exit_with_message( "Saccharide sequence input error: "
						"alpha/beta notation must precede L- or D- prefixes." );
				} else {
					// Set default linkage if missed.
					if ( ! linkage_assigned ) {
						residue_type_name << "->4)-";
						linkage_assigned = true;
					}
					if ( morpheme == "a" ) { morpheme = "alpha"; }
					if ( morpheme == "b" ) { morpheme = "beta"; }
					residue_type_name << morpheme << '-';
					anomer_assigned = true;
					morpheme = "";
				}

				// L/D indication
			} else if ( morpheme == "L" || morpheme == "D" || morpheme == "l" || morpheme == "d" ) {
				if ( L_or_D_assigned ) {
					utility_exit_with_message( "Saccharide sequence input error: "
						"the L or D notation cannot be specified twice!" );
				} else {
					// Set other defaults if missed.
					if ( ! linkage_assigned ) {
						residue_type_name << "->4)-";
						linkage_assigned = true;
					}
					if ( ! anomer_assigned ) {
						residue_type_name << "alpha-";
						anomer_assigned = true;
					}
					if ( morpheme == "l" ) { morpheme = "L"; }
					if ( morpheme == "d" ) { morpheme = "D"; }
					residue_type_name << morpheme << '-';
					L_or_D_assigned = true;
					morpheme = "";
				}

				// 3-Letter code (must be found in map of allowed 3-letter codes) and suffix(es)
			} else if ( morpheme.length() >= 3 ) {
				string const code( morpheme.substr( 0, 3 ) );
				string suffix( morpheme.substr( 3 ) );

				if ( ! CarbohydrateInfoManager::is_valid_sugar_code( code ) ) {
					// TODO: Add code for handling deoxy sugars here.
					// I'll have to check if the 1st character is 'd' and, if not, check if the 1st character is a
					// numeral and the 2nd character is a 'd'.
					// I'll also have to change the scope of the modifications variable and append to it by any
					// modifications coming from the suffix below.
					// What a pain.
					utility_exit_with_message( "Saccharide sequence input error: "
						"Unrecognized sugar 3-letter code: " + code );
				}

				char const sterochem( CarbohydrateInfoManager::default_stereochem_from_code( code ) );

				// Set defaults if missed.
				if ( ! linkage_assigned ) {
					residue_type_name << "->4)-";
				}
				if ( ! anomer_assigned ) {
					residue_type_name << "alpha-";
				}

				if ( ( L_or_D_assigned ) && ( sterochem == '*' ) ) {
					utility_exit_with_message( "Saccharide sequence input error: " +
						code + "cannot have a stereochemical prefix." );
				}
				if ( ( ! L_or_D_assigned ) && ( sterochem != '*' ) ) {
					residue_type_name << sterochem << '-';
				}

				// Assign 3-letter code.
				residue_type_name << code;

				// Check for suffixes.
				if ( suffix != "" ) {
					// Check for 1-letter ring-size-designating suffix.
					char const first_char( suffix[ 0 ] );
					if ( CarbohydrateInfoManager::is_valid_ring_affix( first_char ) ) {
						residue_type_name << first_char;
						suffix = suffix.substr( 1 );  // Strip ring-size designation from suffix.
					} else {
						// If it's not a valid ring affix, we assume that a linear saccharide is intended.
						TR.Debug << "Ring size not specified for residue " << code << ": assuming linear.";
					}

					// Now check for valid sugar modifications.
					if ( suffix != "" ) {
						vector1< pair< uint, string > > const modifications(
							io::carbohydrates::sugar_modifications_from_suffix( suffix ) );
						Size const n_modifications( modifications.size() );
						for ( uint i( 1 ); i <= n_modifications; ++i ) {
							uint position( modifications[ i ].first );
							string const modification( modifications[ i ].second );

							if ( CarbohydrateInfoManager::is_valid_modification_affix( modification ) ) {
								if ( position == 0 ) {  // This signals that a default value should be assigned.
									position = CarbohydrateInfoManager::default_position_from_affix( modification );
								}
								if ( position == 0 ) {  // If there's no default, we don't know what to do here.
									utility_exit_with_message( "Saccharide sequence input error: "
										"The position of the requested sugar modification " + modification +
										" is unclear." );
								}
							} else {
								utility_exit_with_message( "Saccharide sequence input error: "
									"Unrecognized modified sugar indicated by suffix(es)." );
							}
							residue_type_name << ':';  // Patches begin with a colon.
							string const patch_name_without_position(
								CarbohydrateInfoManager::patch_name_from_affix( modification ) );
							// so that I can avoid hardcoding exceptions and string look-ups here.  ~Labonte
							if ( ! CarbohydrateInfoManager::affix_has_inherent_position( modification ) ) {
								residue_type_name << position << '-';
							}
							residue_type_name << patch_name_without_position;
						}
					}
				}

				// Select a matching ResidueType and add to list (or exit without a match).
				ResidueTypeCOP residue( residue_set.name_mapOP( residue_type_name.str() ) );
				if ( residue == nullptr ) {
					utility_exit_with_message("Unrecognized residue " + residue_type_name.str() + " when building saccharide sequence " + sequence );
				}
				Size const n_branches( branch_points.size() );
				for ( uint i( 1 ); i <= n_branches; ++i ) {
					VariantType const variant(
						CarbohydrateInfoManager::branch_variant_type_from_position( branch_points[ i ] ) );
					residue = residue_set.get_residue_type_with_variant_added( *residue, variant ).get_self_ptr();
				}
				main_residue_types.push_back( residue );

				// Reset variables.
				morpheme = "";
				residue_type_name.str( "" );
				branch_points.resize( 0 );
				linkage_assigned = false;
				anomer_assigned = false;
				L_or_D_assigned = false;

			} else {  // Unrecognized morpheme
				utility_exit_with_message( "Saccharide sequence input error: "
					"Unrecognized sugar and/or notation in sequence.  Did you forget to indicate ring size?" );
			}
		}
	}  // next chr_num

	// Set the termini.
	// (Remember, the residues are in reverse order.)
	main_residue_types.front() = residue_set.get_residue_type_with_variant_added(
		*main_residue_types.front(), UPPER_TERMINUS_VARIANT ).get_self_ptr();

	branch_residue_types.append( main_residue_types );
	return branch_residue_types;
}  // residue_types_from_saccharide_sequence()


// TODO: Move to core::pose::carbohydrates?
// Append an empty or current Pose with saccharide residues, building branches as necessary.
/// @details  This function was written primarily as a subroutine for code shared by
/// make_pose_from_saccharide_sequence() and pose::carbohydrates:glycosylate_pose().  You probably do not want to call
/// it directly.
/// @param    <pose>: The Pose must be either empty, in which case a new oligo- or polysaccharide pose will be created,
/// or it must have a BRANCH_POINT variant residue to extend.  This Residue of the Pose must have unsatisfied
/// ResidueConnections.  By default, it is assumed to be the last residue of the Pose.
/// @param    <residue_types>: A list of ResidueTypes from which to build new residues for the Pose.  These must be
/// ordered such that a branch is completely finished before a new one is begun.  Branches are constructed by assuming
/// that branch connections are "satisfied" in the order in which they were created.  ResidueTypes must be of the
/// correct VariantType to be appended properly.  (See pose::carbohydrates:glycosylate_pose().)
/// @param    <resnum_to_be_appended>: The sequence position to which the glycan will be appended.  It must be a
/// BRANCH_POINT with an unsatisfied connection.  By default, it is assumed to be the last residue of the Pose.
void
append_pose_with_glycan_residues(
	pose::Pose & pose,
	chemical::ResidueTypeCOPs residue_types,
	core::uint resnum_to_be_appended /* = 0 */ )
{
	using namespace std;
	using namespace utility;
	using namespace chemical::carbohydrates;

	// Keep track of the branching as we go.
	list< pair< uint, string > > branch_points;

	if ( pose.empty() ) {
		TR.Debug << "Creating new oligosaccharide from provided ResidueTypes..." << endl;
	} else {
		TR.Debug << "Appending Pose with provided ResidueTypes..." << endl;

		// Check the anchor residue of the input Pose to see if it has unsatisfied branch points.
		if ( ! resnum_to_be_appended ) { resnum_to_be_appended = pose.size(); }
		Residue const & residue_to_be_appended( pose.residue( resnum_to_be_appended ) );

		if ( residue_to_be_appended.is_branch_point() ) {
			vector1< string > const branch_atom_names( residue_to_be_appended.type().branch_connect_atom_names() );
			Size const n_branches( branch_atom_names.size() );
			for ( uint i( 1 ); i <= n_branches; ++i ) {
				branch_points.emplace_back( resnum_to_be_appended, branch_atom_names[ i ] );
			}
		} else {
			TR.Warning << "append_pose_with_glycan_residues( " <<
				"pose::Pose & pose, chemical::ResidueTypeCOPs residue_types ): " <<
				"The provided anchor residue of <pose> has no branch point at which to append the glycan chain." << endl;
		}
	}

	// Append the residues.
	uint const n_types( residue_types.size() );
	for ( uint i( 1 ); i <= n_types; ++i ) {
		chemical::ResidueType const & rsd_type( *residue_types[ i ] );
		ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );

		TR.Debug << "Connecting current residue " << new_rsd->name3();
		// The residue_types are in a correct order here, so the lower of residue N is connected to the upper of residue N-1
		// If the upper of N-1 doesn't exist, we're a branch connected elsewhere
		if ( branch_points.size() && ( i == 1 || residue_types[ i - 1 ]->is_upper_terminus() ) ) {
			uint const anchor_residue( branch_points.front().first );  // First branch made is first branch connected.
			string const anchor_atom( branch_points.front().second );
			string const upper_atom( new_rsd->carbohydrate_info()->anomeric_carbon_name() );
			TR.Debug << " as branch from " << anchor_atom << " on residue " << anchor_residue <<
				" to " << upper_atom << "..." << endl;
			pose.append_residue_by_atoms( *new_rsd, true, upper_atom, anchor_residue, anchor_atom, true );
			branch_points.pop_front();
		} else if ( pose.size() != 0 && pose.residue_type( pose.size() ).upper_connect_id() != 0 ) {
			TR.Debug << " by extending the current main chain or branch..." << endl;
			pose.append_residue_by_bond( *new_rsd, true );
		} else {
			TR << " by appending by jump..." << endl;
			pose.append_residue_by_jump( *new_rsd, pose.size() );
		}
		if ( new_rsd->is_branch_point() ) {
			vector1< string > const branch_atom_names( new_rsd->type().branch_connect_atom_names() );
			Size const n_branches( branch_atom_names.size() );
			for ( uint j( 1 ); j <= n_branches; ++j ) {
				branch_points.emplace_back( pose.size(), branch_atom_names[ j ] );
			}
		}
	}
	debug_assert( branch_points.empty() );
}

void make_pose_from_sequence(
	pose::Pose & pose,
	chemical::ResidueTypeCOPs requested_types,
	bool const /* auto_termini=true */)
{
	// clear the pose
	pose.clear();

	// make the pose
	bool jump_to_next = false;
	for ( Size i = 1, ie = requested_types.size(); i <= ie; ++i ) {
		// grab the new residue
		if ( requested_types[ i ] == nullptr ) {
			// Gaps in the sequence can result in this, so we skip them.
			TR.Warning << "When making pose from sequence, position " << i << " can't be handled by Rosetta. Ignoring." << std::endl;
			continue;
		}
		chemical::ResidueType const & rsd_type = *requested_types[ i ];
		core::conformation::ResidueOP new_rsd( nullptr );
		new_rsd = conformation::ResidueFactory::create_residue( rsd_type );

		TR.Trace << "make_pose_from_sequence():  seqpos: " << i << " " << new_rsd->aa() << std::endl;

		// do the actual append
		// AMW: you may want to append a polymer UNK, like a ncaa, by a bond!
		if ( rsd_type.is_lower_terminus(  ) ||
				rsd_type.has_variant_type( chemical::N_ACETYLATION ) ||
				(new_rsd->aa() == chemical::aa_unk && !rsd_type.is_polymer()) ||
				new_rsd->aa() == chemical::aa_vrt ||
				new_rsd->aa() == chemical::aa_h2o ||
				jump_to_next ) {
			if ( ( new_rsd->aa() == chemical::aa_unk && !rsd_type.is_polymer() ) ||
					new_rsd->aa() == chemical::aa_vrt ) {
				//fpd TR.Warning << "found unknown aminoacid or X in sequence at position " << i <<  std::endl;
				//fpd if ( i< ie ) {
				//fpd  utility_exit_with_message(
				//fpd "found unknown aminoacid or X in sequence\n this leads to a seg-fault if we keep going...\n");
				//fpd }

				// if you don't think so ... make the code more stable and remove this
				// but only if this sequence doesn't seg-fault: KPAFGTNQEDYASYIXNGIIK" );

				//fpd ^^^ the problem is that the residue following the X should be connected by a jump as well.
				//     it should be of LOWER_TERMINUS variant type, but if not,
				//     we'll recover & spit out a warning for now.  same thing for ligands???
				jump_to_next = true;
			} else if ( jump_to_next ) {
				jump_to_next = false;
				if ( !rsd_type.is_lower_terminus(  ) ) {
					TR.Debug << "Residue! following X, Z, or an upper terminus is _not_ a lower terminus type!  " <<
						"Continuing ..." << std::endl;
				}
			}
			// each time this happens, a new chain should be started

			pose.append_residue_by_jump( *new_rsd, 1, "", "", true );
		} else {
			pose.append_residue_by_bond( *new_rsd, true );

			//fpd If res i is an upper terminus but (i+1) is not a lower terminus,
			//fpd the code exits on a failed assertion
			//fpd Don't let this happen; always jump in these cases
			if ( rsd_type.is_upper_terminus(  ) ) jump_to_next = true;
		}
	}

	TR.Debug << "sequence in pose: " << pose.sequence() << std::endl;
	TR.Debug << "annotated seq: " << pose.annotated_sequence() << std::endl;

} // core::pose::make_pose_from_sequence

////////////////////////////////////////////////////////////////////////////////
/// @details Given a Pose, a protein sequence where each character represents an
/// amino acid, and a ResidueTypeSet, give the Pose a conformation of covalently
/// linked residues that match the sequence. NOTE: support making pose from a
/// fully annotated sequence now, that is, for each residue variant or ligand
/// which cannot be deduced from one letter code directly, a [] is added
/// directly following the one letter code containing the residue's fullname, e.g.
/// K[lys:NtermProteinFull]ADFGCH[HIS_D]QNVE[glu:CtermProteinFull]Z[ZN].
/// This allows a pose to be constructed with full features from a silent output
/// file, such as with distinguished HIS tautomers, various chain termini and
/// cutpoint variants etc. Currently not working with disulfide variant CYD, but
/// this is on to-do list.
void make_pose_from_sequence(
	pose::Pose & pose,
	std::string const & sequence_in,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini /* true */
)
{
	// grab residue types
	chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( sequence_in, residue_set, auto_termini );
	debug_assert( core::pose::annotated_to_oneletter_sequence( sequence_in ).length() == requested_types.size() );

	make_pose_from_sequence(
		pose, requested_types, auto_termini);
}

////////////////////////////////////////////////////////////////////////////////
/// @details Given a Pose, a protein sequence where each character represents an
/// amino acid, and a ResidueTypeSet, give the Pose a conformation of covalently
/// linked residues that match the sequence. NOTE: support making pose from a
/// fully annotated sequence now, that is, for each residue variant or ligand
/// which cannot be deduced from one letter code directly, a [] is added
/// directly following the one letter code containing the residue's fullname, e.g.
/// K[lys:NtermProteinFull]ADFGCH[HIS_D]QNVE[glu:CtermProteinFull]Z[ZN].
/// This allows a pose to be constructed with full features from a silent output
/// file, such as with distinguished HIS tautomers, various chain termini and
/// cutpoint variants etc. Currently not working with disulfide variant CYD, but
/// this is on to-do list.
void make_pose_from_sequence(
	pose::Pose & pose,
	std::string const & sequence_in,
	chemical::ResidueTypeSetCOP residue_set,
	bool const auto_termini /* true */
)
{
	// grab residue types
	chemical::ResidueTypeCOPs requested_types = core::pose::residue_types_from_sequence( sequence_in, *residue_set, auto_termini );
	debug_assert( core::pose::annotated_to_oneletter_sequence( sequence_in ).length() == requested_types.size() );

	make_pose_from_sequence(
		pose, requested_types, auto_termini);
}

////////////////////////////////////////////////////////////////////////////////
/// @details overloaded version of make_pose_from_sequence, does the same
/// function, but reads in a string of the residue type set instead of a
/// ResidueTypeSet object.
void make_pose_from_sequence(
	pose::Pose & pose,
	std::string const & sequence_in,
	std::string const & type_set_name,
	//chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini /* true */
) {
	chemical::ResidueTypeSetCOP residue_set( chemical::ChemicalManager::get_instance()->residue_type_set( type_set_name ) );
	core::pose::make_pose_from_sequence( pose, sequence_in, *residue_set, auto_termini );
}


// Creates a Pose from an annotated, linear, IUPAC polysaccharide sequence <sequence> with ResidueTypeSet <residue_set>
// and store it in <pose>.
/// @param[in] <pose>: the Pose to fill
/// @param[in] <sequence>: an annotated IUPAC polysaccharide sequence,
/// e.g., "alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp"
/// @param[in] <residue_set>: the desired residue set
/// @param[in] <auto_termini>: if true (default) creates termini variants of terminal residues for the main chain
/// @details   Format for <sequence>:\n
/// Prefixes apply to the residue to which they are attached, below indicated by residue n.\n
/// Residues are listed from N to 1, where N is the total number of residues in the saccharide.\n
/// The sequence is parsed by reading to the next hyphen, so hyphens are crucial.\n
/// Linkage indication: "(a->x)-" specifies the linkage of residue n, where a is the anomeric carbon number of residue
/// (n+1) and x is the oxygen number of residue n.  The first residue listed in the annotated sequence (residue N)
/// need not have the linkage prefix.  A ->4) ResidueType will automatically be assigned by default if not specified.\n
/// Anomer indication: The strings "alpha-" or "beta-" are supplied next, which determines the stereochemistry of the
/// anomeric carbon of the residue to which it is prefixed.  An alpha ResidueType will automatically be assigned by
/// default.\n
/// Stereochemical indication: "L-" or "D-" specifies whether residue n is an L- or D-sugar.  The default is "D-".\n
/// 3-Letter code: A three letter code (in sentence case) MUST be supplied next.  This specifies the "base sugar name",
/// e.g., Glc is for glucose.  (A list of all recognized 3-letter codes for sugars can be found in
/// database/chemical/carbohydrates/codes_to_roots.map.)\n
/// 1-Letter suffix: If no suffix follows, residue n will be linear.  If a letter is present, it indicates the ring
/// size, where "f" is furanose, "p" is puranose, and "s" is septanose.\n
/// Branches are indicated using nested brackets and are best explained by example:\n
/// beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc is:\n
/// beta-D-Galp-(1->4)-D-GlcpNAc\n
///                       |\n
///     alpha-L-Fucp-(1->3)
/// @remarks   The order of residues in the created pose is in the opposite direction as the annotated sequence of
/// saccharide residues, as sugars are named with the 1st residue as the "main chain", with all other residues named as
/// substituents and written as prefixes.  In other words, sugars are usually drawn and named with residue 1 to the
/// right.
/// At present time, param files only exist for a limited number of sugars! ~ Labonte
void
make_pose_from_saccharide_sequence( pose::Pose & pose,
	std::string const & sequence,
	chemical::ResidueTypeSet const & residue_set,
	bool const auto_termini, /*true*/
	bool const idealize_linkages /*true*/ )
{
	using namespace std;
	using namespace utility;
	using namespace chemical;
	using namespace conformation;
	using namespace pose;
	using namespace pose::carbohydrates;
	using namespace conformation::carbohydrates;

	// Get list of carbohydrate ResidueTypes from which to construct the Pose.
	ResidueTypeCOPs residue_types( residue_types_from_saccharide_sequence( sequence, residue_set ) );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Residues to connect: ";
		for ( auto residue_type : residue_types ) {
			TR.Debug << residue_type->name() << ' ' << endl;
		}
	}

	// Fix the termini.
	if ( auto_termini ) {
		residue_types.front() = residue_set.get_residue_type_with_variant_added(
			*residue_types.front(), LOWER_TERMINUS_VARIANT ).get_self_ptr();
	}
	if ( ! auto_termini ) {
		residue_types.front() = residue_set.get_residue_type_with_variant_removed(
			*residue_types.front(), UPPER_TERMINUS_VARIANT ).get_self_ptr();
	}

	// Now we can build the Pose.
	pose.clear();
	append_pose_with_glycan_residues( pose, residue_types );

	// Let the Conformation know that it contains sugars.
	pose.conformation().contains_carbohydrate_residues( true );

	// Finally, set the PDB information.
	PDBInfoOP info( new PDBInfo( pose ) );
	info->name( pose.chain_sequence( 1 ) );  // Use the main-chain sequence as the default name.
	pose.pdb_info( info );

	//Create the GlycanTreeSetObserver if we don't already have it in the pose.
	if ( ! pose.glycan_tree_set() ) {
		pose.conformation().setup_glycan_trees();
	}

	if ( idealize_linkages ) {
		TR.Debug << "Idealizing glycosidic torsions." << endl;
		Size const n_glycans_added( residue_types.size() );
		pose::carbohydrates::idealize_last_n_glycans_in_pose( pose, n_glycans_added );
	}

	TR.Debug << "Created carbohydrate pose with main-chain sequence: " << pose.chain_sequence( 1 ) << endl;
}

// Creates a Pose from an annotated, linear, IUPAC polysaccharide sequence <sequence> with residue type set name
// <type_set_name> and store it in <pose>.
/// @details Overloaded version of make_pose_from_saccharide_sequence() that takes the string name for a residue type
/// set instead of a ResidueTypeSet object.  A convenience method for PyRosetta.
void
make_pose_from_saccharide_sequence( pose::Pose & pose,
	std::string const & sequence,
	std::string const & type_set_name, /*"fa_standard"*/
	bool const auto_termini, /*true*/
	bool const idealize_linkages /*true*/ )
{
	using namespace chemical;

	ResidueTypeSetCOP type_set( ChemicalManager::get_instance()->residue_type_set( type_set_name ) );
	make_pose_from_saccharide_sequence( pose, sequence, *type_set, auto_termini, idealize_linkages );
}

// Return a Pose from an annotated, linear, IUPAC polysaccharide sequence <sequence> with residue type set name
// <type_set_name>.
/// @details A convenience method for PyRosetta.
pose::PoseOP
pose_from_saccharide_sequence( std::string const & sequence,
	std::string const & type_set_name /*"fa_standard"*/,
	bool const auto_termini, /*true*/
	bool const idealize_linkages /*true*/ )
{
	using namespace pose;

	PoseOP pose( new Pose() );
	make_pose_from_saccharide_sequence( *pose, sequence, type_set_name, auto_termini, idealize_linkages );
	return pose;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string annotated_to_oneletter_sequence(
	std::string const & annotated_seq
) {
	bool add( true );
	std::string oneletter_seq;
	for ( char i : annotated_seq ) {
		if ( i == '/' ) continue;
		if ( i == '[' ) add = false;
		if ( add ) oneletter_seq += i;
		if ( i == ']' ) add = true;
	}

	return oneletter_seq;
}

/// @details ResidueTypeFinder finds simplest residue type with this AA & requested termini.
chemical::ResidueTypeCOP
get_rsd_type_from_aa( chemical::ResidueTypeSet const & residue_set,
	chemical::AA const & my_aa, bool const & is_lower_terminus, bool const & is_upper_terminus )
{
	using namespace core::chemical;
	ResidueTypeCOP rsd_type = ResidueTypeFinder( residue_set ).aa( my_aa ).get_representative_type();
	utility::vector1< VariantType > variants;
	if ( rsd_type->is_polymer() ) {
		if ( is_lower_terminus ) variants.push_back( LOWER_TERMINUS_VARIANT );
		if ( is_upper_terminus ) variants.push_back( UPPER_TERMINUS_VARIANT );
	}
	return ResidueTypeFinder( residue_set ).aa( my_aa ).variants( variants ).get_representative_type();
}

} // namespace core
} // namespace pose
