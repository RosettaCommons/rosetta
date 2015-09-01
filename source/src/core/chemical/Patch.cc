// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @file Patch.cc
///
/// @brief Implementation class for abstract class Residue
///
/// @details
/// My understanding of this piece of code comes from various conversations with Oliver Lang, Summer Thyme, and Andrew
/// Leaver-Fay. If the ideas are incorrect, feel free to correct.
///
/// General Overview
///
/// What are patches? Patches are modifications to the original amino acid residue types.  These modifications include
/// capping the N and C terminus, adding phosphorylation to residues, etc, etc.  The actual files that are read and used
/// are found in: database/chemical/residue_type_sets/fa_standard/patches
///
/// All of these patches are read when ResidueTypeSet is instantiated and are applied to all residues that are
/// identified as needing to be patched.  (Identification of residues is explained later.)  Identification of residues
/// that need to be patched varies by what is in the actual patch file.  For example, peptide residues that are C and/or
/// N termini are patched with the N and/or C terminus patches, which are defined in the files NtermProteinFull.txt and
/// CtermProteinFull.txt.
///
/// Overview of a Patch File
///
/// Let's look at a patch file to see what goes on.  The file we are looking at is NtermProteinFull.txt.  The file is
/// shortened to make it easier to read.  Important parts are shown:
///
/// NAME NtermProteinFull  # actual name of the patch
/// TYPES LOWER_TERMINUS   # the type that this patch belongs to.  Types are defined in VariantType.cc/.hh
///
/// # This section is for general selection rules to apply the patch.
/// BEGIN_SELECTOR  # Here is where we define how to select amino acids for the patch.
///                 # Properties/variant types needed for the patch are found between the BEGIN_SELECTOR and
///                 # END_SELECTOR "titles"
/// PROPERTY PROTEIN  # For this patch to apply, the residue must be a protein, as defined in the paramaters file.
/// NOT VARIANT_TYPE LOWER_TERMINUS  # We do not want to patch the variant type LOWER_TERMINUS.  This is because, we do
///                                  # not want to double patch a residue (if the residue is already variant type
///                                  # lower_terminus).
/// END_SELECTOR  # ending the selection process for the patch
///
/// # This section is to modify specific residues that the patch encounters.
/// BEGIN_CASE ## PROLINE  # Within the BEGIN_CASE and END_CASE section is where the residue is modified.
///                        # These are the operations that occur to apply the patch.
/// BEGIN_SELECTOR  # Once again, we are defining what to do in this specific BEGIN_CASE/END_CASE block by using the
///                 # selector.
/// AA PRO        # We only want to modify aa PRO in the following way between the BEGIN_CASE and END_CASE block.
/// END_SELECTOR  # End selection requirements
///
/// ADD_ATOM 1H   Hpol HC 0.24  # straight forward; add atom 1H
/// ADD_ATOM 2H   Hpol HC 0.24
/// ADD_BOND N 1H               # Atoms need bonds; add those bonds.
/// ADD_BOND N 2H
/// SET_POLYMER_CONNECT LOWER NONE # setting a property
///
/// ## totally making these up:
/// SET_ICOOR 1H 120 60 1 N CA C ## like to use CD but CD's parent is CG  # Once new atoms are placed,
/// SET_ICOOR 2H 120 60 1 N CA 1H                                         # new internal cooridantes need to be made.
///
/// ## modify properties of existing atoms
/// SET_ATOM_TYPE N Nlys
/// SET_MM_ATOM_TYPE N NP
/// SET_ATOMIC_CHARGE N -0.07
/// SET_ATOMIC_CHARGE CA 0.16
/// SET_ATOMIC_CHARGE HA 0.09
/// SET_ATOMIC_CHARGE CD 0.16
/// SET_ATOMIC_CHARGE 1HD 0.09
/// SET_ATOMIC_CHARGE 2HD 0.09
/// ADD_PROPERTY LOWER_TERMINUS ## implies terminus
/// END_CASE
///
/// BEGIN_CASE ### THE GENERAL CASE  # Here is the general case, specified at the very beginning of the patch file.
///
/// ## these are the operations involved
/// DELETE_ATOM H ## deletes all bonds to this atom
/// ADD_ATOM 1H   Hpol HC 0.33
/// ADD_ATOM 2H   Hpol HC 0.33
/// ADD_ATOM 3H   Hpol HC 0.33
/// ADD_BOND N 1H
/// ADD_BOND N 2H
/// ADD_BOND N 3H
/// SET_POLYMER_CONNECT LOWER NONE
///
/// ## totally making these up:
/// SET_ICOOR 1H 120 60 1 N CA C
/// SET_ICOOR 2H 120 60 1 N CA 1H
/// SET_ICOOR 3H 120 60 1 N CA 2H
///
/// ## modify properties of existing atoms
/// SET_ATOM_TYPE N Nlys
/// SET_MM_ATOM_TYPE N NH3
/// SET_ATOMIC_CHARGE N -0.3
/// SET_ATOMIC_CHARGE CA 0.21
/// SET_ATOMIC_CHARGE HA 0.10
/// ADD_PROPERTY LOWER_TERMINUS ## implies terminus
/// END_CASE
///
///
/// So, let's go through what was put up there.  In general, you name your patch and assign it a type.  Then, you tell
/// the patch system in the BEGIN_SELECTOR / END_SELECTOR block what type of properties that you are looking for to
/// apply the patch.  You then, in the BEGIN_CASE / END_CASE block, add operations to the residue.  You can specify
/// specific selectors within the BEGIN_CASE / END_CASE by having another BEGIN_SELECTOR / END_SELECTOR block.  People
/// generally specify specific amino acids within the BEGIN_SELECTOR / END_SELECTOR blocks within the BEGIN_CASE /
/// END_CASE blocks.  At the end, there is one last BEGIN_CASE / END_CASE with no BEGIN_SELECTOR / END_SELECTOR block.
/// This block is used to specify what to do with the general case, which was defined at the very beginning of the file.
/// A few more details follow below.
///
/// First, you must name your patch.  Then you have to give it a "type".  This type is defined in VariantType.cc/.hh.
/// All the type does is adds a name that can be used later on in different protocols.  It also helps the patch system
/// keep track of what residues are patched with what type.  The type has no magical meaning; it's a name that you give
/// it that is "defined" in the variantType.cc.  In short, it's a name handler.  Once the name and type have been
/// assigned, you must define a general case in which your patch applies.  Generally, this means that you want it to
/// apply to a specific type of aa, or all residues of a specific type or property.  After this block comes the
/// BEGIN_CASE / END_CASE statements.  In these statements, you specify what modifications you want to apply to the
/// residues.  If you have different requirements for different amino acids, you must specify this in a BEGIN_SELECTOR /
/// END_SELECTOR block within the BEGIN_CASE / END_CASE block.  You must specify specific amino acids first, then in a
/// new BEGIN_CASE / END_CASE block, you specify the general modifications defined by the first selector.  Confusing,
/// huh?  The reason for having specific selectors defined before the general case is because..., I don't know.  I
/// suspect it's because of how the file is read in and applied.  Regardless, if you want to do something to a specific
/// amino acid/residue before you get to the general case, you need to have a BEGIN_CASE / END_CASE block with a
/// selector in there before you get to the general BEGIN_CASE / END_CASE block.
///
/// Using Patch Selector to Specify Application of a Patch to Specific Residues
///
/// Woah there!  Now that you know how to use patches, you don't need to go all gun-ho crazy and create a ton of
/// patches. This is because all patches are read when ResidueTypeSet is created and applied.  This creates overhead if
/// you have a ton of patch files.  To circumvent this, someone wrote an option called -chemical:patch_selectors.  I
/// don't know who wrote it, or where it's located, sorry.  This option allows you to create patches that are only
/// loaded when you use the option -chemical:patch_selectors.  To get this to work, you create a patch like normal.
/// Then, in the BEGIN_SELECTOR / END_SELECTOR line, you add a line CMDLINE_SELECTOR <name of command line>.  This
/// means, if you use the CMDLINE_SELECTOR sc_orbitals, your patch will only be loaded if you use the option
/// -chemical:patch_selectors sc_orbitals.  Cool, huh?
///
///
/// Trouble Shooting!!!
/// So, you added a new patch and used the command line selector and nothing happens.  Did you modify VariantTypes.cc to
/// have your new type?  Did you name things correctly?  Did you edit the patches.txt to include your patch?  Are you
/// trying to use a command that does not exist?  Look at patchoperations.cc for commands that exist.  These are all
/// things to check for.
///
/// So, you have added a new patch, and you get segmentation faults.  Most common seg fault for me is in the stub atom.
/// This is because your patch might add atoms/delete atoms that are used in other patches.  In order to fix this, you
/// must modify other patches to include your modifications.  Remember patches are all loaded whenever ResidueTypeSet is
/// instantiated, which is usually near the beginning of a Rosetta protocol when a pose is 1st loaded!  If you have
/// conflicts with other patches in your patch, even if you are using the command line selector, you must resolve those
/// problems.  Write an integration test so that people won't screw up your patches.
///
///
/// @author
/// Phil Bradely
/// Steven Combs only added comments
/////////////////////////////////////////////////////////////////////////


// Unit headers
#include <core/chemical/Patch.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <fstream>


namespace core {
namespace chemical {

/// @details Auto-generated virtual destructor
Patch::~Patch() {}

/// @details Auto-generated virtual destructor
PatchCase::~PatchCase() {}

static thread_local basic::Tracer tr( "core.chemical" );

/// @brief the string used to generate new residue names
std::string const PATCH_LINKER( ":" );
// Changed from '_p:' to ':' for conciseness. Backward compatibility fixes throughout code.
// std::string const patch_linker( "_p:" );

std::string
residue_type_base_name( ResidueType const & rsd_type )
{
	return rsd_type.name().substr( 0, rsd_type.name().find( PATCH_LINKER ) );
}

std::string
residue_type_all_patches_name( ResidueType const & rsd_type )
{
	Size spos = rsd_type.name().find( PATCH_LINKER );
	if ( spos < rsd_type.name().length() ) return rsd_type.name().substr( spos );
	else return "";
}


/// @brief handy function, return the first word from a line
std::string
tag_from_line( std::string const & line )
{
	std::string tag;
	std::istringstream l( line );
	l >> tag;
	if ( l.fail() ) return "";
	else return tag;
}

/// @brief create a PatchCase from input lines
/// @details add selector_ from lines enclosed by "BEGIN_SELECTOR" and "END_SELECTOR".\n
/// add operations_ from each input line containing a single operation
PatchCaseOP
case_from_lines(
	utility::vector1< std::string > const & lines
)
{
	PatchCaseOP pcase( new PatchCase() );

	bool in_selector( false );
	for ( uint i=1; i<= lines.size(); ++i ) {
		std::string const tag( tag_from_line( lines[i] ) );

		if ( tag == "BEGIN_SELECTOR" ) {
			debug_assert( !in_selector );
			in_selector = true;
		} else if ( tag == "END_SELECTOR" ) {
			in_selector = false;
		} else if ( in_selector ) {
			pcase->selector().add_line( lines[i] );
		} else {
			PatchOperationOP operation( patch_operation_from_patch_file_line( lines[i] ) );
			if ( operation ) pcase->add_operation( operation );
		}
	}

	return pcase;
}


/// @details First clone the base ResidueType.  Then patching for this case is done by applying all the operations.
/// finalize() is called after the VariantTypes and name are set by Patch::apply().
/// @note    If you call this method without calling finalize(), your ResidueType may not have the correct derived data!
ResidueTypeOP
PatchCase::apply( ResidueType const & rsd_in, bool const instantiate /* = true */ ) const
{
	ResidueTypeOP rsd;
	if ( instantiate ) {
		rsd = rsd_in.clone();
	} else {
		rsd = rsd_in.placeholder_clone();
	}

	for ( utility::vector1< PatchOperationOP >::const_iterator iter = operations_.begin(),
			iter_end = operations_.end(); iter != iter_end; ++iter ) {
		if ( !instantiate && !( *iter )->applies_to_placeholder() ) { continue; }
		bool const fail( ( *iter )->apply( *rsd ) );
		if ( fail ) {
			utility_exit_with_message( "Failed to apply a PatchOperation to " + rsd->name() );
			return 0;
		}
	}

	//std::cout << "amw PatchCase::apply Checking on igroup for res " << rsd->name() << ": " << rsd->interchangeability_group() << std::endl;

	return rsd;
}


/// @details Go through patch operations in this PatchCase, and compile list of any atom names that are added.
utility::vector1< std::string >
PatchCase::adds_atoms() const
{
	utility::vector1< std::string > atom_names;
	for ( utility::vector1< PatchOperationOP >::const_iterator iter = operations_.begin(),
			iter_end = operations_.end(); iter != iter_end; ++iter ) {
		std::string const atom_name = ( *iter )->adds_atom();
		if ( atom_name.size() > 0 ) atom_names.push_back( atom_name );
	}
	return atom_names;
}

/// @details Go through patch operations in this PatchCase, and compile list of any atom names that are deleted.
utility::vector1< std::string >
PatchCase::deletes_atoms() const
{
	utility::vector1< std::string > atom_names;
	for ( utility::vector1< PatchOperationOP >::const_iterator iter = operations_.begin(),
			iter_end = operations_.end(); iter != iter_end; ++iter ) {
		std::string const atom_name = ( *iter )->deletes_atom();
		if ( atom_name.size() > 0 ) atom_names.push_back( atom_name );
	}
	return atom_names;
}

/// @details Go through patch operations in this PatchCase, and compile list of any property names that are added.
utility::vector1< std::string >
PatchCase::adds_properties() const
{
	utility::vector1< std::string > property_names;
	for ( utility::vector1< PatchOperationOP >::const_iterator iter = operations_.begin(),
			iter_end = operations_.end(); iter != iter_end; ++iter ) {
		std::string const property_name = ( *iter )->adds_property();
		if ( property_name.size() > 0 ) property_names.push_back( property_name );
	}
	return property_names;
}

/// @details Go through patch operations in this PatchCase, and compile list of any property names that are deleted.
utility::vector1< std::string >
PatchCase::deletes_properties() const
{
	utility::vector1< std::string > property_names;
	for ( utility::vector1< PatchOperationOP >::const_iterator iter = operations_.begin(),
			iter_end = operations_.end(); iter != iter_end; ++iter ) {
		std::string const property_name = ( *iter )->deletes_property();
		if ( property_name.size() > 0 ) property_names.push_back( property_name );
	}
	return property_names;
}

/// @details Go through patch operations in this PatchCase, and compile list of any property names that are deleted.
utility::vector1< std::string >
PatchCase::deletes_variants() const
{
	utility::vector1< std::string > variant_names;
	for ( utility::vector1< PatchOperationOP >::const_iterator iter = operations_.begin(),
			iter_end = operations_.end(); iter != iter_end; ++iter ) {
		std::string const variant_name = ( *iter )->deletes_variant();
		if ( variant_name.size() > 0 ) variant_names.push_back( variant_name );
	}
	return variant_names;
}


/// @details - first read in all lines from the file, discarding # comment lines
/// - parse input lines for Patch name and variant types (NAME, TYPES)
/// - parse input lines for general ResidueTypeSelector defined for this Patch (BEGIN_SELECTOR, END_SELECTOR)
/// - parse input lines to create each case accordingly (BEGIN_CASE, END_CASE)
/// @note keep the order to avoid triggering parsing errors
void
Patch::read_file( std::string const & filename )
{
	// clear old data
	tr.Debug << "Reading patch file: " << filename << std::endl;

	name_ = "";
	types_.clear();
	selector_.clear();
	cases_.clear();
	replaces_residue_type_ = false;

	utility::vector1< std::string > lines;
	{ // read the lines file
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message("Cannot find patch file: "+filename);
		}
		std::string line;
		while ( getline( data,line ) ) {
			std::string const tag( tag_from_line( line ) );
			if ( tag.size() && tag[0] != '#' ) lines.push_back( line );
		}
	}

	// misc parsing
	for ( uint i=1; i<= lines.size(); ++i ) {
		std::istringstream l(lines[i]);
		std::string tag;
		l >> tag;
		if ( tag == "NAME" ) {
			l >> name_;
		} else if ( tag == "TYPES" ) {
			std::string t;
			l >> t;
			while ( !l.fail() ) {
				types_.push_back( t );
				l >> t;
			}
		} else if ( tag == "REPLACE_RES_TYPE" ) {
			replaces_residue_type_ = true;
		}
	}

	// build the residue selector
	{
		bool in_selector( false );
		for ( uint i=1; i<= lines.size(); ++i ) {
			std::string tag( tag_from_line( lines[i] ) );
			if ( tag == "BEGIN_CASE" ) {
				break;
			} else if ( tag == "BEGIN_SELECTOR" ) {
				debug_assert( !in_selector );
				in_selector = true;
			} else if ( tag == "END_SELECTOR" ) {
				debug_assert( in_selector );
				in_selector = false;
			} else if ( in_selector ) {
				selector_.add_line( lines[i] );
			}
		}
	}

	// get the cases
	utility::vector1< std::string > case_lines;
	bool in_case( false );
	while ( !lines.empty() ) {
		// look for a case
		std::string tag( tag_from_line( lines[1] ) );
		if ( tag == "BEGIN_CASE" ) {
			debug_assert( case_lines.empty() );
			debug_assert( !in_case );
			in_case = true;
		} else if ( tag == "END_CASE" ) {
			PatchCaseOP new_case( case_from_lines( case_lines ) );
			if ( new_case ) cases_.push_back( new_case );
			case_lines.clear();
			in_case = false;
		} else if ( in_case ) case_lines.push_back( lines[1] );

		lines.erase( lines.begin() );
	}
}

/// @details loop through the cases in this patch and if it is applicable to this ResidueType, the corresponding patch
/// operations are applied to create a new variant type of the basic ResidueType.  The new types's name and its
/// variant type info are updated together with all other primary and derived ResidueType data.
/// Finally, call finalize() to update all primary and derived data for the new ResidueType.
ResidueTypeOP
Patch::apply( ResidueType const & rsd_type, bool const instantiate /* = true */ ) const
{
	if ( !applies_to( rsd_type ) ) { return 0; }  // I don't know how to patch this residue.

	using namespace basic;

	for ( utility::vector1< PatchCaseOP >::const_iterator iter= cases_.begin(),
			iter_end = cases_.end(); iter != iter_end; ++iter ) {

		if ( ( *iter )->applies_to( rsd_type ) ) {
			// this patch case applies to this rsd_type
			ResidueTypeOP patched_rsd_type( ( *iter )->apply( rsd_type, instantiate ) );
			if ( patched_rsd_type ) {
				// patch succeeded!
				if ( !replaces_residue_type_ ) { // This is bananas. Shouldn't just forget that patch was applied. -- rhiju.
					for ( utility::vector1< std::string >::const_iterator iter=types_.begin(),
							iter_end = types_.end(); iter != iter_end; ++iter ) {
						patched_rsd_type->add_variant_type( *iter );
					}
					// AMW: Special case for the D patch. In ONLY THIS CASE,
					// application of the D patch prepends the letter D. No :
					std::string name_new;
					if ( name_ == "D" ) {
						name_new = "D" + patched_rsd_type->name();
					} else {
						name_new = patched_rsd_type->name() + PATCH_LINKER + name_;
					}
					patched_rsd_type->name( name_new );
					//std::cout << "amw Patch::apply Checking on igroup for res " << patched_rsd_type->name() << ": " << patched_rsd_type->interchangeability_group() << std::endl;
				}

				//std::cout << "amw Patch::apply after renaming Checking on igroup for res " << patched_rsd_type->name() << ": " << patched_rsd_type->interchangeability_group() << std::endl;

				if ( instantiate ) {
					patched_rsd_type->finalize();
					tr.Debug << "successfully patched: " << rsd_type.name() <<
						" to: " << patched_rsd_type->name() << std::endl;
				}

				//std::cout << "amw Patch::apply after finalize Checking on igroup for res " << patched_rsd_type->name() << ": " << patched_rsd_type->interchangeability_group() << std::endl;

				return patched_rsd_type;
			}
		}
	}

	// I am commenting out this exit below, because I have patches with complicated selection criteria.  It can happen
	// that a patch's main selector applies to a ResidueType yet none of the sub-cases apply.  This should not fail;
	// the ResidueType should just be skipped.  If this point is reached, a null pointer is returned, which will be
	// properly handled by ResidueTypeSet, which checks that the new ResidueType exists before adding it to the Set.
	// If there is a genuine error in patching, the PatchOperations (or ResidueType) will error during the actual
	// patching, so I believe this is completely safe. ~Labonte
	//utility_exit_with_message( "no patch applied? " + name_ + " to " + rsd_type.name() );
	return 0;
}

/// @details loop through the cases in this patch and if it is applicable to this ResidueType, compile
/// a list of any atom_names that are added.
utility::vector1< std::string >
Patch::adds_atoms( ResidueType const & rsd_type ) const
{
	utility::vector1< std::string > atom_names;
	if ( !applies_to( rsd_type ) ) return atom_names;  // I don't know how to patch this residue.

	for ( utility::vector1< PatchCaseOP >::const_iterator iter= cases_.begin(),
			iter_end = cases_.end(); iter != iter_end; ++iter ) {

		if ( (*iter)->applies_to( rsd_type ) ) {
			// this patch case applies to this rsd_type
			utility::vector1< std::string > atom_names_for_patch_case = ( *iter )->adds_atoms();
			for ( Size n = 1; n <= atom_names_for_patch_case.size(); n++ ) {
				atom_names.push_back(  atom_names_for_patch_case[ n ] );
			}
		}
	}

	return atom_names;
}

/// @details loop through the cases in this patch and if it is applicable to this ResidueType, compile
/// a list of any atom_names that are deleted.
utility::vector1< std::string >
Patch::deletes_atoms( ResidueType const & rsd_type ) const
{
	utility::vector1< std::string > atom_names;
	if ( !applies_to( rsd_type ) ) return atom_names;  // I don't know how to patch this residue.

	for ( utility::vector1< PatchCaseOP >::const_iterator iter= cases_.begin(),
			iter_end = cases_.end(); iter != iter_end; ++iter ) {

		if ( (*iter)->applies_to( rsd_type ) ) {
			// this patch case applies to this rsd_type
			utility::vector1< std::string > atom_names_for_patch_case = ( *iter )->deletes_atoms();
			for ( Size n = 1; n <= atom_names_for_patch_case.size(); n++ ) {
				atom_names.push_back(  atom_names_for_patch_case[ n ] );
			}
		}
	}

	return atom_names;
}

/// @details loop through the cases in this patch and if it is applicable to this ResidueType, compile
/// a list of any properties that are added.
utility::vector1< std::string >
Patch::adds_properties( ResidueType const & rsd_type ) const
{
	utility::vector1< std::string > properties;
	if ( !applies_to( rsd_type ) ) return properties;  // I don't know how to patch this residue.

	for ( utility::vector1< PatchCaseOP >::const_iterator iter= cases_.begin(),
			iter_end = cases_.end(); iter != iter_end; ++iter ) {

		if ( (*iter)->applies_to( rsd_type ) ) {
			// this patch case applies to this rsd_type
			utility::vector1< std::string > properties_for_patch_case = ( *iter )->adds_properties();
			for ( Size n = 1; n <= properties_for_patch_case.size(); n++ ) {
				properties.push_back(  properties_for_patch_case[ n ] );
			}
		}
	}

	return properties;
}

/// @details loop through the cases in this patch and if it is applicable to this ResidueType, compile
/// a list of any properties that are deleted.
utility::vector1< std::string >
Patch::deletes_properties( ResidueType const & rsd_type ) const
{
	utility::vector1< std::string > properties;
	if ( !applies_to( rsd_type ) ) return properties;  // I don't know how to patch this residue.

	for ( utility::vector1< PatchCaseOP >::const_iterator iter= cases_.begin(),
			iter_end = cases_.end(); iter != iter_end; ++iter ) {

		if ( (*iter)->applies_to( rsd_type ) ) {
			// this patch case applies to this rsd_type
			utility::vector1< std::string > properties_for_patch_case = ( *iter )->deletes_properties();
			for ( Size n = 1; n <= properties_for_patch_case.size(); n++ ) {
				properties.push_back(  properties_for_patch_case[ n ] );
			}
		}
	}

	return properties;
}

/// @details loop through the cases in this patch and if it is applicable to this ResidueType, compile
/// a list of any variants that are deleted.
utility::vector1< std::string >
Patch::deletes_variants( ResidueType const & rsd_type ) const
{
	utility::vector1< std::string > variants;
	if ( !applies_to( rsd_type ) ) return variants;  // I don't know how to patch this residue.

	for ( utility::vector1< PatchCaseOP >::const_iterator iter= cases_.begin(),
			iter_end = cases_.end(); iter != iter_end; ++iter ) {

		if ( (*iter)->applies_to( rsd_type ) ) {
			// this patch case applies to this rsd_type
			utility::vector1< std::string > variants_for_patch_case = ( *iter )->deletes_variants();
			for ( Size n = 1; n <= variants_for_patch_case.size(); n++ ) {
				variants.push_back(  variants_for_patch_case[ n ] );
			}
		}
	}

	return variants;
}


} // chemical
} // core
