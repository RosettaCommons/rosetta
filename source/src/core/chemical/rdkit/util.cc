// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/util.cc
/// @brief Utilities for interacting with the RDKit library.
/// @author Rocco Moretti (rmorettiase@gmail.com)


#include <core/chemical/rdkit/util.hh>

#include <core/chemical/AtomRefMapping.hh>
#include <core/types.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

#include <rdkit/RDGeneral/utils.h>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/ForceField/ForceField.h>
#include <rdkit/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <rdkit/GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <rdkit/GraphMol/ForceFieldHelpers/UFF/AtomTyper.h>
#include <rdkit/GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/FMCS/FMCS.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmartsWrite.h>
#include <rdkit/GraphMol/Descriptors/MolDescriptors.h>

namespace core {
namespace chemical {
namespace rdkit {

static basic::Tracer TR("core.chemical.rdkit.util");
// The following tracer is not meant for actual output, but is instead intended to be used
// to control the RDKit output levels via the Rosetta tracer controls.
static basic::Tracer TR_RDKit("RDKit");

/// @brief Initialize the RDKit random number generator.
/// @details Note that seed is an int to match the seed generated in core/init.cc
void initialize_rdkit_random( int seed ) {
	// Can't pass seed directly to getRandomGenerator(), as negative values are ignored.
	::RDKit::getRandomGenerator().seed(seed);
	TR.Debug << "Initializing RDKit random number generator with seed: " << seed << std::endl;
}

void initialize_rdkit_tracers() {

	::RDLog::InitLogs();

	// This uses internal implementation details of the RDLog mechanism - it may not work with other version of RDKit
	if ( ! TR_RDKit.visible(basic::t_debug) ) {
		rdDebugLog->df_enabled = false;
	} else {
		// Reset to std::cout
		rdDebugLog->dp_dest = &std::cout;
	}
	if ( ! TR_RDKit.visible(basic::t_info) ) {
		rdInfoLog->df_enabled = false;
	} else {
		// Goes to std::cout normally.
	}
	if ( ! TR_RDKit.visible(basic::t_warning) ) {
		rdWarningLog->df_enabled = false;
	} else {
		// Reset to std::cout
		rdDebugLog->dp_dest = &std::cout;
	}
	if ( ! TR_RDKit.visible(basic::t_error) ) {
		rdErrorLog->df_enabled = false;
	} else {
		// Reset to std::cout
		rdDebugLog->dp_dest = &std::cout;
	}

}

::RDKit::Bond::BondType
convert_to_rdkit_bondtype( core::chemical::BondName bondtype, bool aro2double) {
	switch ( bondtype ) {
	case SingleBond :
		return ::RDKit::Bond::SINGLE;
	case DoubleBond :
		return ::RDKit::Bond::DOUBLE;
	case TripleBond :
		return ::RDKit::Bond::TRIPLE;
	case AromaticBond :
		return ( aro2double )? ::RDKit::Bond::DOUBLE : ::RDKit::Bond::AROMATIC;
	default :
		TR.Warning << "Treating Rosetta bond type " << bondtype << " as UNSPECIFIED in RDKit. " << std::endl;
		return ::RDKit::Bond::UNSPECIFIED;
	}
}

core::chemical::BondName
convert_from_rdkit_bondtype( ::RDKit::Bond::BondType bondtype) {
	switch ( bondtype ) {
	case ::RDKit::Bond::SINGLE :
		return SingleBond;
	case ::RDKit::Bond::DOUBLE :
		return DoubleBond;
	case ::RDKit::Bond::TRIPLE :
		return TripleBond;
	case ::RDKit::Bond::AROMATIC :
		return AromaticBond;
	default :
		TR.Warning << "Treating RDKit bond type " << bondtype << " as Unknown." << std::endl;
		return UnknownBond;
	}
}

/// @brief Get the name of the RDMol
std::string
get_name(::RDKit::ROMol const & mol) {
	return mol.getProp< std::string>(::RDKit::common_properties::_Name);
}

bool
has_physical_Hs(::RDKit::ROMol const & mol) {
	for ( core::Size ii(0); ii < mol.getNumAtoms(); ++ii ) {
		if ( mol.getAtomWithIdx(ii)->getAtomicNum() == 1 ) {
			return true;
		}
	}
	return false;
}

bool
has_explicit_Hs(::RDKit::ROMol const & mol) {
	for ( core::Size ii(0); ii < mol.getNumAtoms(); ++ii ) {
		if ( mol.getAtomWithIdx(ii)->getNumExplicitHs() > 0 ) {
			return true;
		}
	}
	return false;
}

bool
has_implicit_Hs(::RDKit::ROMol const & mol) {
	for ( core::Size ii(0); ii < mol.getNumAtoms(); ++ii ) {
		if ( mol.getAtomWithIdx(ii)->getNumImplicitHs() > 0 ) {
			return true;
		}
	}
	return false;
}

::RDKit::ForceFieldOP
get_forcefield(::RDKit::ROMol & mol, int conf_num /*=-1*/) {
	// Try MMFF, if possible, else try UFF
	if ( ::RDKit::MMFF::MMFFMolProperties(mol).isValid() ) {
		TR << "Using MMFF to optimize molecule" << std::endl;
		// 5.0 => Only consider non-bonded within 5 Ang
		::RDKit::ForceFieldOP ff( ::RDKit::MMFF::constructForceField(mol, 5.0, conf_num, true ) );
		ff->initialize();
		return ff;
	} else if ( ::RDKit::UFF::getAtomTypes(mol).second ) {
		TR << "Using UFF to optimize molecule" << std::endl;
		// 1.5 => Only consider non-bonded within ~3.5*1.5  Ang
		::RDKit::ForceFieldOP ff( ::RDKit::UFF::constructForceField(mol, 1.5, conf_num, true ) );
		ff->initialize();
		return ff;
	} else {
		TR.Warning << "Cannot find a suitable forcefield for optimizing molecule: " << ::RDKit::MolToSmiles( mol ) << std::endl;
		return ::RDKit::ForceFieldOP(nullptr);
	}
}

void
softSanitize(::RDKit::RWMol & mol) {
	using namespace ::RDKit::MolOps;
	mol.clearComputedProps();
	cleanUp(mol);
	mol.updatePropertyCache(/*strict*/false);
	::RDKit::VECT_INT_VECT arings;
	symmetrizeSSSR(mol, arings);
	// Don't kekulize -- this is where most of the sanitization errors occur.
	// (I believe the others are kekulization dependent and don't seem to be needed for the current set of search strings.)
}

/// @brief Remove any excess hydrogens, where "excess" is defined as any which contribute to a positive formal charge
/// (Assumes a graph-hydrogen removed form.)
void
remove_excess_protons(::RDKit::RWMol & rdmol) {
	for ( core::Size ii(0); ii < rdmol.getNumAtoms(); ++ii ) {
		::RDKit::Atom & atom( *rdmol.getAtomWithIdx(ii) );
		int const & charge( atom.getFormalCharge() );
		if ( charge > 0 ) {
			unsigned int nhydro( atom.getNumExplicitHs() );
			if ( nhydro > 0 ) {
				unsigned int delta( std::min((unsigned int)(charge), nhydro) );
				atom.setNumExplicitHs( nhydro - delta );
				atom.setFormalCharge( charge - delta );
				atom.updatePropertyCache(false); // Need to reset implicit/explict counts.
			}
			// Now deal with any residual implicits.
			int const & charge2( atom.getFormalCharge() );
			nhydro = atom.getNumImplicitHs();
			if ( charge2 > 0 && nhydro > 0 ) {
				unsigned int delta( std::min((unsigned int)(charge2), nhydro) );
				atom.setFormalCharge( charge2 - delta );
				atom.updatePropertyCache(false); // Need to reset implicit/explict counts.
			}
		}
	}
}

void
remove_most_charges(::RDKit::RWMol & rdmol) {
	for ( core::Size ii(0); ii < rdmol.getNumAtoms(); ++ii ) {
		::RDKit::Atom & atom( *rdmol.getAtomWithIdx(ii) );
		int const & orig_charge( atom.getFormalCharge() );
		if ( orig_charge == 0 ) { continue; }
		atom.setFormalCharge( 0 );
		atom.updatePropertyCache(false); // Need to reset implicit/explict counts.
		unsigned int valence = atom.getTotalValence(); // What's the total valence.

		// NOTE: I (RM) don't claim that the handling here is all that great -- there's probably edge cases I'm not handling properly.
		// But it's a "good enough" right now to get things going.
		switch ( atom.getAtomicNum() ) {
		case 1 : // H
		case 9 : // F
		case 17 : // Cl
		case 35 : // Br
		case 53 : // I
			break; // Should never be charged - keep the zeroed one.
		case 6 : // C
			break; // Charges should never be necessary
		case 7 : // N
			if ( valence > 3 ) {
				atom.setFormalCharge( orig_charge );
				atom.updatePropertyCache(false);
			} else if ( valence == 3 ) {
				// Don't change charge if we're attached to a negatively charged oxygen with a single connection (us)
				for ( auto const & nbri: boost::make_iterator_range(rdmol.getAtomNeighbors(&atom)) ) {
					auto nbr = rdmol[nbri];
					if ( nbr->getAtomicNum() == 8 && nbr->getTotalValence() == 1 && nbr->getFormalCharge() == -1 ) {
						atom.setFormalCharge( orig_charge );
						atom.updatePropertyCache(false);
						break; // Out of the for loop
					}
				}
			}
			break;
		case 8 : // O
			if ( valence > 2 ) {
				atom.setFormalCharge( orig_charge );
				atom.updatePropertyCache(false);
			} else if ( valence == 1 ) {
				// Don't change charge if we're attached to a positively charged nitrogen (ylide/nitro)
				for ( auto const & nbri: boost::make_iterator_range(rdmol.getAtomNeighbors(&atom)) ) {
					auto nbr = rdmol[nbri];
					if ( nbr->getAtomicNum() == 7 && nbr->getFormalCharge() == +1 ) {
						atom.setFormalCharge( orig_charge );
						atom.updatePropertyCache(false);
						break; // Out of the for loop
					}
				}
			}
			break;
		case 15 : // P
			if ( valence > 3 ) {
				atom.setFormalCharge( orig_charge );
				atom.updatePropertyCache(false);
			}
			break;
		case 16 : // S
			if ( valence > 2 ) {
				atom.setFormalCharge( orig_charge );
				atom.updatePropertyCache(false);
			}
			break;
		default :
			atom.setFormalCharge( orig_charge ); // Don't mess with the charge
			atom.updatePropertyCache(false);
			break;
		}
	} // For all atoms
}


void
apply_charge_transforms( ::RDKit::RWMol & rdmol, ChargeTransformList const & transforms ) {
	for ( auto const & transform: transforms ) {
		::RDKit::ROMol const & site( *transform.first );
		int formal_charge = transform.second;

		std::vector< ::RDKit::MatchVectType > matchvect;
		::RDKit::SubstructMatch(rdmol, site, matchvect);
		for ( ::RDKit::MatchVectType const & match: matchvect ) {
			if ( match.empty() ) { continue; }

			// We assume the first atom in the match is the one we want.
			::RDKit::Atom & atom( *rdmol.getAtomWithIdx(match[0].second) );

			// Remove any explicit hydrogens that come along with positive charges
			// (addition of hydrogens happens implicitly)
			if ( formal_charge < atom.getFormalCharge() ) {
				core::Size charge_delta( atom.getFormalCharge() - formal_charge ); // Reduction in charge
				if ( charge_delta < atom.getNumExplicitHs() ) {
					// Subtract one proton for each unit of charge we're removing.
					atom.setNumExplicitHs( atom.getNumExplicitHs() - charge_delta );
				} else {
					atom.setNumExplicitHs( 0 ); // Remove all of them.
				}
			}
			atom.setNoImplicit(false); // Make sure we can add/remove implicit Hs.
			atom.setFormalCharge( formal_charge );
			atom.updatePropertyCache(false); // Need to reset implicit/explict counts.
		}
	}
}

void
reprotonate_rdmol(::RDKit::RWMol & rdmol) {

	// The transforms to use with the respective charge
	static const ChargeTransformList PH_TRANSFORMS {
		// First thing, remove any residual issues that might be found previously but are no longer useful.
		{ ::RDKit::RWMolOP(::RDKit::SmartsToMol("[O-1$([OD1]),$([OD2])]")), 0 }, // Singly connected (e.g. a hydroxyl no longer next to a double-bonded oxygen) or Doubly connected (e.g. a negatively charged oxygen that's now in a ring.)
		{ ::RDKit::RWMolOP(::RDKit::SmartsToMol("[S-1$([SD1]),$([SD2])]")), 0 },
		// Do we need to fix up any issue which may arise?
		{ ::RDKit::RWMolOP(::RDKit::SmartsToMol("[O-0D1$(O-[*]=[O,S,o,s])]")), -1 },
		{ ::RDKit::RWMolOP(::RDKit::SmartsToMol("[S-0D1$(S-[*]=[O,S,o,s])]")), -1 },
		{ ::RDKit::RWMolOP(::RDKit::SmartsToMol("[N+0&!$(N=[*])&!$(N#[*])&!$(N-[*]=[*])&!$(N-[*]#[*])&!$(N-[*]:[*])]")), +1 },
		{ ::RDKit::RWMolOP(::RDKit::SmartsToMol("[N+0D1$(N=C-[N+0])]")), +1 }, // amidine & guanidine
		{ ::RDKit::RWMolOP(::RDKit::SmartsToMol("[O-0D1H1$(O-[N+]=O)]")), -1 } // (e.g. nitro group ylide representation)
		};

	::RDKit::MolOps::removeHs( rdmol, false, false, /*sanitize=*/ false );
	// We need to do a limited amount of sanitization, else the queries will fail on certain cases
	softSanitize(rdmol);

	remove_excess_protons( rdmol );
	apply_charge_transforms( rdmol, PH_TRANSFORMS );

	::RDKit::MolOps::sanitizeMol( rdmol );
	::RDKit::MolOps::addHs( rdmol, false, /*addCoords=*/ true );
}

void
final_neutralize(
	RDKit::RWMOL_SPTR const & rdmol
) {
	bool progress = true;
	while ( progress ) {
		progress = false;
		for ( core::Size ii(0); ii < rdmol->getNumAtoms(); ++ii ) {
			::RDKit::Atom & a( *rdmol->getAtomWithIdx(ii) );
			int formal_charge( a.getFormalCharge() );
			if ( formal_charge == 0 ) {
				continue;
			} else if ( formal_charge > 0 ) {
				// Positive formal charge.
				// Heuristic - if there's an attached hydrogen, remove it and reduce the charge.
				if ( a.getNumExplicitHs() >= 1 ) {
					a.setNumExplicitHs( a.getNumExplicitHs() - 1 );
					a.setFormalCharge( formal_charge - 1 );
					progress = true;
				} else if ( a.getNumImplicitHs() >= 1 ) {
					// We don't need to reset implicit -- will be taken care of automatically
					a.setFormalCharge( formal_charge - 1 );
					progress = true;
				}
			} else {
				// Negative formal charge
				// Heuristic - add a hydrogen to neutralize negatively charged oxygen
				// (Chlorides & other negatively charged atoms aren't covered.)
				// Exception: if we're in a charge separation complex (e.g. Nitro groups)
				// We don't want to adjust the charges.
				if ( a.getAtomicNum() == 8 ) { // Oxygen
					bool ylid( false );
					::RDKit::ROMol::ADJ_ITER itr, itr_end;
					// Yes, address of. RDKit wants a pointer.
					for ( boost::tie(itr, itr_end) = rdmol->getAtomNeighbors(&a); itr != itr_end; ++itr ) {
						if ( rdmol->getAtomWithIdx(*itr)->getFormalCharge() >= 1 ) {
							ylid = true;
						}
					}
					if ( ! ylid ) {
						a.setNumExplicitHs( a.getNumExplicitHs() - 1 );
						a.setFormalCharge( formal_charge - 1 );
						progress = true;
					}
				}
			}
		}
	}
	// Re-sanitize the molecule -- the one from Roset
	try {
		::RDKit::MolOps::sanitizeMol(*rdmol);
	} catch (::RDKit::MolSanitizeException &se){
		TR.Error << "Cannot Sanitize molecule with RDKit after charge neutralization: " << se.what() << std::endl;
		TR.Error << "    molecule: " << ::RDKit::MolToSmiles( *rdmol ) << std::endl;
		utility_exit_with_message("Encountered molecule which cannot properly be represented in RDKit.");
	}
}

void
neutralize_rdmol(::RDKit::RWMol & rdmol, bool addHs) {

	static const ChargeTransformList CHARGE_FIXUPS {
		{ ::RDKit::RWMolOP(::RDKit::SmartsToMol("[O-0D1H1$(O-[N+]=O)]")), -1 } // (e.g. nitro group ylide representation)
		};

	::RDKit::MolOps::removeHs( rdmol, false, false, /*sanitize=*/ false );
	// We need to do a limited amount of sanitization, else the queries will fail on certain cases
	softSanitize(rdmol);

	// Do the "safe" way first, to make sure we remove any explicit H information.
	remove_excess_protons( rdmol );

	// Smash most of the rest of the charges, many of which will be spurious.
	remove_most_charges( rdmol );

	// Add back on some charges which we removed mistakenly
	apply_charge_transforms( rdmol, CHARGE_FIXUPS );

	::RDKit::MolOps::sanitizeMol( rdmol );
	if ( addHs ) {
		::RDKit::MolOps::addHs( rdmol, false, /*addCoords=*/ true );
	}
}

/// @brief Label a molecule with it's index values (for find_mapping, later)
void
label_with_index(  ::RDKit::ROMol & rdmol, std::string const & index_prop /* = "Orig_Index" */ ) {
	for ( core::Size ii(0); ii < rdmol.getNumAtoms(); ++ii ) {
		// Note: This *needs* to be an unsigned int, as otherwise the template instantiation
		// in external/rdkit/RDGeneral/Dict.cpp won't match up and you'll get a linker error
		rdmol.getAtomWithIdx(ii)->setProp< unsigned int >( index_prop, ii );
	}
}

/// @brief Convert the MatchVectType to an IndexIndex map, going query->molecule
/// Use -1 as the invalid value, as zero is a valid one.
core::chemical::IndexIndexMapping
convert_match_vect_to_index_index_map( ::RDKit::MatchVectType const & match_vect ) {
	core::chemical::IndexIndexMapping retval( core::Size(-1), core::Size(-1) ); // Use -1 as invalid
	for ( ::RDKit::MatchVectType::const_iterator iter( match_vect.begin() ), iter_end( match_vect.end() );
			iter != iter_end; ++iter ) {
		retval[ iter->first ] = iter->second;
	}
	return retval;
}

/// @brief Find a mapping from one RDMol to another
core::chemical::IndexIndexMapping
find_mapping( ::RDKit::ROMOL_SPTR from, ::RDKit::ROMOL_SPTR to, std::string const & index_prop /*= ""*/ ) {

	if ( index_prop.size() ) {
		core::chemical::IndexIndexMapping retval( core::Size(-1), core::Size(-1) );

		for ( core::Size ii(0); ii < to->getNumAtoms(); ++ii ) {
			if ( to->getAtomWithIdx(ii)->hasProp(index_prop) ) {
				core::Size orig_index = to->getAtomWithIdx(ii)->getProp< unsigned int >(index_prop);
				if ( from->getAtomWithIdx(orig_index)->hasProp(index_prop) &&
						from->getAtomWithIdx(orig_index)->getProp< unsigned int >(index_prop) == orig_index ) {
					retval[ orig_index ] = ii;
				} else {
					TR.Warning << "Original index mismatch in mapping. Property " << index_prop << " not properly set on 'from' molecule." << std::endl;
				}
			}
		}

		if ( ! retval.empty() ) {
			TR << "Using mapping for property " << index_prop << " to map before/after atoms." << std::endl;
			return retval;
		}
		TR.Warning << "No mapping entries for " << index_prop << " found, attempting to use substructure matching for mapping." << std::endl;
		// Fall through.
	}

	::RDKit::MatchVectType match_vect;

	// Try the easy way first - one molecule is a subset of the other.
	if ( to->getNumHeavyAtoms() >= from->getNumHeavyAtoms() &&
			::RDKit::SubstructMatch( *to, *from, match_vect, true, true /* useChirality */ ) ) {
		// 'from' is a substructure of 'to'.
		return convert_match_vect_to_index_index_map( match_vect );
	}

	if ( from->getNumHeavyAtoms() >= to->getNumHeavyAtoms() &&
			::RDKit::SubstructMatch( *from, *to, match_vect, true, true /* useChirality */ ) ) {
		// 'to' is a substructure of 'from' - need to reverse the match vector orientation
		return convert_match_vect_to_index_index_map( match_vect ).reverse();
	}

	TR.Warning << "Using Slow, MCS-based substructure mapping. (There's probably a better way to do this.)" << std::endl;
	// as in, whatever is calling this should be keeping better track of things.
	// (Also, it looks like there's a more efficient MCS algorithm: Chem.rdFMCS -- http://www.rdkit.org/docs/GettingStartedInPython.html#maximum-common-substructure)

	// Have to do this the hard way - find the maximum common substructure, then find the indices of each
	std::vector< ::RDKit::ROMOL_SPTR > molvector;
	molvector.push_back( from );
	molvector.push_back( to );
	::RDKit::MCSParameters mcsparams;
	mcsparams.AtomCompareParameters.MatchChiralTag = true; // Match atom chirality
	::RDKit::MCSResult mcsresult( ::RDKit::findMCS( molvector, &mcsparams ) ); // Yes, pass a raw pointer
	TR << "MCS Smarts: " << mcsresult.SmartsString << std::endl;
	::RDKit::RWMOL_SPTR mcs( ::RDKit::SmartsToMol( mcsresult.SmartsString ) );
	::RDKit::MatchVectType match_vect_from, match_vect_to;
	bool from_match( ::RDKit::SubstructMatch( *from, *mcs, match_vect_from, true, true ) );
	bool to_match( ::RDKit::SubstructMatch( *to, *mcs, match_vect_to, true, true ) );
	if ( ! from_match || ! to_match ) {
		// The MCS can't be mapped to either the reactant or product, for some reason
		TR.Warning << "Substructure mapping not found!" << std::endl;
		if ( TR.Debug.visible() ) {
			TR.Debug << "     MCS:         " << ::RDKit::MolToSmarts( *mcs ) << std::endl;
			TR.Debug << "     'From' mol:  " << ::RDKit::MolToSmiles( *from ) << std::endl;
			TR.Debug << "     'To' mol:    " << ::RDKit::MolToSmiles( *to ) << std::endl;
		}
		core::chemical::IndexIndexMapping empty_map( core::Size(-1), core::Size(-1) );
		return empty_map;
	}
	// from<-query->to
	return combine( convert_match_vect_to_index_index_map(match_vect_from).reverse(), convert_match_vect_to_index_index_map(match_vect_to) );
}

///@brief Load molecules from file and append them to mol_vector
void load_sdf( std::string const & filename, utility::vector1< ::RDKit::ROMolOP > & mol_vector, bool removeHs) {
	utility::io::izstream data( filename );
	if ( ! data.good() ) {
		TR.Error << "Cannot open sdf file '" << filename << "'" << std::endl;
		utility_exit_with_message("Cannot open sdf file "+filename);
	}

	::RDKit::SDMolSupplier supplier( &data(), /*takeOwnership=*/ false, /*sanitize=*/ true, /*removeHs=*/ false); // Yes, pointer to std::istream - that's what RDKit wants
	while ( ! supplier.atEnd() ) {
		::RDKit::ROMolOP mol( supplier.next() );
		if ( mol ) {
			if ( removeHs ) {
				TR.Trace << "Removing hydrogens from loaded SDF structures." << std::endl;
				::RDKit::ROMolOP newmol;
				// When we remove Hs, make them "explicits"
				try {
					// It's an ROMol, so we can't modify it in place -- it makes a new item which we're then responsible for deleting.
					newmol = ::RDKit::ROMolOP( ::RDKit::MolOps::removeHs(*mol, /*implicitOnly=*/ false, /*updateExplicitCount=*/ true) );
				} catch (::RDKit::MolSanitizeException &se){
					TR.Warning << "Error when removing hydrogens from read-in SDF: " << se.what() << std::endl;
					continue;
				}
				mol = newmol; // Swap.
			}
			mol_vector.push_back( mol );
		}
	}
}

/// @brief Return a set containing all the valid names for the rdkit_metric() function mapped to short descriptions
std::map< std::string, std::string >
get_metric_names() {
	std::map< std::string, std::string > names;
	names["MolWt"] = "Molecular weight";
	names["HeavyAtomCount"] = "Number of heavy atoms in the molecule";
	names["HeavyAtomMolWt"] = "Molecular weight of just the heavy atoms";
	names["NOCount"] = "Number of standard Lipinski hydrogen bond acceptors";
	names["NHOHCount"] = "Number of standard Lipinski hydrogen bond donors";
	names["NumHAcceptors"] = "Number of standard Lipinski hydrogen bond acceptors";
	names["NumHDonors"] = "Number of standard Lipinski hydrogen bond donors";
	names["NumRotatableBonds"] = "Number of rotatable bonds (excluding amides, esters, etc.)";
	names["NumHeteroatoms"] = "Number of heteroatoms";
	names["RingCount"] = "Number of rings";
	names["NumAromaticRings"] = "Number of aromatic rings";
	names["NumAliphaticRings"] = "Number of rings with at least one non-aromatic bond";
	names["NumSaturatedRings"] = "Number of saturated rings";
	names["NumHeterocycles"] = "Number of heterocycle rings";
	names["TPSA"] = "Total polar surface area as in Ertl et al. J.Med.Chem 43:3714";
	names["LabuteASA"] = "Approximate surface area as in Labute J.Mol.Graph 18:464";
	names["MolLogP"] = "clogP value as in Wildman & Crippen JCICS 39:868";
	names["MolMR"] = "Molar refractivity value as in Wildman & Crippen JCICS 39:868";
	return names;
}

/// @brief Return the value of a given RDKit metric for the given mol
core::Real
rdkit_metric(::RDKit::ROMol const & mol, std::string const & metric ) {
	// If you add something to this if/else block, you need to add it above as well.
	if ( metric == "MolWt" ) {
		return ::RDKit::Descriptors::calcAMW(mol, /*onlyHeavy=*/ false);
	} else if ( metric == "HeavyAtomCount" ) {
		return mol.getNumHeavyAtoms();
	} else if ( metric == "HeavyAtomMolWt" ) {
		return ::RDKit::Descriptors::calcAMW(mol, /*onlyHeavy=*/ true);
	} else if ( metric == "NOCount" ) {
		return ::RDKit::Descriptors::calcLipinskiHBA(mol);
	} else if ( metric == "NHOHCount" ) {
		return ::RDKit::Descriptors::calcLipinskiHBD(mol);
	} else if ( metric == "NumHAcceptors" ) {
		return ::RDKit::Descriptors::calcNumHBA(mol);
	} else if ( metric == "NumHDonors" ) {
		return ::RDKit::Descriptors::calcNumHBD(mol);
	} else if ( metric == "NumRotatableBonds" ) {
		return ::RDKit::Descriptors::calcNumRotatableBonds(mol);
	} else if ( metric == "NumHeteroatoms" ) {
		return ::RDKit::Descriptors::calcNumHeteroatoms(mol);
	} else if ( metric == "RingCount" ) {
		return ::RDKit::Descriptors::calcNumRings(mol);
	} else if ( metric == "NumAromaticRings" ) {
		return ::RDKit::Descriptors::calcNumAromaticRings(mol);
	} else if ( metric == "NumAliphaticRings" ) {
		return ::RDKit::Descriptors::calcNumAliphaticRings(mol);
	} else if ( metric == "NumSaturatedRings" ) {
		return ::RDKit::Descriptors::calcNumSaturatedRings(mol);
	} else if ( metric == "NumHeterocycles" ) {
		return ::RDKit::Descriptors::calcNumHeterocycles(mol);
	} else if ( metric == "TPSA" ) {
		return ::RDKit::Descriptors::calcTPSA(mol);
	} else if ( metric == "LabuteASA" ) {
		return ::RDKit::Descriptors::calcLabuteASA(mol);
	} else if ( metric == "MolLogP" || metric == "MolMR" ) {
		double logp, mr;
		::RDKit::Descriptors::calcCrippenDescriptors(mol, logp, mr);
		if (  metric == "MolLogP" ) {
			return logp;
		} else if ( metric == "MolMR" ) {
			return mr;
		}
	}
	utility_exit_with_message("Rosetta does not understand the RDKit metric '"+metric+"'");
	return 0;
}


} // namespace rdkit
} // namespace chemical
} // namespace core
