// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/substitution_support.cc
/// @brief use RDKit to substitute items based on matched templates
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/substitution_support.hh>

#include <core/chemical/rdkit/RDKit.fwd.hh>
#include <core/chemical/rdkit/util.hh>

#include <basic/Tracer.hh>

#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/Geometry/Transform3D.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h> // For MolToSmiles


namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.substructure_support");

::RDKit::Bond::BondType
get_first_bondtype( ::RDKit::ROMol const & mol, ::RDKit::Atom const * atom ) {
	::RDKit::ROMol::OBOND_ITER_PAIR bond_itrs( mol.getAtomBonds( atom ) );
	// We assume that each dummy atom is only bonded to a single other atom.
	debug_assert( bond_itrs.first != bond_itrs.second );
	return mol[ *bond_itrs.first ]->getBondType();
}

::RDKit::Bond::BondType
get_first_bondtype( ::RDKit::ROMol const & mol, unsigned int idx ) {
	return get_first_bondtype( mol, mol.getAtomWithIdx( idx ) );
}

::RDKit::Bond const &
get_first_bond( ::RDKit::ROMol const & mol, unsigned int idx ) {
	::RDKit::Atom const * atom( mol.getAtomWithIdx( idx ) );
	::RDKit::ROMol::OBOND_ITER_PAIR bond_itrs( mol.getAtomBonds( atom ) );
	// We assume that each dummy atom is only bonded to a single other atom.
	debug_assert( bond_itrs.first != bond_itrs.second );
	return *mol[ *bond_itrs.first ];
}

unsigned int
get_bonded_atom( ::RDKit::ROMol const & mol, unsigned int idx ) {
	::RDKit::Bond const & bond( get_first_bond( mol, idx) );
	return bond.getOtherAtomIdx( idx );
}

bool
is_dummy( ::RDKit::ROMol const & mol, unsigned int idx ) {
	debug_assert( idx < mol.getNumAtoms() );
	return mol.getAtomWithIdx(idx)->getAtomicNum() == 0;
}

void
find_dummies( ::RDKit::ROMol const & mol, utility::vector1< unsigned int > & dummy_list ) {
	dummy_list.clear();
	for ( core::Size idx(0); idx < mol.getNumAtoms(); ++idx ) { // std::vector
		if ( is_dummy( mol, idx) ) {
			dummy_list.push_back( idx );
		}
	}
}

bool
bond_is_in_ring( ::RDKit::ROMol const & rdmol, unsigned int atm1, unsigned int atm2 ) {
	if ( rdmol.getRingInfo()->isInitialized() ) {
		::RDKit::MolOps::findSSSR(rdmol);
	}
	debug_assert( rdmol.getBondBetweenAtoms(atm1, atm2) );
	unsigned int bond_idx( rdmol.getBondBetweenAtoms(atm1, atm2)->getIdx() );
	return rdmol.getRingInfo()->numBondRings( bond_idx ) != 0;
}

std::map< unsigned int, unsigned int >
convert_matchvect_to_map( ::RDKit::MatchVectType const & pairings ) {
	std::map< unsigned int, unsigned int > retval;
	for ( core::Size ii(0); ii < pairings.size(); ++ii ) { // std::vector
		int t_idx( pairings[ii].first );
		int m_idx( pairings[ii].second );
		retval[ m_idx ] = t_idx;
	}
	return retval;
}

unsigned int
AtomSubstitution::idx(MoleculeSelection sele) const {
	switch( sele ) {
	case OriginalMol :
		return mdx();
	case TemplateMol :
		return tdx();
	case ReplaceMol :
		return rdx();
	case NewMol :
		return ndx();
	default :
		utility_exit_with_message("Not a valid MoleculeSelection");
		return AtomSubstitution::invalid_index;
	}
}

void
AtomSubstitution::set_idx(MoleculeSelection sele, unsigned int setting) {
	switch( sele ) {
	case OriginalMol :
		set_mdx(setting);
		break;
	case TemplateMol :
		set_tdx(setting);
		break;
	case ReplaceMol :
		set_rdx(setting);
		break;
	case NewMol :
		set_ndx(setting);
		break;
	default :
		utility_exit_with_message("Not a valid MoleculeSelection");
	}
}

MoleculeSubstitution::MoleculeSubstitution( ::RDKit::ROMolOP mol ):
	mol_(mol)
{
	// Add entries for all mol items
	for ( core::Size mdx(0); mdx < mol->getNumAtoms(); ++mdx ) {
		AtomSubstitutionOP atom_sub( new AtomSubstitution( mdx ) );
		atom_substitutions_.push_back( atom_sub );
		by_mdx_[ mdx ] = atom_sub;
	}
}

/// @brief Add a template to this MoleculeSubstitution.
void
MoleculeSubstitution::add_template(::RDKit::ROMolOP templt, ::RDKit::MatchVectType const & pairings )
{
	templt_ = templt;

	// Set pairings (should have all template entries)
	for ( core::Size ii(0); ii < pairings.size(); ++ii ) { // std::vector
		int tdx( pairings[ii].first );
		int mdx( pairings[ii].second );
		// Should have the mdx entry already
		debug_assert( by_mdx_.count( mdx ) == 1 );
		by_tdx_[ tdx ] = by_mdx_[ mdx ];
		by_tdx_[ tdx ]->set_tdx( tdx );
	}

	find_dummies( *templt_, template_dummies_ );
}

::RDKit::MatchVectType
MoleculeSubstitution::make_match_vector() const {
	::RDKit::MatchVectType match_vect;
	for ( core::Size ii(1); ii <= atom_substitutions_.size(); ++ii ) {
		if ( atom_substitutions_[ii]->mdx() != AtomSubstitution::invalid_index &&
				atom_substitutions_[ii]->tdx() != AtomSubstitution::invalid_index ) {
			match_vect.push_back( std::pair<int,int>( atom_substitutions_[ii]->tdx(), atom_substitutions_[ii]->mdx() ) );
		}
	}
	return match_vect;
}

/// @brief Create a mapping from mdx to ndx
core::chemical::IndexIndexMapping
MoleculeSubstitution::find_mdx_to_ndx_mapping() const {
	core::chemical::IndexIndexMapping retval( core::Size(-1), core::Size(-1) ); // Use -1 as invalid as zero is valid
	for ( core::Size ii(1); ii <= atom_substitutions_.size(); ++ii ) {
		if ( atom_substitutions_[ii]->mdx() != AtomSubstitution::invalid_index &&
				atom_substitutions_[ii]->ndx() != AtomSubstitution::invalid_index ) {
			retval[ atom_substitutions_[ii]->mdx() ] = atom_substitutions_[ii]->ndx();
		}
	}
	return retval;
}

/// @brief Make a *new* AtomSubstitution class with the replacement information
MoleculeSubstitutionOP
MoleculeSubstitution::add_replacement(
	::RDKit::ROMolOP replacement,
	std::map< unsigned int, unsigned int > & r_to_t_mapping
) const {
	// First we need to deep copy the MoleculeSubstitution
	MoleculeSubstitutionOP new_molsub( new MoleculeSubstitution );
	new_molsub->mol_ = mol_;
	new_molsub->templt_ = templt_;
	new_molsub->template_dummies_ = template_dummies_;
	for ( core::Size ii(1); ii <= atom_substitutions_.size(); ++ii ) {
		AtomSubstitutionOP atomsub_copy( new AtomSubstitution( *atom_substitutions_[ii] ) );
		new_molsub->atom_substitutions_.push_back( atomsub_copy );
		new_molsub->by_mdx_[ atomsub_copy->mdx() ] = atomsub_copy;
		new_molsub->by_tdx_[ atomsub_copy->tdx() ] = atomsub_copy;
	}

	new_molsub->replace_ = replacement;
	find_dummies( *replacement, new_molsub->replace_dummies_ );
	// Now make entries for each replacement atom, using existing items if possible.
	for ( unsigned int rdx(0); rdx < replacement->getNumAtoms(); ++rdx ) { // std::vector
		if ( r_to_t_mapping.count( rdx ) > 0 ) {
			unsigned int tdx( r_to_t_mapping[ rdx ] );
			debug_assert( new_molsub->by_tdx_.count(tdx) );
			new_molsub->by_tdx_[ tdx ]->set_rdx( rdx );
			new_molsub->by_rdx_[ rdx ] = new_molsub->by_tdx_[ tdx ];
		} else {
			AtomSubstitutionOP rdx_atomsub( new AtomSubstitution );
			rdx_atomsub->set_rdx( rdx );
			new_molsub->atom_substitutions_.push_back(rdx_atomsub);
			new_molsub->by_rdx_[ rdx ] = rdx_atomsub;
		}
	}

	return new_molsub;
}

::RDKit::ROMolOP
MoleculeSubstitution::get_romol( MoleculeSelection sele ) {
	switch( sele ) {
	case OriginalMol :
		return mol();
	case TemplateMol :
		return templt();
	case ReplaceMol :
		return replace();
	case NewMol :
		utility_exit_with_message("NewMol entry not set!");
	default :
		utility_exit_with_message("Not a valid MoleculeSelection");
	}
}

bool
MoleculeSubstitution::tdx_is_dummy( unsigned int tdx ) const {
	return is_dummy( *templt_, tdx );
}

bool
MoleculeSubstitution::rdx_is_dummy( unsigned int rdx ) const {
	return is_dummy( *replace_, rdx );
}

AtomSubstitution &
MoleculeSubstitution::substitution_for_idx( MoleculeSelection sele, unsigned int idx ) {
	switch( sele ) {
	case OriginalMol :
		return substitution_for_mdx(idx);
	case TemplateMol :
		return substitution_for_tdx(idx);
	case ReplaceMol :
		return substitution_for_rdx(idx);
	case NewMol :
		utility_exit_with_message("ndx's aren't indexed!");
	default :
		utility_exit_with_message("Not a valid MoleculeSelection");
	}
}

AtomSubstitution &
MoleculeSubstitution::substitution_for_mdx( unsigned int idx ) {
	debug_assert( by_mdx_.count( idx ) > 0 );
	return *by_mdx_[ idx ];
}
AtomSubstitution &
MoleculeSubstitution::substitution_for_tdx( unsigned int idx ) {
	debug_assert( by_tdx_.count( idx ) > 0 );
	return *by_tdx_[ idx ];
}
AtomSubstitution &
MoleculeSubstitution::substitution_for_rdx( unsigned int idx ) {
	debug_assert( by_rdx_.count( idx ) > 0 );
	return *by_rdx_[ idx ];
}

/// @brief Pick a template to use, return it and the atom pairing of template->rdmol
MoleculeSubstitutionOP
pick_template(
	::RDKit::ROMolOP rdmol,
	utility::vector1< ::RDKit::ROMolOP > & templates,
	bool dummy_only
)
{
	utility::vector1< ::RDKit::ROMolOP > matching_templates;
	utility::vector1< ::RDKit::MatchVectType > pairings;

	TR.Trace << "Finding template for " << ::RDKit::MolToSmiles( *rdmol ) << std::endl;
	for ( core::Size tt(1); tt <= templates.size(); ++tt ) {
		std::vector< ::RDKit::MatchVectType > matches;
		::RDKit::SubstructMatch(*rdmol, *templates[tt], matches, /*uniquify=*/ false); // Want redundancy in orientation
		core::Size nmatches(0);
		for ( core::Size ii(0); ii < matches.size(); ++ii ) { // NOTE: std::vector iteration!
			bool match_good( true );
			// Want to make sure we don't have any heavy atoms connected to the template by non-dummy means.
			if ( dummy_only ) {
				std::map< unsigned int, unsigned int > match_map( convert_matchvect_to_map( matches[ii] ) ); // molecule:template mapping
				for ( unsigned int jj(0); jj < rdmol->getNumAtoms(); ++jj ) { // No guarantees that RDKit atoms are ordered
					if ( match_map.count(jj) == 1 || rdmol->getAtomWithIdx(jj)->getAtomicNum() == 1 ) {
						// Ignore atoms matched to the template, and hydrogens
						continue;
					}
					// We're a heavy atom not matched to the template - we want to make sure we're not bonded to a non-dummy atom in the template.
					for ( ::RDKit::ROMol::OBOND_ITER_PAIR bonds( rdmol->getAtomBonds( rdmol->getAtomWithIdx(jj) ) );
							bonds.first != bonds.second; ++bonds.first ) {
						unsigned int other_jj( (*rdmol)[ *bonds.first ]->getOtherAtomIdx(jj) );
						if ( match_map.count(other_jj) && templates[tt]->getAtomWithIdx( match_map[other_jj])->getAtomicNum() != 0 ) {
							match_good = false;
							break;
						}
					}
					if ( ! match_good ) break;
				}
			}

			if ( match_good ) {
				++nmatches;
				matching_templates.push_back( templates[tt] );
				pairings.push_back( matches[ii] );
			}
		}
		TR.Trace << "Template " << tt << " -- " << ::RDKit::MolToSmiles( *templates[tt] ) << " -- has " << nmatches << " matches." << std::endl;
	}

	TR.Debug << "Found " << matching_templates.size() << " substructures compatible with current molecule." << std::endl;

	if ( matching_templates.size() == 0 ) {
		TR << "[NOTICE]: No matching templates found." << std::endl;
		return MoleculeSubstitutionOP();
	}

	debug_assert( matching_templates.size() == pairings.size() );

	core::Size picked( numeric::random::random_range( 1, matching_templates.size() ) );

	MoleculeSubstitutionOP molsub( new MoleculeSubstitution( rdmol ) );
	molsub->add_template(matching_templates[ picked ], pairings[ picked ]);
	return molsub;
}

MoleculeSubstitutionOP
test_replacement(
	MoleculeSubstitutionOP current_molsub,
	::RDKit::ROMolOP possible_replacement,
	core::Real dist_threshold
) {
	::RDKit::ROMolOP mol( current_molsub->mol() );
	::RDKit::ROMolOP templt( current_molsub->templt() );

	if ( templt.get() == possible_replacement.get() ) {
		// Replacement is the template - don't propose for replacement
		return nullptr;
	}

	::RDKit::Conformer const & replace_conf( possible_replacement->getConformer() );
	::RDKit::Conformer const & templt_conf( templt->getConformer() );

	utility::vector1< unsigned int > const & template_dummies( current_molsub->template_dummies() );
	utility::vector1< unsigned int > replacement_dummies;
	find_dummies( *possible_replacement, replacement_dummies );

	std::map< unsigned int, unsigned int > r_to_t_mapping;

	// For each matched dummy, find the closest replacement dummy.
	for ( core::Size ii(1); ii <= template_dummies.size(); ++ii ) {
		unsigned int tdx( template_dummies[ii] );
		unsigned int mdx( current_molsub->substitution_for_tdx(tdx).mdx() );
		::RDKit::Bond::BondType t_bt( get_first_bondtype( *current_molsub->templt(), tdx ) );
		::RDGeom::Point3D const & t_pos( templt_conf.getAtomPos( tdx ) );
		unsigned int min_rdx = AtomSubstitution::invalid_index;
		core::Real min_dist = dist_threshold;
		for ( core::Size jj(1); jj <= replacement_dummies.size(); ++jj ) {
			::RDKit::Bond::BondType r_bt( get_first_bondtype( *possible_replacement, replacement_dummies[jj] ) );
			// The bond types of the two dummies need to match.
			if ( r_bt != t_bt ) { continue; }
			::RDGeom::Point3D const & r_pos( replace_conf.getAtomPos( replacement_dummies[jj] ) );
			core::Real dist( (r_pos-t_pos).length() );
			if ( dist < min_dist ) {
				min_dist = dist;
				min_rdx = replacement_dummies[jj];
			}
		}
		// See if we have successfully found a match
		if ( min_rdx == AtomSubstitution::invalid_index ) {
			if ( mol->getAtomWithIdx( mdx )->getAtomicNum() <= 1 ) {
				// Stub to remove is just a hydrogen - just ignore
				continue;
			} else {
				// Can't find appropriate stub for heavyatom - the replacement won't work
				return nullptr;
			}
		} else if ( r_to_t_mapping.count( min_rdx ) != 0 ) {
			// We're doubling up on a replacement dummy
			if ( current_molsub->mol()->getAtomWithIdx( mdx )->getAtomicNum() <= 1 ) {
				// Current stub is just a hydrogen - just ignore it
				continue;
			}
			unsigned int old_tdx( r_to_t_mapping[ min_rdx ] );
			if ( mol->getAtomWithIdx( current_molsub->substitution_for_tdx(old_tdx).mdx() )->getAtomicNum() > 1 ) {
				// OOPS - we're trying to place two heavy atoms on this same item.
				// (Theoretically we could match up with second best, but I'm not doing that now.)
				return nullptr;
			} else {
				// We're heavy and the previous one is hydrogen - replace it
				r_to_t_mapping[ min_rdx ] = tdx;
			}
		} else {
			// First assigment to this valid rdx
			r_to_t_mapping[ min_rdx ] = tdx;
		}
	}

	MoleculeSubstitutionOP replacement_molsub( current_molsub->add_replacement( possible_replacement, r_to_t_mapping ) );

	// We should also fail if we have any unpaired replacement dummies which cannot be replaced with hydrogens.
	utility::vector1< unsigned int > const & replace_dummies( replacement_molsub->replace_dummies() );
	for ( core::Size jj(1); jj <= replace_dummies.size(); ++jj ) {
		unsigned int tdx( replacement_molsub->substitution_for_rdx( replace_dummies[jj] ).tdx() );
		if ( tdx == AtomSubstitution::invalid_index ) {
			::RDKit::Bond::BondType r_bt( get_first_bondtype( *possible_replacement, replace_dummies[jj] ) );
			if ( r_bt != ::RDKit::Bond::SINGLE ) {
				// We can't replace a missing non-single bond with a hydrogen
				return nullptr;
			}
		}
	}

	return replacement_molsub;
}

MoleculeSubstitutionOP
pick_replacement(
	MoleculeSubstitutionOP current_molsub,
	utility::vector1< ::RDKit::ROMolOP > & possible_replacements,
	core::Real distance_threshold,
	std::string weighting_property /*= ""*/
) {
	utility::vector1< MoleculeSubstitutionOP > replacements;
	numeric::random::WeightedSampler replacements_sampler;

	for ( core::Size ii(1); ii <= possible_replacements.size(); ++ii ) {
		TR.Trace << "Testing for match: possible replacement #" << ii << " -- " << ::RDKit::MolToSmiles( *possible_replacements[ii] ) << std::endl;
		MoleculeSubstitutionOP replace_molsub( test_replacement( current_molsub, possible_replacements[ii], distance_threshold ) );
		if ( replace_molsub ) {
			core::Real weight( 1.0 );
			if ( weighting_property.size() && possible_replacements[ii]->hasProp(weighting_property) ) {
				weight = possible_replacements[ii]->getProp<core::Real>(weighting_property);
			}
			TR.Trace << "Possible replacement #" << ii << " is a match. " << std::endl;
			replacements.push_back( replace_molsub );
			replacements_sampler.add_weight( weight );
		}
	}

	if ( replacements.size() == 0 ) {
		TR << "[NOTICE]: No replacements found." << std::endl;
		return nullptr;
	}

	// Pick a replacement from the possibilities
	core::Size replacement_num( replacements_sampler.random_sample() );
	return replacements[ replacement_num ];
}

/// @brief Copies all atom and bonds attached to start in the source molecule to the new molecule.
/// If skip is a valid index, it will be supressed in the source molecule, and will not be counted for copying or for attachment purposes
/// The transform is applied to the coordinates of the atom to get the coordinates of the new atom.
void
copy_attached_atoms( MoleculeSubstitution & molsub, MoleculeSelection source, ::RDGeom::Transform3D const & transform, unsigned int start, unsigned int skip ) {
	::RDKit::ROMolOP source_mol( molsub.get_romol(source) );
	::RDKit::RWMolOP new_mol( molsub.newmol() );

	debug_assert( start < source_mol->getNumAtoms() );
	debug_assert( skip != start );

	utility::vector1< unsigned int > atoms_to_copy;
	atoms_to_copy.push_back( start );

	for ( core::Size anum(1); anum <= atoms_to_copy.size(); ++anum ) { // Add to list while iterating - don't cache size!
		// Copy atom.
		unsigned int idx( atoms_to_copy[ anum ] );
		AtomSubstitution & atom_sub( molsub.substitution_for_idx( source, idx ) );
		if ( atom_sub.ndx() != AtomSubstitution::invalid_index ) { // Already made a copy - skip it.
			continue;
		}
		unsigned int ndx( new_mol->addAtom( source_mol->getAtomWithIdx( idx ) ) );
		::RDGeom::Point3D old_pos( source_mol->getConformer().getAtomPos( idx ) );
		new_mol->getConformer().setAtomPos( ndx, transform * old_pos );
		atom_sub.set_ndx( ndx );
		//TR << "Adding atom " << idx << " from molecule " << source << " as atom " << ndx << std::endl;

		// Now add any bonded atoms that aren't supressed.
		// (We're non-ring, so we shouldn't cross that bond by another route.)
		for ( ::RDKit::ROMol::OBOND_ITER_PAIR bonds( source_mol->getAtomBonds( source_mol->getAtomWithIdx( idx ) ) );
				bonds.first != bonds.second; ++bonds.first ) {
			unsigned int idx_bnd( (*source_mol)[ *bonds.first ]->getOtherAtomIdx( idx ) );
			if ( idx_bnd == skip ) { // Don't follow trace through skipped atom. (Should never trigger for invalid entry skip
				continue;
			}
			atoms_to_copy.push_back( idx_bnd ); // Copy the bonded atom
		}
	}

	// Convert to a set to uniquify and for fast lookup
	std::set< unsigned int > copied_atoms(atoms_to_copy.begin(), atoms_to_copy.end());

	for ( std::set< unsigned int>::const_iterator itr( copied_atoms.begin() ), itr_end( copied_atoms.end() );
			itr != itr_end; ++itr ) {
		for ( ::RDKit::ROMol::OBOND_ITER_PAIR bonds( source_mol->getAtomBonds( source_mol->getAtomWithIdx( *itr ) ) );
				bonds.first != bonds.second; ++bonds.first ) {
			unsigned int idx_bnd( (*source_mol)[ *bonds.first ]->getOtherAtomIdx( *itr ) );
			if ( copied_atoms.count( idx_bnd ) != 1 ) { // Don't make bonds to atoms we didn't copy.
				continue;
			}
			unsigned int bgn_ndx( molsub.substitution_for_idx(source, *itr).ndx() );
			unsigned int end_ndx( molsub.substitution_for_idx(source, idx_bnd).ndx() );
			debug_assert( bgn_ndx != AtomSubstitution::invalid_index );
			debug_assert( end_ndx != AtomSubstitution::invalid_index );
			if ( new_mol->getBondBetweenAtoms(bgn_ndx, end_ndx) == nullptr ) { // Bond should not already exist {
				//TR << "Adding bond between new atoms " << bgn_ndx << " and " << end_ndx << " ( was " << *itr << " -- " << idx_bnd << " ) " << std::endl;
				new_mol->addBond( bgn_ndx, end_ndx, source_mol->getBondBetweenAtoms(*itr,idx_bnd)->getBondType() );
			}
		}
	}
}

}
}
