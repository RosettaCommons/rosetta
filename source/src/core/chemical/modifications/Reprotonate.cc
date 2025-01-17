// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/modifications/Reprotonate.cc
/// @brief  Class to add hydrogens to a ResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <core/chemical/modifications/Reprotonate.hh>

#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/Elements.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>
#include <core/chemical/modifications/ValenceHandler.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/atomtype_support.hh>
#include <core/chemical/bond_support.hh>

#include <basic/Tracer.hh>

namespace core {
namespace chemical {
namespace modifications {

static basic::Tracer TR("core.chemical.modifications.Reprotonate");

/// @brief Comes up with a new name for a hydrogen to add. (One that isn't used.)
std::string
get_new_hydrogen_name( MutableResidueType & res ) {
	core::Size index( 0 );
	std::string name;
	do {
		++index;
		name = "H" + std::to_string( index + 1 );
	} while( res.has( name ) );
	// Need to adjust the naming to fit the four-sized value.
	if ( name.size() == 2 ) {
		name = " " + name + " ";
	} else if ( name.size() == 3 ) {
		name = " " + name;
	}
	return name ;
}

/// @brief Add hydrogens to the given Atom, based on existing geometry and Gasteiger atom type
/// Returns the number of hydrogens added.
core::Size
add_hydrogens_to_atom( MutableResidueType & res, VD atom_vd ) {
	utility::vector1<numeric::xyzVector<core::Real> > coords( determine_coordinates(res, atom_vd) );
	TR << "Adding " << coords.size() << " hydrogens to atom " << res.atom_name( atom_vd ) << " of type " << res.atom( atom_vd ).gasteiger_atom_type()->get_name()
		<< " and " << res.number_bonded_hydrogens( atom_vd ) << " existing hydrogens " << std::endl;
	for ( Size new_hydrogens=1; new_hydrogens <= coords.size(); ++new_hydrogens ) {
		std::string name = get_new_hydrogen_name( res );
		TR.Debug << "Added hydrogen name: " << name << std::endl;
		VD vertex_hydrogen = res.add_atom(name);
		res.atom(vertex_hydrogen).ideal_xyz(coords[new_hydrogens]);
		res.atom(vertex_hydrogen).element_type( res.element_set().element( element::H ) ); //dont forget to add the element type!
		res.add_bond(vertex_hydrogen, atom_vd, SingleBond);
	}
	return coords.size();
}

/// @details The current targets are:
/// Deprotonation:
/// * Carboxylates
/// * Thiocarboxylates ( C(=O)-SH as well as C(=S)-OH, if existant)
/// * Sulfates/Phosphates/Nitrates/etc (as well as thio equivalents)
/// * NOT Phenols
///
/// Protonation
/// * (aliphatic) Amines
/// * NOT Amides, anilines, or conjugated amines
///
/// Protonation/deprotonation will not happen on atoms with existing formal charges.
/// The appropriate formal charge will be added to the reprotonated atom.
///
void
Reprotonate::apply( MutableResidueType & res ) {
	// Heavy atoms on which to protonate/deprotonate
	// We notate then alter to avoid modifying the residue type as we're iterating over it.
	utility::vector1< VD > protonate;
	utility::vector1< VD > deprotonate;

	mapping_ = core::chemical::VDVDMapping(true); // Start with the identity mapping.

	for ( VIterPair iters(res.atom_iterators()); iters.first != iters.second; ++iters.first ) {
		VD atom_vd( *iters.first );
		core::chemical::Atom const & atom( res.atom( atom_vd ) );
		// Don't touch things with formal charges already.
		if ( atom.formal_charge() != 0 ) { continue; }

		if ( atom.element() == element::O || atom.element() == element::S ) {

			// Heuristic - deprotonate an oxygen/sulfur if it's bonded to an element which is double bonded
			// to an oxygen/sulfur. This should handle (thio)carboxylates as well as phosphates, nitrates, sulfates, etc.
			if ( res.nbonds( atom_vd ) != 2 ) { continue; } // Only deal with ~X(=O)-O-H motifs
			VD hydrogen( MutableResidueType::null_vertex );
			VD electron_sink( MutableResidueType::null_vertex );
			for ( OutEdgeIterPair eiters( boost::out_edges(atom_vd, res.graph()) ); eiters.first != eiters.second; ++eiters.first ) {
				if ( res.bond( *eiters.first ).order() != SingleBondOrder ) {
					// Oops - the atom of interest needs to be single bonded only
					electron_sink = hydrogen = MutableResidueType::null_vertex;
					break;
				}
				VD attached = boost::target(*eiters.first, res.graph() );
				assert( attached != atom_vd);
				if ( res.atom( attached ).element() == element::H ) {
					hydrogen = attached;
				} else {
					VD putative_sink( find_electron_sink( res, attached ) );
					if ( putative_sink != MutableResidueType::null_vertex && putative_sink != atom_vd ) {
						electron_sink = putative_sink;
					}
				}
			} // For all bonds to atom of interest
			if ( hydrogen != MutableResidueType::null_vertex && electron_sink != MutableResidueType::null_vertex ) {
				deprotonate.push_back( atom_vd );
			}

		} else if ( atom.element() == element::N ) {
			// Heuristic - add a proton to nitrogen if it has three single bonds, and none of the attached atoms are conjugatable.
			if ( res.nbonds( atom_vd ) != 3 ) { continue; } // Only deal with nitrogens with three bonds.
			bool valid( true );
			for ( OutEdgeIterPair eiters( boost::out_edges(atom_vd, res.graph()) ); eiters.first != eiters.second; ++eiters.first ) {
				if ( res.bond( *eiters.first ).order() != SingleBondOrder ) {
					// Oops - the atom of interest needs to be single bonded only
					valid = false;
					break;
				}
				VD attached = boost::target(*eiters.first, res.graph() );
				assert( attached != atom_vd);
				if ( is_conjugatable( res, attached ) ) {
					valid = false;
					break;
				}
			} // For all bonds to atom of interest
			if ( valid ) {
				protonate.push_back( atom_vd );
			}
		} // Element if/else
	} // For all atoms in residue type

	if ( deprotonate.size() == 0 && protonate.size() ==0 ) {
		TR << "No protonation changes for " << res.name() << std::endl;
		return; // Nothing to do, so do nothing.
	}
	// Update charges, delete hydrogens from all modified atoms
	// Then re-do gasteiger atom typing, and add missing hydrogens to the relevant atoms.

	for ( core::Size ii(1); ii <= deprotonate.size(); ++ii ) {
		VD atom( deprotonate[ii] );
		runtime_assert( res.atom( atom ).formal_charge() == 0 );
		TR << "Removing proton from " << res.atom_name( atom ) << " on " << res.name() << std::endl;
		res.atom( atom ).formal_charge( -1 );
		core::Size n_removed = remove_hydrogens( res, atom );
		if ( n_removed != 1 ) {
			TR.Warning << "Attempted to remove single proton from " << res.atom_name( atom ) << " but " << n_removed << " were removed instead." << std::endl;
		}
	}

	for ( core::Size ii(1); ii <= protonate.size(); ++ii ) {
		VD atom( protonate[ii] );
		runtime_assert( res.atom( atom ).formal_charge() == 0 );
		//TR << "Adding proton to " << res.atom_name( atom ) << std::endl;
		res.atom( atom ).formal_charge( +1 );
		// // If we wanted to adjust the positioning of the hydrogens on the atom, we could remove/add
		// // But in all current use cases the Nitrogen should already be tetrahedral.
		//remove_hydrogens( res, atom ); // They'll be added back later
	}

	core::chemical::gasteiger::assign_gasteiger_atom_types( res, /*keep_existing*/ false );

	// // Currently, we only consider di-valent oxygens or sulfurs for deprotonation -
	// // we don't need to worry about rearragment of other existing hydrogens on the atom
	// for( core::Size ii(1); ii <= deprotonate.size(); ++ii ) {
	//    VD atom( deprotonate[ii] );
	//  add_hydrogens_to_atom( res, atom );
	// }

	for ( core::Size ii(1); ii <= protonate.size(); ++ii ) {
		VD atom( protonate[ii] );
		TR << "Adding proton to " << res.atom_name( atom ) << " on " << res.name() << std::endl;
		core::Size n_added = add_hydrogens_to_atom( res, atom );
		if ( ! n_added ) {
			TR.Warning << "Attempted to add single proton to " << res.atom_name( atom ) << " but " << n_added << " were added instead." << std::endl;
		}
	}

	// Do standard semi-derived data updating.

	rosetta_retype_fullatom(res, false); //need to do this, fails currently

	rosetta_recharge_fullatom(res);
	core::chemical::find_bonds_in_rings( res );
	res.assign_internal_coordinates(); // Also sets atom base. Needs nbr atom assignment
	res.autodetermine_chi_bonds(); // Must be after internal coordinate setup
	set_last_status( SUCCESS );
}

/// @brief Find a doubly-bonded oxygen or sulfur attached to neighbor
VD
Reprotonate::find_electron_sink(MutableResidueType const & res, VD neighbor) const {
	// Iterate through all bonds, and see if there are any double bonded oxygens/sulfurs
	for ( OutEdgeIterPair eiters( boost::out_edges(neighbor, res.graph()) ); eiters.first != eiters.second; ++eiters.first ) {
		if ( res.bond( *eiters.first ).order() != DoubleBondOrder ) { continue; }
		VD attached = boost::target(*eiters.first, res.graph() );
		if ( res.atom( attached ).element() == element::O || res.atom( attached ).element() == element::S ) {
			return attached;
		}
	}
	return MutableResidueType::null_vertex;
}

/// @brief Is the provided atom able to participate in conjugation?
bool
Reprotonate::is_conjugatable(MutableResidueType const & res, VD atom) const {
	// Iterate through all bonds, and see if there are any double, triple, atomatic, etc bonds
	for ( OutEdgeIterPair eiters( boost::out_edges(atom, res.graph()) ); eiters.first != eiters.second; ++eiters.first ) {
		Bond const & bond( res.bond( *eiters.first ) );
		core::chemical::BondOrder order( bond.order() );
		if ( order == DoubleBondOrder || order == TripleBondOrder ) { return true; }
		if ( bond.aromaticity() == IsAromaticBond ) { return true; }
	}
	return false;
}

core::Size
Reprotonate::remove_hydrogens(MutableResidueType & restype, VD atom) {
	utility::vector1< VD > to_remove;
	for ( AdjacentIterPair iter( restype.bonded_neighbor_iterators(atom) ); iter.first != iter.second; ++iter.first ) {
		if ( restype.atom( *iter.first ).element() == element::H ) {
			to_remove.push_back( *iter.first );
		}
	}
	for ( core::Size ii(1); ii <= to_remove.size(); ++ii ) {
		mapping_.invalidate( to_remove[ii] );
		restype.delete_atom( to_remove[ii] );
	}
	return to_remove.size();
}

std::string
Reprotonate::class_name() {
	return "Reprotonate";
}

}
}
}


