// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/sdf/MolFileIOData.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/chemical/sdf/MolFileIOData.hh>

// Package headers
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/atomtype_support.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/bond_support.hh>
#include <core/chemical/icoor_support.hh>

#include <numeric/conversions.hh>


// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/numbers.hh>
#include <utility/string_util.hh>
#include <utility/stream_util.hh>

// C++ headers
#include <map>

#include <core/chemical/Element.hh> // AUTO IWYU For Element

namespace core {
namespace chemical {
namespace sdf {

static basic::Tracer TR( "core.chemical.sdf.MolFileIOData" );

void dump_graph( MolFileIOGraph const & graph ) {
	MolFileIOGraph::vertex_iterator aiter, aiter_end;
	for ( boost::tie( aiter, aiter_end ) = boost::vertices( graph ); aiter != aiter_end; ++aiter ) {
		debug_assert( has( graph, *aiter ) );
		MolFileIOAtomCOP atom( graph[*aiter] );
		if ( atom ) {
			TR << "Atom " << *aiter << " " << atom->element() << " '" <<  atom->name() << "'" << std::endl;
		} else {
			TR << "Atom with empty data." << std::endl;
		}
	}

	MolFileIOGraph::edge_iterator eiter, eiter_end;
	for ( boost::tie( eiter, eiter_end ) = boost::edges( graph ); eiter != eiter_end; ++eiter ) {
		MolFileIOBondCOP bond( graph[*eiter] );
		if ( bond ) {
			mioAD source( boost::source(*eiter, graph) );
			mioAD target( boost::target(*eiter, graph) );
			TR << "Bond " << *eiter << " Atoms " << source << " " << target << std::endl;
		} else {
			TR << "Bond with empty data." << std::endl;
		}
	}
}

MolFileIOAtom::MolFileIOAtom() :
	index_( utility::get_undefined_size() ),
	position_( utility::get_undefined_real(), utility::get_undefined_real(), utility::get_undefined_real() ),
	element_(""),
	formal_charge_(0),
	partial_charge_(0)
	//atom_string_data_(),
	//atom_real_data_()
{}

MolFileIOAtom::~MolFileIOAtom() = default;

MolFileIOBond::MolFileIOBond() :
	index_( utility::get_undefined_size() ),
	sdf_type_( utility::get_undefined_size() )
	//order_( utility::get_undefined_size() ),
	//bond_string_data_(),
	//bond_real_data_()
{}

MolFileIOBond::~MolFileIOBond() = default;

MolFileIOMolecule::MolFileIOMolecule() :
	name_(""),
	//name3_(""),
	//name1_(""),
	//nbr_( utility::get_undefined_size() ),
	//nbr_radius_( utility::get_undefined_real() ),
	molecule_string_data_()
	//molecule_real_data_()
{}

MolFileIOMolecule::~MolFileIOMolecule() = default;

/// @brief Retrieve a modifiable atom by index
MolFileIOAtomOP
MolFileIOMolecule::atom_index( core::Size index ) {
	if ( index_atom_map_.count(index) == 0 ) {
		return nullptr; // If not found return default constructor of OP == null pointer
	}
	debug_assert( has(molgraph_, index_atom_map_[index] ) );
	debug_assert( molgraph_[ index_atom_map_[index] ] );
	return molgraph_[ index_atom_map_[index] ];
}

MolFileIOAtomOP
MolFileIOMolecule::atom( std::string const & name ) {
	MolFileIOGraph::vertex_iterator aiter, aiter_end;
	for ( boost::tie( aiter, aiter_end ) = boost::vertices( molgraph_ ); aiter != aiter_end; ++aiter ) {
		debug_assert( has( molgraph_, *aiter ) );
		MolFileIOAtomOP atom( molgraph_[*aiter] );
		debug_assert( atom );
		if (atom->name() == name) {
			return atom;
		}
	}
	return nullptr;
}

void
MolFileIOMolecule::add_atom( MolFileIOAtomOP atom ) {
	// This call has both a "side effect" with the add_vertex, as well as the storage effect
	index_atom_map_[ atom->index() ] = molgraph_.add_vertex( atom );
}

void
MolFileIOMolecule::add_bond( MolFileIOBondOP bond ) {
	mioAD atom1( index_atom_map_[ bond->atom1() ] );
	mioAD atom2( index_atom_map_[ bond->atom2() ] );
	molgraph_.add_edge( atom1, atom2, bond );
}

void
MolFileIOMolecule::add_str_str_data( std::string const & key, std::string const & value ) {
	molecule_string_data_[ key ] = value;
}

bool
MolFileIOMolecule::is_polymeric() const {
	// Probably should be a better way of handling this
	if ( molecule_string_data_.count("Rosetta Properties") ) {
		return molecule_string_data_.at("Rosetta Properties").find("POLYMER") != std::string::npos;
	} else {
		return false;
	}
}

bool
MolFileIOMolecule::is_peptidic() const {
	// Probably should be a better way of handling this
	if ( molecule_string_data_.count("Rosetta Properties") ) {
		return molecule_string_data_.at("Rosetta Properties").find("PROTEIN") != std::string::npos;
	} else {
		return false;
	}
}

bool
MolFileIOMolecule::is_nucleic() const {
	// Probably should be a better way of handling this
	if ( molecule_string_data_.count("Rosetta Properties") ) {
		return (molecule_string_data_.at("Rosetta Properties").find("RNA") != std::string::npos ||
			molecule_string_data_.at("Rosetta Properties").find("DNA") != std::string::npos );
	} else {
		return false;
	}
}

AtomIndex
MolFileIOMolecule::get_free_index() const {
	if ( index_atom_map_.empty() ) {
		return 1;
	} else {
		// std::maps have sorted keys
		// Return one past the end, that should be free.
		return index_atom_map_.rbegin()->first + 1;
	}
}

void
MolFileIOMolecule::normalize() {
	// TODO: Does anything need to go here?
}

MutableResidueTypeOP MolFileIOMolecule::convert_to_ResidueType(
	std::map< AtomIndex, std::string > & index_name_map,
	chemical::AtomTypeSetCOP atom_types,
	chemical::ElementSetCOP elements,
	chemical::MMAtomTypeSetCOP mm_atom_types
) {
	TR.Trace << "Converting " << name_ << " to a ResidueType." << std::endl;

	index_name_map.clear();
	debug_assert( elements );

	// Make sure we're up to date first.
	normalize();

	MutableResidueTypeOP restype( new core::chemical::MutableResidueType( atom_types, elements, mm_atom_types, nullptr ) );

	// Reasonable defaults for:
	// aa_, rotamer_aa_, <properties suite>, variant_types_,

	restype->name( name_ );
	restype->base_name( name_ );
	restype->name3( name_.substr(0,3) );
	restype->interchangeability_group( restype->name3() );
	restype->name1( 'Z' );

	bool uncharged = true; // Have partial charges been set?

	std::map< mioAD, VD > restype_from_mio; // A map of atoms in the input graph to the restype graph
	MolFileIOGraph::vertex_iterator aiter, aiter_end;
	for ( boost::tie( aiter, aiter_end ) = boost::vertices( molgraph_ ); aiter != aiter_end; ++aiter ) {
		debug_assert( has( molgraph_, *aiter ) );
		debug_assert( molgraph_[*aiter] );
		MolFileIOAtom const & atom( *(molgraph_[*aiter]) );
		VD vd = restype->add_atom();
		restype_from_mio[ *aiter ] = vd;
		Atom & restype_atom( restype->atom( vd ) );

		if ( atom.name() == "O1P" ) {
			restype->rename_atom( vd, "OP1" );
		} else if ( atom.name() == "O2P" ) {
			restype->rename_atom( vd, "OP2" );
		} else {
			restype->rename_atom( vd, atom.name() );
		}
		restype_atom.element_type( elements->element( atom.element() ) );
		restype_atom.charge( atom.partial_charge() );
		if ( atom.partial_charge() != 0 ) {
			// Unlike molfile_to_params, any non-zero partial charge turns off recharging.
			uncharged = false;
		}
		restype_atom.formal_charge( atom.formal_charge() );
		restype_atom.ideal_xyz( atom.position() );
		restype_atom.mm_name("VIRT"); // We need to do better on this typing.
	}

	MolFileIOGraph::edge_iterator eiter, eiter_end;
	for ( boost::tie( eiter, eiter_end ) = boost::edges( molgraph_ ); eiter != eiter_end; ++eiter ) {
		MolFileIOBond const & bond( *(molgraph_[*eiter]) );
		mioAD source( boost::source(*eiter, molgraph_) );
		mioAD target( boost::target(*eiter, molgraph_) );
		debug_assert( restype_from_mio.count(source) && restype_from_mio.count(target) );
		core::Size bond_type( bond.sdf_type() );
		if ( bond_type > 4 ) { bond_type = 0; };
		restype->add_bond( restype_from_mio[source], restype_from_mio[target], BondName(bond_type) );
	}

	//fix problem with residues that only have one or two atoms. Rosetta residuetypes need to have at least 3 atoms
	// We're hoping that this isn't going to be colinear
	if ( restype->natoms() == 2 ) {
		VD atom1 = restype->all_atoms()[1];
		VD atom2 = restype->all_atoms()[2];
		//add one extra atom
		create_dummy_atom( restype, "DX1",
			restype->atom(atom2).ideal_xyz()-restype->atom(atom1).ideal_xyz()+Vector(1.0, 0.0, 0.0),
			elements, mm_atom_types);
	} else if ( restype->natoms() == 1 ) {
		//add two extra atoms
		create_dummy_atom( restype, "DX1", core::Vector(1.0, 0.0, 0.0), elements, mm_atom_types);
		create_dummy_atom( restype, "DX2", core::Vector(0.0, 1.0, 0.0), elements, mm_atom_types);
	} else if ( restype->natoms() == 0 ) {
		utility_exit_with_message( "Cannot load in ResidueType for entry with no atoms." );
	}

	// Pull extra data from annotations.
	set_from_extra_data( *restype, restype_from_mio );

	// If coordinates aren't set, we can't generate a restype from this molecule
	// Also, if bonds are messed up (too many bonds to one atom) also kill the molecule
	core::Size n_no_coords(0);
	for ( VD atom_vd: restype->all_atoms() ) {
		if ( restype->atom( atom_vd ).ideal_xyz().is_zero() ) {
			if ( ! restype->is_protein() && restype->atom( atom_vd ).name() != "H" ) {
				// We'll reset the position of the (added) hydrogen later, based on coordinates.
				++n_no_coords;
			}
		}
		if ( restype->nbonds( atom_vd ) > 10 ) {
			TR << "Input molecular structure '" << restype->name() << "' has too many bonds (" << restype->nbonds( atom_vd ) << ") on atom " << restype->atom( atom_vd ).name() << std::endl;
			return MutableResidueTypeOP( nullptr );
		}
		// There's a bunch of assumptions in the code about hydrogen only belonging to a single atom (e.g. in atom ordering).
		if ( restype->atom(atom_vd).element() == core::chemical::element::H && restype->nbonds( atom_vd ) != 1 ) {
			TR << "Input molecular structure '" << restype->name() << "' has " << restype->nbonds( atom_vd ) << " bonds to hydrogen " << restype->atom(atom_vd).name() << " - Hydrogens should only have one bond!" << std::endl;
			return nullptr;
		}
	}
	if ( n_no_coords > 1 ) {
		TR << "Input molecular structure '" << restype->name() << "' has too many zero coordinate atoms (" << n_no_coords << "): cannot convert to ResidueType." << std::endl;
		return nullptr;
	}

	// ///////////////////////
	// Now compute missing data -- be careful not to overwrite things that have been set.

	//Rename atoms early to assist possible debugging
	core::chemical::rename_atoms(*restype, /*preserve=*/true);

	// Now that atom names are assigned, update the return-by-reference mapping
	for ( auto const & elem : index_atom_map_ ) {
		debug_assert( restype_from_mio.find( elem.second ) != restype_from_mio.end() );
		index_name_map[ elem.first ] = restype->atom_name( restype_from_mio[ elem.second ] );
	}

	core::chemical::rosetta_retype_fullatom( *restype, /*preserve=*/true );
	if ( uncharged ) {
		core::chemical::rosetta_recharge_fullatom( *restype );
	}

	// Now annotate which bonds are in rings (not specified by extra info)
	core::chemical::find_bonds_in_rings( *restype );

	// Nbr atom determination needs ring bond settings
	VD nbr_atom = restype->nbr_vertex();
	if ( nbr_atom == MutableResidueType::null_vertex || restype->nbr_radius() == 0 ) {
		if ( restype->nbr_radius() != 0 ) {
			// As a radius without a start point is rather pointless.
			TR.Warning << "neighbor radius specified without neighbor atom specification - ignoring." << std::endl;
		}
		restype->nbr_radius( find_nbr_dist( *restype, nbr_atom ) );
		restype->nbr_atom( nbr_atom );
	}

	// To match molfile_to_params, assume type is LIGAND if not otherwise specified.
	if ( restype->properties().get_list_of_properties().size() == 0 ) {  // TODO: I should add a size() method. ~Labonte
		restype->add_property( "LIGAND" );
	}

	determine_polymeric_connections(restype);

	// TODO: Can't directly specify internal coordinate tree or chi bond info from input files
	// If that changes this needs to be adjusted so as not to overwrite those settings.
	if ( ! lower_atom_.empty()) {
		restype->assign_internal_coordinates( restype->atom_vertex( lower_atom_ ) );
	} else {
		restype->assign_internal_coordinates(); // Also sets atom base. Needs nbr atom assignment
	}
	restype->autodetermine_chi_bonds(); // Must be after internal coordinate setup

	handle_polymeric_assignments(restype);

	handle_rotamers(restype);

	// Find possible disulfide atoms.
	for ( VD atom_vd: restype->all_atoms() ) {
		if ( restype->atom( atom_vd ).element_type()->element() == core::chemical::element::S ) {
			// OK, this is an S. Does it have a bonded H?
			for ( auto const nbr : restype->bonded_neighbors( atom_vd ) ) {
				if ( restype->atom( nbr ).element_type()->element() == core::chemical::element::H ) {
					restype->add_property( "SIDECHAIN_THIOL" );
					restype->set_disulfide_atom_name( restype->atom_name( atom_vd ) );
				}
			}
		}
	}

	//TR.Debug << "First sidechain atom: " << restype->name() << " " << restype->atom_name( restype->first_sidechain_atom() ) << std::endl;
	TR.Debug << "Neighbor atom: " << restype->name() << " " << restype->atom_name( restype->nbr_vertex() ) << std::endl;

	//restype->show( TR, true );

	TR.Trace << "Done Converting " << name_ << " to a ResidueType." << std::endl;
	return restype;
}

void
MolFileIOMolecule::determine_polymeric_connections(MutableResidueTypeOP restype) {
	if ( ! is_polymeric() ) { return; }

	if ( lower_atom_.empty() &&  upper_atom_.empty() ) { // If we've set one, trust the annotation for both.
		TR.Warning << "Upper/Lower atoms not set for nominally polymeric residue " << name_ << "; using potentially inaccurate heuristics" << std::endl;
		if ( is_peptidic() ) {
			if ( atom("N") != nullptr ) {
				lower_atom_ = "N";
				TR.Debug << "Setting lower connect point on protein residue to N, based on atom naming" << std::endl;
			}
			if ( atom("C") != nullptr ) {
				upper_atom_ = "C";
				TR.Debug << "Setting upper connect point on protein residue to C, based on atom naming" << std::endl;
			} else if ( atom("P") != nullptr ) {
				upper_atom_ = "P";
				TR.Debug << "Setting upper connect point on protein residue to P, based on atom naming" << std::endl;
			}
		}
		if ( is_nucleic() ) {
			if ( atom("P") != nullptr ) {
				lower_atom_ = "P";
				TR.Debug << "Setting lower connect point on nucleic residue to P, based on atom naming" << std::endl;
			}
			if ( atom("O3'") != nullptr ) {
				upper_atom_ = "O3'";
				TR.Debug << "Setting upper connect point on nucleic residue to O3', based on atom naming" << std::endl;
			} else if ( atom("C3'") != nullptr ) {
				// Try alternatives -- in particular, any atoms
				// bonded to C3' that aren't hydrogens or C4' or C2'
				MolFileIOGraph::edge_iterator eiter, eiter_end;
				utility::vector1< std::string > possible_upper;
				for ( boost::tie( eiter, eiter_end ) = boost::edges( molgraph_ ); eiter != eiter_end; ++eiter ) {
					mioAD source( boost::source(*eiter, molgraph_) );
					mioAD target( boost::target(*eiter, molgraph_) );
					if ( molgraph_[source]->name() == "C3'" && molgraph_[target]->element() != "H" ) {
						possible_upper.push_back( molgraph_[target]->name() );
					} else if ( molgraph_[target]->name() == "C3'" && molgraph_[source]->element() != "H" ) {
						possible_upper.push_back( molgraph_[source]->name() );
					}
				}
				TR.Trace << "Heavy atoms bonded to C3'" << possible_upper << std::endl;
				for ( std::string const & atm: possible_upper ) {
					if ( atm == "C2'" ) { continue; }
					if ( atm == "C4'" ) { continue; }
					upper_atom_ = atm;
					TR.Debug << "Setting upper connect point on nucleic residue to " << upper_atom_ << " as the (first) heavy atom connected to C3'" << std::endl;
				}
			}
		}
	}

	////////////////// Actually assign the connection points

	// Some RTs -- DOC is an example -- is upper terminal and lacks upper
	if ( lower_atom_.empty() || upper_atom_.empty() ) {
		TR.Warning << "Missing connection point for nominally polymeric residue " << name_ << " LOWER: `" << lower_atom_ << "` UPPER: `" << upper_atom_ << "`" << std::endl;
	}
	if ( !lower_atom_.empty() ) {
		restype->set_lower_connect_atom( lower_atom_ );
	}
	if ( !upper_atom_.empty() ) {
		restype->set_upper_connect_atom( upper_atom_ );
	}

}

void
MolFileIOMolecule::handle_polymeric_assignments(MutableResidueTypeOP restype) {
	// AMW: here is where the party starts:
	// 1. allow chi overrides from extra info in molfile? or alter the
	// autodetermine_chi_bonds function to be smarter for peptidic types
	// 2. definitely provide information like UPPER and LOWER as part of the
	// extra info (what gets input as M lines in molfile to params operations)
	// and export them from CIF somehow.

	if ( ! restype->is_polymer() ) { return; }

	utility::vector1< VD > mainchain_vec = mainchain_path(*restype);
	if ( !mainchain_vec.empty() ) {
		restype->set_mainchain_atoms( mainchain_vec );
	}
	if ( mainchain_vec.size() < 3 ) {
		if ( TR.Trace.visible() ) {
			TR.Trace << "Attempting to grab mainchain equivalents from ICOOR tree." << std::endl;
			pretty_print_atomicoor( TR.Trace, *restype);
		}
		// We're either missing an upper/lower (which results in an empty vector), or we're too short of a backbone.
		// Attempt to find usable atoms such that the ICOOR code below works.
		mainchain_vec.resize(6, INVALID_VD); // 1-3 is N-term entries, 4-6 is C-term. Not 100% correct from an absolute perspective, but works in the context of this function.
		if ( restype->lower_connect_id() != 0 ) {
			mainchain_vec[1] = restype->lower_connect_atom();
			mainchain_vec[2] = restype->atom_base( mainchain_vec[1] );
			mainchain_vec[3] = restype->atom_base( mainchain_vec[2] );

			// Catch case where the base is redundant (likely, if we're the root)
			if ( mainchain_vec[3] == mainchain_vec[1] ) {
				for ( VD bbase: restype->bonded_heavyatoms(mainchain_vec[2]) ) {
					if ( bbase != mainchain_vec[1] ) {
						mainchain_vec[3] = bbase;
						break;
					}
				}
			}

			runtime_assert( mainchain_vec[1] != mainchain_vec[2] );
			runtime_assert( mainchain_vec[2] != mainchain_vec[3] );
			runtime_assert( mainchain_vec[1] != mainchain_vec[3] );

			/// FOLLOWING FOR TEMPORARY TESTING PURPOSES ONLY
			if ( restype->is_RNA() || restype->is_DNA() ) {
				if ( restype->atom_name(mainchain_vec[2]) != "O5'" ) {
					TR << "Found " << restype->atom_name(mainchain_vec[2]) << " instead of O5' for mainchain atom 2" << std::endl;
				}
				if ( restype->atom_name(mainchain_vec[3]) != "C5'" ) {
					TR << "Found " << restype->atom_name(mainchain_vec[3]) << " instead of C5' for mainchain atom 3" << std::endl;
				}
			}
		}
		if ( restype->upper_connect_id() != 0 ) {
			mainchain_vec[6] = restype->upper_connect_atom();
			mainchain_vec[5] = restype->atom_base( restype->upper_connect_atom() );
			mainchain_vec[4] = restype->atom_base( mainchain_vec[5] );

			// Catch case where the base is redundant
			if ( mainchain_vec[4] == mainchain_vec[6] ) {
				for ( VD bbase: restype->bonded_heavyatoms(mainchain_vec[5]) ) {
					if ( bbase != mainchain_vec[6] ) {
						mainchain_vec[4] = bbase;
						break;
					}
				}
			}

			runtime_assert( mainchain_vec[6] != mainchain_vec[5] );
			runtime_assert( mainchain_vec[5] != mainchain_vec[4] );
			runtime_assert( mainchain_vec[6] != mainchain_vec[4] );

			/// FOLLOWING FOR TEMPORARY TESTING PURPOSES ONLY
			if ( restype->is_RNA() || restype->is_DNA() ) {
				if ( restype->atom_name(mainchain_vec[5]) != "C3'" ) {
					TR << "Found " << restype->atom_name(mainchain_vec[5]) << " instead of C3' for mainchain atom 5" << std::endl;
				}
				// Some residues are missing C4' (4JA)
				if ( restype->atom_name(mainchain_vec[4]) != "C4'" && restype->atom_name(mainchain_vec[4]) != "C5'" ) {
					TR << "Found " << restype->atom_name(mainchain_vec[3]) << " instead of C4'/C5' for mainchain atom 3" << std::endl;
				}
			}
		}
	}

	if ( restype->is_protein() ) {

		// If the upper connect point is P, it's a phosphonate
		if ( restype->upper_connect_id() != 0 &&
				restype->atom( restype->upper_connect_atom() ).element() == core::chemical::element::P ) {
			restype->add_property( "PHOSPHONATE" );
		}

		if ( mainchain_vec.size() == 3 ) {
			restype->add_property( "ALPHA_AA" );
		} else if ( mainchain_vec.size() == 4 ) {
			restype->add_property( "BETA_AA" );
		}

		using numeric::conversions::radians;
		core::Real chiral = restype->is_d_aa() ? -1.0 : 1.0;

		if ( restype->lower_connect_id() != 0 ) {
			restype->set_icoor( "LOWER", radians(chiral*149.999985), radians(63.800007), 1.328685,
				restype->atom_name( mainchain_vec[1] ), // the lower atom itself
				restype->atom_name( mainchain_vec[2] ),
				restype->atom_name( mainchain_vec[3] ) );
			// Find the pendant hydrogen to the lower connect, if any.
			auto bonded_hydro = restype->bonded_hydrogens( restype->lower_connect_atom() );
			if ( bonded_hydro.size() == 1 ) {
				// This only works if we have a single n-term hydrogen
				std::string const & hname = restype->atom_name( bonded_hydro[1] );
				restype->set_icoor( hname, radians(chiral*-180.000000), radians(60.849998), 1.010000,
					restype->atom_name( mainchain_vec[1] ), // the lower atom itself
					restype->atom_name( mainchain_vec[2] ),
					"LOWER" );
			} else if ( bonded_hydro.size() == 0 ) {
				TR.Trace << "Protein N-term does not contain an attached hydrogen."  << std::endl;
			} else {
				TR.Trace << "Multiple hydrogens on the protein N-term : ";
				for ( auto vd: bonded_hydro ) { TR.Trace << restype->atom_name(vd) << "  "; }
				TR.Trace << std::endl;
			}
		} else {
			TR.Debug << "Labeling protein residue " << restype->name() << " as lower terminus, due to lack of lower connect." << std::endl;
			restype->add_property( "LOWER_TERMINUS" );
		}
		if ( restype->upper_connect_id() != 0 ) {
			restype->set_icoor( "UPPER", radians(chiral*-150.000000), radians(58.300003), 1.328685,
				restype->atom_name( mainchain_vec[mainchain_vec.size()] ), // the upper atom itself
				restype->atom_name( mainchain_vec[mainchain_vec.size() - 1] ),
				restype->atom_name( mainchain_vec[mainchain_vec.size() - 2] ) );
		} else {
			TR.Debug << "Labeling protein residue " << restype->name() << " as upper terminus, due to lack of upper connect." << std::endl;
			restype->add_property( "UPPER_TERMINUS" );
		}

	} else if ( restype->is_RNA() || restype->is_DNA() ) {

		using numeric::conversions::radians;
		// These magic numbers are the standard upper and lower coordinates used across nucleic acid
		// residue types. This is our best attempt to obtain reasonable polymeric-type behavior for
		// residues where that information is not encoded by default.
		// AMW TODO: in theory we could use the icoor that WOULD be learned from the CIF residue from
		// i.e. the OP3 atom for LOWER, but nothing analogous could be done for UPPER so it may not
		// actually be worth it.
		core::Real chiral = restype->is_l_rna() ? -1.0 : 1.0;
		if ( restype->lower_connect_id() != 0 ) {
			restype->set_icoor( "LOWER", radians(chiral*-60.259000), radians(76.024713), 1.607355,
				restype->atom_name( mainchain_vec[1] ), // the lower atom itself
				restype->atom_name( mainchain_vec[2] ),
				restype->atom_name( mainchain_vec[3] ) );

			//if ( has-pendant-oxygens ) {
			//	find-primary-child-of-root
			//	//std::string OP2_name = restype->has( "OP2" ) ? "OP2" :
			//	//	( restype->has( "O2P" ) ? "O2P" :
			//	//	( restype->has( "S2P" ) ? "S2P" : "N4'" ) );
			//	set-icoor-based-on-lower ("O5'" being mainchain_vec[2]
			//	//restype->set_icoor( OP2_name, radians(chiral*-114.600417), radians(72.020306), 1.484470, "P", "O5'", "LOWER" );
			//}
		} else {
			TR.Debug << "Labeling nucleic residue " << restype->name() << " as lower terminus, due to lack of lower connect." << std::endl;
			restype->add_property( "LOWER_TERMINUS" );
		}
		if ( restype->upper_connect_id() != 0 ) {
			restype->set_icoor( "UPPER", radians(chiral*-139.954848), radians(59.821530), 1.607226,
				restype->atom_name( mainchain_vec[mainchain_vec.size()] ), // the upper atom itself
				restype->atom_name( mainchain_vec[mainchain_vec.size() - 1] ),
				restype->atom_name( mainchain_vec[mainchain_vec.size() - 2] ) );
		} else {
			restype->add_property( "UPPER_TERMINUS" );
			TR.Debug << "Labeling nucleic residue " << restype->name() << " as upper terminus, due to lack of lower connect." << std::endl;
		}

	}
}

void MolFileIOMolecule::handle_rotamers(MutableResidueTypeOP restype) {
	// OK, we now have our chis. Assign chi rotamers if necessary.
	// TODO: There's probably a more general way of doing this (Rarey rotamers?)
	if ( restype->is_protein() ) {
		for ( Size ii = 1; ii <= restype->nchi(); ++ii ) {
			if ( restype->is_proton_chi( ii ) ) continue; // already done.
			restype->add_chi_rotamer( ii,  60, 8 );
			restype->add_chi_rotamer( ii, 180, 8 );
			restype->add_chi_rotamer( ii, 300, 8 );
		}
	}
}

void MolFileIOMolecule::create_dummy_atom(MutableResidueTypeOP restype, std::string atom_name, core::Vector const & xyz_offset, chemical::ElementSetCOP elements, chemical::MMAtomTypeSetCOP ){

	VD vd = restype->add_atom(atom_name, "VIRT", "VIRT", 0.0 );
	Atom & restype_atom( restype->atom( vd ) );

	restype_atom.element_type( elements->element( "X" ) ); // Rosetta-specific non-atom
	restype_atom.formal_charge( 0 );

	VD first_atom = restype->all_atoms()[1];

	restype_atom.ideal_xyz( restype->atom( first_atom ).ideal_xyz() + xyz_offset );

	restype->add_bond( vd, first_atom, UnknownBond );
}

void MolFileIOMolecule::set_from_extra_data(MutableResidueType & restype, std::map< mioAD, VD > & restype_from_mio) {
	//////////// Process Data
	TR.Debug << "Loading " << molecule_string_data_.size() << " extra data entries from sdf file." << std::endl;
	for ( StrStrMap::const_iterator iter( molecule_string_data_.begin() ), iter_end(molecule_string_data_.end() );
			iter != iter_end; ++iter ) {
		std::string const & header( iter->first );
		std::string const & entry( iter->second );
		std::stringstream estream( entry + '\n' ); // Extra '\n' needed for good stream processing.

		// Any entries added here should probably also be added to the ones output in mol_writer.cc
		if ( TR.Debug.visible() ) {
			TR.Debug << "Processing extra data from the tag '" << header << "'" << std::endl;
		}

		if ( header == "Rosetta Name" ) {
			std::string in_name;
			estream >> in_name;
			if ( estream ) {
				if ( restype.name3() == restype.name().substr(0,3) ) {
					restype.name3( in_name.substr(0,3) );
				}
				restype.name( in_name );
				// Also reset the MolFileIO name for consistency ???
				name( in_name );
			} else {
				TR.Warning << "Unable to parse Rosetta Name for " << name() << std::endl;
			}
		} else if ( header == "Rosetta IO_String" || header == "Rosetta IO_string" ) {
			/// Format: 01234
			///         LIG X
			if ( entry.size() >= 5 ) {
				std::string name3( entry.substr(0,3) );
				// Avoid reseting interchangeability group if we explicitly set it.
				if ( restype.name3() == restype.interchangeability_group() ) {
					restype.interchangeability_group( name3 );
				}
				restype.name3( name3 );
				restype.name1( entry[4] );
			} else {
				TR.Warning << "Poorly formatted Rosetta IO_string for " << name();
			}
		} else if ( header == "Rosetta AA" ) {
			AA aa;
			estream >> aa;
			if ( estream ) {
				restype.aa( aa );
			} else {
				TR.Warning << "Unable to parse Rosetta AA for " << name() << std::endl;
			}
		} else if ( header == "Rosetta Interchangeability Group" ) {
			std::string group;
			estream >> group;
			if ( estream ) {
				restype.interchangeability_group( group );
			} else {
				TR.Warning << "Unable to parse Rosetta Interchangeability Group for " << name() << std::endl;
			}
		} else if ( header == "Rosetta nbr_atom" ) {
			core::Size nbr;
			estream >> nbr;
			if ( estream ) {
				if ( ! index_valid(nbr, restype, restype_from_mio) ) {
					TR.Warning << "Cannot find atom listed for Rosetta nbr_atom for "<< name() << std::endl;
					continue;
				}
				restype.nbr_atom( restype_from_mio[ index_atom_map_[nbr] ] );
			} else {
				TR.Warning << "Unable to parse Rosetta nbr_atom for " << name() << std::endl;
			}
		} else if ( header == "Rosetta nbr_radius" ) {
			core::Real radius;
			estream >> radius;
			if ( estream ) {
				restype.nbr_radius( radius );
			} else {
				TR.Warning << "Unable to parse Rosetta nbr_radius for " << name() << std::endl;
			}
		} else if ( header == "Atom Names" ) {
			utility::vector1< std::string > names(index_atom_map_.size(), ""); // Here we assume that all the indicies are consecutive.
			parse_multi( estream, names, "<Atom Names>" );
			for ( core::Size ii(1); ii <= names.size(); ++ii ) {
				if ( names[ii] != "" ) {
					if ( ! index_valid(ii, restype, restype_from_mio) ) {
						TR.Warning << "Atom " << ii << " not found when setting atom name '" << names[ii] << "' on molecule '" << this->name() << "'" << std::endl;
						continue;
					}
					std::string name( names[ii] );
					switch( name.size() ) {
					// Won't be zero, because we checked earlier
					case 1 :
						name = " " + name + "  ";
						break;
					case 2 :
						name = " " + name + " ";
						break;
					case 3 :
						name = " " + name;
						break;
					case 4 :
						break; // Don't change anything
					default :
						TR.Warning << "Atom name '" << name << "' for atom " << ii << " on " << this->name() << " is longer than 4 characters." << std::endl;
						break;
					}
					restype.rename_atom( restype_from_mio[ index_atom_map_[ii] ], name );
				}
			}
		} else if ( header == "Rosetta AtomTypes" ) {
			utility::vector1< std::string > types(index_atom_map_.size(), ""); // Here we assume that all the indicies are consecutive.
			parse_multi( estream, types, "<Rosetta AtomTypes>" );
			for ( core::Size ii(1); ii <= types.size(); ++ii ) {
				std::string const & type( types[ii] );
				if ( type != "" ) {
					if ( ! restype.atom_type_set().has_atom( type ) ) {
						TR.Warning << "Atom type " << type << " not recognized for molecule '" << name() << "' - skipping." << std::endl;
					} else if ( ! index_valid(ii, restype, restype_from_mio) ) {
						TR.Warning << "Atom " << ii << " not found when setting Rosetta type '" << type << "' on molecule '" << name() << "'" << std::endl;
					} else {
						restype.set_atom_type( restype_from_mio[ index_atom_map_[ii] ], type );
					}
				}
			}
		} else if ( header == "Rosetta Properties" ) {
			std::string prop;
			for ( estream >> prop; estream; estream >> prop ) {
				restype.add_property( prop );
			}

			// TODO: should be more here - we should be able to round-trip a customized ResidueType.
			// In particular, I'm skipping explicit Icoor, explicit chis, nu atoms, all rotamer library info
			// mainchain atoms, connections, shadowed atoms, variant types, actcoords, backbone, orbitals
			// string and numeric properties, adducts, force_nbr_atom_orient
			// On Atoms: partial charges, mm types, gasteiger types
			// On Bonds: distances, cut bond status, other annotations not redundant with sdf type.
			// -- cutbond is tricky, as we would want to know that prior to the icoor building.
		} else {
			// Default - just stick it in the ResidueType.
			TR << "Molecule file property " << header << " for " << name() << " not recognized: placing in datamap." << std::endl;
			restype.add_string_property( header, entry );
		}
	}
}

bool MolFileIOMolecule::index_valid(AtomIndex index, MutableResidueType const & restype, std::map< mioAD, core::chemical::VD > & restype_from_mio) {
	return index_atom_map_.count(index) &&
		restype_from_mio.count(index_atom_map_[index]) &&
		restype.has(restype_from_mio[ index_atom_map_[index] ]);
}

/// @brief Parse a multiple value item, say for per-atom or per-bond data
/// This assumes that the passed vector is pre-initialized to the appropriate size.
/// Valid formats are a series of space-separated sequential items, or a space-separated
/// set of items in the form of "(index,value)". Spaces are not valid in either entry.
void MolFileIOMolecule::parse_multi( std::istream & estream, utility::vector1< std::string > & parsed, std::string label) const {
	std::string entry;
	core::Size ii(1);

	for ( ii=1, estream >> entry; estream.good(); ++ii, estream >> entry ) {
		core::Size index( ii );
		if ( entry[0] == '(' && entry[ entry.size()-1 ] == ')' ) {
			// In the form "(index,value)"
			core::Size comma_pos( entry.find( ',' ) );
			if ( comma_pos != std::string::npos ) {
				index = utility::string2Size( entry.substr(1,comma_pos-1) );
				entry = entry.substr(comma_pos+1,entry.size()-comma_pos-2); // 2 = the comma and the end parenthesis
			} // If no comma, assume that the value contains parenthesis
		}
		if ( index == 0 || index > parsed.size() ) {
			TR.Warning << "When parsing " << label << " cannot find position " << index << " - ignoring." << std::endl;
		} else {
			parsed[index] = entry;
		}
	}
}

}
}
}
