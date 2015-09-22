// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/MolFileIOData.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/sdf/MolFileIOData.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/atomtype_support.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/bond_support.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/numbers.hh>

#include <map>

namespace core {
namespace chemical {
namespace sdf {

static THREAD_LOCAL basic::Tracer TR( "core.io.sdf.MolFileIOData" );

void dump_graph(MolFileIOGraph const & graph) {
	MolFileIOGraph::vertex_iterator aiter, aiter_end;
	for ( boost::tie( aiter, aiter_end ) = boost::vertices( graph ); aiter != aiter_end; ++aiter ) {
		debug_assert( has( graph, *aiter ) );
		MolFileIOAtomCOP atom( graph[*aiter] );
		if ( atom ) {
			TR << "Atom " << *aiter << " " << atom->element() << std::endl;
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

MolFileIOAtom::~MolFileIOAtom()
{}

MolFileIOBond::MolFileIOBond() :
	index_( utility::get_undefined_size() ),
	sdf_type_( utility::get_undefined_size() )
	//order_( utility::get_undefined_size() ),
	//bond_string_data_(),
	//bond_real_data_()
{}

MolFileIOBond::~MolFileIOBond()
{}

MolFileIOMolecule::MolFileIOMolecule() :
	name_(""),
	//name3_(""),
	//name1_(""),
	//nbr_( utility::get_undefined_size() ),
	//nbr_radius_( utility::get_undefined_real() ),
	molecule_string_data_()
	//molecule_real_data_()
{}

MolFileIOMolecule::~MolFileIOMolecule()
{}

/// @brief Retrieve a modifiable atom by index
MolFileIOAtomOP
MolFileIOMolecule::atom_index( core::Size index ) {
	if ( index_atom_map_.count(index) == 0 ) {
		return 0; // If not found return default constructor of OP == null pointer
	}
	debug_assert( has(molgraph_, index_atom_map_[index] ) );
	debug_assert( molgraph_[ index_atom_map_[index] ] );
	return molgraph_[ index_atom_map_[index] ];
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

void
MolFileIOMolecule::normalize() {
	// TODO: Does anything need to go here?
}

ResidueTypeOP MolFileIOMolecule::convert_to_ResidueType(
	std::map< AtomIndex, std::string > & index_name_map,
	chemical::AtomTypeSetCOP atom_types,
	chemical::ElementSetCOP elements,
	chemical::MMAtomTypeSetCOP mm_atom_types ) {

	index_name_map.clear();
	debug_assert( elements );

	// Make sure we're up to date first.
	normalize();

	ResidueTypeOP restype( new core::chemical::ResidueType( atom_types, elements, mm_atom_types, NULL ) );

	// Reasonable defaults for:
	// aa_, rotamer_aa_, <properties suite>, variant_types_,

	restype->name( name_ );
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

		restype_atom.name( atom.name() );
		restype_atom.element_type( elements->element( atom.element() ) );
		restype_atom.charge( atom.partial_charge() );
		if ( atom.partial_charge() != 0 ) {
			// Unlike molfile_to_params, any non-zero partial charge turns off recharging.
			uncharged = false;
		}
		restype_atom.formal_charge( atom.formal_charge() );
		restype_atom.ideal_xyz( atom.position() );
		//restype_atom.mm_name("VIRT"); //What is this actually used for? It doesn't look like it gets set.
		restype_atom.mm_atom_type_index( mm_atom_types->atom_type_index("VIRT") ); // We need to do better on this typing.
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

	// Pull extra data from annotations.
	set_from_extra_data(*restype, restype_from_mio);

	// ///////////////////////
	// Now compute missing data -- be careful not to overwrite things that have been set.

	//Rename atoms early to assist possible debugging
	core::chemical::rename_atoms(*restype, /*preserve=*/true);

	// Now that atom names are assigned, update the return-by-reference mapping
	for ( std::map< AtomIndex, mioAD >::const_iterator itr(index_atom_map_.begin()), itr_end(index_atom_map_.end()); itr != itr_end; ++itr ) {
		assert( restype_from_mio.find( itr->second ) != restype_from_mio.end() );
		index_name_map[ itr->first ] = restype->atom_name( restype_from_mio[ itr->second ] );
	}

	core::chemical::rosetta_retype_fullatom(*restype, /*preserve=*/true );
	if ( uncharged ) {
		core::chemical::rosetta_recharge_fullatom(*restype);
	}

	// Now annotate which bonds are in rings (not specified by extra info)
	core::chemical::find_bonds_in_rings( *restype );

	// Nbr atom determination needs ring bond settings
	VD nbr_atom = restype->nbr_vertex();
	if ( nbr_atom == ResidueType::null_vertex || restype->nbr_radius() == 0 ) {
		if ( restype->nbr_radius() != 0 ) {
			// As a radius without a start point is rather pointless.
			TR.Warning << "Warning: neighbor radius specified without neighbor atom specification - ignoring." << std::endl;
		}
		restype->nbr_radius( find_nbr_dist( *restype, nbr_atom ) );
		restype->nbr_atom( nbr_atom );
	}


	// TODO: Can't directly specify internal coordinate tree or chi bond info
	// If that changes this needs to be adjusted so as not to overwrite those settings.

	restype->assign_internal_coordinates(); // Also sets atom base. Needs nbr atom assignment
	restype->autodetermine_chi_bonds(); // Must be after internal coordinate setup

	// To match molfile_to_params, assume type is LIGAND if not otherwise specified.
	if ( restype->properties().get_list_of_properties().size() == 0 ) {  // TODO: I should add a size() method. ~Labonte
		restype->add_property( "LIGAND" );
	}

	restype->finalize();
	return restype;
}

void MolFileIOMolecule::set_from_extra_data(ResidueType & restype, std::map< mioAD, VD > & restype_from_mio) {
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
		} else if ( header == "Rosetta IO_String" ) {
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
			if ( TR.Trace.visible() ) {
				TR.Trace << "Atom Names string : " << entry << std::endl;
			}
			std::string name;
			core::Size ii(1);
			for ( ii=1, estream >> name; estream.good(); ++ii, estream >> name ) {
				if ( ! index_valid(ii, restype, restype_from_mio) ) {
					TR.Warning << "Atom " << ii << " not found when setting atom name '" << name << "' on molecule '" << this->name() << "'" << std::endl;
					continue;
				}
				switch( name.size() ) {
				// Won't be zero, because of >> input
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
					TR.Warning << "Warning: Atom name '" << name << "' for atom " << ii << " on " << this->name() << " is longer than 4 charachters." << std::endl;
					break;
				}
				restype.atom( restype_from_mio[ index_atom_map_[ii] ] ).name( name );
			}
		} else if ( header == "Rosetta AtomTypes" ) {
			std::string type;
			core::Size ii(1);
			for ( ii=1, estream >> type; estream.good(); ++ii, estream >> type ) {
				if ( ! restype.atom_type_set().has_atom( type ) ) {
					TR.Warning << "Atom type " << type << " not recognized for molecule '" << name() << "' - skipping." << std::endl;
					continue;
				}
				if ( ! index_valid(ii, restype, restype_from_mio) ) {
					TR.Warning << "Atom " << ii << " not found when setting Rosetta type '" << type << "' on molecule '" << name() << "'" << std::endl;
					continue;
				}
				restype.set_atom_type( restype_from_mio[ index_atom_map_[ii] ], type );
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

bool MolFileIOMolecule::index_valid(AtomIndex index, ResidueType const & restype, std::map< mioAD, core::chemical::VD > & restype_from_mio) {
	return index_atom_map_.count(index) &&
		restype_from_mio.count(index_atom_map_[index]) &&
		restype.has(restype_from_mio[ index_atom_map_[index] ]);
}

}
}
}
