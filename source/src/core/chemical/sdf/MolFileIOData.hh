// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/sdf/MolFileIOData.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_sdf_MolFileIOData_hh
#define INCLUDED_core_chemical_sdf_MolFileIOData_hh

#include <core/chemical/sdf/MolFileIOData.fwd.hh>

#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/types.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>        // for AtomTypeSetCOP
#include <core/chemical/ElementSet.fwd.hh>         // for ElementSetCOP
#include <core/chemical/MMAtomTypeSet.fwd.hh>      // for MMAtomTypeSetCOP

#include <utility/VirtualBase.hh>

#include <boost/graph/undirected_graph.hpp>

#include <map>
#include <string>

namespace core {
namespace chemical {
namespace sdf {

typedef boost::undirected_graph<
	MolFileIOAtomOP, // struct with properties of a node
	MolFileIOBondOP // struct with properties of an edge
	/*,ResidueType*/
	> MolFileIOGraph;

typedef MolFileIOGraph::vertex_descriptor mioAD; // MolFileIO Atom Descriptor
typedef MolFileIOGraph::edge_descriptor mioBD; // MolFileIO Bond Descriptor

inline bool
has( MolFileIOGraph const & graph, mioAD vd ) {
	typedef boost::graph_traits<MolFileIOGraph>::vertex_iterator viter;
	std::pair<viter, viter> iters( boost::vertices(graph) );
	return std::find( iters.first, iters.second, vd ) != iters.second;
}
inline bool
has( MolFileIOGraph const & graph, mioBD ed ) {
	typedef boost::graph_traits<MolFileIOGraph>::edge_iterator eiter;
	std::pair<eiter,eiter > iters( boost::edges(graph) );
	return std::find( iters.first, iters.second, ed ) != iters.second;
}

void dump_graph(MolFileIOGraph const & graph);

typedef core::Size AtomIndex;
typedef core::Size BondIndex;


typedef std::map<std::string, std::string> StrStrMap;
//typedef std::map<std::string, core::Real> StrRealMap;

class MolFileIOAtom : public utility::VirtualBase
{
public:
	MolFileIOAtom();
	~MolFileIOAtom() override;

	AtomIndex index() const { return index_; }
	std::string const & name() const { return name_; }
	std::string const & element() const { return element_; }
	//std::string const & rosetta_type() const { return rosetta_type_; }
	Vector const & position() const { return position_; }
	int formal_charge() const { return formal_charge_; }
	core::Real partial_charge() const { return partial_charge_; }

	void index( AtomIndex index ) { index_ = index; }
	void name(std::string const & name ) { name_ = name; }
	void element(std::string const & element ) {
		//dirty dirty hack to change deturium elements to hydrogen elements
		if ( element == "D" || element == "T" ) {
			element_ = "H";
		} else {
			element_ = element;
		}
	}
	//void rosetta_type(std::string const & rosetta_type ) { rosetta_type_ = rosetta_type; }
	void position(Vector const & position ) { position_ = position; }
	void formal_charge( int formal_charge ) { formal_charge_ = formal_charge; }
	void partial_charge( core::Real partial_charge ) { partial_charge_ = partial_charge; }

private:
	// Depending on what other file formats we want to use, it may make sense to convert this to a string "designator" instead.
	AtomIndex index_;
	std::string name_; /// Need to pull from file.
	Vector position_;
	std::string element_;
	//std::string rosetta_type_;
	int formal_charge_;
	core::Real partial_charge_; /// Need to pull from file.
	//StrStrMap atom_string_data_;
	//StrRealMap atom_real_data_;
};

class MolFileIOBond : public utility::VirtualBase
{
public:
	MolFileIOBond();
	~MolFileIOBond() override;

	BondIndex index() const { return index_; }
	core::Size atom1() const { return atom1_; }
	core::Size atom2() const { return atom2_; }
	core::Size sdf_type() const { return sdf_type_; }

	void index( core::Size index ) { index_ = index; }
	void atom1(core::Size atom1) { atom1_ = atom1; }
	void atom2(core::Size atom2) { atom2_ = atom2; }
	void sdf_type(core::Size sdf_type) { sdf_type_ = sdf_type; }

private:
	BondIndex index_;
	core::Size atom1_;
	core::Size atom2_;
	core::Size sdf_type_;
	//core::Size order_;
	//StrStrMap bond_string_data_;
	//StrRealMap bond_real_data_;
};


class MolFileIOMolecule : public utility::VirtualBase
{
public:
	MolFileIOMolecule();
	~MolFileIOMolecule() override;

	std::string const & name() const { return name_; }
	std::string const & name3() const { return name3_; }
	std::string const & name1() const { return name1_; }
	//core::Size nbr() const { return nbr_; }
	//core::Real nbr_radius() const { return nbr_radius_; }

	void name(std::string const & name) { name_ = name; }
	void name3(std::string const & name3) { name3_ = name3; }
	void name1(std::string const & name1) { name1_ = name1; }
	//void nbr(core::Size nbr) { nbr_ = nbr; }
	//void nbr_radius(core::Real nbr_radius) { nbr_radius_ = nbr_radius; }

	/// @brief Retrieve a modifiable atom by index
	MolFileIOAtomOP atom_index( AtomIndex index );

	/// @brief Add an atom (takes possession of object)
	void add_atom( MolFileIOAtomOP atom );
	/// @brief Add a bond (takes possession of object)
	void add_bond( MolFileIOBondOP bond );

	void add_str_str_data( std::string const & key, std::string const & value );
	StrStrMap const & get_str_str_data() const {
		return molecule_string_data_;
	};

	/// @brief Generate data for potentially missing fields.
	void normalize();

	/// @brief Make a ResidueType from this object
	/// @details Not const as it calls normalize() to fix up missing data first.
	MutableResidueTypeOP convert_to_ResidueType(
		chemical::AtomTypeSetCOP atom_types,
		chemical::ElementSetCOP elements,
		chemical::MMAtomTypeSetCOP mm_atom_types
	) {
		std::map< AtomIndex, std::string > index_name_map;
		return convert_to_ResidueType( index_name_map, atom_types, elements, mm_atom_types );
	}

	/// @brief Make a ResidueType from this object
	/// Index_name_map will contain a mapping from this object's AtomIndexes
	/// to the atom names in the returned ResidueType
	/// @details Not const as it calls normalize() to fix up missing data first.
	MutableResidueTypeOP convert_to_ResidueType(
		std::map< AtomIndex, std::string > & index_name_map,
		chemical::AtomTypeSetCOP atom_types,
		chemical::ElementSetCOP elements,
		chemical::MMAtomTypeSetCOP mm_atom_types
	);

private:

	/// @brief From the str/str data in MolFile, additional information
	/// Specifically, it's the overall atom and bond information.
	void set_from_extra_data(MutableResidueType & restype, std::map< mioAD, core::chemical::VD > & restype_from_mio);

	bool index_valid(AtomIndex index, MutableResidueType const & restype, std::map< mioAD, core::chemical::VD > & restype_from_mio);

	/// @brief Parse a multiple value item, say for per-atom or per-bond data
	/// This assumes that the passed vector is pre-initialized to the appropriate size.
	/// Valid formats are a series of space-separated sequential items, or a space-separated
	/// set of items in the form of "(index,value)". Spaces are not valid in either entry.
	void parse_multi( std::istream & estream, utility::vector1< std::string > & parsed, std::string label) const;

	void create_dummy_atom(MutableResidueTypeOP restype, std::string atom_name, core::Vector const & xyz_offset, chemical::ElementSetCOP elements, chemical::MMAtomTypeSetCOP mm_atom_types);

private:
	std::string name_;
	std::string name3_;
	std::string name1_;
	MolFileIOGraph molgraph_;
	std::map< AtomIndex, mioAD > index_atom_map_;
	StrStrMap molecule_string_data_;
	//StrRealMap molecule_real_data_;
};

}
}
}

#endif
