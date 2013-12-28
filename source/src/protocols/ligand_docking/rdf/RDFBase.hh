// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/rdf/RDFBase.hh
///
/// @brief A base class that defines an RDF function.  These functions computes one or more RDF values for a pair of atoms
/// @author Sam DeLuca


#ifndef INCLUDED_protocols_ligand_docking_rdf_RDFBase_hh
#define INCLUDED_protocols_ligand_docking_rdf_RDFBase_hh

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/ligand_docking/rdf/RDFBase.fwd.hh>

#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Atom.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <string>
#include <map>

namespace protocols {
namespace ligand_docking {
namespace rdf {

/// @brief a vector of pairs, each pair stores the name of a score produced by the RDF and the value. Used for RDF output
typedef utility::vector1<std::pair<std::string,core::Real> > RDFResultList;
	
/// @brief a simple struct for storing information needed to compute an RDF interaction
struct AtomPairData
{
	
	AtomPairData(core::chemical::AtomType const & protein_atom,core::chemical::AtomType const & ligand_atom	) :
	protein_atom_id(),
	ligand_atom_id(),
	protein_atom_coords(),
	ligand_atom_coords(),
	atom_atom_distance(0.0),
	protein_atom(),
	ligand_atom(),
	protein_atom_charge(0.0),
	ligand_atom_charge(0.0),
	protein_atom_type(protein_atom),
	ligand_atom_type(ligand_atom)
	{
		
	}
	
	core::id::AtomID protein_atom_id;
	core::id::AtomID ligand_atom_id;
	
	core::Vector protein_atom_coords;
	core::Vector ligand_atom_coords;
	
	core::Real atom_atom_distance;
	
	core::conformation::Atom protein_atom;
	core::conformation::Atom ligand_atom;
	
	core::Real protein_atom_charge;
	core::Real ligand_atom_charge;
	
	core::chemical::AtomType protein_atom_type;
	core::chemical::AtomType ligand_atom_type;
};
	
class RDFBase : public utility::pointer::ReferenceCount
{
public:
	
	
	RDFBase(std::string const & name) : name_(name)
	{
		
	}
	virtual ~RDFBase() {};
	
	/// @brief given an AtomPairData object, return a map
	virtual RDFResultList operator()(AtomPairData const &  ) = 0;
	
	/// @brief parse tags
	virtual void parse_my_tag(
		utility::tag::TagCOP ,
		basic::datacache::DataMap & ) = 0;
	
	/// @brief If you have code that needs to be run once per pose (for data caching, etc), put it here.
	virtual void preamble(core::pose::Pose & ) {}
	
	/// @brief add a function name to the list
	void add_function_name(std::string const & name)
	{
		function_names_.push_back(name);
	}
	
	/// @brief get a list of all the function names which will be present in RDFResultList after operator() is called
	utility::vector1<std::string> get_function_names()
	{
		return function_names_;
	}
	
private:
	std::string name_;
	utility::vector1<std::string> function_names_;
};

	
}
}
}

#endif
