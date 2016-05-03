// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/InterGroupByVectorSelector.hh
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_InterGroupInterfaceByVectorSelector_HH
#define INCLUDED_core_select_residue_selector_InterGroupInterfaceByVectorSelector_HH

// Unit headers
#include <core/select/residue_selector/InterGroupInterfaceByVectorSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief  The %InterGroupInterfaceByVectorSelector identifies residues at the interface between
/// two groups of residues using two distance cutoffs: the first, shorter cutoff is based on
/// atom-atom distances, the second, longer cutoff is used if the c-alpha/c-beta vector of the
/// residues in the first group point toward the c-alphas of the other group.  Each group may be
/// specified either explicitly through a std::set of residue indices, through a string of
/// comma-separated values giving either residue indices or PDB ids (e.g. 13B for residue
/// 13 on chain B) or through a ResidueSelector.  Clears the passed ReisueSubset.
class InterGroupInterfaceByVectorSelector : public ResidueSelector {
public:
	// derived from base class
	InterGroupInterfaceByVectorSelector();

	/// @brief Copy constructor
	///
	InterGroupInterfaceByVectorSelector( InterGroupInterfaceByVectorSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	virtual ~InterGroupInterfaceByVectorSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	ResidueSelectorCOP group1_selector() const;
	void group1_selector( ResidueSelectorCOP setting );
	std::set< Size > const & group1_set() const;
	void group1_set( std::set< Size > const & setting );
	std::string const & group1_resstring() const;
	void group1_resstring( std::string const & setting );

	ResidueSelectorCOP group2_selector() const;
	void group2_selector( ResidueSelectorCOP setting );
	std::set< Size > const & group2_set() const;
	void group2_set( std::set< Size > const & setting );
	std::string const & group2_resstring() const;
	void group2_resstring( std::string const & setting );

	Real cb_dist_cut() const;
	Real nearby_atom_cut() const;
	Real vector_angle_cut() const;
	Real vector_dist_cut() const;

	void cb_dist_cut( Real setting );
	void nearby_atom_cut( Real setting );
	void vector_angle_cut( Real setting );
	void vector_dist_cut( Real setting );

	virtual
	std::string
	get_name() const;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	void
	set_from_residue_selector(
		core::pose::Pose const & pose,
		ResidueSelector const & selector,
		std::set< Size > & subset
	) const;

	void
	set_from_residue_list_string(
		core::pose::Pose const & pose,
		std::string const & res_list_string,
		std::set< Size > & subset
	) const;

private:

	ResidueSelectorCOP group1_selector_;
	std::set< core::Size > group1_set_;
	std::string group1_resstring_;

	ResidueSelectorCOP group2_selector_;
	std::set< core::Size > group2_set_;
	std::string group2_resstring_;

	Real cb_dist_cut_;
	Real nearby_atom_cut_;
	Real vector_angle_cut_;
	Real vector_dist_cut_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_InterGroupInterfaceByVectorSelector )
#endif // SERIALIZATION


#endif
