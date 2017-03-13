// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SecondaryStructureSelector.hh
/// @brief  The SecondaryStructureSelector selects residues based on their secondary structure
/// @author Tom Linsky (tlinsky at uw dot edu)

#ifndef INCLUDED_core_select_residue_selector_SecondaryStructureSelector_HH
#define INCLUDED_core_select_residue_selector_SecondaryStructureSelector_HH

// Unit headers
#include <core/select/residue_selector/SecondaryStructureSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

class SecondaryStructureSelector : public ResidueSelector {
public:
	SecondaryStructureSelector();
	SecondaryStructureSelector( std::string const & selected );

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual ResidueSelectorOP clone() const;

	virtual ~SecondaryStructureSelector();

	virtual ResidueSubset apply( core::pose::Pose const & pose ) const;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap );

	virtual std::string get_name() const;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public:
	// mutators

	/// @brief sets number of residues around the given SS elements to select (default=0)
	void set_overlap( core::Size const overlapval );
	/// @brief sets the minimal number of consecutive H residues to keep the secondary structure assignation (default=1)
	void set_minH( core::Size const minHval );
	/// @brief sets the minimal number of consecutive E residues to keep the secondary structure assignation (default=1)
	void set_minE( core::Size const minEval );
	/// @brief sets ss characters to select -- must be set prior to apply()
	void set_selected_ss( std::string const & selected );
	/// @brief if true, one-residue terminal "loops" will be included (default=false)
	void set_include_terminal_loops( bool const inc_term );

	/// @brief Override pose secondary structure. The secondary structure set by this
	///        method will always be used if it is non-empty.
	void
	set_pose_secstruct( std::string const & ss );

	/// @brief If set, dssp will be used to determine secondary structure. Has no effect if pose_secstruct_
	///        is set
	/// @details Determines secondary structure by the following rules:
	///          1. If pose_secstruct_ is user-specified, use that
	///          2. If use_dssp_ is true, run DSSP and use DSSP secstruct
	///          3. If use_dssp_ is false, return pose.secstruct()
	void
	set_use_dssp( bool const use_dssp );

private:
	/// @brief gets the secondary structure to be used by the selector
	/// @param[in] pose  Input pose
	/// @details Determines secondary structure by the following rules:
	///          1. If pose_secstruct_ is user-specified, use that
	///          2. If use_dssp_ is true, run DSSP and use DSSP secstruct
	///          3. If use_dssp_ is false, return pose.secstruct()
	std::string
	get_secstruct( core::pose::Pose const & pose ) const;

	/// @brief fixes the secondary structure according to the expected minimal length of each secondary structure.
	/// @param[in] ss string with the secondary structure definition
	void fix_secstruct_definition( std::string & ss ) const;

	void add_overlap(
		ResidueSubset & matching_ss,
		pose::Pose const & pose,
		std::string const & ss ) const;

	bool check_ss( std::string const & ss ) const;

private:
	std::string pose_secstruct_;
	core::Size overlap_;
	core::Size minH_;
	core::Size minE_;
	bool include_terminal_loops_;
	bool use_dssp_;
	std::set< char > selected_ss_;
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
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_SecondaryStructureSelector )
#endif // SERIALIZATION


#endif
