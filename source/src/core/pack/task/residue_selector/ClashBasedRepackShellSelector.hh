// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ClashBasedRepackShellSelector.hh
/// @brief  The ClashBasedRepackShellSelector identifies all residues that clash with at least one rotamer of a design position
/// @author Noah Ollikainen (nollikai@gmail.com)
/// @author Roland A. Pache, PhD

#ifndef INCLUDED_core_pack_task_residue_selector_ClashBasedRepackShellSelector_HH
#define INCLUDED_core_pack_task_residue_selector_ClashBasedRepackShellSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/ClashBasedRepackShellSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

class ClashBasedRepackShellSelector : public core::select::residue_selector::ResidueSelector {
public:
	ClashBasedRepackShellSelector();

	/// @brief Copy constructor
	///
	ClashBasedRepackShellSelector( ClashBasedRepackShellSelector const &src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	virtual core::select::residue_selector::ResidueSelectorOP clone() const;

	ClashBasedRepackShellSelector( core::pack::task::PackerTaskOP packer_task, core::scoring::ScoreFunctionOP score_fxn );

	virtual ~ClashBasedRepackShellSelector();

	virtual core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const;
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	virtual
	std::string
	get_name() const;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	// helpers
	bool is_sc_sc_clash(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2
	) const;

	bool is_sc_bb_clash(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2
	) const;

	utility::vector1<core::Size> get_clashing_positions(
		core::pose::Pose const & pose,
		core::conformation::Residue const & rsd1,
		core::Size const resnum
	) const;

	// setters
	void set_packer_task( core::pack::task::PackerTaskOP packer_task );
	void set_score_fxn( core::scoring::ScoreFunctionOP score_fxn );
	void set_bump_overlap_factor( core::Real set_bump_overlap_factor );

	// getters
	core::pack::task::PackerTaskOP get_packer_task() const;
	core::scoring::ScoreFunctionOP get_score_fxn() const;
	core::Real get_bump_overlap_factor() const;

private:

	core::pack::task::PackerTaskOP packer_task_;
	core::scoring::ScoreFunctionOP score_fxn_;
	core::Real bump_overlap_factor_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_ClashBasedRepackShellSelector )
#endif // SERIALIZATION


#endif
