// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/Picker.hh
/// @brief Auto-caching fragment picker
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_components_Picker_hh
#define INCLUDED_protocols_denovo_design_components_Picker_hh

// Unit headers
#include <protocols/forge/components/VarLengthBuild.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Package headers

// Core headers

// Basic/Numeric/Utility Headers

// C++ Headers

namespace protocols {
namespace denovo_design {
namespace components {

using namespace protocols::denovo_design::components;

/// @brief class used for picking/caching fragments
class Picker : public protocols::forge::components::VarLengthBuild {
public:
	Picker();
	Picker( core::Size const n_frags );
	virtual ~Picker();

	/// @brief picks and caches fragments for the listed components with size frag_size
	core::fragment::ConstantLengthFragSetOP
	fragments_for_permutation(
		StructureData const & perm,
		utility::vector1< std::string > const & comp_ids,
		core::Size const frag_size );

	/// @brief picks and caches fragments for the listed components with size frag_size
	core::fragment::ConstantLengthFragSetOP
	fragments_for_permutation_take_X_from_pose(
		StructureData const & perm,
		utility::vector1< std::string > const & comp_ids,
		core::Size const frag_size );

	/// @brief pick and cache fragments without considering primary sequence
	core::fragment::ConstantLengthFragSetOP
	pick_and_cache_fragments(
		std::string const & complete_ss,
		utility::vector1< std::string > const & complete_abego,
		core::Size const start_res,
		core::Size const end_res,
		core::Size const frag_length );

	/// @brief pick and cache fragments with a primary sequence in mind
	core::fragment::ConstantLengthFragSetOP
	pick_and_cache_fragments(
		std::string const & complete_aa,
		std::string const & complete_ss,
		utility::vector1< std::string > const & complete_abego,
		core::Size const start_res,
		core::Size const end_res,
		core::Size const frag_length );

protected:
	/// @brief generates a key based on secondary structure to be used in fragcache
	std::string ss_key(
		std::string const & aa,
		std::string const & ss,
		utility::vector1< std::string > const & complete_abego,
		core::Size const start,
		core::Size const end,
		core::Size const fragsize ) const;

	/// @brief pick fragments of a given length, padding when necessary -- DOES NOT CACHE
	/// @param[in] complete_aa The complete amino acid string, typically from a Pose;
	///            can be empty.  If empty, sequence bias is not used to pick fragments.
	/// @param[in] complete_ss The complete secondary structure string, typically from a Pose.
	/// @param[in] complete_abego The complete abego string, typically from a setter, set_abego
	/// @param[in] interval The interval [left, right] to pick fragments from; Pose
	///  numbering (i.e. 1-based indexing).
	/// @param[in] frag_length The desired length of the fragments
	/// @param[in] n_frags The number of fragments to pick per position.
	core::fragment::FrameList
	get_framelist(
		std::string const & complete_aa,
		std::string const & complete_ss,
		utility::vector1< std::string > const & complete_abego,
		core::Size const start_res,
		core::Size const end_res,
		core::Size const frag_length );

private:
	/// @brief number of fragments to pick for each position (default = 200)
	core::Size n_frags_;
	std::map< std::string, core::fragment::FrameList > fragcache_;

};

} // components
} // denovo_design
} // protocols

#endif
