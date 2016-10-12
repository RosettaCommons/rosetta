// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_submotif_SubMotifLibrary_HH
#define INCLUDED_protocols_stepwise_modeler_submotif_SubMotifLibrary_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <set>
#include <map>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace submotif {

typedef utility::vector1< std::string > SubMotifSequenceSet;
typedef std::string PoseTag;
typedef utility::vector1< core::Size > SequenceMapping;

class SubMotifLibrary: public utility::pointer::ReferenceCount {

public:

	//constructor
	SubMotifLibrary( core::chemical::ResidueTypeSetCAP rsd_set,
		bool const include_submotifs_from_jump_library = false,
		bool const use_first_jump_for_submotif = false );

	//destructor
	~SubMotifLibrary();

public:

	core::pose::PoseOP
	create_new_submotif( SequenceMapping const & sequence_mapping,
		PoseTag const & submotif_tag,
		core::pose::Pose const & pose,
		bool const & seed = false ) const;

	utility::vector1< monte_carlo::mover::StepWiseMove >
	get_submotif_moves( core::pose::Pose const & pose ) const;

private:

	void
	initialize();

	void
	initialize_from_directory( std::string const & dir_name );

	void
	initialize_from_jump_library();

	SubMotifSequenceSet
	get_submotif_sequence_set( core::pose::Pose const & pose, bool sort_sequences = true ) const;

	void
	save_pose_as_submotif( core::pose::PoseOP pose, std::string const & tag );

	utility::vector1< SequenceMapping >
	get_matches_for_one_submotif_sequence_set( SubMotifSequenceSet const & submotif_sequence_set,
		core::pose::Pose const & pose,
		bool const use_full_model_info = true ) const;

	void
	get_moves_for_one_submotif( core::pose::PoseCOP submotif_pose, core::pose::Pose const & pose ) const;

	void
	get_matches( utility::vector1< SequenceMapping > & all_matches /* stores matches */,
		SequenceMapping const & matching_residues /* working mapping */,
		std::string const & submotif_sequence,
		utility::vector1< Size > const & submotif_cutpoints,
		std::string const & pose_sequence,
		utility::vector1< Size > const & pose_cutpoints,
		utility::vector1< Size > const & pose_domain_map ) const;

	void
	output_tags() const;

private:

	core::chemical::ResidueTypeSetCAP rsd_set_;
	bool const include_submotifs_from_jump_library_;
	bool const use_first_jump_for_submotif_;

	// central list of 'submotif sets', grouped by chain sequences. So [aa, uu] will not show up twice as [uu, aa].
	std::set< SubMotifSequenceSet > submotif_sequence_sets_;

	// each submotif set could map to a lot of different poses.
	// In fact a submotif set like [g, g] will map into the same g/g mismatch pose in two different orders.
	// Those mappings are all kept separately.
	std::map< SubMotifSequenceSet, utility::vector1< std::pair< PoseTag, SequenceMapping > > > submotif_mappings_by_sequence_set_;
	std::map< PoseTag, core::pose::PoseCOP > submotif_poses_by_tag_;

};

} //submotif
} //monte_carlo
} //stepwise
} //protocols

#endif
