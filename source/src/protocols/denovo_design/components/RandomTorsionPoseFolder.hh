// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/RandomTorsionPoseFolder.hh
/// @brief Folds a pose using random phi/psi torsions
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_components_RandomTorsionPoseFolder_hh
#define INCLUDED_protocols_denovo_design_components_RandomTorsionPoseFolder_hh

// Unit headers
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.fwd.hh>
#include <protocols/denovo_design/components/PoseFolder.hh>

// Protocol headers

// Core headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace denovo_design {
namespace components {

///@brief Folds a pose using random phi/psi torsions
class RandomTorsionPoseFolder : public protocols::denovo_design::components::PoseFolder {
public:
	typedef protocols::denovo_design::components::PoseFolder PoseFolder;
	typedef protocols::denovo_design::components::PoseFolderOP PoseFolderOP;

public:
	RandomTorsionPoseFolder();

	virtual ~RandomTorsionPoseFolder();

	static std::string
	class_name();

	PoseFolderOP
	clone() const;

	virtual void
	apply(
		core::pose::Pose & pose,
		core::select::residue_selector::ResidueSubset const & movable,
		protocols::loops::Loops const & loops ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

private:

};

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_RandomTorsionPoseFolder_hh
