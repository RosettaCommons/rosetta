// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/ExtraThreadingMover.hh
/// @brief class for copying extra stuff (ligands,dna,peptides,etc) from a template Pose
/// @author James Thompson
/// @author TJ Brunette

// libRosetta headers

#ifndef INCLUDED_protocols_comparative_modeling_ExtraThreadingMover_HH
#define INCLUDED_protocols_comparative_modeling_ExtraThreadingMover_HH

#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

class ExtraThreadingMover : public protocols::moves::Mover {

public:
	/// @brief align describes the association between the query and template
	/// sequences, template_pose is the conformation from which to build a
	/// threading model.
	ExtraThreadingMover(
		core::sequence::SequenceAlignment const & align,
		core::pose::Pose const & template_pose,
		utility::vector1< core::Size > const & residue_selection
	);

	virtual ~ExtraThreadingMover() {}

	/// @brief Returns the index of the query sequence in SequenceAlignment
	/// object.
	core::Size query_index() const;

	/// @brief Returns the index of the template sequence in SequenceAlignment
	/// object.
	core::Size template_index() const;

	/// @brief Returns the SequenceAlignment object used in threading.
	core::sequence::SequenceAlignment alignment();

	/// @brief Returns the SequenceMapping between query and template.
	core::id::SequenceMapping get_qt_mapping(
		core::pose::Pose const & query_pose
	) const;

	/// @brief Threads the given Pose onto the template_pose with the
	/// SequenceAlignment provided.
	virtual void apply( core::pose::Pose & query_pose );

	virtual std::string get_name() const;

private: // data members
	core::Size query_index_, template_index_;
	core::pose::Pose template_pose_;
	core::sequence::SequenceAlignment align_;
	utility::vector1< core::Size > residue_selection_;
}; // ExtraThreadingMover

} // comparative_modeling
} // protocols

#endif
