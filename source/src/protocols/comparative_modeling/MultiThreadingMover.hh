// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/MultiThreadingMover.hh
/// @brief
/// @author James Thompson

// libRosetta headers

#ifndef INCLUDED_protocols_comparative_modeling_MultiThreadingMover_HH
#define INCLUDED_protocols_comparative_modeling_MultiThreadingMover_HH

#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/sequence/SequenceAlignment.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>

// C++ headers
#include <string>
#include <utility/vector1.hh>

namespace protocols {
namespace comparative_modeling {

class MultiThreadingMover : public protocols::moves::Mover {

public:
	typedef core::sequence::SequenceAlignment Alignment;
	typedef utility::vector1< Alignment > Alignments;
	typedef utility::vector1< core::pose::Pose > Poses;

public:

	/// @brief align describes the association between the query and template
	/// sequences, template_pose is the conformation from which to build a
	/// threading model.
	MultiThreadingMover(
		Alignments const & aligns,
		Poses const & template_poses
	);

	~MultiThreadingMover() override = default;

	/// @brief Returns the SequenceAlignments used in threading.
	//Alignments alignment();

	/// @brief Manipulate the SequenceAlignments associated with this object.
	//void alignments( Alignments alns );
	//void add_alignment( Alignment const & aln );
	//void clear_alignments();

	//void template_poses( Poses poses );
	//void add_template_pose( core::pose::Pose const & template_pose );
	//void clear_template_poses();

	void build_loops( bool setting );

	void randomize_loop_coords( bool setting );

	void repack_query( bool setting );

	bool build_loops() const;

	bool repack_query() const;

	bool randomize_loop_coords();

	void min_loop_size( core::Size const new_size );

	Size min_loop_size() const;

	utility::vector1< core::fragment::FragSetOP > frag_libs() const;

	void frag_libs(
		utility::vector1< core::fragment::FragSetOP > new_libs
	);

	/// @brief Threads the given Pose onto the template_pose with the
	/// SequenceAlignment provided.
	void apply( core::pose::Pose & query_pose ) override;

	std::string get_name() const override;

private:
	void check_internals() const;

private: // data members
	Poses template_poses_;
	Alignments alignments_;
	bool build_query_loops_;
	bool repack_query_;
	bool randomize_loop_coords_;
	core::Size min_loop_size_;

	utility::vector1< core::fragment::FragSetOP > frag_libs_;
}; // MultiThreadingMover

} // comparative_modeling
} // protocols

#endif
