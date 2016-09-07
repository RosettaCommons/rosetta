// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/comparative_modeling/ThreadingMover.hh
/// @brief
/// @author James Thompson

// libRosetta headers

#ifndef INCLUDED_protocols_comparative_modeling_ThreadingMover_HH
#define INCLUDED_protocols_comparative_modeling_ThreadingMover_HH

#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>

// C++ headers
#include <set>
#include <string>

#include <utility/vector1.hh>

namespace protocols {
namespace comparative_modeling {

class ThreadingMover : public protocols::moves::Mover {

public:
	/// @brief Empty constructor
	ThreadingMover();

	/// @brief  Copy constructor
	ThreadingMover(ThreadingMover const & object_to_copy);

	/// @brief align describes the association between the query and template
	/// sequences, template_pose is the conformation from which to build a
	/// threading model.
	ThreadingMover(
		core::sequence::SequenceAlignment const & align,
		core::pose::Pose const & template_pose
	);

	~ThreadingMover() override = default;

	/// @brief Returns the index of the query sequence in SequenceAlignment
	/// object.
	core::Size query_index() const;

	/// @brief Returns the index of the template sequence in SequenceAlignment
	/// object.
	core::Size template_index() const;

	/// @brief Returns the SequenceAlignment object used in threading.
	core::sequence::SequenceAlignment alignment();

	/// @brief Sets the index of the query sequence in SequenceAlignment object.
	void query_index( core::Size new_index );

	/// @brief Sets the index of the template sequence in SequenceAlignment
	/// object.
	void template_index( core::Size new_index );

	/// @brief Sets the SequenceAlignment associated with this object.
	void alignment( core::sequence::SequenceAlignment new_align );

	void template_pose( core::pose::Pose template_pose );

	//boolean setters
	void build_loops( bool setting );

	void randomize_loop_coords( bool setting );

	void repack_query( bool setting );

	//boolean getters
	bool build_loops() const;

	bool repack_query() const;

	bool randomize_loop_coords();

	void min_loop_size( core::Size const new_size );

	Size min_loop_size() const;

	utility::vector1< core::fragment::FragSetOP > frag_libs() const;

	void frag_libs(
		utility::vector1< core::fragment::FragSetOP > new_libs
	);

	/// @brief Returns the SequenceMapping between query and template.
	core::id::SequenceMapping get_qt_mapping(
		core::pose::Pose const & query_pose
	) const;

	/// @brief Threads the given Pose onto the template_pose with the
	/// SequenceAlignment provided.
	void apply( core::pose::Pose & query_pose ) override;

	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;

private: // methods


	bool atoms_are_equiv( std::string const & a1, std::string const & a2 );

	/// @brief add pair of equivalent atoms to atom_equiv_ table.
	void add_equiv_atoms( std::string const & a1, std::string const & a2 );

	/// @brief initialize table of equivalent atoms between all residues.
	void init_atom_equiv();

private: // data members

	core::Size query_index_, template_index_;
	core::pose::Pose template_pose_;
	core::sequence::SequenceAlignment align_;
	bool build_query_loops_;
	bool repack_query_;
	bool randomize_loop_coords_;
	core::Size min_loop_size_;

	//std::map< std::string, utility::vector1< std::string > > atom_equiv_;
	std::map< std::string, std::set< std::string > > atom_equiv_;
	utility::vector1< core::fragment::FragSetOP > frag_libs_;
}; // ThreadingMover

} // comparative_modeling
} // protocols

#endif
