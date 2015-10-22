// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/environment/SequenceAnnotation.hh
/// @brief An object to store, mainipulate and manage the different labels that might be attached to a ProtectedConformation.
///
/// @author Justin Porter

#ifndef INCLUDED_core_environment_SequenceAnnotation_hh
#define INCLUDED_core_environment_SequenceAnnotation_hh

// Unit Headers
#include <core/environment/SequenceAnnotation.fwd.hh>

// Package headers
#include <core/environment/LocalPosition.hh>

// Project headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

// C++ Headers
#include <map>
#include <string>
#include <vector>

// ObjexxFCL Headers

namespace core {
namespace environment {

typedef core::environment::LocalPosition LocalPosition;

class SequenceAnnotation : public utility::pointer::ReferenceCount {
	typedef utility::vector1< std::map< std::string, core::Size> > NumMap;
	typedef std::map< std::string, utility::vector1< core::Size > > SeqLabelMap;
	typedef std::map< std::string, core::Size > JumpLabelMap;
	//typedef std::map< core::Size, std::string > JumpNumMap;

public:
	SequenceAnnotation( core::Size length );

	virtual ~SequenceAnnotation() {}

	void add_seq_label( std::string const&, std::vector< core::Size > const& );

	void add_seq_label( std::string const&, utility::vector1< core::Size > const& );

	void add_jump_label( std::string const&, core::Size );

	void rm_seq_label( std::string const& );

	// Undefined, commenting out to fix PyRosetta build  void rm_jump_label( std::string const& );

	/// @brief append a single residue with the given label to the end of the Annotation
	void append_seq( std::string const& label );

	// Undefined, commenting out to fix PyRosetta build  void shift_notify( core::Size const seqpos, core::Size const shift_size );

	core::Size resolve_seq( LocalPosition const& ) const;

	utility::vector1< core::Size > const&
	resolve_seq( std::string const& label ) const;

	core::Size resolve_jump( std::string const& label ) const;

	core::Size const& length() const { return length_; }

	core::Size length( std::string const& label ) const;

	bool has_seq_label( std::string const& ) const;

private:
	void _add_seq_label( std::string const&, utility::vector1< core::Size > );

	core::Size length_;

	// For residue positions
	NumMap pose_to_local_numbers_;
	SeqLabelMap label_to_pose_numbers_;

	// For Jumps
	//JumpNumMap jump_number_to_label_;
	JumpLabelMap jump_label_to_number_;

}; // end SequenceAnnotation base class

} // environment
} // core

#endif //INCLUDED_core_environment_SequenceAnnotation_hh
