// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSizedComb.hh
/// @brief Aggregate of multiple rotamer samplers for sampling combinatorially.
/// @detailed
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerSizedComb_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerSizedComb_HH

#include <protocols/rotamer_sampler/RotamerSizedComb.fwd.hh>

#include <protocols/rotamer_sampler/RotamerSized.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerSizedComb : public RotamerSized {
public:
	RotamerSizedComb();

	RotamerSizedComb( RotamerSizedComb const & other );

	RotamerSizedComb & operator=( RotamerSizedComb const & rhs );

	virtual ~RotamerSizedComb();

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if is_random()) rotamer
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Check if reach the end of rotamer list
	virtual bool not_end() const;

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose );

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose & pose, Size const id );

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const {
		runtime_assert( is_init() );
		return size_;
	}

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_rotamer( RotamerSizedOP const rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		size_list_.clear();
		id_list_.clear();
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerSizedComb"; }
private:
	/// @brief Convert input id number to the individual id_list for each
	/// stored rotamer sampler
	utility::vector1<core::Size> id2list( core::Size const id ) const;

	core::Size size_, id_;

	utility::vector1<core::Size> id_list_, size_list_;

	utility::vector1<RotamerSizedOP> rotamer_list_;
};
}
}

#endif

