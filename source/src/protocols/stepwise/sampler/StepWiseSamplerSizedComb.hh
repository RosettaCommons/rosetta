// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerSizedComb.hh
/// @brief Aggregate of multiple rotamer samplers for modeler combinatorially.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_stepwise_sampler_StepWiseSamplerSizedComb_HH
#define INCLUDED_stepwise_sampler_StepWiseSamplerSizedComb_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>

// Note from rhiju -- perhaps this should have multiple inheritance from StepWiseSamplerComb
//  and StepWiseSamplerSized -- would prevent copying code.

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerSizedComb : public StepWiseSamplerSized {

public:

	StepWiseSamplerSizedComb();

	StepWiseSamplerSizedComb( StepWiseSamplerSizedOP outer_loop_rotamer, StepWiseSamplerSizedOP inner_loop_rotamer );

	virtual ~StepWiseSamplerSizedComb();

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if random()) rotamer
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose );

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose & pose, Size const id );

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const {
		runtime_assert( is_init() );
		return size_;
	}

	/// @brief Add one more rotamer sampler to this sampler. This one is 'external' to
	///    any previously existing rotamers, i.e., this new only gets incremented when
	///    all previous existing rotamers complete their cycle.
	virtual void add_external_loop_rotamer( StepWiseSamplerSizedOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Add one more rotamer sampler to this sampler. This one is 'internal' to
	///    existing rotamers,  i.e., when this one completes its cycle, the existing
	//     rotamers are incremented.
	virtual void add_internal_loop_rotamer( StepWiseSamplerSizedOP const & rotamer ) {
		rotamer_list_.insert( rotamer_list_.begin(), rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		size_list_.clear();
		id_list_.clear();
		rotamer_list_.clear();
		set_init( false );
	}

	using StepWiseSamplerSized::fast_forward; // fast forward to very end.

	/// @brief Move sampler to end.
	void fast_forward( Size const sampler_number );

	/// @brief Name of the class
	virtual std::string get_name() const { return "StepWiseSamplerSizedComb"; }

	/// @brief Type of class (see enum in SamplerPlusPlusTypes.hh)
	virtual toolbox::SamplerPlusPlusType type() const { return toolbox::SIZED_COMB; }

	/// @brief Set the random modeler state
	virtual void set_random( bool const setting );

	/// @brief output summary of class
	virtual
	void show( std::ostream & out, Size const indent = 0 ) const;

private:

	void update_rotamer_ids();

protected:

	/// @brief Convert input id number to the individual id_list for each
	/// stored rotamer sampler
	utility::vector1<core::Size> id2list( core::Size const id ) const;

	/// @brief Convert id_list for each stored rotamer sampler
	/// stored to the global id number
	core::Size list2id( utility::vector1<core::Size> const & id_list ) const;

private:

	core::Size size_;
	utility::vector1<core::Size> size_list_;

protected: // can be read out by derived classes.

	utility::vector1<core::Size> id_list_;
	utility::vector1<StepWiseSamplerSizedOP> rotamer_list_;

};
} //sampler
} //stepwise
} //protocols

#endif

