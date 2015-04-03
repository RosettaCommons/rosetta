// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/StepWiseScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_StepWiseScreener_HH
#define INCLUDED_protocols_stepwise_StepWiseScreener_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreenerType.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerBase.fwd.hh>
#include <protocols/moves/CompositionMover.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <string>


#ifdef WIN32
	#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>
	#include <protocols/moves/CompositionMover.hh>
#endif


namespace protocols {
namespace stepwise {
namespace screener {

	class StepWiseScreener: public utility::pointer::ReferenceCount {

	public:

		//constructor
		StepWiseScreener();

		//destructor
		~StepWiseScreener();

	public:

		virtual
		void
		get_update( sampler::StepWiseSamplerBaseOP ){}

		virtual
		void
		apply_mover( moves::CompositionMoverOP, Size const, Size const ){}

		virtual
		void
		add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover );

		virtual
		bool
		check_screen(){ return true;} // = 0;

		virtual
		std::string
		name() const = 0;

		virtual
		StepWiseScreenerType
		type() const = 0;

		virtual
		void
		fast_forward( sampler::StepWiseSamplerBaseOP ) {}

		Size const &
		count() const { return count_; }

		void
		increment_count();

		void
		reset(){ count_ = 0;	}

		void
		set_ok_to_increment( bool const setting ){ ok_to_increment_ = setting; }

	private:

		Size count_;
		bool ok_to_increment_;
	};

} //screener
} //stepwise
} //protocols

#endif
