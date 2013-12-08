// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/screener/ChainBreakScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_screener_ChainBreakScreener_HH
#define INCLUDED_protocols_swa_rna_screener_ChainBreakScreener_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/rna/screener/ChainBreakScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {
namespace screener {

	class ChainBreakScreener: public utility::pointer::ReferenceCount {

	public:

		//constructor
		ChainBreakScreener( pose::Pose const & pose, Size const five_prime_res );


		//destructor
		~ChainBreakScreener();

	public:

		void
		set_reinitialize_CCD_torsions( bool const & setting ){ reinitialize_CCD_torsions_ = setting; };

		void
		add_harmonic_chain_break_constraint( Size const five_prime_res );

		void
		copy_CCD_torsions( pose::Pose & pose ) const;

		void
		copy_CCD_torsions_general( pose::Pose & pose, Size const five_prime_res, Size const three_prime_res ) const;

		bool
		check_loop_closed( pose::Pose const & pose );

		bool
		chain_break_screening_general( pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & chainbreak_scorefxn, Size const five_ );

		bool
		chain_break_screening( pose::Pose & chain_break_screening_pose, core::scoring::ScoreFunctionOP const & chainbreak_scorefxn );

		bool
		check_screen();

		bool
		check_screen( pose::Pose & pose );

		pose::Pose & pose(){ return chain_break_screening_pose_; }

	private:

		pose::Pose chain_break_screening_pose_;
		Size const five_prime_res_;
		bool reinitialize_CCD_torsions_;
		bool verbose_;

		core::scoring::ScoreFunctionOP chain_break_scorefxn_;

		StepWiseRNA_CountStruct count_data_;

	};

} //screener
} //rna
} //swa
} //protocols

#endif
