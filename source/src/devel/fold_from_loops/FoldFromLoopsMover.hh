// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.


/// Where is the canonical header for this file?
/// Who wrote this file?
/// What is this file for?

#ifndef INCLUDED_devel_fold_from_loops_FoldFromLoopsMover_hh
#define INCLUDED_devel_fold_from_loops_FoldFromLoopsMover_hh


#include <core/pose/Pose.hh> /// Replace Pose objects with PoseOP objects to remove this header.
#include <core/fragment/FragSet.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>
#include <devel/fold_from_loops/FoldFromLoopsMover.fwd.hh>

#include <utility/vector1.hh>



namespace devel {
namespace fold_from_loops {

class FoldFromLoopsMover : public  protocols::moves::Mover {


public:

	//default constructor no args

	FoldFromLoopsMover();

	virtual ~FoldFromLoopsMover();




	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;


	void loops( protocols::loops::Loops const & loops_in ) {

			loops_ = loops_in;

		}

	void loops_pdb( core::pose::Pose const & loops_pdb ) {

			loops_pdb_ = loops_pdb;

		}

	void set_frag_3mers( core::fragment::FragSetOP const & frag_3mers){
		frag_small_ = frag_3mers;
	}

	void set_frag_9mers( core::fragment::FragSetOP const &  frag_9mers){

		frag_large_ = frag_9mers;

	}


	void set_ca_csts( core::scoring::constraints::ConstraintSetOP const & ca_csts){

			ca_csts_ = ca_csts;

	}



private:


	protocols::loops::Loops loops_;
	core::pose::Pose loops_pdb_;
	core::fragment::FragSetOP frag_small_;
	core::fragment::FragSetOP frag_large_;
	core::scoring::constraints::ConstraintSetOP ca_csts_;

};

}
}

#endif
