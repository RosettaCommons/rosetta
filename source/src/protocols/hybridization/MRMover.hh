// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Frank DiMaio


#ifndef INCLUDED_protocols_electron_density_MRMover_hh
#define INCLUDED_protocols_electron_density_MRMover_hh

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/Mover.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>


//// C++ headers
#include <string>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace hybridization {

class MRMover : public protocols::moves::Mover {
public:
	MRMover();

	void init();

	void apply( Pose & pose );

	virtual std::string get_name() const { return "MRMover"; }
	protocols::moves::MoverOP clone() const { return new MRMover( *this ); }
	protocols::moves::MoverOP fresh_instance() const { return new MRMover; }

	// density weight
	void set_centroid_density_weight( core::Real newval ) {
		cen1_scorefxn_->set_weight( core::scoring::elec_dens_fast, newval );
		cen2_scorefxn_->set_weight( core::scoring::elec_dens_fast, newval );
	}

	void set_fullatom_density_weight( core::Real newval, bool fast ) {
		if (fast) {
			fa_scorefxn_->set_weight( core::scoring::elec_dens_fast, 10*newval );
			fa_scorefxn_->set_weight( core::scoring::elec_dens_window, 0 );
		} else {
			fa_scorefxn_->set_weight( core::scoring::elec_dens_fast, 0 );
			fa_scorefxn_->set_weight( core::scoring::elec_dens_window, newval );
		}
	}

	// set options
	void set_max_gaplength_to_model( core::Size newval ) { max_gaplength_to_model_=newval; }
	void set_symmdef_file( std::string newval ) { symm_def_file_=newval; }
	void set_censcale( core::Real newval ) { censcale_=newval; }
	void set_disulf( utility::vector1<std::string> newval ) { disulfs_=newval; };

	void set_big_fragments( core::fragment::FragSetOP newval ) { fragments_big_trim_ = fragments_big_ = newval; }
	void set_small_fragments( core::fragment::FragSetOP newval ) { fragments_small_trim_ = fragments_small_ = newval; }

private:
	// remove long gaps from threaded pose
	//    updates loop definitions and fragments
	void trim_target_pose( Pose & query_pose, protocols::loops::Loops &loops , core::Size max_gaplength );

	// repack any missing sidechains after threading (threading does not do this)
	void pack_missing_sidechains( Pose & pose );

private:
	// symmetry
	std::string symm_def_file_;

	// fragments
	core::fragment::FragSetOP fragments_big_, fragments_small_;
	core::fragment::FragSetOP fragments_big_trim_, fragments_small_trim_;

	// scorefunctions
	core::scoring::ScoreFunctionOP cen1_scorefxn_;
	core::scoring::ScoreFunctionOP cen2_scorefxn_;
	core::scoring::ScoreFunctionOP fa_scorefxn_;

	// forced disulfides
	utility::vector1<std::string> disulfs_;

	// other parameters
	core::Size max_gaplength_to_model_, relax_max_iter_, relax_cycles_;
	core::Real censcale_;
};

typedef utility::pointer::owning_ptr< MRMover > MRMoverOP;

}
}

#endif
