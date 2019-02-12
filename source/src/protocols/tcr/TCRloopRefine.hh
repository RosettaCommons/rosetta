// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/TCRloopRefine.hh
/// @brief class for tcr CDR loop modeling/refinement
/// @author Ragul Gowthaman (ragul@umd.edu)

#ifndef INCLUDED_protocols_tcr_TCRloopRefine_hh
#define INCLUDED_protocols_tcr_TCRloopRefine_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/tcr/TCRloopRefine.fwd.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRmodel.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace tcr {

/// @brief class for modeling/refinement of tcr CDR loops
/// @details calls other loop modeling classes/ util functions to do loop rebuilding and refinement
/// The refinement level was choosen based on user provide optional flags
/// TCRloopRefine takes in TCRseqInfo or the parsed tcr segments for corresponding cdr positions for a/b chains
class TCRloopRefine : public utility::pointer::ReferenceCount {

public:
	/// @brief constructor with arguments
	TCRloopRefine( TCRseqInfo const &tcrseq_info );
	TCRloopRefine( TCRseqInfo::tcrsegs const &atcrseqs,
		TCRseqInfo::tcrsegs const &btcrseqs,
		TCRseqInfo::tcrposi const &acdrposi,
		TCRseqInfo::tcrposi const &bcdrposi );
	/// @brief default constructor
	//TCRloopRefine();

	TCRloopRefineOP clone() const;

	/// @brief default destructor
	~TCRloopRefine() override;

	void apply( core::pose::Pose &pose_in );

	//Setters
	void set_tcr_loop_model(bool tcrloopmodel) {
		tcr_loop_model_ = tcrloopmodel;
	}

	/// @brief return true if loop modeling/refinement is applied to the pose
	bool tcr_loop_model() const { return tcr_loop_model_; }

private:
	/// @brief initizialize from options
	void init();

	/// @brief set default values
	void set_default();

private:

	/// @brief TCRseqInfo holds the parsed tcr segments and corresponding seq numbering
	//tcr::TCRseqInfo ts_info_;

	/// @brief parsed tcr segments for a/b chains
	TCRseqInfo::tcrsegs aseqs_;
	TCRseqInfo::tcrsegs bseqs_;

	/// @brief sequence numbering start/end CDR positions for a/b chains
	TCRseqInfo::tcrposi aposi_;
	TCRseqInfo::tcrposi bposi_;

	/// @brief loop modeling/refinement is applied to the pose
	bool tcr_loop_model_;

	/// @brief remodel with IndependentLoopMover (default remodel mover : "perturb_kic")
	bool ind_remodel_cdr3a_loop_;
	bool ind_remodel_cdr3b_loop_;

	///  @brief remodel cdr3 loops for a/b chains
	bool remodel_cdr3a_loop_;
	bool remodel_cdr3b_loop_;
	bool remodel_cdr3_loops_;

	///  @brief refine cdr3 loops for a/b chains
	bool refine_cdr3a_loop_;
	bool refine_cdr3b_loop_;
	bool refine_cdr3_loops_;

	/// @brief refine all cdr loops (CDR1, CDR2, CDR3)
	bool refine_all_cdr_loops_;

	core::scoring::ScoreFunctionOP scorefxn_;

	std::string outname_prefix_;

}; // class TCRloopRefine

} // tcr
} // protocols

#endif // __ANTIBODY_GRAFTING__

#endif
