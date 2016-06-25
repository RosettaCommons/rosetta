// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/normalmode/NormalModeRelaxMover.hh
/// @brief   initialization for NormalMode
/// @details
/// @author  Hahnbeom Park

#ifndef INCLUDED_protocols_normalmode_NormalModeRelaxMover_hh
#define INCLUDED_protocols_normalmode_NormalModeRelaxMover_hh

// Unit headers
#include <protocols/normalmode/NormalMode.hh>
#include <protocols/normalmode/NormalModeRelaxMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/MinimizerOptions.hh>

// Project headers
#include <core/pose/Pose.hh>

//// C++ headers
#include <string>

#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace normalmode {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

class NormalModeRelaxMover : public protocols::moves::Mover
{
public:
	NormalModeRelaxMover();

	NormalModeRelaxMover( core::scoring::ScoreFunctionCOP sfxn,
		bool const cartesian );

	NormalModeRelaxMover( core::scoring::ScoreFunctionCOP sfxn,
		bool cartesian,
		core::kinematics::MoveMapCOP mm,
		std::string const relaxmode,
		Real const distcut
	);

	~NormalModeRelaxMover();

	virtual
	protocols::moves::MoverOP
	clone() const { return protocols::moves::MoverOP( new NormalModeRelaxMover(*this) ); }

	std::string get_name() const; //{ return NormalModeMinimizerCreator::mover_name(); }

	void apply( core::pose::Pose & pose );

	void apply_on_pose( core::pose::Pose & pose );

	/// Virtual functions from Mover
	virtual protocols::moves::MoverOP fresh_instance() const { return clone(); };

	virtual
	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & pose
	);

	void set_default();

	/// Setters
	virtual
	void set_movemap(
		core::pose::Pose const & pose,
		core::kinematics::MoveMapCOP movemap );

	void set_harmonic_constants( Real const k_uniform );

	void set_harmonic_constants( Real const k_connected,
		Real const k_segment,
		Real const k_long );

	void invert_direction(){ direction_ *= -1; }

	void set_extrapolate_scale( Real const scale ){ moving_distance_ = scale; }

	void set_cst_sdev( Real const cst_sdev ){ cst_sdev_ = cst_sdev; }

	// Use this to modify minimizer; for example to run AtomTreeMin in CartRelaxMover, set as false
	void set_cartesian_minimize( bool const value ){ cartesian_minimize_ = value; }

	void set_minoption( optimization::MinimizerOptionsCOP minoption ){ minoption_ = minoption->clone(); }

	void set_relaxmode( std::string const mode ){ relaxmode_=mode; }

	// Mode setup;
	// user-defined mode setup
	void set_mode( utility::vector1< Size > const mode_using,
		utility::vector1< Real > const mode_scales );

	void set_mode( Size const i_mode );

	void set_random_mode( std::string const select_option = "probabilistic",
		Real const importance_portion = 1.0 );

	void refresh_normalmode() { refresh_normalmode_ = true; }

	// Accessors
	protocols::normalmode::NormalMode const NM() const { return NM_; }
	Real cst_sdev() const { return cst_sdev_; }
	Real get_dynamic_scale() const { return scale_dynamic_; }
	std::string relaxmode() const { return relaxmode_; }

private:
	void
	gen_coord_constraint( pose::Pose &pose,
		utility::vector1< Vector > const &excrd ) const;

	void set_default_minoption();

	utility::vector1< Vector >
	extrapolate_mode_on_crd( pose::Pose const &pose ) const;

	pose::Pose
	extrapolate_mode_on_pose( pose::Pose const &pose ) const;

	Real get_RMSD( utility::vector1< Vector > const excrd,
		pose::Pose const &pose ) const;


private:
	Size nmodes_;
	bool mix_modes_;
	bool cartesian_;
	Real pertscale_;
	bool randomselect_;
	Real selection_kT_;
	bool centroid_;
	Size nsample_;
	Real moving_distance_;
	bool refresh_normalmode_;
	Real direction_;
	Real cst_sdev_;
	bool cartesian_minimize_;
	bool dump_silent_;
	std::string outsilent_;

	std::string relaxmode_;
	optimization::MinimizerOptionsCOP minoption_;

	utility::vector1< Size > mode_using_;
	utility::vector1< Real > mode_scale_;

	core::kinematics::MoveMapOP mm_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::scoring::ScoreFunctionOP sfxn_cen_;
	protocols::normalmode::NormalMode NM_;

	// dynamic ones
	mutable Real scale_dynamic_;

	// torsion stuffs
	mutable utility::vector1< Real > dtor_;

}; // end NormalModeRelaxMover

} // normalmode
} // protocols

#endif
