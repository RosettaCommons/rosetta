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
/// @detailed
/// @author  Hahnbeom Park

#ifndef INCLUDED_protocols_normalmode_NormalModeRelaxMover_hh
#define INCLUDED_protocols_normalmode_NormalModeRelaxMover_hh

// Unit headers
#include <protocols/normalmode/NormalMode.hh>

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

using namespace core;

class NormalModeRelaxMover : public protocols::moves::Mover
{
public:

  NormalModeRelaxMover(){}

  ~NormalModeRelaxMover(){}

	/// Virtual functions from Mover
	virtual	protocols::moves::MoverOP	fresh_instance() const { return clone(); };

  virtual
  void set_movemap(
     core::pose::Pose const & pose,
     core::kinematics::MoveMapCOP movemap );

	/// Pure virtual functions inside NMmover
	virtual
	void gen_coord_constraint( pose::Pose &pose,
														 utility::vector1< Vector > const &excrd ) = 0;

	virtual 
	void set_default_minoption() = 0;

	/// Common options
  void set_harmonic_constants( Real const k_uniform );

  void set_harmonic_constants( Real const k_connected,
															 Real const k_segment,
															 Real const k_long );

	void invert_direction(){ direction_ *= -1; }
	
	void set_extrapolate_scale( Real const scale ){ moving_distance_ = scale; }

	void set_cst_sdev( Real const cst_sdev ){ cst_sdev_ = cst_sdev; }

	// Mode setup
	void set_mode( utility::vector1< Size > const mode_using,
								 utility::vector1< Real > const mode_scales );
	
	void set_mode( Size const i_mode );

	// Use this to modify minimizer; for example to run AtomTreeMin in CartRelaxMover, set as false
	void set_cartesian_minimize( bool const value ){ cartesian_minimize_ = value;	}

	void set_random_mode( Size const nmode,
												std::string const select_option = "probabilistic",
												Real const importance_portion = 1.0 );

	void set_minoption( optimization::MinimizerOptionsCOP minoption ){ minoption_ = minoption->clone(); }

	void set_relaxmode( std::string const mode ){ relaxmode_=mode; }

	void refresh_normalmode() { refresh_normalmode_ = true; }

	// Accessors
	protocols::normalmode::NormalMode NM(){ return NM_; }
	Real cst_sdev(){ return cst_sdev_; }
	Real get_dynamic_scale(){ return scale_dynamic_; }
	std::string relaxmode() const { return relaxmode_; }

protected:
	Real direction_;
	Real moving_distance_;
	Real cst_sdev_;
	bool refresh_normalmode_;
	bool cartesian_minimize_;
	//bool mode_changed_;

	optimization::MinimizerOptionsCOP minoption_;

	Real scale_dynamic_;

	std::string relaxmode_;
	utility::vector1< Size > mode_using_;
	utility::vector1< Real > mode_scale_;

	core::kinematics::MoveMapOP mm_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::scoring::ScoreFunctionOP sfxn_cen_;
	protocols::normalmode::NormalMode NM_;

}; // end NormalModeRelaxMover

//// Cartesian
class CartesianNormalModeMover : public protocols::normalmode::NormalModeRelaxMover
{
public:

	CartesianNormalModeMover( core::pose::Pose const &, //pose
														core::scoring::ScoreFunctionCOP sfxn, 
														std::string const relaxmode );

	CartesianNormalModeMover( core::pose::Pose const & pose,
														core::scoring::ScoreFunctionCOP sfxn,
														core::kinematics::MoveMapCOP mm,
														std::string const mode = "CA",	
														Real const distcut = 10.0,
														std::string const relaxmode = "min" );

	~CartesianNormalModeMover();

	virtual
	protocols::moves::MoverOP
	clone() const {	return protocols::moves::MoverOP( new CartesianNormalModeMover(*this) ); }

	virtual 
	void apply( core::pose::Pose & pose );

	virtual
	std::string get_name() const {	return "CartesianNormalModeMover";}

	virtual
	void gen_coord_constraint( pose::Pose &pose,
														 utility::vector1< Vector > const &excrd );

	virtual 
	void set_default_minoption();

private:
	utility::vector1< Vector >
	extrapolate_mode( pose::Pose const &pose );

	Real get_RMSD( utility::vector1< Vector > const excrd,
								 pose::Pose const &pose );

};

//// Torsion
class TorsionNormalModeMover : public protocols::normalmode::NormalModeRelaxMover
{
public:

	TorsionNormalModeMover( core::pose::Pose const &, //pose
													core::scoring::ScoreFunctionCOP sfxn,
													std::string const relaxmode );

	TorsionNormalModeMover( core::pose::Pose const & pose,
													core::scoring::ScoreFunctionCOP sfxn, 
													core::kinematics::MoveMapCOP mm,
													std::string const mode = "CA",	
													Real const distcut = 10.0,
													std::string const relaxmode = "min" );

	~TorsionNormalModeMover();

	virtual
	protocols::moves::MoverOP
	clone() const {	return protocols::moves::MoverOP( new TorsionNormalModeMover(*this) ); }

	virtual 
	void apply( core::pose::Pose & pose );

	virtual
	std::string get_name() const {	return "TorsionNormalModeMover";}

	virtual
	void gen_coord_constraint( pose::Pose &pose,
														 utility::vector1< Vector > const & );

	virtual 
	void set_default_minoption();

private:

	//void set_dynamic_scale( pose::Pose const &pose );
	pose::Pose extrapolate_mode( pose::Pose const &pose );

private:
	utility::vector1< Real > dtor_;

};

} // normalmode
} // protocols

#endif
