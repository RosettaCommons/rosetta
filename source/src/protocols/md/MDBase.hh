// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/md/MDBase.hh
/// @brief   initialization for MD
/// @detailed
/// @author  Hahnbeom Park

#ifndef INCLUDED_protocols_md_MDBase_hh
#define INCLUDED_protocols_md_MDBase_hh

// Unit headers
#include <protocols/md/MDBase.fwd.hh>
#include <protocols/md/MDConstants.hh>
#include <protocols/md/thermostat.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

//// C++ headers
#include <string>

#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace md {

using namespace core;

struct MDscheduleData
{
	Size nstep;
	Real temp0;
};

class MDBase : public protocols::moves::Mover 
{
public:

	typedef protocols::moves::Mover parent;
	typedef utility::vector1<Real> Multivec;

	MDBase(){}

	~MDBase(){}

	//virtual protocols::moves::MoverOP clone() const;

	virtual void apply( core::pose::Pose & pose ) = 0;
	virtual std::string get_name() const = 0;

	//virtual void initialize_movemap( core::kinematics::MoveMap & movemap ) = 0;
	// Undefinded commenting out to fix PyRostta build. static void register_options();

	// Default options -------------------------------------
	void set_defaults(){
		dt_ = 0.001;
	}

	// Undefinded commenting out to fix PyRostta build. void set_scorefxn( core::scoring::ScoreFunctionOP score );
	
	core::scoring::ScoreFunctionOP scorefxn(){ return scorefxn_; }
	
	core::scoring::ScoreFunctionCOP scorefxn() const { return scorefxn_; }
	
	virtual
	void set_movemap(
		core::pose::Pose const & pose,
		core::kinematics::MoveMapCOP movemap) = 0;

	virtual
	core::kinematics::MoveMapOP movemap() const = 0;

	// Accessors
	Real dt(){ return dt_; }
	Size n_dof(){ return n_dof_; }
	Real cummulative_time(){ return cummulative_time_; }

	void set_nstep( core::Size const nstep ){ nstep_ = nstep;	}
	void set_temperature( core::Real const temp0 ){ temp0_ = temp0;	}
	void set_reportstep( core::Size const nstep ){ md_report_stepsize_ = nstep; }
	void set_scorefxn_obj( core::scoring::ScoreFunctionCOP sfxn ){ scorefxn_obj_ = sfxn->clone(); }
	void set_selectmode( std::string const selectmode_in ){ selectmode_ = selectmode_in; }
	void set_context_update_step( Size const value ){ context_update_step_ = value; }
	void set_premin( Size const value ){ ncyc_premin_ = value; }
	void set_report_scorecomp( bool const value ){ report_scorecomp_ = value; }

	void parse_schfile( std::string const schfile );

	core::Size nstep(){ return nstep_; }
	core::Real temp0(){ return temp0_; }

	void set_constraint( Real const sdev );
	void cst_on_pose( pose::Pose &pose );

	void set_store_trj( bool const value ){ store_trj_ = value; }
	bool store_trj() const { return store_trj_; }
	
protected:
	// The movemap
	core::kinematics::MoveMapOP movemap_;

	// Fullatom scoring function used
	core::scoring::ScoreFunctionOP scorefxn_;

	// Setup for final model selection
	std::string selectmode_;
	core::scoring::ScoreFunctionOP scorefxn_obj_;
	core::pose::Pose pose_minobj_;
	Real Emin_obj_;
	Real time_minobj_;

	// Persistent
	Size n_dof_;
	Size n_dof_temp_;
	Real dt_;
	Multivec mass_;
	Real cummulative_time_;
	Size md_report_stepsize_;
	core::pose::Pose pose0_;
	bool report_scorecomp_;
	core::Size context_update_step_;
	Size ncyc_premin_, ncyc_postmin_;

	// Schedule
	bool scheduled_;
	utility::vector1< MDscheduleData > mdsch_;

	// Constrained MD
	bool constrained_;
	Real cst_sdev_;

	// Assigned
	Real temp0_;
	Size nstep_;

	// Dynamic
	Real temperature_;
	Real kinetic_energy_;
	Real potential_energy_;

  Multivec xyz_;
  Multivec vel_;
  Multivec acc_;

	bool store_trj_;
	utility::vector1< Multivec > trj_;

};

}
} // protocols

#endif
