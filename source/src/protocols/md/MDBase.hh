// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/md/MDBase.hh
/// @brief   initialization for MD
/// @details
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
#include <core/optimization/types.hh>

//// C++ headers
#include <string>

#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace md {

struct MDscheduleData
{
	std::string type;
	core::Size nstep;
	core::Real temp0;
};

class MDBase : public protocols::moves::Mover
{
public:

	typedef protocols::moves::Mover parent;

	MDBase();

	~MDBase() override;

	//virtual protocols::moves::MoverOP clone() const;

	void apply( core::pose::Pose & pose ) override = 0;
	std::string get_name() const override = 0;

	// Default options -------------------------------------
	void set_defaults(){
		dt_ = 0.001;
	}

	// Undefinded commenting out to fix PyRostta build. void set_scorefxn( core::scoring::ScoreFunctionOP score );

	virtual
	void set_movemap(
		core::pose::Pose const & pose,
		core::kinematics::MoveMapCOP movemap);

	// Accessors
	core::Real dt() const { return dt_; }
	core::Size n_dof() const { return n_dof_; }
	core::Size n_dof_temp() const { return n_dof_temp_; }
	core::optimization::Multivec mass() const { return mass_; }
	core::Real mass( core::Size iatm ) const { return mass_[iatm]; }
	core::Real cummulative_time() const { return cummulative_time_; }
	core::Size nstep() const { return nstep_; }
	core::Real temp0() const { return temp0_; }

	core::kinematics::MoveMapOP movemap() const { return movemap_; }
	core::scoring::ScoreFunctionOP scorefxn() const { return scorefxn_; }
	core::scoring::ScoreFunctionOP scorefxn_obj() const { return scorefxn_obj_; }

	std::string selectmode() const { return selectmode_; }
	core::pose::Pose pose_minobj() const { return pose_minobj_; }
	core::Real Emin_obj() const { return Emin_obj_; }
	core::Real time_minobj() const { return time_minobj_; }

	core::Size md_report_stepsize() const { return md_report_stepsize_; }
	core::Size md_energy_report_stepsize() const { return md_energy_report_stepsize_; }
	core::Size md_rsr_update_stepsize() const { return md_rsr_update_stepsize_; }
	core::pose::Pose pose0() const { return pose0_; }
	core::Size context_update_step() const { return context_update_step_; }

	core::Size ncyc_premin() const { return ncyc_premin_; }
	core::Size ncyc_postmin() const { return ncyc_postmin_; }
	bool scheduled() const { return scheduled_; }
	utility::vector1< MDscheduleData > mdsch() const { return mdsch_; }
	MDscheduleData mdsch( core::Size istep ) const { return mdsch_[istep]; }
	bool uniform_coord_constrained() const { return uniform_coord_constrained_; }
	core::Real cst_sdev() const { return cst_sdev_; }

	core::Real temperature() const { return temperature_; }
	core::Real kinetic_energy() const { return kinetic_energy_; }
	core::Real potential_energy() const { return potential_energy_; }
	core::optimization::Multivec xyz() const { return xyz_; }
	core::optimization::Multivec vel() const { return vel_; }
	core::optimization::Multivec acc() const { return acc_; }
	core::Real &xyz( core::Size idof ) { return xyz_[idof]; }
	core::Real &vel( core::Size idof ) { return vel_[idof]; }
	core::Real &acc( core::Size idof ) { return acc_[idof]; }
	core::optimization::Multivec const &ref_xyz() const { return ref_xyz_; }
	core::optimization::Multivec const &prv_eqxyz() const { return prv_eqxyz_; }
	core::Real Kappa() const { return Kappa_; }
	core::Real Gamma() const { return Gamma_; }

	bool report_scorecomp() const { return report_scorecomp_; }
	bool report_as_silent() const { return report_as_silent_; }
	bool trj_score_only() const { return trj_score_only_; }
	std::string silentname() const { return silentname_; }
	bool store_trj() const { return store_trj_; }
	utility::vector1< core::optimization::Multivec > trj() const { return trj_; }
	core::optimization::Multivec trj( core::Size itrj ) const { return trj_[itrj]; }
	utility::vector1< core::optimization::Multivec > trj_scratch() const { return trj_scratch_; }
	bool write_dynamic_rsr() const { return write_dynamic_rsr_; }
	std::string rsrfilename() const { return rsrfilename_; }

	// Setters
	void set_dt( core::Real setting ) { dt_ = setting; }
	void set_n_dof( core::Size setting ) { n_dof_ = setting; }
	void set_n_dof_temp( core::Size setting ) { n_dof_temp_ = setting; }
	void set_mass( core::Size iatm, core::Real setting ) { mass_[iatm] = setting; }
	void set_cummulative_time( core::Real setting ) { cummulative_time_ = setting; }
	void set_temp0( core::Real setting ) { temp0_ = setting; }
	void set_nstep( core::Size setting ) { nstep_ = setting; }

	void set_scorefxn( core::scoring::ScoreFunctionCOP setting ) { scorefxn_ = setting->clone(); }
	void set_scorefxn_obj( core::scoring::ScoreFunctionCOP setting ) { scorefxn_obj_ = setting->clone(); }
	void set_scorefxn( core::scoring::ScoreFunction const &setting ) { scorefxn_ = setting.clone(); }
	void set_scorefxn_obj( core::scoring::ScoreFunction const &setting ) { scorefxn_obj_ = setting.clone(); }
	//void set_scorefxn( core::scoring::ScoreFunctionOP setting ) { scorefxn_ = setting; }
	//void set_scorefxn_obj( core::scoring::ScoreFunctionOP setting ) { scorefxn_obj_ = setting; }

	void set_selectmode( std::string setting ) { selectmode_ = setting; }
	void set_pose_minobj( core::pose::Pose setting ) { pose_minobj_ = setting; }
	void set_Emin_obj( core::Real setting ) { Emin_obj_ = setting; }
	void set_time_minobj( core::Real setting ) { time_minobj_ = setting; }

	void set_md_report_stepsize( core::Size setting ) { md_report_stepsize_ = setting; }
	void set_md_energy_report_stepsize( core::Size setting ) { md_energy_report_stepsize_ = setting; }
	void set_md_rsr_update_stepsize( core::Size setting ) { md_rsr_update_stepsize_ = setting; }
	void set_pose0( core::pose::Pose setting ) { pose0_ = setting; }
	void set_context_update_step( core::Size setting ) { context_update_step_ = setting; }

	void set_ncyc_premin( core::Size setting ) { ncyc_premin_ = setting; }
	void set_ncyc_postmin( core::Size setting ) { ncyc_postmin_ = setting; }
	void set_scheduled( bool setting ) { scheduled_ = setting; }
	void set_mdsch( utility::vector1< MDscheduleData > setting ) { mdsch_ = setting; }
	void set_uniform_coord_constrained( bool setting ) { uniform_coord_constrained_ = setting; }
	void set_cst_sdev( core::Real setting ) { cst_sdev_ = setting; }

	void set_temperature( core::Real setting ) { temperature_ = setting; }
	void set_kinetic_energy( core::Real setting ) { kinetic_energy_ = setting; }
	void set_potential_energy( core::Real setting ) { potential_energy_ = setting; }
	void set_xyz( core::optimization::Multivec setting ) { xyz_ = setting; }
	void set_vel( core::optimization::Multivec setting ) { vel_ = setting; }
	void set_acc( core::optimization::Multivec setting ) { acc_ = setting; }
	void set_ref_xyz( core::optimization::Multivec setting ) { ref_xyz_ = setting; }
	void set_prv_eqxyz( core::optimization::Multivec setting ) { prv_eqxyz_ = setting; }
	void set_Kappa( core::Real setting ) { Kappa_ = setting; }
	void set_Gamma( core::Real setting ) { Gamma_ = setting; }

	void set_report_scorecomp( bool setting ) { report_scorecomp_ = setting; }
	void set_report_as_silent( bool setting ) { report_as_silent_ = setting; }
	void set_trj_score_only( bool setting ) { trj_score_only_ = setting; }
	void set_silentname( std::string setting ) { silentname_ = setting; }
	void set_store_trj( bool setting ) { store_trj_ = setting; }
	void set_trj( utility::vector1< core::optimization::Multivec > setting ) { trj_ = setting; }
	void set_trj_scratch( utility::vector1< core::optimization::Multivec > setting ) { trj_scratch_ = setting; }
	void set_write_dynamic_rsr( bool setting ) { write_dynamic_rsr_ = setting; }
	void set_rsrfilename( std::string setting ) { rsrfilename_ = setting; }

	// others
	void init();

	void set_constraint( core::Real const sdev );

	void add_trj( core::optimization::Multivec xyz ) { trj_.push_back( xyz ); }
	void add_trj_scratch( core::optimization::Multivec xyz ) { trj_scratch_.push_back( xyz ); }
	void renew_trj_scratch(){ trj_scratch_.resize( 0 ); }
	void resize_natm_variables();

	void parse_schfile( std::string const schfile );

	void report_silent( core::pose::Pose &pose,
		core::Real rmsd = -1.0, core::Real gdttm = -1.0, core::Real gdtha = -1.0 );

	void set_write_dynamic_rsr( std::string const filename )
	{
		write_dynamic_rsr_ = true;
		rsrfilename_ = filename;
	}

private:
	// The movemap
	core::kinematics::MoveMapOP movemap_;

	// Fullatom scoring function used
	core::scoring::ScoreFunctionOP scorefxn_;

	// Setup for final model selection
	std::string selectmode_;
	core::scoring::ScoreFunctionOP scorefxn_obj_;
	core::pose::Pose pose_minobj_;
	core::Real Emin_obj_;
	core::Real time_minobj_;

	// Persistent
	core::Size n_dof_;
	core::Size n_dof_temp_;
	core::Real dt_;
	core::optimization::Multivec mass_;
	core::Real cummulative_time_;
	core::Size md_report_stepsize_; // for trajectory
	core::Size md_energy_report_stepsize_; // for energy
	core::Size md_rsr_update_stepsize_; // for adaptive restraint update
	core::pose::Pose pose0_;
	bool report_scorecomp_;
	core::Size context_update_step_;
	core::Size ncyc_premin_, ncyc_postmin_;

	// Schedule
	bool scheduled_;
	utility::vector1< MDscheduleData > mdsch_;

	// Constrained MD
	bool uniform_coord_constrained_;
	core::Real cst_sdev_;

	// Assigned variables
	core::Real temp0_;
	core::Size nstep_;

	// Dynamic variables
	core::Real temperature_;
	core::Real kinetic_energy_;
	core::Real potential_energy_;

	core::optimization::Multivec xyz_;
	core::optimization::Multivec vel_;
	core::optimization::Multivec acc_;

	// Trj
	bool report_as_silent_;
	bool trj_score_only_;
	std::string silentname_;
	bool store_trj_;
	utility::vector1< core::optimization::Multivec > trj_;

	// Adaptive restraints
	utility::vector1< core::optimization::Multivec > trj_scratch_;
	core::optimization::Multivec ref_xyz_, prv_eqxyz_; // Init xyz, avrg over prv trj_scratch_
	core::Real Kappa_, Gamma_;
	bool write_dynamic_rsr_;
	std::string rsrfilename_;
};

}
} // protocols

#endif
