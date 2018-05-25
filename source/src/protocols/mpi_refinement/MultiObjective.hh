// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mpi_refinement/MultiObjective.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_mpi_refinement_MultiObjective_hh
#define INCLUDED_protocols_mpi_refinement_MultiObjective_hh

#include <protocols/mpi_refinement/MultiObjective.fwd.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <string>
#include <vector>

#include <utility/pointer/ReferenceCount.hh>

#include <core/kinematics/Jump.hh>

namespace protocols {
namespace mpi_refinement {

class MultiObjective : public utility::pointer::ReferenceCount {

public:
	MultiObjective();
	~MultiObjective() override;

	bool
	update_library_seeds(protocols::wum::SilentStructStore &structs,
		protocols::wum::SilentStructStore &new_structs,
		core::Real const dcut,
		utility::vector1< core::Size > const seeds, // should be 'poolid'
		std::string const prefix_in,
		std::string const objname = "",
		core::Size const maxreplace = 100000
	);

	bool
	update_library_NSGAII(protocols::wum::SilentStructStore &structs,
		protocols::wum::SilentStructStore &new_structs,
		core::Size const nmax,
		bool const update_obj_cut = false
	);

	void
	succeed_substitute_info( core::io::silent::SilentStructOP ss_sub,
		core::io::silent::SilentStructCOP ss_ref,
		bool const reset
	) const;

	void
	filter_similar( protocols::wum::SilentStructStore &structs,
		std::string const measure,
		core::Real const criteria,
		std::string const score_for_priority,
		core::Size const nmax = 0
	);


	void
	add_objective_function_info( core::io::silent::SilentStructOP ss,
		protocols::wum::SilentStructStore & sstore
	) const;

	void
	add_objective_function_info( protocols::wum::SilentStructStore & sstore ) const;

	std::string formatted_objs_values( core::io::silent::SilentStruct const &ss ) const;

	std::string formatted_objs_names() const;

	// Accessors
	core::Size nobjs() const { return fobjnames_.size();}
	std::string fobjnames( core::Size i ) const { return fobjnames_[i]; }
	bool has_score( std::string value ) const { return fobjnames_.contains( value ); }

	core::Real get_fobj( core::io::silent::SilentStruct const &ss, core::Size i ) const
	{ return ss.get_energy( fobjnames_[i] ); }

	core::Real get_fobj( core::io::silent::SilentStruct const &ss, std::string name ) const
	{ return ss.get_energy( name ); }

	core::Real
	get_fobj( core::io::silent::SilentStructCOP pss, core::Size i ) const
	{ return pss->get_energy( fobjnames_[i] ); }

	void set_init_pose( core::pose::Pose inpose ) { init_pose_ = inpose; }
	core::pose::Pose init_pose() const { return init_pose_; }

	void set_iha_cut( core::Real value ) { iha_cut_ = value; }
	core::Real iha_cut() const { return iha_cut_; }

	void set_iha_penalty_slope( core::Real value ) { iha_penalty_slope_ = value; }
	core::Real iha_penalty_slope() const { return iha_penalty_slope_; }

	void set_iha_penalty_mode( std::string value ) { iha_penalty_mode_ = value; }
	std::string iha_penalty_mode() const { return iha_penalty_mode_; }

	void set_nremain_reset( core::Size value ) { nremain_reset_ = value; }
	core::Size nremain_reset() const { return nremain_reset_; }

	void
	calculate_pool_diversity( protocols::wum::SilentStructStore &structs ) const;

	void
	calculate_pool_diversity( protocols::wum::SilentStructStore &structs1,
		protocols::wum::SilentStructStore &structs2 ) const;

	core::scoring::ScoreFunctionOP
	get_scorefxn( core::Size const i ) const { return objsfxnOPs_[i]->clone(); }

private:

	void
	set_defaults();

	bool
	is_dominant( core::io::silent::SilentStructCOP ss1,
		core::io::silent::SilentStructCOP ss2 );

	void
	calculate_structure_diversity( core::io::silent::SilentStructOP ss1,
		protocols::wum::SilentStructStore &structs ) const;


private:
	utility::vector1< core::Real > obj_dominant_cut_;
	utility::vector1< core::Real > obj_cut_increment_;
	utility::vector1< std::string > fobjnames_;
	utility::vector1< core::scoring::ScoreFunctionCOP > objsfxnOPs_;
	core::pose::Pose init_pose_;
	core::Real iha_cut_;
	core::Real iha_penalty_slope_;
	std::string iha_penalty_mode_;
	core::Size nremain_reset_;
};

}
}

#endif // include guard
