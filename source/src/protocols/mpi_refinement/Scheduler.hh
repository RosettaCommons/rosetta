// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mpi_refinement/Scheduler.hh
/// @brief
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_mpi_refinement_Scheduler_hh
#define INCLUDED_protocols_mpi_refinement_Scheduler_hh

#include <core/types.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace mpi_refinement {

struct MethodParams
{
	std::string name;
	std::string movertype; // Should use registered name: md/nm/relax/bbgauss/loophash
	std::string roundtype;
	core::Size nrun;
	core::Size istart; // Starting index
	core::Size irun;   // Current struct index

	// Store no more than below because WU can only handle limited number of args
	core::Real cstw;
	core::Size nperrun;   // for MD, bbgauss
	core::Size score_type;
	core::Size relax_type;
	core::Size rerelax_type;
	core::Size index;
	core::Real fshave1; // shave fraction, before rerelax
	core::Real fshave2; // shave fraction, after rerelax
	core::Real fshave3; // shave fraction, after collecting whole pool
	//bool filter_before_avrg;
};

class Scheduler
{

public:

	Scheduler();
	~Scheduler();

	void
	prepare_search_stage( core::Size const mpi_rank );

	void
	prepare_enrich_stage( protocols::wum::SilentStructStore const &decoys,
		std::string const & scorename );

	void proceed(); // Proceed iter/round

	void clear();

	// Setters
	void set_random( bool const value ) { is_random_ = value; }

	// Accessors
	core::Size n_to_gen() const;
	core::Size nparams() const { return params_.size(); }
	core::Size iter() const { return iter_; }
	core::Size n_to_rerelax() const { return n_to_rerelax_; }
	core::Size n_rerelaxed() const { return n_rerelaxed_; }
	core::Size npick_per_iter() const { return npick_per_iter_; }

	std::string const pick_strategy(){ return pick_strategy_; }
	std::string const pick_objfunction(){ return pick_objfunction_; }

	std::string const methodname( core::Size value ) const {
		if ( methodname_.find( value ) == methodname_.end() ) {
			return "";
		} else {
			return methodname_.find( value )->second;
		}
	}

	void add_rerelaxed( core::Size const value ){
		n_rerelaxed_ += value;
		if ( n_to_rerelax_ > value ) {
			n_to_rerelax_ -= value;
		} else {
			n_to_rerelax_ = 0;
		}
	}
	void add_torerelax( core::Size const value ) { n_to_rerelax_ += value; }

	bool final_iter() const {
		if ( iter_+1 < niter_ ) {
			return true;
		} else {
			return false;
		}
	}

	utility::vector1< core::Size > methods_picked() const { return methods_picked_; }

	MethodParams const get_params() const { return params_[isch_]; }
	MethodParams const get_params( core::Size const imethod ) const;

	std::string roundtype() const { return params_[isch_].roundtype; }
	bool round_expired() const { return isch_ >= params_.size(); }

	//check exclusion rule
	bool is_excluded( utility::vector1< std::string > picked, std::string name2 ) const;

private:

	void set_default();

	void add_fresh_param( std::string const & name );

	void
	read_cmd( std::string const & cmdfile,
		core::Size const mpi_rank,
		core::Size const stage_to_run = 1 );

	utility::vector1< core::Size >
	pick_enrich_methods( protocols::wum::SilentStructStore const &decoys,
		std::string const & scorename ) const;

private:

	bool is_random_;
	utility::vector1< MethodParams > params_;
	utility::vector1< core::Size > methods_picked_;
	std::map< core::Size, std::string > methodname_;
	core::Size isch_;
	core::Size niter_;
	core::Size iter_;
	//core::Size stage_;

	// pick for next round
	core::Size npick_per_iter_;
	std::string pick_strategy_;
	std::string pick_objfunction_;

	// rerelax
	core::Size n_to_rerelax_;
	core::Size n_rerelaxed_;

	// parameters
	core::Size nmethods_enrich_max_;
	core::Real enrich_Zscore_cut_;
	core::Real Zdiff_outstand_;

	// exclusion rules
	utility::vector1< std::pair< std::string, std::string > > exclusions_;

};

}
}

#endif
