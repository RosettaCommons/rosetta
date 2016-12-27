// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/RECCES_Mover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_RECCES_Mover_HH
#define INCLUDED_protocols_recces_RECCES_Mover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/recces/RECCES_Mover.fwd.hh>
#include <protocols/recces/util.hh> // for Histogram
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh> // should really make this a generic sampler.
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace recces {

class RECCES_Mover: public protocols::moves::Mover {

public:

	//constructor
	RECCES_Mover( utility::vector1< core::Real > const & temps,
		utility::vector1< core::Real > const & st_weights );

	//destructor
	~RECCES_Mover();

public:

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const{ return "RECCES_Mover"; }

	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ) { scorefxn_ = scorefxn ; }

	void set_dump_pdb( bool const & setting ){ dump_pdb_ = setting; }
	bool dump_pdb() const { return dump_pdb_; }

	void set_n_cycle( core::Size const & setting ){ n_cycle_ = setting; }
	core::Size n_cycle() const { return n_cycle_; }

	void set_n_dump( core::Size const & setting ){ n_dump_ = setting; }
	core::Size n_dump() const { return n_dump_; }

	void set_save_scores( bool const & setting ){ save_scores_ = setting; }
	bool save_scores() const { return save_scores_; }

	void set_a_form_range( core::Real const & setting ){ a_form_range_ = setting; }
	core::Real a_form_range() const { return a_form_range_; }

	void set_out_prefix( std::string const & setting ){ out_prefix_ = setting; }
	std::string out_prefix() const { return out_prefix_; }

private:

	void
	initialize( utility::vector1< core::Real > const & orig_weights );

	void
	initialize_sampler( core::pose::Pose const & pose );

	void
	run_sampler( core::pose::Pose & pose );

	void
	update_dump_pdb( core::Size const & n,
		core::pose::Pose const & pose,
		utility::vector1< core::Real > const & scores,
		core::pose::Pose & min_pose,
		core::Real & min_score,
		core::Size & curr_dump ) const;

	void
	save_data_to_disk() const;

	void
	set_sampler_gaussian_stdev( core::Real const & temperature, core::pose::Pose const & pose );

private:

	utility::vector1< core::Real > temps_;
	utility::vector1< core::Real > weights_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	protocols::recces::sampler::rna::MC_RNA_MultiSuiteOP sampler_; // ugh should generalize

	core::Size n_cycle_;
	core::Size n_dump_;
	bool dump_pdb_;
	bool save_scores_;
	std::string out_prefix_;
	core::Real a_form_range_;

	utility::vector1<utility::vector1< float > > data_;
	utility::vector1< Histogram > hist_list_;
};

} //recces
} //protocols

#endif
