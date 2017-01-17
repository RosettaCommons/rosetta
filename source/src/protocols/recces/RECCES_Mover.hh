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
#include <protocols/recces/options/RECCES_Options.fwd.hh>
#include <protocols/recces/params/RECCES_Parameters.fwd.hh>
#include <protocols/recces/Histogram.hh>
#include <protocols/recces/sampler/MC_Comb.hh>
#include <protocols/recces/sampler/MC_Loop.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

namespace protocols {
namespace recces {

class RECCES_Mover: public protocols::moves::Mover {

public:

	//constructor
	RECCES_Mover( options::RECCES_OptionsCOP const & options );

	//destructor
	~RECCES_Mover();

public:

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const{ return "RECCES_Mover"; }

	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ) { scorefxn_ = scorefxn ; }

	void set_parameters( params::RECCES_ParametersCOP params ) { params_ = params; }

private:

	void
	initialize();

	void
	run_sampler( core::pose::Pose & pose );

	void
	increment_accepts( Size & n_accept_total,
		std::map< std::string, Size > & num_accepts,
		sampler::MC_LoopCOP loop_sampler ) const;

	void
	save_history(
		core::Size const & curr_counts,
		utility::vector1< float > const & scores,
		core::Size const & temp_id );

	void
	dump_stuff( core::Size const & n,
		core::pose::Pose const & pose,
		utility::vector1< core::Real > const & scores,
		core::pose::Pose & min_pose,
		core::Real & min_score,
		core::Size & curr_dump ) const;


	void
	more_dump_stuff( core::pose::Pose const & pose, Size const & n, core::Real const & temperature ) const;

	void
	final_dump_stuff( core::pose::Pose & pose,
		core::pose::Pose & min_pose ) const;

	void
	save_data_to_disk() const;

	void
	set_sampler_gaussian_stdev( core::Real const & temperature, core::pose::Pose const & pose );

	void
	prepare_output_torsion_ids();

	std::string
	output_pdb_name( std::string const & tag ) const;

	void
	output_torsions( core::pose::Pose const & pose, core::Size const & curr_counts ) const;

private:

	core::scoring::ScoreFunctionCOP scorefxn_;
	options::RECCES_OptionsCOP options_;
	params::RECCES_ParametersCOP params_;
	sampler::MC_CombOP sampler_;

	utility::vector1< core::Real > weights_;
	utility::vector1<utility::vector1< float > > data_;
	utility::vector1< Histogram > hist_list_;
	utility::vector1< core::id::TorsionID > bb_torsion_ids_, chi_torsion_ids_;
	mutable std::ofstream bb_tor_out_, chi_tor_out_;
};

} //recces
} //protocols

#endif
