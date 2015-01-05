// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/ReplicaExchangeMC.hh
/// @brief  implementing a Recplica Exchange Monte Carlo Mover
/// @author Yuan Liu (wendao@u.washington.edu)

#ifndef INCLUDED_protocols_moves_ReplicaExchangeMC_hh
#define INCLUDED_protocols_moves_ReplicaExchangeMC_hh

#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/ReplicaExchangeMC.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

class ReplicaExchangeMC : public  MonteCarlo
{
public:
	typedef MonteCarlo Parent;
	typedef core::Size Size;

	ReplicaExchangeMC(
		Pose const & init_pose, // PoseCOP init_pose,
		ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
		utility::vector1<core::Real> const &tlist,
		core::Size nint
	);

	ReplicaExchangeMC(
		ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
		utility::vector1<core::Real> const &tlist,
		core::Size nint
	);

	void init();

	~ReplicaExchangeMC();

	void build_temperature_list(double *elist);

	using Parent::boltzmann;

    virtual bool
    boltzmann(
        Pose & pose,//PoseOP pose,
        std::string const & move_type = "unk",
        core::Real const proposal_density_ratio = 1,
				core::Real const inner_score_delta_over_temperature = 0
    );

    //void set_noutput(core::Size n){noutput_=n;}

private:
    int rank_;
    int size_;
    core::Size nreplica_frequency_;
    //core::Size noutput_;
    core::Size ntrials_;
    utility::vector1<core::Real> Tlist_;
    utility::vector1< utility::vector1<std::pair<int, int> > > exchange_schedule;
    double *last_energylist;
    int *T_tag;
    int *T_rev;
    int T_ndx;
};

} // moves
} // prot

#endif

