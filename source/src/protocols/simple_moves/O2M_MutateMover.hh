// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/O2M_MutateMover.hh
/// @brief Mover that takes a starting structure and some task ops, and produces every single possible point mutation
/// DOES NOT WORK WITH ROSETTASCRIPTS
/// @author Ken Jung

#ifndef INCLUDED_protocols_simple_moves_O2M_MutateMover_hh
#define INCLUDED_protocols_simple_moves_O2M_MutateMover_hh
#include <protocols/simple_moves/O2M_MutateMover.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/exit.hh>

namespace protocols{
namespace simple_moves{

class O2M_MutateMover : public protocols::moves::Mover
{
public:
	typedef core::pose::PoseSP PoseSP;

	O2M_MutateMover(){}
	~O2M_MutateMover(){}

	void apply( core::io::serialization::PipeMap & pmap);
	
	// stupid mover base class crap
	void apply( Pose & ){
		utility_exit_with_message( "DOESNT WORK WITH SINGLE POSE" );
	}
	std::string get_name() const { return "stupid jd2"; }
	void parse_my_tag( utility::tag::TagCOP const , basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) {
		utility_exit_with_message( "DOESNT WORK WITH JD2" );
	}

	void parse_def( utility::lua::LuaObject const & def,
					utility::lua::LuaObject const & score_fxns,
					utility::lua::LuaObject const & tasks,
					protocols::moves::MoverCacheSP cache );

	protocols::moves::MoverSP create() { return protocols::moves::MoverSP( new O2M_MutateMover ); }
	static std::string name() {
		return "O2M_MutateMover";
	}
private:
	core::pack::task::TaskFactorySP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
};


}//simple_moves
}//protocols

#endif
