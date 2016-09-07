// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/DumpPdb.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_moves_DumpPdb_hh
#define INCLUDED_protocols_simple_moves_DumpPdb_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @brief what you think
/// this can now be assimilated into DumpPdbMover
class DumpPdb : public protocols::moves::Mover
{
public:
	DumpPdb();
	DumpPdb( std::string fname ); // argument is moved
	~DumpPdb() override;
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	protocols::moves::MoverOP clone() const override {
		return( protocols::moves::MoverOP( new DumpPdb( *this ) ) );
	}
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new DumpPdb ); }
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn);
	void tag_time(bool setting) { addtime_ = setting; }
	bool tag_time() const { return addtime_; }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
private:
	std::string fname_;
	/// @brief Dump a scored pdb?
	core::scoring::ScoreFunctionOP scorefxn_;
	/// @brief Add timing information to filename?
	bool addtime_;
	/// @brief Add timing information to filename?
};


} // simple_moves
} // protocols


#endif /*INCLUDED_protocols_moves_DumpPdb_HH*/
