// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/RepackMinimize.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_RepackMinimize_hh
#define INCLUDED_protocols_protein_interface_design_movers_RepackMinimize_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief One round of design/repacking followed by interface sc/bb and rigid-body minimization
class RepackMinimize : public simple_moves::DesignRepackMover
{
public:
	typedef core::pose::Pose Pose;
	typedef core::Real Real;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
public:
	RepackMinimize();
	RepackMinimize(
		ScoreFunctionCOP scorefxn_repack,
		ScoreFunctionCOP scorefxn_minimize,
		utility::vector1< core::Size > const target_residues,
		bool const repack_partner1=false,
		bool const repack_partner2=true,
		core::Real const interface_distance_cutoff=8.0,
		bool const repack_non_ala=true
	);
	virtual ~RepackMinimize();

	void apply( Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new RepackMinimize ); }
	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_RepackMinimize_HH*/
