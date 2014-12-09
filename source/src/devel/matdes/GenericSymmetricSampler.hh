// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief 
/// @author Neil King ( neilking@uw.edu )
/// @author Javier Castellanos ( javiercv@uw.edu )
/// @author Frank DiMaio

#ifndef INCLUDED_devel_matdes_GenericSymmetricSampler_HH
#define INCLUDED_devel_matdes_GenericSymmetricSampler_HH

// Unit headers
#include <devel/matdes/GenericSymmetricSampler.fwd.hh>

// Package headers

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace devel {
namespace matdes {

class GenericSymmetricSampler : public protocols::moves::Mover {
	typedef core::pose::Pose Pose;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::filters::FilterOP FilterOP;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;
	typedef utility::tag::TagCOP TagCOP;

public:
  GenericSymmetricSampler();

  // --- virtual functions from mover ---
  virtual std::string get_name() const { return "GenericSymmetricSampler"; }
  virtual void apply(Pose& pose);

	// --- virtual copy constructors 
	virtual MoverOP clone() const;


	/// @brief create this type of object
	virtual MoverOP fresh_instance() const;


	virtual void parse_my_tag( TagCOP tag,
														 basic::datacache::DataMap & data,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );

private:
	std::string dof_id_;
	core::Real angle_min_,angle_max_,angle_step_;
	core::Real radial_disp_min_,radial_disp_max_,radial_disp_step_;

	bool usescore_, maxscore_;
	ScoreFunctionOP scorefxn_;
	FilterOP filter_;
	MoverOP mover_;
};

} // matdes
} // devel
#endif
