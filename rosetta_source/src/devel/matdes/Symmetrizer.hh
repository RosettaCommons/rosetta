
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

#ifndef INCLUDED_devel_matdes_Symmetrizer_HH
#define INCLUDED_devel_matdes_Symmetrizer_HH

// Unit headers
#include <devel/matdes/Symmetrizer.fwd.hh>

// Package headers

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace devel {
namespace matdes {

class Symmetrizer : public protocols::moves::Mover {
	typedef core::pose::Pose Pose;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::moves::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;
	typedef utility::tag::TagPtr TagPtr;

public:
  Symmetrizer();
	Symmetrizer(const Symmetrizer& rval);

  // --- virtual functions from mover ---
  virtual std::string get_name() const { return "Symmetrizer"; }
  virtual void apply(Pose& pose);

	// --- virtual copy constructors 
	virtual MoverOP clone() const;


	/// @brief create this type of object
	virtual MoverOP fresh_instance() const;


	virtual void parse_my_tag( TagPtr const tag,
														 DataMap & data,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );
private:
	Real get_radial_disp();
	Real get_angle();
	
private:
	std::string symm_file_;
	Real radial_disp_;
	Real angle_;
	char symmetry_axis_;
	
	bool explore_grid_;
};

} // matdes
} // devel
#endif
