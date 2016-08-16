// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  A FilterMover that also calls report() on apply()
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @date   Jul. 30, 2014

#ifndef INCLUDED_protocols_peptide_deriver_FilterReporterMover_hh
#define INCLUDED_protocols_peptide_deriver_FilterReporterMover_hh

// Unit headers
#include <protocols/moves/FilterReporterMover.fwd.hh>

// Project headers
#include <protocols/filters/Filter.hh>
#include <protocols/moves/FilterMover.hh>
#include <protocols/moves/Mover.hh>

// C++ header
#include <iostream>
#include <string>

namespace protocols {
namespace moves {

/// A FilterMover that also calls report() on apply()
// TODO : perhaps this should be added to FilterMover? FilterMover is only being used by DME_FilterMover, currently.
// TODO : is placement in protocols/moves (protocols.1) correct? Or should it be in protocols/simple_moves (protocols.3)?
class FilterReporterMover: public protocols::moves::Mover {
private:
	protocols::moves::MoverOP mover_;
	protocols::filters::FilterOP filter_;

	/// mixin inheritence of FilterMover
	protocols::moves::FilterMoverOP filter_mover_;

	std::ostream out_;

	/// used by cctor and operator= to assign one such object to the other
	static void assign(FilterReporterMover & lhs,
		FilterReporterMover const & rhs);

public:
	// ctors and dtors
	FilterReporterMover();
	FilterReporterMover(protocols::moves::MoverOP const &,
		protocols::filters::FilterOP const &, core::Size const,
		std::ostream &, protocols::moves::MoverStatus const =
		protocols::moves::FAIL_DO_NOT_RETRY);
	FilterReporterMover(FilterReporterMover const &);
	~FilterReporterMover();

	// assignment operator
	FilterReporterMover & operator=(FilterReporterMover const & rhs);

	// pure virtual overrides
	virtual void apply(core::pose::Pose & pose);
	virtual protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	// accessors
	protocols::moves::MoverOP get_mover() const;
	void set_mover(protocols::moves::MoverOP const & mover);

	protocols::filters::FilterOP get_filter() const;
	void set_filter(protocols::filters::FilterOP const & filter);

	std::ostream & get_out();
	void set_out(std::ostream const & out);
};

}  // namespace peptide_deriver
}  // namespace protocols

#endif
