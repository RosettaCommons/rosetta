// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/FormFactor.hh
/// @brief Represents an atomic scattering form factor
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_core_scoring_saxs_FormFactor_hh
#define INCLUDED_core_scoring_saxs_FormFactor_hh

#include <core/scoring/saxs/FormFactor.fwd.hh>

// utility headers
#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/interpolation/spline/Interpolator.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace saxs {

/// @brief
class FormFactor: public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FormFactor();

	/// @ A unique identifier of a form factor object, set by the manager
	Size id_;

	/// @brief  Constructor reads a file with a spline function
	FormFactor(std::string const &,std::string const &);

	/// @brief  evaluates an atomic form factor for a given scattering angle (defined in 1/A)
	Real ff(Real q) const {

		Real y = 0.0;
		Real dy = 0.0;
		spline_interpolator_->interpolate(q, y, dy);

		return y;
	}

	/// @brief Returns tabulated ff-value (computed for i-th value of q-argument)
	Real get(Size q_index) { return ff_values_[q_index]; }

	std::string & name() { return name_; }

	void is_glob(bool flag) { glob_flag_ = flag; }

	bool is_glob() const { return glob_flag_; }

	void tabulate(const utility::vector1<Real> & q);
private:
	bool glob_flag_;
	utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator > spline_interpolator_;
	std::string name_;
	mutable utility::vector1<Real> ff_values_;
};

} // core
} // scoring
} // saxs

#endif /* INCLUDED_core_scoring_saxs_FormFactor_HH */
