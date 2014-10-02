// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/saxs/SAXSEnergyCreatorFA.hh
/// @brief  Declaration for the class that connects SAXSEnergyCreator with the ScoringManager
/// @author Dominik Gront dgront@chem.uw.edu.pl

#ifndef INCLUDED_core_scoring_saxs_SinXOverX_hh
#define INCLUDED_core_scoring_saxs_SinXOverX_hh

#include <core/types.hh>

// utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace saxs {

class SinXOverX : public utility::SingletonBase< SinXOverX >
{
public:
	friend class utility::SingletonBase< SinXOverX >;

public:

	core::Real evaluate(core::Real x) const {
		core::Size tmp_i = ((Size) (x * 100) + 1);
		return sin_x_over_x_[tmp_i];
	}

private:

	SinXOverX();
	void fill_sin_x_over_x_table();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static SinXOverX * create_singleton_instance();

private:

	utility::vector1<Real> sin_x_over_x_;

};

}
}
}

#endif
