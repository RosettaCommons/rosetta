// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCSingleSet.cc
/// @brief   Implementation of class RDCSingleSet
/// @details last Modified: 07/27/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/rdc/RDCSingleSet.hh>

// Package headers
#include <core/scoring/nmr/rdc/RDCSingle.hh>
#include <core/io/nmr/AtomSelection.hh>
#include <core/io/nmr/util.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// C++ headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

// Boost headers
#include <boost/algorithm/string.hpp>

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

static basic::Tracer TR( "core.scoring.nmr.rdc.RDCSingleSet" );

/// @brief construct from filedata
///        set default values for RDCSingleSet weight and single_rdc_weighting_scheme
RDCSingleSet::RDCSingleSet(
	std::string const & filename,
	pose::Pose const & pose
) :
	dataset_name_(utility::file_basename(filename)),
	weight_(1.0),
	single_rdc_weighting_scheme_(CONST)
{
	init_from_rdc_filedata(filename, pose);
}

/// @brief constructor with full argument list
RDCSingleSet::RDCSingleSet(
	std::string const & filename,
	pose::Pose const & pose,
	Real const weight,
	std::string single_rdc_weighting
) :
	dataset_name_(utility::file_basename(filename)),
	weight_(weight)
{
	single_rdc_weighting_scheme_ = convert_string_to_weighting_scheme(single_rdc_weighting);
	init_from_rdc_filedata(filename, pose);
}

/// @brief copy constructor
RDCSingleSet::RDCSingleSet(RDCSingleSet const & other) :
	dataset_name_(other.dataset_name_),
	rdc_single_vec_(other.rdc_single_vec_),
	weight_(other.weight_),
	number_rdc_(other.number_rdc_),
	single_rdc_weighting_scheme_(other.single_rdc_weighting_scheme_),
	rdc_type_(other.rdc_type_)
{ }

/// @brief assignment operator
RDCSingleSet&
RDCSingleSet::operator=(RDCSingleSet const & rhs) {
	if ( this != &rhs ) {
		dataset_name_ = rhs.dataset_name_;
		rdc_single_vec_ = rhs.rdc_single_vec_;
		weight_ = rhs.weight_;
		number_rdc_ = rhs.number_rdc_;
		single_rdc_weighting_scheme_ = rhs.single_rdc_weighting_scheme_;
		rdc_type_ = rhs.rdc_type_;
	}
	return *this;
}

/// @brief destructor
RDCSingleSet::~RDCSingleSet() { }

void
RDCSingleSet::set_single_rdc_weighting_scheme(std::string const & weighting_scheme) {
	single_rdc_weighting_scheme_ = convert_string_to_weighting_scheme(weighting_scheme);
}

void
RDCSingleSet::show(std::ostream & tracer) const {
	tracer << "   * * * RDCSingleSet Summary Report * * *   " << std::endl;
	tracer << " RDC dataset " << dataset_name_ << " contains " << number_rdc_ << " rdc values and has weight " << weight_ << "." << std::endl;
	tracer << " * * * * * * RDC values * * * * * * " << std::endl;
	for ( Size i = 1; i <= number_rdc_; ++i ) {
		rdc_single_vec_[i].show(tracer);
	}
}

/// @brief utility function used in constructor to initialize RDCSingelSet object from data file.
void
RDCSingleSet::init_from_rdc_filedata(
	std::string const & filename,
	pose::Pose const & pose
)
{
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;
	using namespace basic::options;
	using core::io::nmr::AtomSelection;

	conformation::symmetry::SymmetryInfoCOP syminfo_ptr;
	if ( is_symmetric( pose ) ) {
		syminfo_ptr = symmetry_info( pose );
	}
	dataset_name_ = utility::file_basename(filename);

	utility::vector1< utility::vector1< AtomSelection > > spinsA, spinsB;
	utility::vector1<Real> values, errors;
	core::io::nmr::read_rdc_datafile(filename, spinsA, spinsB, values, errors);
	runtime_assert( (spinsA.size() == spinsB.size()) && (spinsA.size() == values.size()) && (spinsA.size() == errors.size()) );
	number_rdc_ = values.size();

	TR.Info << "Creating RDCSingleSet from file " << dataset_name_ << ". Number of RDCs: " << number_rdc_ << std::endl;

	for ( Size i = 1; i <= number_rdc_; ++i ) {
		runtime_assert(spinsA[i].size() == spinsB[i].size());
		utility::vector1<std::pair<AtomSelection, AtomSelection> > spinsAB_per_rdc;
		spinsAB_per_rdc.reserve(spinsA[i].size());
		for ( Size j = 1, j_end = spinsA[i].size(); j <= j_end; ++j ) {
			spinsAB_per_rdc.push_back(std::make_pair(spinsA[i][j], spinsB[i][j]));
		}
		RDCSingle single_rdc(spinsAB_per_rdc, pose, values[i], errors[i]);

		// Make sure that RDC values are all of the same type
		if ( i == 1 ) { rdc_type_ = single_rdc.get_rdc_type(); }
		if ( rdc_type_ != single_rdc.get_rdc_type() ) {
			utility_exit_with_message("ERROR in creation of RDCSingleSet for dataset " + dataset_name_ + ". RDC input datafile must contain only RDCs of the same type.");
		}
		// If we perform RDC calculation with automatic symmetry deduction, make sure that only RDCs for the ASU are provided
		if ( basic::options::option[ basic::options::OptionKeys::nmr::rdc::use_symmetry_calc ]() && is_symmetric(pose) ) {
			for ( Size j = 1, j_end = single_rdc.get_spinsAB().size(); j <= j_end; ++j ) {
				// We test that the spins belong to subunit with index 1
				if ( syminfo_ptr->subunit_index(single_rdc.get_spinsAB()[j].first.rsd())  != 1 ||
						syminfo_ptr->subunit_index(single_rdc.get_spinsAB()[j].second.rsd()) != 1 ) {
					utility_exit_with_message("ERROR in creation of RDCSingleSet for dataset " + dataset_name_
						+ ". For RDC calculation with automatic symmetry deduction the input datafile must contain only residue selections for the asymmetric subunit.");
				}
			}
		}
		rdc_single_vec_.push_back(single_rdc);
	}
}

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core
