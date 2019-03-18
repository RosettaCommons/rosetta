// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PRESingleSet.cc
/// @brief   Implementation of class PRESingleSet
/// @details last Modified: 08/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pre/PRESingleSet.hh>

// Package headers
#include <core/scoring/nmr/pre/PRESingle.hh>
#include <core/io/nmr/AtomSelection.hh>
#include <core/io/nmr/util.hh>

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
#include <iostream>
#include <iomanip>
#include <string>
#include <numeric>

namespace core {
namespace scoring {
namespace nmr {
namespace pre {

static basic::Tracer TR( "core.scoring.nmr.pre.PRESingleSet" );

/// @brief constructor with arguments
///        initialize PRESingleSet from data file, use mutator methods to set experimental conditions
PRESingleSet::PRESingleSet(
	std::string const & filename,
	pose::Pose const & pose
) :
	dataset_name_(utility::file_basename(filename)),
	weight_(1.0),
	scaling_factor_(1.0),
	rate_type_(R2_PARA),
	single_pre_weighting_scheme_(CONST),
	field_strength_(600.0),
	normalized_data_(false)
{
	register_options();
	init_from_cml();
	// read in the pre values and errors from file and fill the vector of PRESingle objects
	init_from_filedata(filename, pose);
}

/// @brief constructor with full argument list
///        construct PRESingleSet from data file and set computation, weighting and rate type
PRESingleSet::PRESingleSet(
	std::string const & filename,
	pose::Pose const & pose,
	Real const weight,
	std::string rate,
	std::string single_pre_weigting
) :
	dataset_name_(utility::file_basename(filename)),
	weight_(weight),
	scaling_factor_(1.0),
	field_strength_(600.0),
	normalized_data_(false)
{
	register_options();
	init_from_cml();
	single_pre_weighting_scheme_ = convert_string_to_weighting_scheme(single_pre_weigting);
	rate_type_ = convert_string_to_rate_type(rate);
	// read in the pre values and errors from file and fill the vector of PRESingle objects
	init_from_filedata(filename, pose);
}

/// @brief copy constructor
PRESingleSet::PRESingleSet(PRESingleSet const & other) :
	dataset_name_(other.dataset_name_),
	pre_single_vec_(other.pre_single_vec_),
	weight_(other.weight_),
	scaling_factor_(other.scaling_factor_),
	number_pre_(other.number_pre_),
	rate_type_(other.rate_type_),
	single_pre_weighting_scheme_(other.single_pre_weighting_scheme_),
	field_strength_(other.field_strength_),
	normalized_data_(other.normalized_data_)
{ }

/// @brief assignment operator
PRESingleSet&
PRESingleSet::operator=(PRESingleSet const & rhs) {
	if ( this != &rhs ) {
		dataset_name_ = rhs.dataset_name_;
		pre_single_vec_ = rhs.pre_single_vec_;
		weight_ = rhs.weight_;
		scaling_factor_ = rhs.scaling_factor_;
		number_pre_ = rhs.number_pre_;
		rate_type_ = rhs.rate_type_;
		single_pre_weighting_scheme_ = rhs.single_pre_weighting_scheme_;
		field_strength_ = rhs.field_strength_;
		normalized_data_ = rhs.normalized_data_;
	}
	return *this;
}

/// @brief destructor
PRESingleSet::~PRESingleSet() { }

void
PRESingleSet::init_from_filedata(
	std::string const & filename,
	pose::Pose const & pose
)
{
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	utility::vector1<Real> errors, values;
	utility::vector1< utility::vector1< core::io::nmr::AtomSelection > > protein_spins_all;
	io::nmr::read_pre_datafile(filename, protein_spins_all, values, errors);
	dataset_name_ = utility::file_basename(filename);
	number_pre_ = values.size();

	TR.Info << "Creating PRESingleSet from file " << dataset_name_ << ". Number of PREs: " << number_pre_ << std::endl;

	conformation::symmetry::SymmetryInfoCOP syminfo_ptr;
	if ( basic::options::option[ basic::options::OptionKeys::nmr::pre::use_symmetry_calc ]() && is_symmetric( pose ) ) {
		syminfo_ptr = symmetry_info( pose );
	}

	if ( normalized_data_ ) {
		// Normalize PRE data if needed and set the scaling factor
		Real mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
		Real sq_sum = 0.0;
		for ( const Real & v : values ) {
			sq_sum += (v - mean) * (v - mean);
		}
		scaling_factor_ = std::sqrt(sq_sum / values.size());
	} else {
		scaling_factor_ = 1.0;
	}

	for ( Size i = 1; i <= number_pre_; ++i ) {
		Real norm_pre = values[i] / scaling_factor_;
		Real norm_err = errors[i] / scaling_factor_;
		PRESingle single_pre(protein_spins_all[i], pose, norm_pre, norm_err);
		// If we perform PRE calculation with automatic deduction of protein symmetry,
		// make sure that only PREs for the ASU are provided
		if ( basic::options::option[ basic::options::OptionKeys::nmr::pre::use_symmetry_calc ]() && is_symmetric( pose ) ) {
			for ( Size j = 1, j_end = single_pre.get_protein_spins().size(); j <= j_end; ++j ) {
				// We test that the spins belong to subunit with index 1
				if ( syminfo_ptr->subunit_index(single_pre.get_protein_spins()[j].rsd()) != 1 ) {
					utility_exit_with_message("ERROR in creation of PRESingleSet for dataset " + dataset_name_
						+ ". For PRE calculation with automatic symmetry deduction the input datafile must contain only residue selections for the asymmetric subunit.");
				}
			}
		}
		pre_single_vec_.push_back(single_pre);
	}
}

/// @brief register options
void
PRESingleSet::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::pre::normalize_data);
}

void
PRESingleSet::init_from_cml() {
	using namespace basic::options;
	normalized_data_ = option[ basic::options::OptionKeys::nmr::pre::normalize_data ]();
}

/// @brief return gyromagnetic ratio of the nuclear spin in rad/(s*T) (dimension is 10^6)
Real
PRESingleSet::calc_gamma_I() const {
	PRE_SPIN_TYPE nuclear_spin = pre_single_vec_[1].get_pre_spin_type();
	Real gI(0.0);
	switch(nuclear_spin) {
	case PRE_SPIN_TYPE_H :
		//gI = 42.576e+6; // Hz/T
		gI = 267.513e+6; // rad/(s*T)
		break;
	case PRE_SPIN_TYPE_C :
		//gI = 10.705e+6; // Hz/T
		gI = 67.262e+6; // rad/(s*T)
		break;
	case PRE_SPIN_TYPE_N :
		//gI = -4.316e+6; // Hz/T
		gI = -27.116e+6; // rad/(s*T)
		break;
	default :
		utility_exit_with_message("ERROR while calculating the nuclear spin gyromagnetic ratio. No valid PRE spin type is set.");
	}
	return gI;
}

/// @brief calculate nuclear spin frequency at given field strength in rad/s
Real
PRESingleSet::calc_omega_I() const {
	runtime_assert_msg(field_strength_ > 0.0, "ERROR in calculating nuclear spin frequency. Magnetic field strength must be a positive value.");
	// Convert field strength from MHz into Tesla
	Real B0 = field_strength_ / 42.576;
	Real g = calc_gamma_I(); // g is in 10^6 rad/(s*T)
	return g * B0;
}

void
PRESingleSet::show(std::ostream & tracer) const {
	// Store old iostream manipulator flags
	std::ios oldState(nullptr);
	oldState.copyfmt(tracer);
	tracer << "   * * * PRESingleSet Summary Report * * *   " << std::endl;
	tracer << " PRE dataset " << dataset_name_ << " contains " << number_pre_ << " PRE values and has weight " << weight_ << "." << std::endl;
	tracer << " * * * * * * PRE values  * * * * * * " << std::endl;
	for ( Size i = 1; i <= number_pre_; ++i ) {
		pre_single_vec_[i].show(tracer);
	}
	tracer << " * * * Experimental conditions * * * " << std::endl;
	tracer << std::setprecision(1) << std::fixed;
	tracer << "Paramagnetic rate    = " << (rate_type_ == R2_PARA ? "R2" : "R1") << std::endl;
	tracer << "Field strength (MHz) = " << field_strength_ << std::endl;
	tracer << "Scaling factor       = " << std::scientific << scaling_factor_ << std::endl;
	tracer.copyfmt(oldState);
}

void
PRESingleSet::set_single_pre_weighting_scheme(std::string const & weighting_scheme) {
	single_pre_weighting_scheme_ = convert_string_to_weighting_scheme(weighting_scheme);
}

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core
