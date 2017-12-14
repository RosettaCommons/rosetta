// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/NumericalDerivCheckResult.hh
/// @brief  Declaration for nuerical derivative check results classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_optimization_NumericalDerivCheckResult_hh
#define INCLUDED_core_optimization_NumericalDerivCheckResult_hh

// Unit headers
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/DOF_ID.hh>

// Utility headers
#include <utility>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace optimization {

class DerivCheckDataPoint {
public:
	DerivCheckDataPoint() :
		num_deriv_( 0.0 ),
		ana_deriv_( 0.0 ),
		ratio_( 0.0 ),
		f11_( 0.0 ),
		f00_( 0.0 ),
		f22_( 0.0 ),
		dof_val_( 0.0 )
	{}

	DerivCheckDataPoint(
		Real   num_deriv,
		Real   ana_deriv,
		Real   ratio,
		Real   f11,
		Real   f00,
		Real   f22,
		Real   dof_val
	) :
		num_deriv_( num_deriv ),
		ana_deriv_( ana_deriv ),
		ratio_( ratio ),
		f11_( f11 ),
		f00_( f00 ),
		f22_( f22 ),
		dof_val_( dof_val )
	{}

	Real num_deriv() const { return num_deriv_; } // numerically computed derivative
	Real ana_deriv() const { return ana_deriv_; } // analytically computed derivative
	Real ratio() const { return ratio_; }
	Real f11() const { return f11_; }
	Real f00() const { return f00_; }
	Real f22() const { return f22_; }
	Real dof_val() const { return dof_val_; }

private:
	Real       num_deriv_; // numerically computed derivative
	Real       ana_deriv_; // analytically computed derivative
	Real       ratio_;
	Real       f11_;
	Real       f00_;
	Real       f22_;
	Real       dof_val_;
};

class SimpleDerivCheckResult : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~SimpleDerivCheckResult() override;
	SimpleDerivCheckResult( Size const nangles, Size const nsteps ) :
		step_data_( nangles ),
		abs_deriv_dev_( nangles, 0.0 ),
		rel_deriv_dev_( nangles )
	{
		for ( Size ii = 1; ii <= nangles; ++ii ) {
			step_data_[ ii ].resize( nsteps );
		}
	}

	Size nangles() const { return step_data_.size(); }
	Size nsteps() const { if ( step_data_.size() > 0 ) { return step_data_[ 1 ].size(); } return 0; }

	void abs_deriv_dev( Size dof_ind, Real val ) { abs_deriv_dev_[ dof_ind ] = val; }
	void rel_deriv_dev( Size dof_ind, Real val ) { rel_deriv_dev_[ dof_ind ] = val; }

	void best_cos_theta( Real val ) { best_cos_theta_ = val; }
	void best_abs_log_norm_ratio( Real val ) { best_abs_log_norm_ratio_ = val; }
	void best_norm_analytic( Real val ) { best_norm_analytic_ = val; }
	void best_norm_numeric( Real val ) { best_norm_numeric_ = val; }

	void step_data( Size dof_ind, Size step_ind, DerivCheckDataPoint const & dp ) { step_data_[ dof_ind ][ step_ind ] = dp; }

	Real abs_deriv_dev( Size dof_ind ) const { return abs_deriv_dev_[ dof_ind ]; }
	Real rel_deriv_dev( Size dof_ind ) const { return rel_deriv_dev_[ dof_ind ]; }

	Real best_cos_theta() const { return best_cos_theta_; }
	Real best_abs_log_norm_ratio() const { return best_abs_log_norm_ratio_; }
	Real best_norm_analytic() const { return best_norm_analytic_; }
	Real best_norm_numeric() const { return best_norm_numeric_; }

	DerivCheckDataPoint const & step_data( Size dof_ind, Size step_ind ) const { return step_data_[ dof_ind ][ step_ind ]; }

private:

	utility::vector1< utility::vector1< DerivCheckDataPoint > > step_data_;

	utility::vector1< Real > abs_deriv_dev_;
	utility::vector1< Real > rel_deriv_dev_;

	Real best_cos_theta_;
	Real best_abs_log_norm_ratio_;
	Real best_norm_analytic_;
	Real best_norm_numeric_;

};

class DOF_DataPoint {
public:
	DOF_DataPoint() :
		natoms_( 0 ),
		num_deriv_( 0.0 ),
		ana_deriv_( 0.0 ),
		ratio_( 0.0 ),
		f11_( 0.0 ),
		f00_( 0.0 ),
		f22_( 0.0 ),
		dof_val_( 0.0 )
	{}

	DOF_DataPoint(
		id::DOF_ID const & dof_id,
		id::DOF_ID const & parent_id,
		Size   natoms,
		Real   num_deriv,
		Real   ana_deriv,
		Real   ratio,
		Real   f11,
		Real   f00,
		Real   f22,
		Real   dof_val
	) :
		dof_id_( dof_id ),
		parent_id_( parent_id ),
		natoms_( natoms ),
		num_deriv_( num_deriv ),
		ana_deriv_( ana_deriv ),
		ratio_( ratio ),
		f11_( f11 ),
		f00_( f00 ),
		f22_( f22 ),
		dof_val_( dof_val )
	{}

	id::DOF_ID const & dof_id() const { return dof_id_; }
	id::DOF_ID const & parent_id() const { return parent_id_; }
	Size natoms() const { return natoms_; }
	Real num_deriv() const { return num_deriv_; } // numerically computed derivative
	Real ana_deriv() const { return ana_deriv_; } // analytically computed derivative
	Real ratio() const { return ratio_; }
	Real f11() const { return f11_; }
	Real f00() const { return f00_; }
	Real f22() const { return f22_; }
	Real dof_val() const { return dof_val_; }

private:
	id::DOF_ID dof_id_;
	id::DOF_ID parent_id_;
	Size       natoms_;
	Real       num_deriv_; // numerically computed derivative
	Real       ana_deriv_; // analytically computed derivative
	Real       ratio_;
	Real       f11_;
	Real       f00_;
	Real       f22_;
	Real       dof_val_;
};


class NumDerivCheckData : public SimpleDerivCheckResult {
public:
	typedef SimpleDerivCheckResult parent;

public:
	NumDerivCheckData( Size const nangles, Size const nsteps ):
		parent( nangles, nsteps ),
		dof_step_data_( nangles )
	{
		for ( Size ii = 1; ii <= nangles; ++ii ) {
			dof_step_data_[ ii ].resize( nsteps );
		}
	}

	~NumDerivCheckData() override = default;

	void dof_step_data( Size dof_ind, Size step_ind, DOF_DataPoint const & dofdp ) {
		dof_step_data_[ dof_ind ][ step_ind ] = dofdp;
		DerivCheckDataPoint dp( dofdp.num_deriv(), dofdp.ana_deriv(), dofdp.ratio(), dofdp.f11(), dofdp.f00(), dofdp.f22(), dofdp.dof_val() );
		parent::step_data( dof_ind, step_ind, dp );
	}
	DOF_DataPoint const & dof_step_data( Size dof_ind, Size step_ind ) const { return dof_step_data_[ dof_ind ][ step_ind ]; }


private:
	utility::vector1< utility::vector1< DOF_DataPoint > > dof_step_data_;

};

class NumericalDerivCheckResult : public utility::pointer::ReferenceCount {
public:
	typedef utility::pointer::ReferenceCount parent;

public:
	NumericalDerivCheckResult();
	~NumericalDerivCheckResult() override;

	void send_to_stdout( bool setting ) { send_to_stdout_ = setting; }
	bool send_to_stdout() const { return send_to_stdout_; }

	void add_deriv_data(
		NumDerivCheckDataOP deriv_check_data
	);

	Size n_deriv_check_results() const;

	NumDerivCheckData const &
	deriv_check_result( Size ind ) const;

private:
	bool send_to_stdout_;
	utility::vector1< NumDerivCheckDataOP > deriv_check_results_;
};


} // namespace optimization
} // namespace core


#endif
