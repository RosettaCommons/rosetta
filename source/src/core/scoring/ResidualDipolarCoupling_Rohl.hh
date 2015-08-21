// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ResidualDipolarCoupling.hh
/// @brief  Uses NMR RDC for scoring
/// @author Srivatsan Raman

#ifndef INCLUDED_core_scoring_ResidualDipolarCoupling_Rohl_hh
#define INCLUDED_core_scoring_ResidualDipolarCoupling_Rohl_hh

#include <core/scoring/ResidualDipolarCoupling_Rohl.fwd.hh>

#include <core/types.hh>
#include <basic/datacache/CacheableData.hh>
#include <numeric/numeric.functions.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

void store_RDC_ROHL_in_pose( ResidualDipolarCoupling_RohlOP, core::pose::Pose& );
ResidualDipolarCoupling_RohlOP retrieve_RDC_ROHL_from_pose(  core::pose::Pose& );
ResidualDipolarCoupling_RohlCOP retrieve_RDC_ROHL_from_pose(  core::pose::Pose const& );

class ResidualDipolarCoupling_Rohl : public basic::datacache::CacheableData {

public: // typedefs
	typedef core::Real Real;
	typedef core::Size Size;
	typedef utility::vector1< core::scoring::RDC_Rohl > RDC_lines;

public:

	ResidualDipolarCoupling_Rohl(){
		read_RDC_file();
	}

	ResidualDipolarCoupling_Rohl( RDC_lines data_in ) : All_RDC_lines( data_in ) {};

	//  ResidualDipolarCoupling_Rohl( ResidualDipolarCoupling_Rohl const & src ){}

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new ResidualDipolarCoupling_Rohl( *this ) );
	}


	void read_RDC_file();

	Size get_RDC_data_type(
		std::string const & atom1,
		std::string const & atom2
	);

	inline RDC_lines get_RDC_data() const
	{
		return All_RDC_lines;
	}

private:
	RDC_lines All_RDC_lines;

};

/////////////////////////////////////////////////
//@brief short class that stores the RDC data lines
/////////////////////////////////////////////////
class RDC_Rohl {

public:
	RDC_Rohl(){}

	RDC_Rohl(
		Size type,
		Size res,
		Real Jdipolar,
		Real weight = 1.0//for alignment calculation
		//  core::Real Reduced_Jdipolar
	) :
		type_( type ),
		res_( res ),
		Jdipolar_( Jdipolar ),
		weight_( weight )
		//  Reduced_Jdipolar_( Reduced_Jdipolar )
	{}

	inline Size type() const
	{
		return type_;
	}

	inline Size res() const
	{
		return res_;
	}

	inline Real Jdipolar() const
	{
		return Jdipolar_;
	}

	/* inline core::Real Reduced_Jdipolar() const
	{
	return Reduced_Jdipolar_;
	}
	*/

	inline Real fixed_dist() const
	{
		Real fixed_dist(0.0);
		if ( type_ ==  1 ) {
			fixed_dist = 1.01;
		} else if ( type_ == 2 ) {
			fixed_dist = 1.08;
		} else if ( type_ == 3 ) {
			fixed_dist = 1.52325877;
		}
		//************* ADD MORE TYPES LATER !! ***************
		return fixed_dist;

	}

	inline Real Reduced_Jdipolar() const
	{
		using namespace numeric;

		//   Real invDcnst(0.0);
		//     if ( type_ == 1 )
		//       invDcnst = 0.0000821215;
		//     else if ( type_ == 2 )
		//       invDcnst = -0.0000331025;
		//     else if ( type_ == 3 )
		//       invDcnst = 0.000326533;
		return Jdipolar_*invDcnst()*numeric::cube( fixed_dist() );

	}

	Real weight() const {
		return weight_;
	}

	void weight( Real w_in ) {
		weight_ = w_in;
	}

	inline Real invDcnst() const
	{

		core::Real invDcnst(0.0);
		if ( type_ == 1 ) {
			invDcnst = 0.0000821215;
		} else if ( type_ == 2 ) {
			invDcnst = -0.0000331025;
		} else if ( type_ == 3 ) {
			invDcnst = 0.000326533;
		}
		return invDcnst;
	}

private:
	Size type_, res_;
	Real Jdipolar_, Reduced_Jdipolar_;
	Real weight_;
};


} //scoring
} //core

#endif
