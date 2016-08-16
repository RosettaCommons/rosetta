// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: LoopModeler
///
///
/// @author Vatsan Raman


#ifndef INCLUDED_devel_residual_dipolar_coupling_RDC_main_hh
#define INCLUDED_devel_residual_dipolar_coupling_RDC_main_hh


//Core
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers
#include <numeric/numeric.functions.hh>

//// C++ headers
#include <map>

#include <utility/vector1_bool.hh>
#include <iostream>


namespace devel {
namespace residual_dipolar_coupling {


class RDC {

public:
	RDC(){}

	RDC(
		core::Size type,
		core::Size res1,
		core::Size res2,
		core::Real Jdipolar
		//  core::Real Reduced_Jdipolar
	) :
		type_( type ),
		res1_( res1 ),
		res2_( res2 ),
		Jdipolar_( Jdipolar )
		//  Reduced_Jdipolar_( Reduced_Jdipolar )
	{}

	inline core::Size type() const
	{
		return type_;
	}

	inline core::Size res1() const
	{
		return res1_;
	}

	inline core::Size res2() const
	{
		return res2_;
	}

	inline core::Real Jdipolar() const
	{
		return Jdipolar_;
	}

	/* inline core::Real Reduced_Jdipolar() const
	{
	return Reduced_Jdipolar_;
	}
	*/

	inline core::Real fixed_dist() const
	{
		core::Real fixed_dist(0.0);
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

	inline core::Real Reduced_Jdipolar() const
	{
		using namespace numeric;

		core::Real invDcnst(0.0);
		if ( type_ == 1 ) {
			invDcnst = 0.0000821215;
		} else if ( type_ == 2 ) {
			invDcnst = -0.0000331025;
		} else if ( type_ == 3 ) {
			invDcnst = 0.000326533;
		}
		return Jdipolar_*invDcnst*cube( fixed_dist() );

	}


	inline core::Real invDcnst() const
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
	core::Size type_, res1_, res2_;
	core::Real Jdipolar_, Reduced_Jdipolar_;

};


/////////////////////////////////////////////////////////////////

class RDC_data {

public:
	RDC_data(){}

	RDC_data(
		std::string filename
	) :
		filename_( filename )
	{
		read_RDC_file( filename_ );
	}

	void read_RDC_file(
		std::string const & filename
	);

	core::Size get_RDC_data_type(
		std::string const & atom1,
		std::string const & atom2
	);

	core::Real get_invDmax(
		std::string const & atom1,
		std::string const & atom2
	);

	inline std::map< core::Size, utility::vector1<devel::residual_dipolar_coupling::RDC> > get_RDC_data()
	{
		return RDC_data_lines;
	};

private:
	std::string filename_;
	std::map< core::Size, utility::vector1<devel::residual_dipolar_coupling::RDC> > RDC_data_lines;


};

/////////////////////////////////////////////////////////////////

class RDC_data_set {

public:
	RDC_data_set(){}

	RDC_data_set(
		core::Size data_type
	) :
		data_type_( data_type )
	{}


	inline void add_data_line(
		core::Size type,
		core::Size res1,
		core::Size res2,
		core::Real Jdipolar
	)
	{
		RDC_lines_.push_back( RDC( type, res1, res2, Jdipolar ) );
	}

	//private:
	utility::vector1< devel::residual_dipolar_coupling::RDC> RDC_lines_;
	core::Size data_type_;

};


utility::vector1< devel::residual_dipolar_coupling::RDC > read_RDC_file(
	std::string const & filename
);

core::Size get_RDC_data_type(
	std::string const & atom1,
	std::string const & atom2
);


} //ResidualDipolarCoupling
} //devel

#endif
