// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.




//#include <core/scoring/rna/RNA_Util.hh>

#include <utility/vector1.hh>
#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/MathVector_operations.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility/io/izstream.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <basic/database/open.hh>

#include <utility/exit.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>


namespace core{
namespace scoring{
namespace orbitals{

//this function reads in the potentials from database. The potentials are already in the form of Energy.
OrbitalsLookup::OrbitalsLookup( utility::vector1< std::string > const & DHO_energies, utility::vector1< std::string > const & AOH_energies) :
	number_stats_(6),
	number_elements_(600)
{

	if(basic::options::option[ basic::options::OptionKeys::in::add_orbitals] != 1 ){
		utility_exit_with_message( "Trying to run orbitals score without orbitals! Pass the flag -add_orbitals!" );
	}

	//set up initial vectors. Data read from the database will be placed into these vectors
	utility::vector1< core::Real > E_vector(number_elements_, 0);
	utility::vector1< utility::vector1< core::Real > > DHO_Hpol_scOrbH_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > DHO_Hpol_bbOrbH_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > DHO_Haro_scOrbH_vector(static_cast <core::Size> (number_stats_), E_vector);

	std::string line;
	utility::io::izstream stream;
	basic::database::open( stream, DHO_energies[1] );

	//sc hpol energy file processing
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> DHO_Hpol_scOrbH_vector[1][count] >> DHO_Hpol_scOrbH_vector[2][count] >> DHO_Hpol_scOrbH_vector[3][count] >> DHO_Hpol_scOrbH_vector[4][count] >> DHO_Hpol_scOrbH_vector[5][count] >> DHO_Hpol_scOrbH_vector[6][count];
	}
	stream.close();

	//backbone hpol energy file processing
	basic::database::open( stream, DHO_energies[2] );
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> DHO_Hpol_bbOrbH_vector[1][count] >> DHO_Hpol_bbOrbH_vector[2][count] >> DHO_Hpol_bbOrbH_vector[3][count] >> DHO_Hpol_bbOrbH_vector[4][count] >> DHO_Hpol_bbOrbH_vector[5][count] >> DHO_Hpol_bbOrbH_vector[6][count];
	}
	stream.close();

	//haro energy file processing
	basic::database::open( stream, DHO_energies[3] );
	//read in database file. This is the haro file, put files into the array constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> DHO_Haro_scOrbH_vector[1][count] >> DHO_Haro_scOrbH_vector[2][count] >> DHO_Haro_scOrbH_vector[3][count] >> DHO_Haro_scOrbH_vector[4][count] >> DHO_Haro_scOrbH_vector[5][count] >> DHO_Haro_scOrbH_vector[6][count];
	}
	stream.close();

	//initial construction of a vector of MathMatrix. We will be pushing back matrixes into this
	utility::vector1< numeric::MathMatrix<core::Real> > DHO_Hpol_scOrbH_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > DHO_Hpol_bbOrbH_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > DHO_Haro_scOrbH_vector_matrix;

	for( core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count ) {
		//MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
		numeric::MathMatrix<core::Real> Hpol_scOrbH_matrix(30, 20, & DHO_Hpol_scOrbH_vector[count][1] );
		numeric::MathMatrix<core::Real> Hpol_bbOrbH_matrix(30, 20, & DHO_Hpol_bbOrbH_vector[count][1] );
		numeric::MathMatrix<core::Real> Haro_scOrbH_matrix(30, 20, & DHO_Haro_scOrbH_vector[count][1]);

		DHO_Hpol_scOrbH_vector_matrix.push_back(Hpol_scOrbH_matrix);
		DHO_Hpol_bbOrbH_vector_matrix.push_back(Hpol_bbOrbH_matrix);
		DHO_Haro_scOrbH_vector_matrix.push_back(Haro_scOrbH_matrix);

	}



	//the spline will behave as a natural spline. This is controlled through an enum and put into an array. Using the natural spline
	//allows for the border conditions to actuallay be met, unlike derivatives, where border conditions are never met
	numeric::interpolation::spline::BorderFlag behavior[2] = {
		numeric::interpolation::spline::e_Natural,
		numeric::interpolation::spline::e_Natural };

	core::Real const start[2] = {0.00, -1.0}; //starting position of values in the matrix.0
	core::Real const delta[2] = {0.1,  0.05  }; //change in values. Needed to increment the values
	bool const linear_cont[2] = {false, false }; //if outside of range continue linearly

	//what is the characteristic of the spline at 0 values. Should be smoothed out at edges. Powerful if used with e_FristDer
	std::pair< core::Real, core::Real > const first_deriv[2] = {
		std::pair< core::Real, core::Real >(0.0,0.0),
		std::pair< core::Real, core::Real >(0.0,0.0) };


	numeric::interpolation::spline::BicubicSpline bicubic_spline;// apl -- fixing memory leak = new numeric::interpolation::spline::BicubicSpline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> DHO_Hpol_scOrbH_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> DHO_Hpol_bbOrbH_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> DHO_Haro_scOrbH_vector_spline;

	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		bicubic_spline.train(behavior, start, delta, DHO_Hpol_scOrbH_vector_matrix[count], linear_cont, first_deriv);
		DHO_Hpol_scOrbH_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, DHO_Hpol_bbOrbH_vector_matrix[count], linear_cont, first_deriv);
		DHO_Hpol_bbOrbH_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, DHO_Haro_scOrbH_vector_matrix[count], linear_cont, first_deriv);
		DHO_Haro_scOrbH_vector_spline.push_back( bicubic_spline );

	}

	//store all splines in a private member look up table. This allows for looking up the values in the actual score and provides an
	//interface for the OrbitalsScore class.
	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		DHO_Hpol_scOrbH_splines_.push_back(DHO_Hpol_scOrbH_vector_spline[count]);
		DHO_Hpol_bbOrbH_splines_.push_back(DHO_Hpol_bbOrbH_vector_spline[count]);
		DHO_Haro_scOrbH_splines_.push_back(DHO_Haro_scOrbH_vector_spline[count]);
	}


 //some verbose checking. Not needed

/*

	for(core::Real i=0.00; i <= 30; i+=0.1){
	core::Size number=0;
			for(core::Real j=-1.0; j<=0; j+=0.05){
				++number;
				if( number == 20){
					std::cout << DHO_Hpol_scOrbH_vector_spline[5].F((numeric::MakeVector(i,j))) << std::endl;

				}else{
					std::cout << DHO_Hpol_scOrbH_vector_spline[5].F((numeric::MakeVector(i,j))) << " ";
					//std::cout << i << " " << j << std::endl;
				}
			}
		}
	std::cout << "###################################################################################\n#####################3" << std::endl;

	for(core::Real i=0.00; i <= 30; i+=0.1){
	core::Size number=0;
			for(core::Real j=-1.0; j<=0; j+=0.05){
				++number;
				if( number == 20){
					std::cout << i << ", " << j << std::endl;

				}else{
					std::cout << i << ", " << j << " ";
					//std::cout << i << " " << j << std::endl;
				}
			}
		}
*/

//std::cout << "###################################################################################\n#####################3" << std::endl;
	utility::vector1< utility::vector1< core::Real > > AOH_Hpol_scOrbH_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > AOH_Hpol_bbOrbH_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > AOH_Haro_scOrbH_vector(static_cast <core::Size> (number_stats_), E_vector);

	basic::database::open( stream, AOH_energies[1] );
	//sc hpol energy file processing
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> AOH_Hpol_scOrbH_vector[1][count] >> AOH_Hpol_scOrbH_vector[2][count] >> AOH_Hpol_scOrbH_vector[3][count] >> AOH_Hpol_scOrbH_vector[4][count] >> AOH_Hpol_scOrbH_vector[5][count] >> AOH_Hpol_scOrbH_vector[6][count];
	}
	stream.close();

	//backbone hpol energy file processing
	basic::database::open( stream, AOH_energies[2] );
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> AOH_Hpol_bbOrbH_vector[1][count] >> AOH_Hpol_bbOrbH_vector[2][count] >> AOH_Hpol_bbOrbH_vector[3][count] >> AOH_Hpol_bbOrbH_vector[4][count] >> AOH_Hpol_bbOrbH_vector[5][count] >> AOH_Hpol_bbOrbH_vector[6][count];
	}
	stream.close();

	//haro energy file processing
	basic::database::open( stream, AOH_energies[3] );
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> AOH_Haro_scOrbH_vector[1][count] >> AOH_Haro_scOrbH_vector[2][count] >> AOH_Haro_scOrbH_vector[3][count] >> AOH_Haro_scOrbH_vector[4][count] >> AOH_Haro_scOrbH_vector[5][count] >> AOH_Haro_scOrbH_vector[6][count];
	}
	stream.close();



	//initial construction of a vector of MathMatrix. We will be pushing back matrixes into this
	utility::vector1< numeric::MathMatrix<core::Real> > AOH_Hpol_scOrbH_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > AOH_Hpol_bbOrbH_vector_matrix;


	utility::vector1< numeric::MathMatrix<core::Real> > AOH_Haro_scOrbH_vector_matrix;

	for( core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count ) {
		//MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
		numeric::MathMatrix<core::Real> Hpol_scOrbH_matrix(30, 20, & AOH_Hpol_scOrbH_vector[count][1] );
		numeric::MathMatrix<core::Real> Hpol_bbOrbH_matrix(30, 20, & AOH_Hpol_bbOrbH_vector[count][1] );
		numeric::MathMatrix<core::Real> Haro_scOrbH_matrix(30, 20, & AOH_Haro_scOrbH_vector[count][1]);

		AOH_Hpol_scOrbH_vector_matrix.push_back(Hpol_scOrbH_matrix);
		AOH_Hpol_bbOrbH_vector_matrix.push_back(Hpol_bbOrbH_matrix);
		AOH_Haro_scOrbH_vector_matrix.push_back(Haro_scOrbH_matrix);

	}

	utility::vector1<numeric::interpolation::spline::BicubicSpline> AOH_Hpol_scOrbH_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> AOH_Hpol_bbOrbH_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> AOH_Haro_scOrbH_vector_spline;

	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		bicubic_spline.train(behavior, start, delta, AOH_Hpol_scOrbH_vector_matrix[count], linear_cont, first_deriv);
		AOH_Hpol_scOrbH_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, AOH_Hpol_bbOrbH_vector_matrix[count], linear_cont, first_deriv);
		AOH_Hpol_bbOrbH_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, AOH_Haro_scOrbH_vector_matrix[count], linear_cont, first_deriv);
		AOH_Haro_scOrbH_vector_spline.push_back( bicubic_spline );

	}

	//store all splines in a private member look up table. This allows for looking up the values in the actual score and provides an
	//interface for the OrbitalsScore class.
	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		AOH_Hpol_scOrbH_splines_.push_back(AOH_Hpol_scOrbH_vector_spline[count]);
		AOH_Hpol_bbOrbH_splines_.push_back(AOH_Hpol_bbOrbH_vector_spline[count]);
		AOH_Haro_scOrbH_splines_.push_back(AOH_Haro_scOrbH_vector_spline[count]);
	}



}//end orbitalsLookup

//get the energy from the bicubic spline. This is done by using the mm atom name as a key in the map. We pass in the distance and angle.
void OrbitalsLookup::OrbHdist_cosDHO_energy
(
	const h_type h_enum,
	const core::chemical::orbitals::orbital_type_enum orb_type_name,
	const core::Real distance,
	const core::Real DHO_angle,
	core::Real & energy,
	core::Real & distance_derivative,
	core::Real & angle_derivative,
	bool check_derivative
) const
{
	if ( distance > 3.0 ||  DHO_angle > 0.0 ) { energy = distance_derivative = angle_derivative = 0.0; return; }

	numeric::MathVector<core::Real> dist_angle_pair(numeric::MakeVector(distance, DHO_angle));
	if(check_derivative){
		if(h_enum == Hpol_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline( DHO_Hpol_scOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);

		}
		if(h_enum == Haro_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Haro_scOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
		if(h_enum == Hpol_bbOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Hpol_bbOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
	}else{
		if(h_enum == Hpol_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Hpol_scOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(dist_angle_pair);
		}
		if(h_enum == Haro_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Haro_scOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(dist_angle_pair);
		}
		if(h_enum == Hpol_bbOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Hpol_bbOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(dist_angle_pair);
		}
	}
}


void OrbitalsLookup::OrbHdist_cosAOH_energy
(
	const h_type h_enum,
	const core::chemical::orbitals::orbital_type_enum orb_type_name,
	const core::Real distance,
	const core::Real AOH_angle,
	core::Real & energy,
	core::Real & distance_derivative,
	core::Real & angle_derivative,
	bool check_derivative
) const
{
	if ( distance > 3.0 ||  AOH_angle > 0.0 ) { energy = distance_derivative = angle_derivative = 0.0; return; }

	numeric::MathVector<core::Real> dist_angle_pair(numeric::MakeVector(distance, AOH_angle));
	if(check_derivative){
		if(h_enum == Hpol_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline( AOH_Hpol_scOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);

		}
		if(h_enum == Haro_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(AOH_Haro_scOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
		if(h_enum == Hpol_bbOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(AOH_Hpol_bbOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
	}else{
		if(h_enum == Hpol_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(AOH_Hpol_scOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(dist_angle_pair);
		}
		if(h_enum == Haro_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(AOH_Haro_scOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(dist_angle_pair);
		}
		if(h_enum == Hpol_bbOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(AOH_Hpol_bbOrbH_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(dist_angle_pair);
		}
	}
}



}//namespace orbitals

}//namespace scoring

}//namespace core
