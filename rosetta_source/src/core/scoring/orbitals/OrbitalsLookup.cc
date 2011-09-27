// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.




#include <core/scoring/rna/RNA_Util.hh>

#include <utility/vector1.hh>
#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/MathVector_operations.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility/io/izstream.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/scoring/orbitals/OrbitalsScore.hh>
#include <core/scoring/orbitals/OrbitalsScoreCreator.hh>
#include <core/chemical/orbitals/OrbitalTypeMapper.hh>
#include <basic/database/open.hh>
#include <map>
#include <basic/options/option.hh>
#include <basic/options/keys/orbitals.OptionKeys.gen.hh>



namespace core{
namespace scoring{
namespace orbitals{

//this function reads in the potentials from database. The potentials are already in the form of Energy.
OrbitalsLookup::OrbitalsLookup( utility::vector1< std::string > const & filename ) :
	number_stats_(core::chemical::orbitals::num_orbital_types),
	number_elements_(1000)
{

	//set up initial vectors. Data read from the database will be placed into these vectors
	utility::vector1< core::Real > E_vector(number_elements_, 0);
	utility::vector1< utility::vector1< core::Real > > hpol_E_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > sc_orb_bb_H_E_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > haro_E_vector(static_cast <core::Size> (number_stats_), E_vector);
	//utility::vector1< utility::vector1< core::Real > > orb_orb_E_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > orb_bb_E_vector(static_cast <core::Size> (number_stats_), E_vector);



	std::string line;

	utility::io::izstream stream;
	basic::database::open( stream, filename[1] );


	//read in database file. This is the hpol file, put files into the vectors constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> hpol_E_vector[1][count] >> hpol_E_vector[2][count] >> hpol_E_vector[3][count] >> hpol_E_vector[4][count] >> hpol_E_vector[5][count] >> hpol_E_vector[6][count] >> hpol_E_vector[7][count];
		//std::cout << "printing out numbers " << hpol_E_vector[1][count] << hpol_E_vector[2][count] << hpol_E_vector[3][count] << hpol_E_vector[4][count] << hpol_E_vector[5][count] << hpol_E_vector[6][count] << hpol_E_vector[7][count];
	}
	stream.close();

	//haro data
	basic::database::open( stream, filename[2] );
	//read in database file. This is the haro file, put files into the array constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> haro_E_vector[1][count] >> haro_E_vector[2][count] >> haro_E_vector[3][count] >> haro_E_vector[4][count] >> haro_E_vector[5][count] >> haro_E_vector[6][count] >> haro_E_vector[7][count];
	}
	stream.close();

	//orb_bb
	basic::database::open( stream, filename[3] );
	//read in database file. This is the haro file, put files into the array constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> orb_bb_E_vector[1][count] >> orb_bb_E_vector[2][count] >> orb_bb_E_vector[3][count] >> orb_bb_E_vector[4][count] >> orb_bb_E_vector[5][count] >> orb_bb_E_vector[6][count] >> orb_bb_E_vector[7][count];
		//std::cout << orb_e1[count] << " <- orb_e1, orb_e2->" << orb_e2[count] << std::endl;
	}
	stream.close();

	//orb_bb
	basic::database::open( stream, filename[4] );
	//read in database file. This is the haro file, put files into the array constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> sc_orb_bb_H_E_vector[1][count] >> sc_orb_bb_H_E_vector[2][count] >> sc_orb_bb_H_E_vector[3][count] >> sc_orb_bb_H_E_vector[4][count] >> sc_orb_bb_H_E_vector[5][count] >> sc_orb_bb_H_E_vector[6][count] >> sc_orb_bb_H_E_vector[7][count];
		//std::cout << orb_e1[count] << " <- orb_e1, orb_e2->" << orb_e2[count] << std::endl;
	}
	stream.close();


	//initial construction of a vector of MathMatrix. We will be pushing back matrixes into this
	utility::vector1< numeric::MathMatrix<core::Real> > hpol_E_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > sc_orb_bb_H_E_vector_matrix;

	utility::vector1< numeric::MathMatrix<core::Real> > haro_E_vector_matrix;
	//utility::vector1< numeric::MathMatrix<core::Real> > orb_orb_E_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > orb_bb_E_vector_matrix;

	for( core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count ) {
		//MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
		numeric::MathMatrix<core::Real> hpol_E_matrix(50, 20, & hpol_E_vector[count][1] );
		numeric::MathMatrix<core::Real> sc_orb_bb_H_E_matrix(50, 20, & sc_orb_bb_H_E_vector[count][1] );

		numeric::MathMatrix<core::Real> haro_E_matrix(50, 20, & haro_E_vector[count][1]);
		numeric::MathMatrix<core::Real> orb_bb_E_matrix(50, 20, & orb_bb_E_vector[count][1]);

		hpol_E_vector_matrix.push_back(hpol_E_matrix);
		sc_orb_bb_H_E_vector_matrix.push_back(sc_orb_bb_H_E_matrix);

		haro_E_vector_matrix.push_back(haro_E_matrix);
		orb_bb_E_vector_matrix.push_back(orb_bb_E_matrix);
	}

	//the spline will behave as a first derivative. This means that the borders of the spline will be smoothed out if the first_derivative
	//is set to 0. This is controlled through an enum and put into an array.
	numeric::interpolation::spline::BorderFlag behavior[2] = {
		numeric::interpolation::spline::e_FirstDer,
		numeric::interpolation::spline::e_FirstDer };

	core::Real const start[2] = {0.05, -0.95}; //starting position of values in the matrix.0
	core::Real const delta[2] = {0.1,  0.1  }; //change in values. Needed to increment the values
	bool const linear_cont[2] = {true, true }; //if outside of range continue linearly


	//slls

	//what is the characteristic of the spline at 0 values. Should be smoothed out at edges. Powerful if used with e_FristDer
	std::pair< core::Real, core::Real > const first_deriv[2] = {
		std::pair< core::Real, core::Real >(0.0,0.0),
		std::pair< core::Real, core::Real >(0.0,0.0) };


	numeric::interpolation::spline::BicubicSpline bicubic_spline;// apl -- fixing memory leak = new numeric::interpolation::spline::BicubicSpline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> hpol_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> sc_orb_bb_H_vector_spline;

	utility::vector1<numeric::interpolation::spline::BicubicSpline> haro_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> orb_bb_vector_spline;

	hpol_vector_spline.reserve( number_stats_ );
	sc_orb_bb_H_vector_spline.reserve( number_stats_ );
	haro_vector_spline.reserve( number_stats_ );
	orb_bb_vector_spline.reserve( number_stats_ );

	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		bicubic_spline.train(behavior, start, delta, hpol_E_vector_matrix[count], linear_cont, first_deriv);
		hpol_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, sc_orb_bb_H_E_vector_matrix[count], linear_cont, first_deriv);
		sc_orb_bb_H_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, haro_E_vector_matrix[count], linear_cont, first_deriv);
		haro_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, orb_bb_E_vector_matrix[count], linear_cont, first_deriv);
		orb_bb_vector_spline.push_back( bicubic_spline );

	}

	//store all splines in a private member look up table. This allows for looking up the values in the actual score and provides an
	//interface for the OrbitalsScore class.
	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		hpol_interpolation_splines_.push_back(hpol_vector_spline[count]);
		hpol_sc_orb_bb_Hinterpolation_splines_.push_back(sc_orb_bb_H_vector_spline[count]);
		haro_interpolation_splines_.push_back(haro_vector_spline[count]);
		orb_bb_interpolation_splines_.push_back(orb_bb_vector_spline[count]);

	}


	for(int count=core::chemical::orbitals::C_pi_sp2; count <= static_cast <core::Size> (number_stats_); ++count){
		hpol_splines_.insert(std::pair<core::Size, numeric::interpolation::spline::BicubicSpline>(count, hpol_interpolation_splines_[count] ) );
		hpol_sc_orb_bb_splines_.insert(std::pair<core::Size, numeric::interpolation::spline::BicubicSpline>(count, hpol_sc_orb_bb_Hinterpolation_splines_[count] ) );
		haro_splines_.insert(std::pair<core::Size, numeric::interpolation::spline::BicubicSpline>(count, haro_interpolation_splines_[count] ) );
		orb_bb_splines_.insert(std::pair<core::Size, numeric::interpolation::spline::BicubicSpline>(count, orb_bb_interpolation_splines_[count] ) );
	}



 //some verbose checking. Not needed
/*
	for(core::Real i=0.05; i <= 50; i+=0.1){
	core::Size number=0;
			for(core::Real j=-0.95; j<=1; j+=0.1){
				++number;
				if( number == 20){
					std::cout << hpol_vector_spline[1].F((numeric::MakeVector(i,j))) << std::endl;

				}else{
					std::cout << hpol_vector_spline[1].F((numeric::MakeVector(i,j))) << " ";
					//std::cout << i << " " << j << std::endl;
				}
			}
		}
*/







}//end orbitalsLookup

//get the energy from the bicubic spline. This is done by using the mm atom name as a key in the map. We pass in the distance and angle.
void OrbitalsLookup::get_energy
(
	const h_type h_enum,
	const core::chemical::orbitals::orbital_type_enum orb_type_name,
	const core::Real distance,
	const core::Real angle,
	core::Real & energy,
	core::Real & distance_derivative,
	core::Real & angle_derivative,
	bool check_derivative
) const
{
	if ( distance > 5.05 ) { energy = distance_derivative = angle_derivative = 0.0; return; }

	numeric::interpolation::spline::BicubicSpline spline;
	if(check_derivative){
		if(h_enum == hpol){
			spline = hpol_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			distance_derivative = spline.dFdx(numeric::MakeVector(distance, angle));
			angle_derivative = spline.dFdy(numeric::MakeVector(distance, angle));

		}
		if(h_enum == sc_orb_bb_H){
			spline = hpol_sc_orb_bb_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			distance_derivative = spline.dFdx(numeric::MakeVector(distance, angle));
			angle_derivative = spline.dFdy(numeric::MakeVector(distance, angle));
		}
		if(h_enum == haro){
			spline = haro_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			distance_derivative = spline.dFdx(numeric::MakeVector(distance, angle));
			angle_derivative = spline.dFdy(numeric::MakeVector(distance, angle));
		}
		if(h_enum == bb){
			spline = orb_bb_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			distance_derivative = spline.dFdx(numeric::MakeVector(distance, angle));
			angle_derivative = spline.dFdy(numeric::MakeVector(distance, angle));
		}
		energy = spline.F(numeric::MakeVector(distance, angle));
	}else{
		if(h_enum == hpol){
			spline = hpol_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			energy = spline.F(numeric::MakeVector(distance, angle));
		}
		if(h_enum == sc_orb_bb_H){
			spline = hpol_sc_orb_bb_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			energy = spline.F(numeric::MakeVector(distance, angle));
		}
		if(h_enum == haro){
			spline = haro_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			energy = spline.F(numeric::MakeVector(distance, angle));
		}
		if(h_enum == bb){
			spline = orb_bb_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			energy = spline.F(numeric::MakeVector(distance, angle));
		}
	}
	/*if(energy <= 5.0 && energy >= -0.0001){
		energy=0;
	}
	if(distance_derivative <= 5.0 && distance_derivative >= -0.0001){
		distance_derivative=0;
	}
	if(angle_derivative <= 5.0 && angle_derivative >= -0.0001){
		angle_derivative=0;
	}*/

//	std::cout<< "ENERGY: " << energy << std::endl;

}





}//namespace orbitals

}//namespace scoring

}//namespace core
