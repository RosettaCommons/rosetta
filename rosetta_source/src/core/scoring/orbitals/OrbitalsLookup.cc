// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.




// AUTO-REMOVED #include <core/scoring/rna/RNA_Util.hh>

#include <utility/vector1.hh>
#include <numeric/interpolation/spline/Bicubic_spline.hh>
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/MathVector_operations.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility/io/izstream.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
// AUTO-REMOVED #include <core/scoring/orbitals/OrbitalsScore.hh>
// AUTO-REMOVED #include <core/scoring/orbitals/OrbitalsScoreCreator.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/OrbitalTypeMapper.hh>
#include <basic/database/open.hh>


namespace core{
namespace scoring{
namespace orbitals{

//this function reads in the potentials from database. The potentials are already in the form of Energy.
OrbitalsLookup::OrbitalsLookup( utility::vector1< std::string > const & filename, utility::vector1< std::string > const & cubic_files) :
	number_stats_(core::chemical::orbitals::num_orbital_types),
	number_elements_(1000)
{

	//set up initial vectors. Data read from the database will be placed into these vectors
	utility::vector1< core::Real > E_vector(number_elements_, 0);
	utility::vector1< utility::vector1< core::Real > > HPOL_sc_H_sc_orb_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > HPOL_bb_H_sc_orb_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > HARO_sc_H_sc_orb_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > HPOL_sc_H_bb_orb_vector(static_cast <core::Size> (number_stats_), E_vector);






	std::string line;

	utility::io::izstream stream;
	basic::database::open( stream, filename[1] );


	//read in database file. This is the hpol file, put files into the vectors constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> HPOL_sc_H_sc_orb_vector[1][count] >> HPOL_sc_H_sc_orb_vector[2][count] >> HPOL_sc_H_sc_orb_vector[3][count] >> HPOL_sc_H_sc_orb_vector[4][count] >> HPOL_sc_H_sc_orb_vector[5][count] >> HPOL_sc_H_sc_orb_vector[6][count] >> HPOL_sc_H_sc_orb_vector[7][count];
		//std::cout << "printing out numbers " << HPOL_sc_H_sc_orb_vector[1][count] << HPOL_sc_H_sc_orb_vector[2][count] << HPOL_sc_H_sc_orb_vector[3][count] << HPOL_sc_H_sc_orb_vector[4][count] << HPOL_sc_H_sc_orb_vector[5][count] << HPOL_sc_H_sc_orb_vector[6][count] << HPOL_sc_H_sc_orb_vector[7][count];
	}
	stream.close();

	//orb_bb
	basic::database::open( stream, filename[3] );
	//read in database file. This is the haro file, put files into the array constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> HPOL_sc_H_bb_orb_vector[1][count] >> HPOL_sc_H_bb_orb_vector[2][count] >> HPOL_sc_H_bb_orb_vector[3][count] >> HPOL_sc_H_bb_orb_vector[4][count] >> HPOL_sc_H_bb_orb_vector[5][count] >> HPOL_sc_H_bb_orb_vector[6][count] >> HPOL_sc_H_bb_orb_vector[7][count];
		//std::cout << orb_e1[count] << " <- orb_e1, orb_e2->" << orb_e2[count] << std::endl;
	}
	stream.close();

	//orb_bb
	basic::database::open( stream, filename[4] );
	//read in database file. This is the haro file, put files into the array constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> HPOL_bb_H_sc_orb_vector[1][count] >> HPOL_bb_H_sc_orb_vector[2][count] >> HPOL_bb_H_sc_orb_vector[3][count] >> HPOL_bb_H_sc_orb_vector[4][count] >> HPOL_bb_H_sc_orb_vector[5][count] >> HPOL_bb_H_sc_orb_vector[6][count] >> HPOL_bb_H_sc_orb_vector[7][count];
		//std::cout << orb_e1[count] << " <- orb_e1, orb_e2->" << orb_e2[count] << std::endl;
	}
	stream.close();


	//haro data
	basic::database::open( stream, filename[2] );
	//read in database file. This is the haro file, put files into the array constructed above
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> HARO_sc_H_sc_orb_vector[1][count] >> HARO_sc_H_sc_orb_vector[2][count] >> HARO_sc_H_sc_orb_vector[3][count] >> HARO_sc_H_sc_orb_vector[4][count] >> HARO_sc_H_sc_orb_vector[5][count] >> HARO_sc_H_sc_orb_vector[6][count] >> HARO_sc_H_sc_orb_vector[7][count];
	}
	stream.close();



	//initial construction of a vector of MathMatrix. We will be pushing back matrixes into this
	utility::vector1< numeric::MathMatrix<core::Real> > HPOL_sc_H_sc_orb_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > HPOL_bb_H_sc_orb_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > HPOL_sc_H_bb_orb_vector_matrix;


	utility::vector1< numeric::MathMatrix<core::Real> > HARO_sc_H_sc_orb_vector_matrix;

	for( core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count ) {
		//MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
		numeric::MathMatrix<core::Real> HPOL_sc_H_sc_orb_matrix(50, 20, & HPOL_sc_H_sc_orb_vector[count][1] );
		numeric::MathMatrix<core::Real> HPOL_bb_H_sc_orb_matrix(50, 20, & HPOL_bb_H_sc_orb_vector[count][1] );

		numeric::MathMatrix<core::Real> HARO_sc_H_sc_orb_matrix(50, 20, & HARO_sc_H_sc_orb_vector[count][1]);
		numeric::MathMatrix<core::Real> HPOL_sc_H_bb_orb_matrix(50, 20, & HPOL_sc_H_bb_orb_vector[count][1]);

		HPOL_sc_H_sc_orb_vector_matrix.push_back(HPOL_sc_H_sc_orb_matrix);
		HPOL_bb_H_sc_orb_vector_matrix.push_back(HPOL_bb_H_sc_orb_matrix);

		HARO_sc_H_sc_orb_vector_matrix.push_back(HARO_sc_H_sc_orb_matrix);
		HPOL_sc_H_bb_orb_vector_matrix.push_back(HPOL_sc_H_bb_orb_matrix);
	}



	//the spline will behave as a first derivative. This means that the borders of the spline will be smoothed out if the first_derivative
	//is set to 0. This is controlled through an enum and put into an array.
	numeric::interpolation::spline::BorderFlag behavior[2] = {
		numeric::interpolation::spline::e_FirstDer, //
		numeric::interpolation::spline::e_Periodic };//periodic because it's an angle

	core::Real const start[2] = {0.05, -0.95}; //starting position of values in the matrix.0
	core::Real const delta[2] = {0.1,  0.1  }; //change in values. Needed to increment the values
	bool const linear_cont[2] = {true, false }; //if outside of range continue linearly


	////////////////////////////////////////////
	////////////////////////////////////////////
	////////////////////////////////////////////



	//what is the characteristic of the spline at 0 values. Should be smoothed out at edges. Powerful if used with e_FristDer
	std::pair< core::Real, core::Real > const first_deriv[2] = {
		std::pair< core::Real, core::Real >(0.0,0.0),
		std::pair< core::Real, core::Real >(0.0,0.0) };


	numeric::interpolation::spline::BicubicSpline bicubic_spline;// apl -- fixing memory leak = new numeric::interpolation::spline::BicubicSpline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> HPOL_sc_H_sc_orb_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> HPOL_bb_H_sc_orb_vector_spline;

	utility::vector1<numeric::interpolation::spline::BicubicSpline> HARO_sc_H_sc_orb_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> HPOL_sc_H_bb_orb_vector_spline;

	//HPOL_sc_H_sc_orb_vector_spline.reserve( number_stats_ );
	//HPOL_bb_H_sc_orb_vector_spline.reserve( number_stats_ );
	//HARO_sc_H_sc_orb_vector_spline.reserve( number_stats_ );
	//HPOL_sc_H_bb_orb_vector_spline.reserve( number_stats_ );

	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		bicubic_spline.train(behavior, start, delta, HPOL_sc_H_sc_orb_vector_matrix[count], linear_cont, first_deriv);
		HPOL_sc_H_sc_orb_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, HPOL_bb_H_sc_orb_vector_matrix[count], linear_cont, first_deriv);
		HPOL_bb_H_sc_orb_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, HARO_sc_H_sc_orb_vector_matrix[count], linear_cont, first_deriv);
		HARO_sc_H_sc_orb_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, HPOL_sc_H_bb_orb_vector_matrix[count], linear_cont, first_deriv);
		HPOL_sc_H_bb_orb_vector_spline.push_back( bicubic_spline );

	}

	//store all splines in a private member look up table. This allows for looking up the values in the actual score and provides an
	//interface for the OrbitalsScore class.
	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		HPOL_sc_H_sc_orb_splines_.push_back(HPOL_sc_H_sc_orb_vector_spline[count]);
		HPOL_bb_H_sc_orb_splines_.push_back(HPOL_bb_H_sc_orb_vector_spline[count]);
		HPOL_sc_H_bb_orb_splines_.push_back(HPOL_sc_H_bb_orb_vector_spline[count]);
		HARO_sc_H_sc_orb_splines_.push_back(HARO_sc_H_sc_orb_vector_spline[count]);
	}

/*	HPOL_sc_H_sc_orb_splines_.resize(number_stats_);
	haro_spline_vector_.resize(number_stats_);
	for(int count=core::chemical::orbitals::C_pi_sp2; count <= static_cast <core::Size> (number_stats_); ++count){
		hpol_splines_.insert(std::pair<core::Size, numeric::interpolation::spline::BicubicSpline>(count, HPOL_sc_H_sc_orb_splines_[count] ) );
		hpol_sc_orb_bb_splines_.insert(std::pair<core::Size, numeric::interpolation::spline::BicubicSpline>(count, HPOL_bb_H_sc_orb_splines_[count] ) );
		haro_splines_.insert(std::pair<core::Size, numeric::interpolation::spline::BicubicSpline>(count, HARO_sc_H_sc_orb_splines_[count] ) );
		orb_bb_splines_.insert(std::pair<core::Size, numeric::interpolation::spline::BicubicSpline>(count, HPOL_sc_H_bb_orb_splines_[count] ) );



		HPOL_sc_H_sc_orb_splines_[count] = HPOL_sc_H_sc_orb_splines_[count];
		haro_spline_vector_[count] = HARO_sc_H_sc_orb_splines_[count];


	}*/



 //some verbose checking. Not needed

/*
	for(core::Real i=0.05; i <= 50; i+=0.1){
	core::Size number=0;
			for(core::Real j=-0.95; j<=1; j+=0.1){
				++number;
				if( number == 20){
					std::cout << HPOL_sc_H_sc_orb_vector_spline[2].F((numeric::MakeVector(i,j))) << std::endl;

				}else{
					std::cout << HPOL_sc_H_sc_orb_vector_spline[2].F((numeric::MakeVector(i,j))) << " ";
					//std::cout << i << " " << j << std::endl;
				}
			}
		}

*/


	////////////////////////////////////////
	////////////////////////////////////////
	////////////////////////////////////////
	////////////////////////////////////////
	////////////////////////////////////////
	////////////////////////////////////////
	////////////////////////////////////////



	//start the cubic splines
	utility::vector1< utility::vector1< core::Real > > cubic_HPOL_sc_H_sc_orb_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > cubic_HPOL_sc_H_bb_orb_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > cubic_HPOL_bb_H_sc_orb_vector(static_cast <core::Size> (number_stats_), E_vector);
	utility::vector1< utility::vector1< core::Real > > cubic_HARO_sc_H_sc_orb_vector(static_cast <core::Size> (number_stats_), E_vector);



	//cubic hpol sc_orb_sc_h
	basic::database::open( stream, cubic_files[1] );
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> cubic_HPOL_sc_H_sc_orb_vector[1][count] >> cubic_HPOL_sc_H_sc_orb_vector[2][count] >> cubic_HPOL_sc_H_sc_orb_vector[3][count] >> cubic_HPOL_sc_H_sc_orb_vector[4][count] >> cubic_HPOL_sc_H_sc_orb_vector[5][count] >> cubic_HPOL_sc_H_sc_orb_vector[6][count] >> cubic_HPOL_sc_H_sc_orb_vector[7][count];
		//std::cout << orb_e1[count] << " <- orb_e1, orb_e2->" << orb_e2[count] << std::endl;
	}
	stream.close();

	basic::database::open( stream, cubic_files[2] );
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> cubic_HPOL_sc_H_bb_orb_vector[1][count] >> cubic_HPOL_sc_H_bb_orb_vector[2][count] >> cubic_HPOL_sc_H_bb_orb_vector[3][count] >> cubic_HPOL_sc_H_bb_orb_vector[4][count] >> cubic_HPOL_sc_H_bb_orb_vector[5][count] >> cubic_HPOL_sc_H_bb_orb_vector[6][count] >> cubic_HPOL_sc_H_bb_orb_vector[7][count];
		//std::cout << orb_e1[count] << " <- orb_e1, orb_e2->" << orb_e2[count] << std::endl;
	}
	stream.close();


	basic::database::open( stream, cubic_files[3] );
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> cubic_HPOL_bb_H_sc_orb_vector[1][count] >> cubic_HPOL_bb_H_sc_orb_vector[2][count] >> cubic_HPOL_bb_H_sc_orb_vector[3][count] >> cubic_HPOL_bb_H_sc_orb_vector[4][count] >> cubic_HPOL_bb_H_sc_orb_vector[5][count] >> cubic_HPOL_bb_H_sc_orb_vector[6][count] >> cubic_HPOL_bb_H_sc_orb_vector[7][count];
		//std::cout << orb_e1[count] << " <- orb_e1, orb_e2->" << orb_e2[count] << std::endl;
	}
	stream.close();

	basic::database::open( stream, cubic_files[4] );
	for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
		std::istringstream l( line );
		l >> cubic_HARO_sc_H_sc_orb_vector[1][count] >> cubic_HARO_sc_H_sc_orb_vector[2][count] >> cubic_HARO_sc_H_sc_orb_vector[3][count] >> cubic_HARO_sc_H_sc_orb_vector[4][count] >> cubic_HARO_sc_H_sc_orb_vector[5][count] >> cubic_HARO_sc_H_sc_orb_vector[6][count] >> cubic_HARO_sc_H_sc_orb_vector[7][count];
		//std::cout << orb_e1[count] << " <- orb_e1, orb_e2->" << orb_e2[count] << std::endl;
	}
	stream.close();

	//initial construction of a vector of MathMatrix. We will be pushing back matrixes into this
	utility::vector1< numeric::MathVector<core::Real> > cubic_HPOL_sc_H_sc_orb_vector_matrix;
	utility::vector1< numeric::MathVector<core::Real> > cubic_HPOL_sc_H_bb_orb_vector_matrix;
	utility::vector1< numeric::MathVector<core::Real> > cubic_HPOL_bb_H_sc_orb_vector_matrix;
	utility::vector1< numeric::MathVector<core::Real> > cubic_HARO_sc_H_sc_orb_vector_matrix;

	for( core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count ) {
		//MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
		numeric::MathVector<core::Real> cubic_HPOL_sc_H_sc_orb_matrix(20, & cubic_HPOL_sc_H_sc_orb_vector[count][1] );
		numeric::MathVector<core::Real> cubic_HPOL_sc_H_bb_orb_matrix(20, & cubic_HPOL_sc_H_bb_orb_vector[count][1] );
		numeric::MathVector<core::Real> cubic_HPOL_bb_H_sc_orb_matrix(20, & cubic_HPOL_bb_H_sc_orb_vector[count][1] );
		numeric::MathVector<core::Real> cubic_HARO_sc_H_sc_orb_matrix(20, & cubic_HARO_sc_H_sc_orb_vector[count][1] );

		cubic_HPOL_sc_H_sc_orb_vector_matrix.push_back(cubic_HPOL_sc_H_sc_orb_matrix);
		cubic_HPOL_sc_H_bb_orb_vector_matrix.push_back(cubic_HPOL_sc_H_bb_orb_matrix);
		cubic_HPOL_bb_H_sc_orb_vector_matrix.push_back(cubic_HPOL_bb_H_sc_orb_matrix);
		cubic_HARO_sc_H_sc_orb_vector_matrix.push_back(cubic_HARO_sc_H_sc_orb_matrix);
	}


	//the spline will behave as a first derivative. This means that the borders of the spline will be smoothed out if the first_derivative
	//is set to 0. This is controlled through an enum and put into an array.
	numeric::interpolation::spline::BorderFlag cubic_behavior = numeric::interpolation::spline::e_Periodic; //periodic because it's an angle

	core::Real const cubic_start( -0.95); //starting position of values in the matrix.0
	core::Real const cubic_delta(0.1  ); //change in values. Needed to increment the values
//	bool const linear_cont[2] = {true, false }; //if outside of range continue linearly

	std::pair< core::Real, core::Real > const cubic_first_deriv(std::pair< core::Real, core::Real >(0.0,0.0));






	numeric::interpolation::spline::CubicSpline cubic_spline;
	utility::vector1<numeric::interpolation::spline::CubicSpline> cubic_HPOL_sc_H_sc_orb_vector_spline;
	utility::vector1<numeric::interpolation::spline::CubicSpline> cubic_HPOL_sc_H_bb_orb_vector_spline;
	utility::vector1<numeric::interpolation::spline::CubicSpline> cubic_HPOL_bb_H_sc_orb_vector_spline;
	utility::vector1<numeric::interpolation::spline::CubicSpline> cubic_HARO_sc_H_sc_orb_vector_spline;


	cubic_HPOL_sc_H_sc_orb_vector_spline.reserve( number_stats_ );

	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		cubic_spline.train(cubic_behavior, cubic_start, cubic_delta, cubic_HPOL_sc_H_sc_orb_vector_matrix[count], cubic_first_deriv);
		cubic_HPOL_sc_H_sc_orb_vector_spline.push_back( cubic_spline );

		cubic_spline.train(cubic_behavior, cubic_start, cubic_delta, cubic_HPOL_sc_H_bb_orb_vector_matrix[count], cubic_first_deriv);
		cubic_HPOL_sc_H_bb_orb_vector_spline.push_back( cubic_spline );

		cubic_spline.train(cubic_behavior, cubic_start, cubic_delta, cubic_HPOL_bb_H_sc_orb_vector_matrix[count], cubic_first_deriv);
		cubic_HPOL_bb_H_sc_orb_vector_spline.push_back( cubic_spline );

		cubic_spline.train(cubic_behavior, cubic_start, cubic_delta, cubic_HARO_sc_H_sc_orb_vector_matrix[count], cubic_first_deriv);
		cubic_HARO_sc_H_sc_orb_vector_spline.push_back( cubic_spline );
	}



	for(core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count){
		cubic_HPOL_sc_H_sc_orb_splines_.push_back(cubic_HPOL_sc_H_sc_orb_vector_spline[count]);
		cubic_HPOL_sc_H_bb_orb_splines_.push_back(cubic_HPOL_sc_H_bb_orb_vector_spline[count]);
		cubic_HPOL_bb_H_sc_orb_splines_.push_back(cubic_HPOL_bb_H_sc_orb_vector_spline[count]);
		cubic_HARO_sc_H_sc_orb_splines_.push_back(cubic_HARO_sc_H_sc_orb_vector_spline[count]);
	}




/*	core::Size number=0;
			for(core::Real j=-0.95; j<=1; j+=0.1){
				++number;
				std::cout << cubic_HARO_sc_H_sc_orb_splines_[1].F(j) << " ";
					//std::cout << i << " " << j << std::endl;

			}*/


}//end orbitalsLookup

//get the energy from the bicubic spline. This is done by using the mm atom name as a key in the map. We pass in the distance and angle.
void OrbitalsLookup::get_dist_AOH_angle_energy
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
	if ( distance > 5.05 ) { energy = distance_derivative = angle_derivative = 0.0; return; }

	numeric::MathVector<core::Real> dist_angle_pair(numeric::MakeVector(distance, AOH_angle));
	if(check_derivative){
		if(h_enum == hpol){
			const numeric::interpolation::spline::BicubicSpline &spline( HPOL_sc_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			//spline = hpol_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);

		}
		if(h_enum == sc_orb_bb_H){
			const numeric::interpolation::spline::BicubicSpline &spline(HPOL_bb_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
		if(h_enum == bb){
			const numeric::interpolation::spline::BicubicSpline &spline(HPOL_sc_H_bb_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
		if(h_enum == haro){
			//spline = haro_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			const numeric::interpolation::spline::BicubicSpline &spline(HARO_sc_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
		//energy = spline.F(numeric::MakeVector(distance, AOH_angle));
	}else{
		if(h_enum == hpol){
			const numeric::interpolation::spline::BicubicSpline &spline(HPOL_sc_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			//spline = hpol_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			energy = spline.F(dist_angle_pair);
		}
		if(h_enum == sc_orb_bb_H){
			const numeric::interpolation::spline::BicubicSpline &spline(HPOL_bb_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(dist_angle_pair);
		}
		if(h_enum == haro){
			//spline = haro_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			const numeric::interpolation::spline::BicubicSpline &spline(HARO_sc_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);

			energy = spline.F(dist_angle_pair);
		}
		if(h_enum == bb){
			const numeric::interpolation::spline::BicubicSpline &spline( HPOL_sc_H_bb_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(dist_angle_pair);
		}
	}
}

void OrbitalsLookup::get_DHO_angle_energy(
		const h_type h_enum,
		const core::chemical::orbitals::orbital_type_enum orb_type_name,
		const core::Real DHO_angle,
		core::Real & energy,
		core::Real & angle_derivative,
		bool check_derivative
) const{
	if(h_enum==haro){
		if(DHO_angle > -0.05){
			energy=0.0; return;
		}
	}
	if(h_enum==hpol ){
			if(orb_type_name == core::chemical::orbitals::C_pi_sp2){
				if(DHO_angle > 0.50){
					energy=0.0; return;
				}
			}
			if(orb_type_name == core::chemical::orbitals::N_pi_sp2){
				if(DHO_angle > 0.05){
					energy=0.0; return;
				}
			}
			if(orb_type_name == core::chemical::orbitals::N_p_sp2){
				energy=0.0; return;
			}
			if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
				if(DHO_angle > -0.15){
					energy=0.0; return;
				}
			}
			if(orb_type_name == core::chemical::orbitals::O_p_sp2){
				if(DHO_angle > -0.25){
					energy=0.0; return;
				}
			}
			if(orb_type_name == core::chemical::orbitals::O_p_sp3){
				if(DHO_angle > -.05){
					energy=0.0; return;
				}
			}
			if(orb_type_name == core::chemical::orbitals::S_p_sp3){
				energy=0.0; return;
			}
		}
	if(h_enum==sc_orb_bb_H){
		if(orb_type_name == core::chemical::orbitals::C_pi_sp2){
			if(DHO_angle > 0.65 || DHO_angle < -0.35){
				energy=0.0; return;
			}
		}
		if(orb_type_name == core::chemical::orbitals::N_pi_sp2){
			energy=0.0; return;
		}
		if(orb_type_name == core::chemical::orbitals::N_p_sp2){
			energy=0.0; return;
		}
		if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
			if(DHO_angle > -0.45){
				energy=0.0; return;
			}
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp2){
			if(DHO_angle > -0.55){
				energy=0.0; return;
			}
		}
/*		if(orb_type_name == core::chemical::orbitals::O_p_sp3){
			if(AOH_angle > .25){
				return false;
			}else return true;
		}*/
		if(orb_type_name == core::chemical::orbitals::S_p_sp3){
			energy=0.0; return;
		}
	}
	if(h_enum==bb){
		if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
			if(DHO_angle > -0.05){
				energy=0.0; return;
			}
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp2){
			if(DHO_angle > -0.35){
				energy=0.0; return;
			}
		}
	}





//	numeric::interpolation::spline::CubicSpline spline;
	if(check_derivative){
		if(h_enum == hpol){
			const numeric::interpolation::spline::CubicSpline &spline(cubic_HPOL_sc_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			angle_derivative = spline.dF(DHO_angle);

		}
		if(h_enum == sc_orb_bb_H){
			const numeric::interpolation::spline::CubicSpline &spline(cubic_HPOL_bb_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			angle_derivative = spline.dF(DHO_angle);
		}
		if(h_enum == bb){
			const numeric::interpolation::spline::CubicSpline &spline(cubic_HPOL_sc_H_bb_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			angle_derivative = spline.dF(DHO_angle);
		}
		if(h_enum == haro){
			//spline = haro_splines_.find(static_cast <core::Size>(orb_type_name))->second;
			const numeric::interpolation::spline::CubicSpline &spline(cubic_HARO_sc_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			angle_derivative = spline.dF(DHO_angle);
		}
		//energy = spline.F(numeric::MakeVector(distance, AOH_angle));
	}else{
		if(h_enum == hpol){
			const numeric::interpolation::spline::CubicSpline &spline(cubic_HPOL_sc_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(DHO_angle);
		}
		if(h_enum == sc_orb_bb_H){
			const numeric::interpolation::spline::CubicSpline &spline(cubic_HPOL_bb_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(DHO_angle);
		}
		if(h_enum == haro){
			const numeric::interpolation::spline::CubicSpline &spline(cubic_HARO_sc_H_sc_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(DHO_angle);
		}
		if(h_enum == bb){
			const numeric::interpolation::spline::CubicSpline &spline(cubic_HPOL_sc_H_bb_orb_splines_[static_cast <core::Size>(orb_type_name)]);
			energy = spline.F(DHO_angle);
		}
	}

}

bool OrbitalsLookup::check_distance(
	const core::Real distance,
	const h_type h_enum,
	const core::chemical::orbitals::orbital_type_enum orb_type_name
) const
{

	//a set of rules to stop hitting the bicubic spline so hard. This set of rules speeds up the overall scoring function
	if(h_enum==haro){
		if(distance > 10.89){
			return false;
		}else return true;
	}
	//sc orb to sc h rules!
	if(h_enum==hpol){
		if(orb_type_name == core::chemical::orbitals::C_pi_sp2){
			if(distance > 10.89 || distance < 2.4025){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::N_pi_sp2){
			if(distance > 9.61 || distance < 4.2025){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::N_p_sp2){
			return false;
		}
		if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
			if(distance > 3.0625 || distance < 1.5625){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp2){
			if(distance > 3.0625 || distance < 0.7225){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp3){
			if(distance > 2.25 || distance < 0.9025){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::S_p_sp3){
			return false;
		}
	}
	if(h_enum==sc_orb_bb_H){
		if(orb_type_name == core::chemical::orbitals::C_pi_sp2){
			if(distance > 12.6025 || distance < 3.0625){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::N_pi_sp2){
			return false;
		}
		if(orb_type_name == core::chemical::orbitals::N_p_sp2){
			return false;
		}
		if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
			if(distance > 3.4225 || distance < 1.5625){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp2){
			if(distance > 3.0625 || distance < 1.1025){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp3){
			if(distance > 5.0625 || distance < 1.3225){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::S_p_sp3){
			return false;
		}
	}
	if(h_enum==bb){
		if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
			if(distance > 3.8025 || distance < 1.8225){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp2){
			if(distance > 3.4225 || distance < 1.1025){
				return false;
			}else return true;
		}
	}
}




bool OrbitalsLookup::check_AOH_angle(
	const core::Real AOH_angle,
	const h_type h_enum,
	const core::chemical::orbitals::orbital_type_enum orb_type_name
) const{
	if(h_enum==haro){
		if(AOH_angle > -0.40){
			return false;
		}else return true;
	}
	//sc orb to sc h rules!
	if(h_enum==hpol ){
		if(orb_type_name == core::chemical::orbitals::C_pi_sp2){
			if(AOH_angle > -0.70){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::N_pi_sp2){
			if(AOH_angle > -0.85){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::N_p_sp2){
			return false;
		}
		if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
			if(AOH_angle > -0.35){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp2){
			if(AOH_angle > -0.45){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp3){
			if(AOH_angle > -.35){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::S_p_sp3){
			return false;
		}
	}
	if(h_enum==sc_orb_bb_H){
		if(orb_type_name == core::chemical::orbitals::C_pi_sp2){
			if(AOH_angle > -0.70){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::N_pi_sp2){
			return false;
		}
		if(orb_type_name == core::chemical::orbitals::N_p_sp2){
			return false;
		}
		if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
			if(AOH_angle > -0.35){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp2){
			if(AOH_angle > -0.35){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp3){
			if(AOH_angle > .25){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::S_p_sp3){
			return false;
		}
	}
	if(h_enum==bb){
		if(orb_type_name == core::chemical::orbitals::O_pi_sp2){
			if(AOH_angle > -0.35){
				return false;
			}else return true;
		}
		if(orb_type_name == core::chemical::orbitals::O_p_sp2){
			if(AOH_angle > -0.45){
				return false;
			}else return true;
		}
	}
}


}//namespace orbitals

}//namespace scoring

}//namespace core
