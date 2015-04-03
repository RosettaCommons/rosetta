// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//#include <core/chemical/rna/util.hh>

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
#include <map>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/chemical/orbitals/OrbitalTypeMapper.hh>


namespace core{
namespace scoring{
namespace orbitals{

//this function reads in the potentials from database. The potentials are already in the form of Energy.
OrbitalsLookup::OrbitalsLookup(
		utility::vector1< std::string > const & DHO_energies,
		utility::vector1< std::string > const & AOH_energies,
		utility::vector1< std::string > const & AOD_orb_orb_energies,
		utility::vector1< std::string > const & DOA_orb_orb_energies,
		utility::vector1< std::string > const & ACO_AOH_orb_Hpol_energies
) :
	number_stats_(6),
	number_elements_(600)
{

	if(basic::options::option[ basic::options::OptionKeys::in::add_orbitals] != true ){
		utility_exit_with_message( "Trying to run orbitals score without orbitals! Pass the flag -add_orbitals!" );
	}


	std::map<core::Size, std::pair<core::Size, core::Size> > DHO_Hpol_scOrbscH_map;
	utility::vector1< utility::vector1< core::Real > > DHO_Hpol_scOrbH_vector = parse_files(DHO_energies[1], DHO_Hpol_scOrbscH_map);
	std::map<core::Size, std::pair<core::Size, core::Size> > DHO_bbOrbscH_map;
	utility::vector1< utility::vector1< core::Real > > DHO_Hpol_bbOrbH_vector = parse_files(DHO_energies[2], DHO_bbOrbscH_map);
	std::map<core::Size, std::pair<core::Size, core::Size> > DHO_HARO_scOrbscH_map;
	utility::vector1< utility::vector1< core::Real > > DHO_Haro_scOrbH_vector = parse_files(DHO_energies[3], DHO_HARO_scOrbscH_map);

	//initial construction of a vector of MathMatrix. We will be pushing back matrixes into this
	utility::vector1< numeric::MathMatrix<core::Real> > DHO_Hpol_scOrbH_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > DHO_Hpol_bbOrbH_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > DHO_Haro_scOrbH_vector_matrix;

	for( core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count ) {
		//MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
		numeric::MathMatrix<core::Real> Hpol_scOrbH_matrix(DHO_Hpol_scOrbscH_map[count].second, DHO_Hpol_scOrbscH_map[count].first, & DHO_Hpol_scOrbH_vector[count][1] );
		//std::cout << "second " << DHO_scOrbscH_map[count].second << "first " << DHO_scOrbscH_map[count].first << std::endl;
		numeric::MathMatrix<core::Real> Hpol_bbOrbH_matrix(DHO_bbOrbscH_map[count].second, DHO_bbOrbscH_map[count].first, & DHO_Hpol_bbOrbH_vector[count][1] );
		numeric::MathMatrix<core::Real> Haro_scOrbH_matrix(DHO_HARO_scOrbscH_map[count].second, DHO_HARO_scOrbscH_map[count].first, & DHO_Haro_scOrbH_vector[count][1]);

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
	for(core::Real i=0.00; i <= 3; i+=0.1){
		core::Size number=0;
		for(core::Real j=-1.0; j<=0; j+=0.05){
			++number;
			if( number == 20){
				std::cout << DHO_Hpol_scOrbH_vector_spline[1].F((numeric::MakeVector(i,j))) << std::endl;

			}else{
				std::cout << DHO_Hpol_scOrbH_vector_spline[1].F((numeric::MakeVector(i,j))) << " ";
				//std::cout << i << " " << j << std::endl;
			}
		}
	}
	std::cout << "###################################################################################\n#####################3" << std::endl;
*/


/*

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


	std::map<core::Size, std::pair<core::Size, core::Size> > AOH_Hpol_scOrbscH_map;
	utility::vector1< utility::vector1< core::Real > > AOH_Hpol_scOrbH_vector = parse_files(AOH_energies[1], AOH_Hpol_scOrbscH_map);
	std::map<core::Size, std::pair<core::Size, core::Size> > AOH_bbOrbscH_map;
	utility::vector1< utility::vector1< core::Real > > AOH_Hpol_bbOrbH_vector = parse_files(AOH_energies[2], AOH_bbOrbscH_map);
	std::map<core::Size, std::pair<core::Size, core::Size> > AOH_HARO_scOrbscH_map;
	utility::vector1< utility::vector1< core::Real > > AOH_Haro_scOrbH_vector = parse_files(AOH_energies[3], AOH_HARO_scOrbscH_map);


	//initial construction of a vector of MathMatrix. We will be pushing back matrixes into this
	utility::vector1< numeric::MathMatrix<core::Real> > AOH_Hpol_scOrbH_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > AOH_Hpol_bbOrbH_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > AOH_Haro_scOrbH_vector_matrix;

	for( core::Size count=1; count <= static_cast< core::Size > (number_stats_); ++count ) {
		//MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
		numeric::MathMatrix<core::Real> Hpol_scOrbH_matrix(AOH_Hpol_scOrbscH_map[count].second, AOH_Hpol_scOrbscH_map[count].first, & AOH_Hpol_scOrbH_vector[count][1] );
		numeric::MathMatrix<core::Real> Hpol_bbOrbH_matrix(AOH_bbOrbscH_map[count].second, AOH_bbOrbscH_map[count].first, & AOH_Hpol_bbOrbH_vector[count][1] );
		numeric::MathMatrix<core::Real> Haro_scOrbH_matrix(AOH_HARO_scOrbscH_map[count].second, AOH_HARO_scOrbscH_map[count].first, & AOH_Haro_scOrbH_vector[count][1]);

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


	//std::cout << "###################################################################################\n#####################3" << std::endl;
/*	for(core::Real i=0.00; i <= 4; i+=0.1){
			core::Size number=0;
			for(core::Real j=-1.0; j<=1; j+=0.05){
				++number;
				if( number == 80){
					std::cout << i << " " << j << ", " << std::endl;

				}else{
					std::cout <<  i << " " << j << ", ";
					//std::cout << i << " " << j << std::endl;
				}
			}
		}


	std::cout << "########################################################################################################3" << std::endl;

	for(core::Real i=0.00; i <= 4; i+=0.01){
		core::Size number=0;
		for(core::Real j=-1.0; j<=1; j+=0.025){
			++number;
			if( number == 80){
				std::cout << AOH_Hpol_scOrbH_splines_[6].F((numeric::MakeVector(i,j))) <<  "\n" << std::endl;
			}else{
				std::cout << AOH_Hpol_scOrbH_splines_[6].F((numeric::MakeVector(i,j))) << " ";
				//std::cout << i << " " << j << std::endl;
			}
		}
	}*/
/*
	std::cout << "########################################################################################################3" << std::endl;
	std::cout << "########################################################################################################3" << std::endl;

	std::cout << "table ";
	core::Size number=0;
	for(core::Real j=-1.0; j<=1; j+=0.025){
		++number;
		if( number == 80){
			std::cout << j << std::endl;
		}else{
			std::cout << j << " ";
		}
	}
	for(core::Real i=0.00; i <= 4; i+=0.05){
		core::Size number=0;
		std::cout << i << " ";
		for(core::Real j=-1.0; j<=1; j+=0.025){
			++number;
			if( number == 80){
				std::cout << AOH_Hpol_scOrbH_splines_[1].F((numeric::MakeVector(i,j))) << std::endl;
			}else{
				std::cout << AOH_Hpol_scOrbH_splines_[1].F((numeric::MakeVector(i,j))) << " ";
				//std::cout << i << " " << j << std::endl;
			}
		}
	}


*/


	std::map<core::Size, std::pair<core::Size, core::Size> > AOD_orb_orb_map;
	utility::vector1< utility::vector1< core::Real > > AOD_orb_orb_vector = parse_files(AOD_orb_orb_energies[1], AOD_orb_orb_map);

	std::map<core::Size, std::pair<core::Size, core::Size> > DOA_orb_orb_map;
	utility::vector1< utility::vector1< core::Real > > DOA_orb_orb_vector = parse_files(DOA_orb_orb_energies[1], DOA_orb_orb_map);

	utility::vector1< numeric::MathMatrix<core::Real> > AOD_orb_orb_vector_matrix;
	utility::vector1< numeric::MathMatrix<core::Real> > DOA_orb_orb_vector_matrix;

	for( core::Size count=1; count <= 5; ++count ) {
		//MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
		numeric::MathMatrix<core::Real> AOD_orb_orb_matrix(AOD_orb_orb_map[count].second, AOD_orb_orb_map[count].first, & AOD_orb_orb_vector[count][1] );
		numeric::MathMatrix<core::Real> DOA_orb_orb_matrix(DOA_orb_orb_map[count].second, DOA_orb_orb_map[count].first, & DOA_orb_orb_vector[count][1] );

		AOD_orb_orb_vector_matrix.push_back(AOD_orb_orb_matrix);
		DOA_orb_orb_vector_matrix.push_back(DOA_orb_orb_matrix);
	}

	utility::vector1<numeric::interpolation::spline::BicubicSpline> AOD_orb_orb_vector_spline;
	utility::vector1<numeric::interpolation::spline::BicubicSpline> DOA_orb_orb_vector_spline;

	for(core::Size count=1; count <= 5; ++count){
		bicubic_spline.train(behavior, start, delta, AOD_orb_orb_vector_matrix[count], linear_cont, first_deriv);
		AOD_orb_orb_vector_spline.push_back( bicubic_spline );

		bicubic_spline.train(behavior, start, delta, DOA_orb_orb_vector_matrix[count], linear_cont, first_deriv);
		DOA_orb_orb_vector_spline.push_back( bicubic_spline );
	}

	//store all splines in a private member look up table. This allows for looking up the values in the actual score and provides an
	//interface for the OrbitalsScore class.
	for(core::Size count=1; count <= 5; ++count){
		AOD_orb_orb_splines_.push_back(AOD_orb_orb_vector_spline[count]);
		DOA_orb_orb_splines_.push_back(DOA_orb_orb_vector_spline[count]);
	}


/*
	for(core::Real i=0.00; i <= 4; i+=0.1){
		core::Size number=0;
		for(core::Real j=-1.0; j<=0; j+=0.05){
			++number;
			if( number == 40){
				if(AOD_orb_orb_splines_[1].F((numeric::MakeVector(i,j))) > 0 || AOD_orb_orb_splines_[1].F((numeric::MakeVector(i,j))) > -1e-5){
					std::cout << "0" << std::endl;
				}else {
					std::cout << AOD_orb_orb_splines_[1].F((numeric::MakeVector(i,j))) << std::endl;
				}
			}else{
				if(AOD_orb_orb_splines_[1].F((numeric::MakeVector(i,j))) > 0 || AOD_orb_orb_splines_[1].F((numeric::MakeVector(i,j))) > -1e-5){
					std::cout << "0" << " ";
				}else {
					std::cout << AOD_orb_orb_splines_[1].F((numeric::MakeVector(i,j))) << " ";
				}
				//std::cout << i << " " << j << std::endl;
			}
		}
	}
	std::cout << "###################################################################################\n#####################3" << std::endl;
*/


	std::map<core::Size, std::pair<core::Size, core::Size> > ACO_AOH_orb_Hpol_map;
	utility::vector1< utility::vector1< core::Real > > ACO_AOH_orb_Hpol_vector = parse_files(ACO_AOH_orb_Hpol_energies[1], ACO_AOH_orb_Hpol_map);
	utility::vector1< numeric::MathMatrix<core::Real> > ACO_AOH_orb_Hpol_vector_matrix;


 //MathMatrix requires an array, not a vector. To get an array from a vector, we can use the & vector[1]. See wikipedia!
	numeric::MathMatrix<core::Real> ACO_AOH_orb_Hpol_matrix(ACO_AOH_orb_Hpol_map[1].second, ACO_AOH_orb_Hpol_map[1].first, & ACO_AOH_orb_Hpol_vector[1][1] );
	ACO_AOH_orb_Hpol_vector_matrix.push_back(ACO_AOH_orb_Hpol_matrix);

	utility::vector1<numeric::interpolation::spline::BicubicSpline> ACO_AOH_orb_Hpol_vector_spline;


		bicubic_spline.train(behavior, start, delta, ACO_AOH_orb_Hpol_vector_matrix[1], linear_cont, first_deriv);
		ACO_AOH_orb_Hpol_vector_spline.push_back( bicubic_spline );


	//store all splines in a private member look up table. This allows for looking up the values in the actual score and provides an
	//interface for the OrbitalsScore class.
		ACO_AOH_orb_Hpol_splines_.push_back(ACO_AOH_orb_Hpol_vector_spline[1]);


}//end orbitalsLookup


utility::vector1< utility::vector1< core::Real > > OrbitalsLookup::parse_files(
    std::string const & file,
    std::map<core::Size, std::pair<core::Size, core::Size> > & orbital_angle_dist_map
)const
{
  utility::vector1< core::Real > E_vector(number_elements_, 0);
  utility::vector1< utility::vector1< core::Real > > energy_vector(static_cast <core::Size> (number_stats_), E_vector); //600 default value for KBP, resized later
  std::string line;
  utility::io::izstream stream;
  basic::database::open( stream, file );
  core::Size orbital_type(0);
  core::Size overall_count(1);
  for( core::Size count=1; utility::io::getline(stream, line); ++count ) {
    utility::vector1< std::string > split_string = utility::string_split(line, '\t'); //file is tab-delimenated
    if(split_string[1]=="Orbital"){
      orbital_type = static_cast< core::Size > (core::chemical::orbitals::OrbitalTypeMapper::get_instance()->get_orbital_enum(split_string[2]));
      overall_count=1;
    }else if(split_string[1] == "Size"){
      core::Size angle_bins = utility::string2int(split_string[2]);
      core::Size dist_bins = utility::string2int(split_string[3]);
      std::pair<core::Size, core::Size> angle_dist(angle_bins, dist_bins);
      orbital_angle_dist_map.insert(std::pair<core::Size, std::pair<core::Size, core::Size> >(orbital_type, angle_dist));
      energy_vector[orbital_type].resize((angle_bins*dist_bins), 0);//need to resize. KBP have differing amount of elements
    }else{
      for(core::Size x=1; x<= orbital_angle_dist_map[orbital_type].first; ++x){//iterate through angles and put into vector
        energy_vector[orbital_type][overall_count] = utility::string2float(split_string[x]);
        ++overall_count;
      }
    }
  }
  stream.close();
  return energy_vector;
}


//get the energy from the bicubic spline. This is done by using the mm atom name as a key in the map. We pass in the distance and angle.
void OrbitalsLookup::OrbHdist_cosDHO_energy
(
	const h_type h_enum,
	const core::Size orb_type_name,
	const core::Real distance,
	const core::Real DHO_angle,
	core::Real & energy,
	core::Real & distance_derivative,
	core::Real & angle_derivative,
	bool check_derivative
) const
{

	//if ( (distance > 4.0 ) ) { energy = distance_derivative = angle_derivative = 0.0; return; }

	//if(orb_type_name == core::chemical::orbitals::N_pi_sp2 && h_enum==Hpol_scOrbH){energy = distance_derivative = angle_derivative = 0.0; return;}

	numeric::MathVector<core::Real> dist_angle_pair(numeric::MakeVector(distance, DHO_angle));
	if(check_derivative){
		if(h_enum == Hpol_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline( DHO_Hpol_scOrbH_splines_[orb_type_name]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
			if(orb_type_name==chemical::orbitals::C_pi_sp2){
				//this is to match the orbital orbital and orbital hpol weight. Set the weight of cation_pi interactions to that of orb orb
				distance_derivative = (distance_derivative*0.01);//legacy multiplication
				angle_derivative = (angle_derivative*.01); //legacy multiplication. happened when I converted orbs_orbs to pci_cation_pi and pci_pi_pi
			}

		}
		if(h_enum == Haro_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Haro_scOrbH_splines_[orb_type_name]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
		if(h_enum == Hpol_bbOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Hpol_bbOrbH_splines_[orb_type_name]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
	}else{
		if(h_enum == Hpol_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Hpol_scOrbH_splines_[orb_type_name]);
			energy = spline.F(dist_angle_pair);
			if(orb_type_name==chemical::orbitals::C_pi_sp2){
				//this is to match the orbital orbital and orbital hpol weight. Set the weight of cation_pi interactions to that of orb orb
				energy = (energy*0.01);
			}
		}
		if(h_enum == Haro_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Haro_scOrbH_splines_[orb_type_name]);
			energy = spline.F(dist_angle_pair);
		}
		if(h_enum == Hpol_bbOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(DHO_Hpol_bbOrbH_splines_[orb_type_name]);
			energy = spline.F(dist_angle_pair);
		}
	}
}


void OrbitalsLookup::OrbHdist_cosAOH_energy
(
	const h_type h_enum,
	const core::Size orb_type_name,
	const core::Real distance,
	const core::Real AOH_angle,
	core::Real & energy,
	core::Real & distance_derivative,
	core::Real & angle_derivative,
	bool check_derivative,
	bool ACO //action centor orbital?
) const
{
	//if ( distance > 4.0 ) { energy = distance_derivative = angle_derivative = 0.0; return; }
	//if(orb_type_name == core::chemical::orbitals::N_pi_sp2 && h_enum==Hpol_scOrbH){energy = distance_derivative = angle_derivative = 0.0; return;}
//	if(orb_type_name == core::chemical::orbitals::C_pi_sp2 && h_enum==Hpol_bbOrbH){energy = distance_derivative = angle_derivative = 0.0; return;}

	numeric::MathVector<core::Real> dist_angle_pair(numeric::MakeVector(distance, AOH_angle));
	if(check_derivative){
		if(h_enum == Hpol_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline( AOH_Hpol_scOrbH_splines_[orb_type_name]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
			if(orb_type_name==chemical::orbitals::C_pi_sp2){
				//this is to match the orbital orbital and orbital hpol weight. Set the weight of cation_pi interactions to that of orb orb
				distance_derivative = (distance_derivative*0.01);
				angle_derivative = (angle_derivative*0.01);
			}

		}
		if(h_enum == Haro_scOrbH){
			if(ACO){
				const numeric::interpolation::spline::BicubicSpline &spline(ACO_AOH_orb_Hpol_splines_[orb_type_name]);
				distance_derivative = spline.dFdx(dist_angle_pair);
				angle_derivative = spline.dFdy(dist_angle_pair);
			}else{
				const numeric::interpolation::spline::BicubicSpline &spline(AOH_Haro_scOrbH_splines_[orb_type_name]);
				distance_derivative = spline.dFdx(dist_angle_pair);
				angle_derivative = spline.dFdy(dist_angle_pair);
			}

		}
		if(h_enum == Hpol_bbOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(AOH_Hpol_bbOrbH_splines_[orb_type_name]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}
	}else{
		if(h_enum == Hpol_scOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(AOH_Hpol_scOrbH_splines_[orb_type_name]);
			energy = spline.F(dist_angle_pair);
			if(orb_type_name==chemical::orbitals::C_pi_sp2){
				//this is to match the orbital orbital and orbital hpol weight. Set the weight of cation_pi interactions to that of orb orb
				energy = (energy*0.01); //legacy. comes from splitting terms into individual pci_terms
			}
		}
		if(h_enum == Haro_scOrbH){
			if(ACO){
				const numeric::interpolation::spline::BicubicSpline &spline(ACO_AOH_orb_Hpol_splines_[orb_type_name]);
				energy = spline.F(dist_angle_pair);
			}else{
				const numeric::interpolation::spline::BicubicSpline &spline(AOH_Haro_scOrbH_splines_[orb_type_name]);
				energy = spline.F(dist_angle_pair);
			}
		}
		if(h_enum == Hpol_bbOrbH){
			const numeric::interpolation::spline::BicubicSpline &spline(AOH_Hpol_bbOrbH_splines_[orb_type_name]);
			energy = spline.F(dist_angle_pair);
		}
	}
}

void OrbitalsLookup::OrbOrbDist_cosAOD_energy(
		const core::Size orb_type_name1,
		const core::Size orb_type_name2,
		const core::Real distance,
		const core::Real AOO_angle,
		core::Real & energy,
		core::Real & distance_derivative,
		core::Real & angle_derivative,
		bool check_derivative
)const{

	//if(AOO_angle > 0) energy=distance_derivative=angle_derivative=0.0;
	numeric::MathVector<core::Real> dist_angle_pair(numeric::MakeVector(distance, AOO_angle));
	if(check_derivative){
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::N_pi_sp2))
			{
				const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[2]);
				distance_derivative = spline.dFdx(dist_angle_pair) * 0.033;
				angle_derivative = spline.dFdy(dist_angle_pair) * 0.033;
			}
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::O_pi_sp2)){
				const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[3]);
				distance_derivative = spline.dFdx(dist_angle_pair);
				angle_derivative = spline.dFdy(dist_angle_pair);
			}
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
				const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[1]);
				distance_derivative = (spline.dFdx(dist_angle_pair)/10); //legacy from splitting terms individually up
				angle_derivative = (spline.dFdy(dist_angle_pair)/10);
			}

		}
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::N_pi_sp2) && orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[4]);
			distance_derivative = spline.dFdx(dist_angle_pair) *.033;
			angle_derivative = spline.dFdy(dist_angle_pair) * 0.033;
		}
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::O_pi_sp2) && orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[5]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}

	}else{
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			if(orb_type_name2==core::chemical::orbitals::N_pi_sp2)
			{
				const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[2]);
				energy = spline.F(dist_angle_pair) * 0.033;
			}
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::O_pi_sp2)){
				const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[3]);
				energy = spline.F(dist_angle_pair);
			}
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
				const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[1]);
				energy = (spline.F(dist_angle_pair)/10);
			}

		}
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::N_pi_sp2) && orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[4]);
			energy = spline.F(dist_angle_pair) *0.033;
		}
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::O_pi_sp2) && orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			const numeric::interpolation::spline::BicubicSpline &spline(AOD_orb_orb_splines_[5]);
			energy = spline.F(dist_angle_pair);
		}
	}
}
void OrbitalsLookup::OrbOrbDist_cosDOA_energy(
		const core::Size orb_type_name1,
		const core::Size orb_type_name2,
		const core::Real distance,
		const core::Real DOO_angle,
		core::Real & energy,
		core::Real & distance_derivative,
		core::Real & angle_derivative,
		bool check_derivative

)const{
	//if(DOO_angle > 0) energy=distance_derivative=angle_derivative=0.0;
	numeric::MathVector<core::Real> dist_angle_pair(numeric::MakeVector(distance, DOO_angle));

	if(check_derivative){
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			if(orb_type_name2==core::chemical::orbitals::N_pi_sp2)
			{
				const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[2]);
				distance_derivative = (spline.dFdx(dist_angle_pair) *0.033);
				angle_derivative = (spline.dFdy(dist_angle_pair) *0.033);
			}
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::O_pi_sp2)){
				const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[3]);
				distance_derivative = spline.dFdx(dist_angle_pair);
				angle_derivative = spline.dFdy(dist_angle_pair);
			}
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
				const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[1]);
				distance_derivative = (spline.dFdx(dist_angle_pair)/10);
				angle_derivative = (spline.dFdy(dist_angle_pair)/10);
			}

		}
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::N_pi_sp2) && orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[4]);
			distance_derivative = (spline.dFdx(dist_angle_pair) *0.033);
			angle_derivative = (spline.dFdy(dist_angle_pair) *0.033);
		}
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::O_pi_sp2) && orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[5]);
			distance_derivative = spline.dFdx(dist_angle_pair);
			angle_derivative = spline.dFdy(dist_angle_pair);
		}

	}else{
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			if(orb_type_name2==core::chemical::orbitals::N_pi_sp2)
			{
				const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[2]);
				energy = (spline.F(dist_angle_pair) *0.033);
			}
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::O_pi_sp2)){
				const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[3]);
				energy = spline.F(dist_angle_pair);
			}
			if(orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
				const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[1]);
				energy = (spline.F(dist_angle_pair)/10);
			}

		}
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::N_pi_sp2) && orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[4]);
			energy = (spline.F(dist_angle_pair) *0.033);
		}
		if(orb_type_name1==static_cast <core::Size>(core::chemical::orbitals::O_pi_sp2) && orb_type_name2==static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			const numeric::interpolation::spline::BicubicSpline &spline(DOA_orb_orb_splines_[5]);
			energy = spline.F(dist_angle_pair);
		}
	}
}


}//namespace orbitals

}//namespace scoring

}//namespace core
