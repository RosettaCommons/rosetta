// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin OrbitalStatistics
///
/// @brief
/// A class for generating statistics from orbitals
///
/// @detailed
/// This is an attempt to be transparent in how the orbital-hydrogen interactions were converted into
/// a KBP. In general terms, a set of proteins from PDB were taken and the shortest distance from an
/// orbital to a polar or aromatic hydrogen was recorded and binned. In addition to the shortest distance,
/// an angle was taken. The angle considered is the angle between the base atom that contains the orbital, the
/// orbital, and the hydrogen. For polar hydrogens, only sidechain polar hydrogens were considered.
///
/// For protein interactions, there are 7 classes (orbital types) of orbitals that statistics are generated.
/// These 7 types are mapped using a map to enum data structure. For each sidechain interaction, only the shortest
/// distance between any given orbital type to a hydrogen is calculated. That means, that for each sidechain interaction
/// only 1 distance, 1 angle, and 1 class is recorded.
///
/// Bin sizes were calculated by .1A for distance and cos of 1 for angles. See below:
///
///           angle
///       -1 -.9 -.8 -.7..........
/// d .1 | 0   0   0   0
/// i .2 | 500 0   0   0
/// s .3 | 25  0   0   0
/// t .4 |  0  30  5   0
///
/// This is not the original code that was used to generate the statistics. The original code was much more
/// convoluted than this because I had no idea how to program. I wrote this piece for clarity. I have tested
/// it and it produces the same results.
///
/// @authors
/// Steven Combs
///
///
/////////////////////////////////////////////////////////////////////////

#include <core/scoring/orbitals/OrbitalsStatistics.hh>
#include <core/pose/Pose.hh>
#include <numeric/histograms/TwoDHistogram.hh>
#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>
#include <core/chemical/orbitals/OrbitalType.hh>
#include <core/pose/PDB_Info.hh>


//utility headers
#include <utility/vector1.hh>


// AUTO-REMOVED #include <math.h> // REQUIRED FOR WINDOWS

#include <map>


//option headers
#include <basic/options/option.hh>
#include <basic/options/keys/orbitals.OptionKeys.gen.hh>


namespace core{
namespace scoring{
namespace orbitals{

// REQUIRED FOR WINDOWS
#ifdef _WIN32                       // REQUIRED FOR WINDOWS
	double round( double d ) {      // REQUIRED FOR WINDOWS
		return (d > 0.0) ? floor(d + 0.5) : ceil(d - 0.5);      // REQUIRED FOR WINDOWS
	}                               // REQUIRED FOR WINDOWS
#endif                              // REQUIRED FOR WINDOWS
// REQUIRED FOR WINDOWS

OrbitalsStatistics::OrbitalsStatistics()
{
	if(
		basic::options::option[basic::options::OptionKeys::orbitals::Haro] ||
		basic::options::option[basic::options::OptionKeys::orbitals::Hpol])
	{
		number_of_histograms_=7;
	}
	for(core::Size i=0; i <= 100; ++i){
		for(core::SSize j=-10; j<= 10; ++j){
			std::pair<core::Size, core::SSize> pair(i,j);
			twoD_histogram_.insert_data(pair, 0);

		}
	}

	//map orbial type to enum
	orbital_type_2_enum_["C.pi.sp2"]=C_pi_sp2;
	orbital_type_2_enum_["N.pi.sp2"]=N_pi_sp2;
	orbital_type_2_enum_["N.p.sp2"]=N_p_sp2;
	orbital_type_2_enum_["O.pi.sp2"]=O_pi_sp2;
	orbital_type_2_enum_["O.p.sp2"]=O_p_sp2;
	orbital_type_2_enum_["O.p.sp3"]=O_p_sp3;
	orbital_type_2_enum_["S.p.sp3"]=S_p_sp3;



}




numeric::histograms::TwoDHistogram<core::Size, core::SSize> OrbitalsStatistics::get_2D_histogram()
{
	return twoD_histogram_;
}

utility::vector1< numeric::histograms::TwoDHistogram<core::Size, core::SSize> >  OrbitalsStatistics::get_histogram_vector()
{
	return histogram_vector_;
}

core::Size OrbitalsStatistics::get_number_of_histograms()
{
	return number_of_histograms_;
}


///@brief increment the orbital histogram bin based upon a distance and angle.
/// Currently the statistics are best with .1 incremented bins
void OrbitalsStatistics::increment_histogram_bin(
		core::Real & distance,
		core::Real & angle,
		numeric::histograms::TwoDHistogram<core::Size, core::SSize> & histogram
)
{

	if(angle >= -1 && angle <= -.9){
		angle = -.95;
	}
	if(angle >= -.89 && angle <= -.8){
		angle = -.85;
	}
	if(angle >= -.79 && angle <= -.7){
		angle = -.75;
	}
	if(angle >= -.69 && angle <= -.6){
		angle = -.65;
	}
	if(angle >= -.59 && angle <= -.5){
		angle = -.55;
	}
	if(angle >= -.49 && angle <= -.4){
		angle = -.45;
	}
	if(angle >= -.39 && angle <= -.3){
		angle = -.35;
	}
	if(angle >= -.29 && angle <= -.2){
		angle = -.25;
	}
	if(angle >= -.19 && angle <= -.1){
		angle = -.15;
	}
	if(angle >= -.09 && angle <= 0){
		angle = -.05;
	}
	if(angle >= 0.01 && angle <= .1){
		angle = .05;
	}
	if(angle >= .11 && angle <= .2){
		angle = .15;
	}
	if(angle >= .21 && angle <= .3){
		angle = .25;
	}
	if(angle >= .31 && angle <= .4){
		angle = .35;
	}
	if(angle >= .41 && angle <= .5){
		angle = .45;
	}
	if(angle >= .51 && angle <= .6){
		angle = .55;
	}
	if(angle >= .61 && angle <= .7){
		angle = .65;
	}
	if(angle >= .71 && angle <= .8){
		angle = .75;
	}
	if(angle >= .81 && angle <= .9){
		angle = .85;
	}
	if(angle >= .91 && angle <= 1){
		angle = .95;
	}
	core::SSize new_angle= static_cast<core::SSize> (round(angle/.1));
	core::Size new_dist= static_cast<core::Size> (round(distance/.1));
	std::pair< core::Size, core::SSize> new_pair(new_dist, new_angle);
	histogram.increase_count(new_pair);
	//++histogram[new_pair];
//	std::cout << " outputing paired distance that was put in increment_histogram_bind "<<histogram.lookup_counts(new_pair) << std::endl;;
}


///@brief get statistics based upon hydrogen to orbital distance/angle
void
OrbitalsStatistics::sc_H_orbital( core::pose::Pose & pdb )
{
	//std::cout << "Made it to the sc_H_orbital function!!!!!" << std::endl;
	numeric::histograms::TwoDHistogram<core::Size, core::SSize> histogram=twoD_histogram_;

	//a vector that will contain all the histograms.
	utility::vector1<numeric::histograms::TwoDHistogram<core::Size, core::SSize> > histogram_vector;

	//add histograms to the previously created vector. Notice that the number created is based
	//upon number_of_histograms_ which is constructed in the
	for(core::Size x=1; x<= number_of_histograms_; ++x){
		histogram_vector.push_back(histogram);
	}

	//Let the magic begin!
	for (core::Size res_num1 = 1; res_num1 <= pdb.n_residue(); ++res_num1) {
		core::conformation::Residue resid1 = pdb.residue(res_num1);
		if (resid1.has_sc_orbitals()) {
			numeric::xyzVector<core::Real> final_orb_xyz(0.0);
			core::Real low_D = 11;
			numeric::xyzVector<core::Real> res2_H_xyz;
			std::string orbital_type("");
			core::Real angle=2;
			core::Size atom_index_min_dist(0);
			core::Size res2seqpos(0);
			std::string res2name;
			core::Real angle2=2;


			//std::cout << "inside the resid1.has_sc_orbitals()!" << std::endl;
			for (core::Size res_num2 = 1; res_num2 <= pdb.n_residue(); ++res_num2) {
				core::conformation::Residue resid2 = pdb.residue(res_num2);
				if (resid1.seqpos() != resid2.seqpos()) {
					if (resid1.nbr_atom_xyz().distance(resid2.nbr_atom_xyz()) <= 10.0) {
						numeric::xyzVector<core::Real> res2_H_xyz;
						//iterate through atoms with orbitals
						for ( core::chemical::AtomIndices::const_iterator
								atom_index  =  resid1.atoms_with_orb_index().begin(),
								atom_end = resid1.atoms_with_orb_index().end();
								atom_index != atom_end; ++atom_index
						){
							if(!resid1.atom_is_backbone(*atom_index)){
								utility::vector1<core::Size> const & orbital_indices(resid1.bonded_orbitals(*atom_index));
								//iterate through the orbitals
								for(
										utility::vector1<core::Size>::const_iterator
										orbital_index = orbital_indices.begin(),
										orbital_end = orbital_indices.end();
										orbital_index != orbital_end; ++orbital_index
								){
									numeric::xyzVector<core::Real> orb_xyz = resid1.orbital_xyz(*orbital_index);
									numeric::xyzVector<core::Real> bonded_atom_xyz(resid1.atom(*atom_index).xyz());
									if ( basic::options::option[basic::options::OptionKeys::orbitals::Hpol] ) {
										//iterate only throught the sidechain polar hydrogens
										for(
												core::chemical::AtomIndices::const_iterator
												hpol_index = resid2.Hpol_index().begin(),
												hpol_end = resid2.Hpol_index().end(); hpol_index != hpol_end; ++hpol_index
										) {
											res2_H_xyz = resid2.atom( *hpol_index ).xyz();
											numeric::xyzVector<core::Real> DHO_atom_xyz(resid2.atom(resid2.bonded_neighbor(*hpol_index)[1]).xyz() );
											core::Real container = orb_xyz.distance( res2_H_xyz );
											if(container <= low_D){
												final_orb_xyz = orb_xyz;
												angle = cos_of(bonded_atom_xyz, orb_xyz, res2_H_xyz );
												angle2 = cos_of(DHO_atom_xyz, res2_H_xyz, orb_xyz);
												low_D = container;
												orbital_type = resid1.orbital_type(*orbital_index).name();
												res2name= resid2.name3();
												res2seqpos=resid2.seqpos();
												atom_index_min_dist=*hpol_index;
											}
										}
									}  // end if Hpol
									//haro
									if(basic::options::option[basic::options::OptionKeys::orbitals::Haro]){
										for(core::chemical::AtomIndices::const_iterator
												haro_index = resid2.Haro_index().begin(),
												haro_end = resid2.Haro_index().end(); haro_index != haro_end; ++haro_index) {
											res2_H_xyz = resid2.atom( *haro_index ).xyz();
											numeric::xyzVector<core::Real> DHO_atom_xyz(resid2.atom(resid2.bonded_neighbor(*haro_index)[1]).xyz() );
											core::Real container = orb_xyz.distance( res2_H_xyz  );
											if (container <= low_D) {
												angle = cos_of(bonded_atom_xyz, orb_xyz, res2_H_xyz );
												angle2 = cos_of(DHO_atom_xyz, res2_H_xyz, orb_xyz);
												orbital_type = resid1.orbital_type(*orbital_index).name();
												low_D = container;
												res2name= resid2.name3();
												res2seqpos=resid2.seqpos();
												atom_index_min_dist=*haro_index;
											}
										}
									}
								}  // end for orbital_index != orbital_end
							}  // end if !resid1.atom_is_backbone(*atom_index)
						}  // end for atom_index != atom_end
					}  // end if resid1.nbr_atom_xyz().distance(resid2.nbr_atom_xyz()) <= 10.0
				}  // end if resid1.seqpos() != resid2.seqpos()
			}  // end for res_num2 <= pdb.n_residue()

			if ( low_D <= 10 ) {
				core::conformation::Residue res2 = pdb.residue(res2seqpos);
				//std::cout << low_D << " angle: " << angle << std::endl;
				//decompose_score(resid1, resid2, pdb, low_D, angle, electron_name);
				if(orbital_type=="O.p.sp3" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])){
					statistics_output_.open("O.p.sp3.backbone", std::ios_base::app);

					statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " "<< pdb.pdb_info()->name() << " " << final_orb_xyz.x() << " " << final_orb_xyz.y() << " " << final_orb_xyz.z() << std::endl;
					statistics_output_.close();
				}else if(orbital_type=="S.p.sp3" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
					statistics_output_.open("S.p.sp3.backbone", std::ios_base::app);
					statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
					statistics_output_.close();
				}else if(orbital_type=="C.pi.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
					statistics_output_.open("C.pi.sp2.backbone", std::ios_base::app);
					statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
					statistics_output_.close();
				}
				else if(orbital_type=="N.pi.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
					statistics_output_.open("N.pi.sp2.backbone", std::ios_base::app);
					statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
					statistics_output_.close();
				}else if(orbital_type=="N.p.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
					statistics_output_.open("N.p.sp2.backbone", std::ios_base::app);
					statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
					statistics_output_.close();
				}else if(orbital_type=="O.pi.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
					statistics_output_.open("O.pi.sp2.backbone", std::ios_base::app);
					statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
					statistics_output_.close();
				}else if(orbital_type=="O.p.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
					statistics_output_.open("O.p.sp2.backbone", std::ios_base::app);
					statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
					statistics_output_.close();
				}else {
					statistics_output_.open(orbital_type.c_str(), std::ios_base::app);
					statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
					statistics_output_.close();
					//core::Size index_of_histogram(orbital_type_2_enum_.find(orbital_type)->second);
					//std::cout << "index of histogram " << index_of_histogram << std::endl;
					//increment_histogram_bin(low_D, angle, histogram_vector[ index_of_histogram ]);
				}
			}  // end if low_D <= 10
		}  //end if resid1.has_sc_orbitals()
	}  // end for res_num1 <= pdb.n_residue()
	histogram_vector_=histogram_vector;
}


void OrbitalsStatistics::bb_stats(
		core::pose::Pose & pdb
){
	//std::cout << "Made it to the sc_H_orbital function!!!!!" << std::endl;
		numeric::histograms::TwoDHistogram<core::Size, core::SSize> histogram=twoD_histogram_;

		//a vector that will contain all the histograms.
		utility::vector1<numeric::histograms::TwoDHistogram<core::Size, core::SSize> > histogram_vector;

		//add histograms to the previously created vector. Notice that the number created is based
		//upon number_of_histograms_ which is constructed in the
		for(core::Size x=1; x<= number_of_histograms_; ++x){
			histogram_vector.push_back(histogram);
		}

		//Let the magic begin!
		for(core::Size res_num1 = 1; res_num1 <= pdb.n_residue(); ++res_num1){
			core::conformation::Residue resid1 = pdb.residue(res_num1);

				core::Real low_D = 11;
				numeric::xyzVector<core::Real> res2_H_xyz;
				std::string orbital_type("");
				core::Real angle=2;
				core::Size atom_index_min_dist(0);
				core::Size res2seqpos(0);
				std::string res2name;
				core::Real angle2=2;

				//std::cout << "inside the resid1.has_sc_orbitals()!" << std::endl;
				for(core::Size res_num2 = 1; res_num2 <= pdb.n_residue(); ++res_num2){
					core::conformation::Residue resid2 = pdb.residue(res_num2);
					if(resid1.seqpos() != resid2.seqpos()){
						if(resid1.nbr_atom_xyz().distance(resid2.nbr_atom_xyz()) <= 10.0)
						{
						numeric::xyzVector<core::Real> res2_H_xyz;
						//iterate through atoms with orbitals
						for ( core::chemical::AtomIndices::const_iterator
								atom_index  =  resid1.atoms_with_orb_index().begin(),
								atom_end = resid1.atoms_with_orb_index().end();
								atom_index != atom_end; ++atom_index
						){
							if(resid1.atom_is_backbone(*atom_index)){

									utility::vector1<core::Size> const & orbital_indices(resid1.bonded_orbitals(*atom_index));
									//iterate through the orbitals
									for(
											utility::vector1<core::Size>::const_iterator
											orbital_index = orbital_indices.begin(),
											orbital_end = orbital_indices.end();
											orbital_index != orbital_end; ++orbital_index
									){

										numeric::xyzVector<core::Real> orb_xyz = resid1.orbital_xyz(*orbital_index);
										numeric::xyzVector<core::Real> bonded_atom_xyz(resid1.atom(*atom_index).xyz());
										if(
												basic::options::option[basic::options::OptionKeys::orbitals::Hpol]
										){

											//iterate only throught the sidechain polar hydrogens
											for(
													core::chemical::AtomIndices::const_iterator
													hpol_index = resid2.Hpol_index().begin(),
													hpol_end = resid2.Hpol_index().end(); hpol_index != hpol_end; ++hpol_index
											)
											{
												res2_H_xyz = resid2.atom( *hpol_index ).xyz();
												core::Real container = orb_xyz.distance( res2_H_xyz );
												numeric::xyzVector<core::Real> DHO_atom_xyz(resid2.atom(resid2.bonded_neighbor(*hpol_index)[1]).xyz() );
												if(container <= low_D){
													angle = cos_of(bonded_atom_xyz, orb_xyz, res2_H_xyz );
													angle2 = cos_of(DHO_atom_xyz, res2_H_xyz, orb_xyz);
													low_D = container;
													orbital_type = resid1.orbital_type(*orbital_index).name();
													res2name= resid2.name3();
													res2seqpos=resid2.seqpos();
													atom_index_min_dist=*hpol_index;
												}
											}
										}
										//haro
										if(basic::options::option[basic::options::OptionKeys::orbitals::Haro]){
											for(core::chemical::AtomIndices::const_iterator
													haro_index = resid2.Haro_index().begin(),
													haro_end = resid2.Haro_index().end(); haro_index != haro_end; ++haro_index)
											{
												res2_H_xyz = resid2.atom( *haro_index ).xyz();
												core::Real container = orb_xyz.distance( res2_H_xyz  );
												numeric::xyzVector<core::Real> DHO_atom_xyz(resid2.atom(resid2.bonded_neighbor(*haro_index)[1]).xyz() );
												if(container <= low_D){
													angle = cos_of(bonded_atom_xyz, orb_xyz, res2_H_xyz );
													angle2 = cos_of(DHO_atom_xyz, res2_H_xyz, orb_xyz);
													orbital_type = resid1.orbital_type(*orbital_index).name();
													low_D = container;
													res2name= resid2.name3();
													res2seqpos=resid2.seqpos();
													atom_index_min_dist=*haro_index;
												}
											}
										}
									}
								}
							}
						}
					}
				}// end for resid2

				if(low_D <= 10){
					core::conformation::Residue res2 = pdb.residue(res2seqpos);
					//std::cout << low_D << " angle: " << angle << std::endl;
					//decompose_score(resid1, resid2, pdb, low_D, angle, electron_name);
					if(orbital_type=="O.p.sp3" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])){
						statistics_output_.open("O.p.sp3.backbone", std::ios_base::app);
						statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
						statistics_output_.close();
					}else if(orbital_type=="S.p.sp3" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
						statistics_output_.open("S.p.sp3.backbone", std::ios_base::app);
						statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
						statistics_output_.close();
					}else if(orbital_type=="C.pi.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
							statistics_output_.open("C.pi.sp2.backbone", std::ios_base::app);
							statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
							statistics_output_.close();
					}
					else if(orbital_type=="N.pi.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
							statistics_output_.open("N.pi.sp2.backbone", std::ios_base::app);
							statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
							statistics_output_.close();
					}else if(orbital_type=="N.p.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
							statistics_output_.open("N.p.sp2.backbone", std::ios_base::app);
							statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
							statistics_output_.close();
					}else if(orbital_type=="O.pi.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
							statistics_output_.open("O.pi.sp2.backbone", std::ios_base::app);
							statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
							statistics_output_.close();
					}else if(orbital_type=="O.p.sp2" && res2.atom_is_backbone(res2.bonded_neighbor(atom_index_min_dist)[1])) {
							statistics_output_.open("O.p.sp2.backbone", std::ios_base::app);
							statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
							statistics_output_.close();
					}else {
						statistics_output_.open(orbital_type.c_str(), std::ios_base::app);
						statistics_output_ << resid1.name3() << resid1.seqpos() << " " << res2name << res2seqpos << " "<< low_D << " " << angle << " " << angle2 << " " << pdb.pdb_info()->name() << std::endl;
						statistics_output_.close();
					}
				}
		}
		histogram_vector_=histogram_vector;
}


}
}
}
