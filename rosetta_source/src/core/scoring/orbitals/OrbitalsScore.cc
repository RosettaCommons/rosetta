// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#include <core/pose/Pose.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/scoring/orbitals/OrbitalsScore.hh>
#include <core/scoring/orbitals/OrbitalsScoreCreator.hh>
#include <core/scoring/orbitals/OrbitalsAssigned.hh>
#include <core/chemical/orbitals/OrbitalTypeMapper.fwd.hh>
#include <map>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/distance_deriv.hh>


#include <core/conformation/Residue.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/orbitals/OrbitalType.hh>
#include <core/scoring/ContextGraphTypes.hh>

#include <core/scoring/EnergyMap.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <numeric/xyzTriple.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/orbitals.OptionKeys.gen.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/pose/util.hh>

#include <core/pose/Pose.hh>
#include <list>

namespace core{
namespace scoring{
namespace orbitals{

methods::EnergyMethodOP
OrbitalsScoreCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const
{
	return new OrbitalsScore(options);
}

ScoreTypes
OrbitalsScoreCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( orbitals_hpol );
	sts.push_back( orbitals_haro );
	sts.push_back(orbitals_orbitals);
	sts.push_back(orbitals_hpol_bb);
	return sts;
}

static basic::Tracer TR("core.scoring.orbitals_hpol");

OrbitalsScore::OrbitalsScore(methods::EnergyMethodOptions const &) :
	parent( new OrbitalsScoreCreator ),
	lookup_table_(core::scoring::ScoringManager::get_instance()->get_OrbitalsLookupTable()),
	//lookup_Etable_(core::scoring::etable::EtableEnergy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options)),
	max_dist_squared_(100.0),
	max_orbital_dist_squared_(100.0),
	nbr_distance_squared_(basic::options::option[basic::options::OptionKeys::orbitals::nbr_distance_squared])
	//etable_( *ScoringManager::get_instance()->etable( options.etable_type() ))

{
	min_orb_dist_enum_map_[ core::chemical::orbitals::C_pi_sp2 ] = 100.0;
	min_orb_dist_enum_map_[ core::chemical::orbitals::N_pi_sp2 ] = 100.0;
	min_orb_dist_enum_map_[ core::chemical::orbitals::N_p_sp2  ] = 100.0;
	min_orb_dist_enum_map_[ core::chemical::orbitals::O_pi_sp2 ] = 100.0;
	min_orb_dist_enum_map_[ core::chemical::orbitals::O_p_sp2  ] = 100.0;
	min_orb_dist_enum_map_[ core::chemical::orbitals::O_p_sp3  ] = 100.0;
	min_orb_dist_enum_map_[ core::chemical::orbitals::S_p_sp3  ] = 100.0;
}

methods::EnergyMethodOP
OrbitalsScore::clone() const
{
	return new OrbitalsScore(*this);
}

void OrbitalsScore::setup_for_scoring(pose::Pose & pose, ScoreFunction const &) const
{
	pose.update_residue_neighbors();
}

void OrbitalsScore::setup_for_derivatives(pose::Pose &, ScoreFunction const &) const
{
	// do nothing for now, otherwise do: pose.update_residue_neighbors();
}



bool OrbitalsScore::defines_intrares_energy(core::scoring::EnergyMap const &) const{
	return false;
}




void OrbitalsScore::eval_intrares_energy(
	core::conformation::Residue const &,
	core::pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const {}

/// This is very wrong.
core::Distance OrbitalsScore::atomic_interaction_cutoff() const{
	return  6.0;
}


void OrbitalsScore::indicate_required_context_graphs(utility::vector1< bool > & ) const
{}


void
OrbitalsScore::residue_pair_energy(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	core::pose::Pose const &,
	core::scoring::ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( res1.has_sc_orbitals() || res2.has_sc_orbitals() ) {
		//if ( res1.nbr_atom_xyz().distance_squared(res2.nbr_atom_xyz()) <= nbr_distance_squared_ ) {

			if ( res1.is_aromatic() && res2.is_aromatic() ) {
				std::list< orbital_info_for_score > all_orbital_haro_pair_energies;
				core::Real haro_E=0.0;
				bool was_scored=false;
				bool check_derivative=false;
				all_orbital_haro_pair_energies = get_E_haro_one_way(res1, res2, haro_E, check_derivative, was_scored);
				add_energies_to_emap(emap, all_orbital_haro_pair_energies);
				all_orbital_haro_pair_energies = get_E_haro_one_way(res2, res1, haro_E, check_derivative, was_scored);
				add_energies_to_emap(emap, all_orbital_haro_pair_energies);
			} else {
				//PROF_START( basic::ORBITALS_E_1_WAY );
				std::list< orbital_info_for_score > all_orbital_hpol_pair_energies;
				if ( res1.has_sc_orbitals() ) {
					core::Real hpol_E(0.0);
					core::Real bb_hpol_E(0.0);
					bool check_derivative=false;
					bool was_scored=false;
					all_orbital_hpol_pair_energies = get_E_hpol_one_way(res1, res2, hpol_E, bb_hpol_E, check_derivative,was_scored);
					add_energies_to_emap(emap, all_orbital_hpol_pair_energies);
				}
				if ( res2.has_sc_orbitals() ) {
					core::Real hpol_E(0.0);
					core::Real bb_hpol_E(0.0);
					bool check_derivative=false;
					bool was_scored=false;
					all_orbital_hpol_pair_energies = get_E_hpol_one_way(res2, res1, hpol_E, bb_hpol_E, check_derivative, was_scored);
					add_energies_to_emap(emap, all_orbital_hpol_pair_energies);
				}
			}


			//we only want to consider those residues that have polar side chain hydrogens. The only atoms on the bb that have
			//orbitals are the carbonyl oxygens. The only orbitals considered are the sp2 orbitals. These are hardcorded in the
			//params file. Maybe one day the ligand code will be put into svn and the params files wont be hardcoded.
/*			if ( res2.is_polar() || res2.name3() == "TYR" || res2.name3() == "TRP" ) {
				core::Real sc_bb_E(0.0);
				bool check_derivative=false;
				bool was_scored=false;
				std::list< orbital_info_for_score > not_needed_struct;
				not_needed_struct= get_sc_bb_E_one_way(res1, res2, sc_bb_E, check_derivative, was_scored);
				emap[orbitals_hpol_bb] += sc_bb_E;
			}
			if ( res1.is_polar() || res1.name3() == "TYR" || res1.name3() == "TRP" ) {
				core::Real sc_bb_E(0.0);
				bool check_derivative=false;
				bool was_scored=false;
				std::list< orbital_info_for_score > not_needed_struct;
				not_needed_struct= get_sc_bb_E_one_way(res2, res1, sc_bb_E, check_derivative, was_scored);
				emap[orbitals_hpol_bb] += sc_bb_E;
			}*/
			//PROF_STOP( basic::ORBITALS_E_1_WAY );

		//}


	}



}



void OrbitalsScore::add_energies_to_emap(
		EnergyMap & emap,
		std::list< orbital_info_for_score > & orbital_hydrogen_energies
) const
{
	for(
			std::list< orbital_info_for_score >::iterator
			energy = orbital_hydrogen_energies.begin(),
			energy_end = orbital_hydrogen_energies.end();
			energy != energy_end; ++energy
	)
	{

		//std::cout << "ENRGY_IN_EMAP" << energy->energy << std::endl;
		if(energy->htype == lookup_table_.haro){
			emap[orbitals_haro] += energy->energy;
		}
		if(energy->htype == lookup_table_.hpol){
			emap[orbitals_hpol] += energy->energy;
		}
	}
}


std::list< orbital_info_for_score >
OrbitalsScore::get_E_haro_one_way(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	core::Real & haro_E,
	bool check_derivative,
	bool & was_scored
) const
{
	bool check_orb_distance=false;


	std::list<std::pair< core::Size, utility::vector1< core::Size > >  > atoms_with_orbs_and_hydrogen_pairs;

	get_bonded_orb_atom_and_Haro_atom(res1, res2, atoms_with_orbs_and_hydrogen_pairs, check_orb_distance);


	bool get_haro_score=false;
	core::Real min_dist_squared(max_orbital_dist_squared_);
	core::chemical::orbitals::orbital_type_enum orbital_type;
	std::list< orbital_info_for_score > min_orbital_dist_tau;

	if ( check_orb_distance ) {
		min_orbital_dist_tau = get_orb_H_distance_and_angle(res1, res2, atoms_with_orbs_and_hydrogen_pairs, min_dist_squared, orbital_type, get_haro_score);
	}

	if ( get_haro_score ) {
		was_scored=true;
		core::Real distance_derivative(0.0);
		core::Real angle_derivative(0.0);
		for( std::list< orbital_info_for_score >::iterator
				orbital_min_dist_tau = min_orbital_dist_tau.begin(),
				orbital_min_dist_tau_end = min_orbital_dist_tau.end();
				orbital_min_dist_tau != orbital_min_dist_tau_end; ++orbital_min_dist_tau ) {
			//core::Size struct_index=*orbital_min_dist_tau;
			get_orbital_H_score(
				lookup_table_.haro,
				orbital_min_dist_tau->orbital_type,
				orbital_min_dist_tau->distance,
				orbital_min_dist_tau->tau,
				haro_E,
				distance_derivative,
				angle_derivative,
				check_derivative
			);
			orbital_min_dist_tau->distance_derivative=distance_derivative;
			orbital_min_dist_tau->angle_derivative=angle_derivative;
			orbital_min_dist_tau->energy=haro_E;
			orbital_min_dist_tau->htype=lookup_table_.haro;
			//if ( check_derivative ) {
			if ( false ) {
				std::cout << "Haro " << res1.seqpos() << " " << res2.seqpos() << " a_index: " << orbital_min_dist_tau->atom_index
					<< " o_index: " << orbital_min_dist_tau->orbital_index << " h_index: "<<  orbital_min_dist_tau->h_index
					<< " E: " << haro_E << " d_der: " << distance_derivative << " " <<  orbital_min_dist_tau->distance_derivative
					<< " a_der: " << angle_derivative << " " << orbital_min_dist_tau->angle_derivative << " dist: "
					<< orbital_min_dist_tau->distance << " angle: " << orbital_min_dist_tau->tau << " enum: " << orbital_min_dist_tau->orbital_type << std::endl;

					Real delta = 1e-4;
					Real hpol_E_dpdelta, distderiv_dpdelta, angderiv_dpdelta;
					get_orbital_H_score(
						lookup_table_.haro,
						orbital_min_dist_tau->orbital_type,
						orbital_min_dist_tau->distance + delta,
						orbital_min_dist_tau->tau,
						hpol_E_dpdelta,
						distderiv_dpdelta,
						angderiv_dpdelta,
						check_derivative );
					Real hpol_E_dmdelta, distderiv_dmdelta, angderiv_dmdelta;
					get_orbital_H_score(
						lookup_table_.haro,
						orbital_min_dist_tau->orbital_type,
						orbital_min_dist_tau->distance - delta,
						orbital_min_dist_tau->tau,
						hpol_E_dmdelta,
						distderiv_dmdelta,
						angderiv_dmdelta,
						check_derivative );
					Real num_dderiv( (hpol_E_dpdelta - hpol_E_dmdelta ) / (2*delta) );
					std::cout << "distance derivative: analytic " << distance_derivative << " numeric deriv: " << num_dderiv << std::endl;

					Real hpol_E_apdelta, distderiv_apdelta, angderiv_apdelta;
					get_orbital_H_score(
						lookup_table_.haro,
						orbital_min_dist_tau->orbital_type,
						orbital_min_dist_tau->distance,
						orbital_min_dist_tau->tau + delta,
						hpol_E_apdelta,
						distderiv_apdelta,
						angderiv_apdelta,
						check_derivative );
					Real hpol_E_amdelta, distderiv_amdelta, angderiv_amdelta;
					get_orbital_H_score(
						lookup_table_.haro,
						orbital_min_dist_tau->orbital_type,
						orbital_min_dist_tau->distance,
						orbital_min_dist_tau->tau - delta,
						hpol_E_amdelta,
						angderiv_apdelta,
						angderiv_amdelta,
						check_derivative );
					Real num_aderiv( (hpol_E_apdelta - hpol_E_amdelta ) / (2*delta) );
					std::cout << "angle derivative: analytic " << angle_derivative << " numeric deriv: " << num_aderiv << std::endl;

			}
		}

	}
	return min_orbital_dist_tau;
}


std::list< orbital_info_for_score >
OrbitalsScore::get_E_hpol_one_way(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	core::Real & hpol_E,
	core::Real & sc_bb_E,
	bool check_derivative,
	bool & was_scored
) const
{

	bool check_orb_distance=false;

	std::list<std::pair< core::Size, utility::vector1< core::Size > >  > atoms_with_orbs_and_hydrogen_pairs;

	//this function checks to see if an atom with orbitals and a Hydrogen pair fall within 4A. If it does, it returns
	//a vector with the paired atoms and hydrogens
	get_bonded_orb_atom_and_Hpol_atom(res1, res2, atoms_with_orbs_and_hydrogen_pairs, check_orb_distance);

	bool get_hpol_score=false;
	core::Real min_dist_squared(max_orbital_dist_squared_);
	core::chemical::orbitals::orbital_type_enum orbital_type;
	std::list< orbital_info_for_score > min_orbital_dist_tau;

	if ( check_orb_distance ) {
		min_orbital_dist_tau = get_orb_H_distance_and_angle(res1, res2, atoms_with_orbs_and_hydrogen_pairs, min_dist_squared, orbital_type, get_hpol_score);
	}


	if ( get_hpol_score ) {
		was_scored=true;
		for( std::list< orbital_info_for_score >::iterator
				orbital_min_dist_tau = min_orbital_dist_tau.begin(),
				orbital_min_dist_tau_end = min_orbital_dist_tau.end();
				orbital_min_dist_tau != orbital_min_dist_tau_end; ++orbital_min_dist_tau ) {
			if ( orbital_min_dist_tau->lookup_backbone_energy ) {
				core::Real angle_derivative(0.0);
				core::Real distance_derivative(0.0);
				get_orbital_H_score(
					lookup_table_.sc_orb_bb_H,
					orbital_min_dist_tau->orbital_type,
					orbital_min_dist_tau->distance,
					orbital_min_dist_tau->tau,
					sc_bb_E,
					angle_derivative,
					distance_derivative,
					check_derivative );
				orbital_min_dist_tau->energy=sc_bb_E;
				orbital_min_dist_tau->distance_derivative=distance_derivative;
				orbital_min_dist_tau->angle_derivative=angle_derivative;

				//if ( check_derivative ) {
				if ( false ) {
					std::cout << "bb_H "<< res1.seqpos() << " " << res2.seqpos() << " a_index: " << orbital_min_dist_tau->atom_index
						<< " o_index: " << orbital_min_dist_tau->orbital_index << " h_index: "<<  orbital_min_dist_tau->h_index
						<< " E: " << sc_bb_E << " d_der: " << distance_derivative << " " <<  orbital_min_dist_tau->distance_derivative
						<< " a_der: " << angle_derivative << " " << orbital_min_dist_tau->angle_derivative << std::endl;
				}

			} else {
				core::Real angle_derivative(0.0);
				core::Real distance_derivative(0.0);

				get_orbital_H_score(
					lookup_table_.hpol,
					orbital_min_dist_tau->orbital_type,
					orbital_min_dist_tau->distance,
					orbital_min_dist_tau->tau,
					hpol_E,
					distance_derivative,
					angle_derivative,
					check_derivative );
				orbital_min_dist_tau->energy=hpol_E;
				orbital_min_dist_tau->distance_derivative=distance_derivative;
				orbital_min_dist_tau->angle_derivative=angle_derivative;
				orbital_min_dist_tau->htype=lookup_table_.hpol;

				if ( false ) {
				//if ( check_derivative ) {
					std::cout << "Hpol "<< res1.seqpos() << " " << res2.seqpos() << " a_index: " << orbital_min_dist_tau->atom_index
						<< " o_index: " << orbital_min_dist_tau->orbital_index << " h_index: "<<  orbital_min_dist_tau->h_index
						<< " E: " << hpol_E << " d_der: " << distance_derivative << " " <<  orbital_min_dist_tau->distance_derivative
						<< " a_der: " << angle_derivative << " " << orbital_min_dist_tau->angle_derivative << " dist: "
						<< orbital_min_dist_tau->distance << " angle: " << orbital_min_dist_tau->tau << " enum: " << orbital_min_dist_tau->orbital_type << std::endl;

					Real delta = 1e-4;
					Real hpol_E_dpdelta, distderiv_dpdelta, angderiv_dpdelta;
					get_orbital_H_score(
						lookup_table_.hpol,
						orbital_min_dist_tau->orbital_type,
						orbital_min_dist_tau->distance + delta,
						orbital_min_dist_tau->tau,
						hpol_E_dpdelta,
						distderiv_dpdelta,
						angderiv_dpdelta,
						check_derivative );
					Real hpol_E_dmdelta, distderiv_dmdelta, angderiv_dmdelta;
					get_orbital_H_score(
						lookup_table_.hpol,
						orbital_min_dist_tau->orbital_type,
						orbital_min_dist_tau->distance - delta,
						orbital_min_dist_tau->tau,
						hpol_E_dmdelta,
						distderiv_dmdelta,
						angderiv_dmdelta,
						check_derivative );
					Real num_dderiv( (hpol_E_dpdelta - hpol_E_dmdelta ) / (2*delta) );
					std::cout << "distance derivative: analytic " << distance_derivative << " numeric deriv: " << num_dderiv << std::endl;

					Real hpol_E_apdelta, distderiv_apdelta, angderiv_apdelta;
					get_orbital_H_score(
						lookup_table_.hpol,
						orbital_min_dist_tau->orbital_type,
						orbital_min_dist_tau->distance,
						orbital_min_dist_tau->tau + delta,
						hpol_E_apdelta,
						distderiv_apdelta,
						angderiv_apdelta,
						check_derivative );
					Real hpol_E_amdelta, distderiv_amdelta, angderiv_amdelta;
					get_orbital_H_score(
						lookup_table_.hpol,
						orbital_min_dist_tau->orbital_type,
						orbital_min_dist_tau->distance,
						orbital_min_dist_tau->tau - delta,
						hpol_E_amdelta,
						angderiv_apdelta,
						angderiv_amdelta,
						check_derivative );
					Real num_aderiv( (hpol_E_apdelta - hpol_E_amdelta ) / (2*delta) );
					std::cout << "angle derivative: analytic " << angle_derivative << " numeric deriv: " << num_aderiv << std::endl;

				}
			}
		}
	}
	return min_orbital_dist_tau;

}




std::list< orbital_info_for_score >
OrbitalsScore::get_sc_bb_E_one_way(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	core::Real & sc_bb_E,
	bool check_derivative,
	bool & was_scored
) const
{

	bool check_orb_distance=false;

	std::list<std::pair< core::Size, utility::vector1< core::Size > >  > atoms_with_orbs_and_hydrogen_pairs;

	//this function checks to see if an atom with orbitals and a Hydrogen pair fall within 4A. If it does, it returns
	//a vector with the paired atoms and hydrogens
	get_bonded_bb_orb_atom_and_Hpol_atom(res1, res2, atoms_with_orbs_and_hydrogen_pairs, check_orb_distance);

	bool get_hpol_score=false;
	core::Real min_dist_squared(max_orbital_dist_squared_);
	core::chemical::orbitals::orbital_type_enum orbital_type;
	std::list< orbital_info_for_score > min_orbital_dist_tau;

	if ( check_orb_distance ) {
		min_orbital_dist_tau = get_orb_H_distance_and_angle(res1, res2, atoms_with_orbs_and_hydrogen_pairs, min_dist_squared, orbital_type, get_hpol_score);
	}

	if ( get_hpol_score ) {
		was_scored=true;
		for ( std::list< orbital_info_for_score >::iterator
				orbital_min_dist_tau = min_orbital_dist_tau.begin(),
				orbital_min_dist_tau_end = min_orbital_dist_tau.end();
				orbital_min_dist_tau != orbital_min_dist_tau_end; ++orbital_min_dist_tau ) {
			core::Real distance_derivative(0.0);
			core::Real angle_derivative(0.0);
			get_orbital_H_score(
				lookup_table_.bb,
				orbital_min_dist_tau->orbital_type,
				orbital_min_dist_tau->distance,
				orbital_min_dist_tau->tau,
				sc_bb_E,
				distance_derivative,
				angle_derivative,
				check_derivative );
			orbital_min_dist_tau->distance_derivative=distance_derivative;
			orbital_min_dist_tau->angle_derivative=angle_derivative;
			//if ( check_derivative ) {
			if ( false ) {
				std::cout << "bb_hpol "<< res1.seqpos() << " " << res2.seqpos() << " a_index: " << orbital_min_dist_tau->atom_index
					<< " orb_index: " << orbital_min_dist_tau->orbital_index << " h_index: "<<  orbital_min_dist_tau->h_index
					<< " E: " << sc_bb_E << " d_deriv: " << distance_derivative << " " <<  orbital_min_dist_tau->distance_derivative
					<< " a_deriv: " << angle_derivative << " " << orbital_min_dist_tau->angle_derivative << " dist: "
					<< orbital_min_dist_tau->distance << " angle: " << orbital_min_dist_tau->tau << " enum: " << orbital_min_dist_tau->orbital_type << std::endl;
			}
		}

	}
	return min_orbital_dist_tau;
}





void OrbitalsScore::get_bonded_orb_atom_and_Hpol_atom(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	std::list<std::pair< core::Size, utility::vector1< core::Size > >  > & atoms_with_orbs_and_hydrogen_pairs,
	bool & check_orbital_distance
) const
{
	//std::cout << "ENTERING HPOL" << std::endl;

	core::Real min_dist_squared(max_dist_squared_);
	for ( chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.atoms_with_orb_index().begin(),
			atoms_with_orb_index_end = res1.atoms_with_orb_index().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index ) {
		utility::vector1< core::Size > hydrogens;
		bool add_to_vector=false;
		if ( !res1.atom_is_backbone(*atoms_with_orb_index) ) {
			for ( chemical::AtomIndices::const_iterator
					hpol_index = res2.Hpos_polar_sc().begin(),
					hpol_end = res2.Hpos_polar_sc().end();
					hpol_index != hpol_end; ++hpol_index ) {
				numeric::xyzVector<core::Real> const & atom_xyz = res1.atom(*atoms_with_orb_index).xyz();
				numeric::xyzVector<core::Real> const & H_xyz = res2.atom(*hpol_index).xyz();
				core::Real temp_dist = atom_xyz.distance_squared(H_xyz);
				if ( temp_dist < min_dist_squared ) {
					//min_dist_squared = temp_dist;
					hydrogens.push_back(*hpol_index);
					add_to_vector=true;
				}

			}
		}
		if ( add_to_vector ) {
/*			for(core::Size x=1; x<= hydrogens.size(); ++x ) {
				std::cout << res1.seqpos() << " " << res2.seqpos() << " a_index " << *atoms_with_orb_index << " hpol-index " << hydrogens[x] << std::endl;
			}*/
			std::pair< core::Size, utility::vector1< core::Size > > atom_index_H_index(std::make_pair(*atoms_with_orb_index, hydrogens));
			atoms_with_orbs_and_hydrogen_pairs.push_back(atom_index_H_index);
			check_orbital_distance=true;
		}

	}




}

void
OrbitalsScore::get_bonded_orb_atom_and_Haro_atom(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	std::list<std::pair< core::Size, utility::vector1< core::Size  > >  > & atoms_with_orbs_and_hydrogen_pairs,
	bool & check_orbital_distance
) const
{
	core::Real min_dist_squared(max_dist_squared_);
//	std::cout << "ENTERING HARO" << std::endl;
	for ( chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.atoms_with_orb_index().begin(),
			atoms_with_orb_index_end = res1.atoms_with_orb_index().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index ) {
		utility::vector1< core::Size > hydrogens;
		//std::cout << "NEW_ATOM " << std::endl;
		bool add_to_vector=false;
		if ( !res1.atom_is_backbone(*atoms_with_orb_index) ) {
			for ( chemical::AtomIndices::const_iterator
					haro_index = res2.Haro_index().begin(),
					haro_end = res2.Haro_index().end();
					haro_index != haro_end; ++haro_index ) {
				//std::cout << "h-aro index! " << *haro_index << std::endl;
				numeric::xyzVector<core::Real> const & atom_xyz = res1.atom(*atoms_with_orb_index).xyz();
				numeric::xyzVector<core::Real> const & H_xyz = res2.atom(*haro_index).xyz();
				core::Real temp_dist = atom_xyz.distance_squared(H_xyz);
				if ( temp_dist < min_dist_squared ) {
					hydrogens.push_back(*haro_index);
					add_to_vector=true;
				}

			}
		}
		if ( add_to_vector ) {
/*			for(core::Size x=1; x<= hydrogens.size(); ++x ) {
				std::cout << res1.seqpos() << " " << res2.seqpos() << " a_index " << *atoms_with_orb_index << " haro-index " << hydrogens[x] << std::endl;
			}*/
			std::pair< core::Size, utility::vector1< core::Size > > atom_index_H_index(std::make_pair(*atoms_with_orb_index, hydrogens));
			atoms_with_orbs_and_hydrogen_pairs.push_back(atom_index_H_index);
			check_orbital_distance=true;
		}

	}
}

void OrbitalsScore::get_bonded_bb_orb_atom_and_Hpol_atom(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	std::list<std::pair< core::Size, utility::vector1< core::Size > >  > & atoms_with_orbs_and_hydrogen_pairs,
	bool & check_orbital_distance
) const {
	core::Real min_dist_squared(max_dist_squared_);
	//std::cout << "ENTERING BB-HPOL" << std::endl;

	for ( chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.accpt_pos().begin(),
			atoms_with_orb_index_end = res1.accpt_pos().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index ) {
		bool add_to_vector=false;
		utility::vector1< core::Size > hydrogens;
		if ( res1.atom_is_backbone(*atoms_with_orb_index) ) {
			for ( chemical::AtomIndices::const_iterator
					hpol_index = res2.Hpos_polar_sc().begin(),
					hpol_end = res2.Hpos_polar_sc().end();
					hpol_index != hpol_end; ++hpol_index ) {
				numeric::xyzVector<core::Real> const & atom_xyz = res1.atom(*atoms_with_orb_index).xyz();
				numeric::xyzVector<core::Real> const & H_xyz = res2.atom(*hpol_index).xyz();
				core::Real temp_dist = atom_xyz.distance_squared(H_xyz);
				if ( temp_dist < min_dist_squared ) {
					//min_dist_squared = temp_dist;
					hydrogens.push_back(*hpol_index);
					add_to_vector=true;
				}

			}
		}
		if ( add_to_vector ) {
/*			for(core::Size x=1; x<= hydrogens.size(); ++x ) {
				std::cout << res1.seqpos() << " " << res2.seqpos() << " bb_a_index " << *atoms_with_orb_index << " hpol-index " << hydrogens[x] << std::endl;
			}*/
			std::pair< core::Size, utility::vector1< core::Size > > atom_index_H_index(std::make_pair(*atoms_with_orb_index, hydrogens));
			atoms_with_orbs_and_hydrogen_pairs.push_back(atom_index_H_index);
			check_orbital_distance=true;
		}

	}



}



std::list< orbital_info_for_score >
OrbitalsScore::get_orb_H_distance_and_angle(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	std::list< std::pair< core::Size, utility::vector1< core::Size > > > const & atoms_with_orbs_and_hydrogen_pairs,
	core::Real & min_dist_squared,
	core::chemical::orbitals::orbital_type_enum & orbital_type,
	bool & get_score
) const
{

	min_dist_squared=100.0;
	std::list< orbital_info_for_score > orbital_dist_tau;

	for ( std::list<std::pair< core::Size, utility::vector1< core::Size > >  >::const_iterator
			atoms_hydrogen_pair = atoms_with_orbs_and_hydrogen_pairs.begin(),
			atoms_hydrogen_pair_end = atoms_with_orbs_and_hydrogen_pairs.end();
			atoms_hydrogen_pair != atoms_hydrogen_pair_end; ++atoms_hydrogen_pair ) {
		core::Size atom_index(atoms_hydrogen_pair->first);
		utility::vector1< core::Size > const & H_indices(atoms_hydrogen_pair->second);
		utility::vector1< core::Size > const & orbital_indices(res1.bonded_orbitals(atom_index));

		for ( utility::vector1< core::Size >::const_iterator
				orbital_index = orbital_indices.begin(),
				orbital_index_end = orbital_indices.end();
				orbital_index != orbital_index_end; ++orbital_index ) {
			numeric::xyzVector< core::Real > const orbital_xyz(res1.orbital_xyz(*orbital_index) );
			//std::cout << orbital_xyz.x() << " " <<  orbital_xyz.y() << " " << orbital_xyz.z() << std::endl;
			for ( utility::vector1< core::Size >::const_iterator
					H_index = H_indices.begin(),
					H_index_end = H_indices.end();
					H_index != H_index_end; ++H_index ) {
				numeric::xyzVector< core::Real > const H_xyz =res2.atom( *H_index ).xyz();
				core::Real temp_dist_squared = orbital_xyz.distance_squared( H_xyz );
				if ( temp_dist_squared < min_dist_squared ) {
					//std::cout << "resid: " << res1.seqpos() << " " << res2.seqpos() << " orb_index: " << *orbital_index << std::endl;
					//min_dist_squared = temp_dist_squared;
					orbital_type = res1.orbital_type(*orbital_index).orbital_enum();
					orbital_info_for_score orbital_info_struct;
					get_score=true;
					if ( res2.atom_is_backbone(*H_index) ) {
						orbital_info_struct.lookup_backbone_energy=true;
						//std::cout << "lookup_backbone_energy=true!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
					} else {
						orbital_info_struct.lookup_backbone_energy=false;
					}
					orbital_info_struct.atom_index=atom_index;
					orbital_info_struct.h_index=*H_index;
					orbital_info_struct.orbital_index=*orbital_index;
					orbital_info_struct.distance = std::sqrt(temp_dist_squared);
					orbital_info_struct.orbital_type = orbital_type;
					orbital_info_struct.tau = cos_of(res1.atom(atom_index).xyz(), orbital_xyz, H_xyz);
					orbital_dist_tau.push_back(orbital_info_struct);

				}
			}
		}
	}
	return orbital_dist_tau;

}

//make a lot of checks here. Taking the spline actually takes some time, so its better to avoid checking it if need possible.
void OrbitalsScore::get_orbital_H_score(
	OrbitalsLookup::h_type htype,
	core::chemical::orbitals::orbital_type_enum orbital_type,
	core::Real const & min_dist,
	core::Real const & tau,
	core::Real & energy,
	core::Real & distance_derivative,
	core::Real & angle_derivative,
	bool check_derivative
) const
{
	lookup_table_.get_energy( htype, orbital_type, min_dist, tau, energy,
		distance_derivative, angle_derivative, check_derivative );
}

void
OrbitalsScore::eval_residue_pair_derivatives(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const &, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {

	if ( res1.has_sc_orbitals() || res2.has_sc_orbitals() ) {
		//if ( res1.nbr_atom_xyz().distance_squared(res2.nbr_atom_xyz()) <= nbr_distance_squared_ ) {

		if ( res1.is_aromatic() && res2.is_aromatic() ) {
			std::list< orbital_info_for_score > min_orbital_dist_tau;
			core::Real haro_E=0.0;
			bool check_derivative=true;
			bool was_scored=false;
			min_orbital_dist_tau = get_E_haro_one_way( res1, res2, haro_E, check_derivative, was_scored );
			if ( was_scored ) {
				assign_atom_derivatves( res1, res2, weights, min_orbital_dist_tau, r1_atom_derivs,r2_atom_derivs );
			}
			was_scored=false;

			min_orbital_dist_tau = get_E_haro_one_way( res2, res1, haro_E, check_derivative, was_scored);
			if ( was_scored ) {
				assign_atom_derivatves( res2, res1, weights, min_orbital_dist_tau, r2_atom_derivs,r1_atom_derivs );
			}
		} else {
			//PROF_START( basic::ORBITALS_E_1_WAY );
			std::list< orbital_info_for_score > min_orbital_dist_tau;
			if ( res1.has_sc_orbitals() ) {
				core::Real hpol_E(0.0);
				core::Real bb_hpol_E(0.0);
				bool check_derivative=true;
				bool was_scored=false;
				min_orbital_dist_tau = get_E_hpol_one_way( res1, res2, hpol_E, bb_hpol_E, check_derivative,was_scored);
				if ( was_scored ) {
					assign_atom_derivatves( res1, res2, weights, min_orbital_dist_tau, r1_atom_derivs,r2_atom_derivs );
				}
			}
			if ( res2.has_sc_orbitals() ) {
				core::Real hpol_E(0.0);
				core::Real bb_hpol_E(0.0);
				bool check_derivative=true;
				bool was_scored=false;
				min_orbital_dist_tau = get_E_hpol_one_way( res2, res1, hpol_E, bb_hpol_E, check_derivative, was_scored);
				if ( was_scored ) {
					assign_atom_derivatves( res2, res1, weights, min_orbital_dist_tau, r2_atom_derivs,r1_atom_derivs );
				}
			}
		}

		//commented out until I can get sidechain minimization to work properly
		//we only want to consider those residues that have polar side chain hydrogens. The only atoms on the bb that have
		//orbitals are the carbonyl oxygens. The only orbitals considered are the sp2 orbitals. These are hardcorded in the
		//params file. Maybe one day the ligand code will be put into svn and the params files wont be hardcoded.
/*			if ( res2.is_polar() || res2.name3() == "TYR" || res2.name3() == "TRP" ) {
			std::list< orbital_info_for_score > min_orbital_dist_tau;
			core::Real sc_bb_E(0.0);
			bool check_derivative=true;
			bool was_scored=false;
			min_orbital_dist_tau=get_sc_bb_E_one_way(res1, res2, sc_bb_E, check_derivative, was_scored);
			if ( was_scored ) {
				assign_atom_derivatves(res1, res2, min_orbital_dist_tau, r1_atom_derivs,r2_atom_derivs );
			}
		}
		if ( res1.is_polar() || res1.name3() == "TYR" || res1.name3() == "TRP" ) {
			std::list< orbital_info_for_score > min_orbital_dist_tau;
			core::Real sc_bb_E(0.0);
			bool check_derivative=true;
			bool was_scored=false;
			min_orbital_dist_tau=get_sc_bb_E_one_way(res2, res1, sc_bb_E, check_derivative, was_scored);
			if ( was_scored ) {
				assign_atom_derivatves(res2, res1, min_orbital_dist_tau, r2_atom_derivs,r1_atom_derivs );
			}
		}*/
		//PROF_STOP( basic::ORBITALS_E_1_WAY );

	//}
	}




}


void OrbitalsScore::assign_atom_derivatves(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	EnergyMap const & weights,
	std::list< orbital_info_for_score > min_orbital_dist_tau,
	utility::vector1< DerivVectorPair > & res1_atom_derivs,
	utility::vector1< DerivVectorPair > & res2_atom_derivs
) const
{
	for (
			std::list< orbital_info_for_score >::const_iterator
			orbital_min_dist_tau = min_orbital_dist_tau.begin(),
			orbital_min_dist_tau_end = min_orbital_dist_tau.end();
			orbital_min_dist_tau != orbital_min_dist_tau_end;
			++orbital_min_dist_tau
	)
	{
		Vector pAB = res1.xyz(orbital_min_dist_tau->atom_index);
		Vector pAO = res1.orbital_xyz(orbital_min_dist_tau->orbital_index);
		Vector pH  = res2.xyz(orbital_min_dist_tau->h_index);

		//core::Size orbital_surrogate_atom_index = res1.bonded_neighbor(orbital_min_dist_tau->atom_index)[1];
		core::Size orbital_surrogate_atom_index = orbital_min_dist_tau->atom_index;



		Vector f1ab,f2ab;
		Vector f1ao,f2ao;
		Vector f1h,f2h;
		core::Real tau(0.0);

		numeric::deriv::angle_p1_deriv( pAB, pAO, pH, tau, f1ab, f2ab  );
		numeric::deriv::angle_p2_deriv( pAB, pAO, pH, tau, f1ao, f2ao  );
		numeric::deriv::angle_p1_deriv( pH, pAO, pAB, tau, f1h, f2h  );

		//std::cout << "ANGLE/DISTANCE DERIVATIVES " << orbital_min_dist_tau->angle_derivative << " " << orbital_min_dist_tau->distance_derivative<< std::endl;

		//std::cout << "1) comp_tau: " << tau
		//		<< " struct_tau: " << orbital_min_dist_tau->tau
		//		<< " cos_of: " << cos_of(pAB, pAO, pH)
		//		<<  " angle_of: " << angle_of(pAB, pAO, pH)
		//		<< " dist: "  << orbital_min_dist_tau->distance << " " << res1.seqpos() << " " << res2.seqpos() <<  " a_index " << orbital_min_dist_tau->atom_index
		//		<< " orb_ind " << orbital_min_dist_tau->orbital_index << " h_index " << orbital_min_dist_tau->h_index << std::endl;



		//std::cout << "deriv "<< res1.seqpos() << " " << res2.seqpos() << " a_index: " << orbital_min_dist_tau->atom_index
		//		<< " o_index: " << orbital_min_dist_tau->orbital_index << " h_index: "<<  orbital_min_dist_tau->h_index
		//		<< " E: " << "0" << " d_der: " << "0" << " " <<  orbital_min_dist_tau->distance_derivative
		//		<< " a_der: " << "0" << " " << orbital_min_dist_tau->angle_derivative << " dist: "
		//		<< orbital_min_dist_tau->distance << " angle: " << orbital_min_dist_tau->tau << " enum: " << orbital_min_dist_tau->orbital_type << std::endl;

		Real weight(0.0);
		if (orbital_min_dist_tau->htype == lookup_table_.haro) {
			weight = weights[orbitals_haro];
		} else if (orbital_min_dist_tau->htype == lookup_table_.hpol) {
			weight = weights[orbitals_hpol];
		}

		Real neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;
		//Real neg_sine_tau = 1.0;

		//std::cout << "computed_tau: " << tau << " struct_tau: " << orbital_min_dist_tau->tau << std::endl;
		//std::cout << "angle derivative " << orbital_min_dist_tau->angle_derivative << " dist deriv " << orbital_min_dist_tau->distance_derivative << std::endl;
		res1_atom_derivs[ orbital_min_dist_tau->atom_index ].f1() += weight * neg_sine_tau * orbital_min_dist_tau->angle_derivative * f1ab;
		res1_atom_derivs[ orbital_min_dist_tau->atom_index ].f2() += weight * neg_sine_tau * orbital_min_dist_tau->angle_derivative * f2ab;

		res1_atom_derivs[ orbital_surrogate_atom_index ].f1() += weight * neg_sine_tau * orbital_min_dist_tau->angle_derivative * f1ao;
		res1_atom_derivs[ orbital_surrogate_atom_index ].f2() += weight * neg_sine_tau * orbital_min_dist_tau->angle_derivative * f2ao;

		res2_atom_derivs[ orbital_min_dist_tau->h_index ].f1() += weight * neg_sine_tau * orbital_min_dist_tau->angle_derivative * f1h;
		res2_atom_derivs[ orbital_min_dist_tau->h_index ].f2() += weight * neg_sine_tau * orbital_min_dist_tau->angle_derivative * f2h;


		f1ao = f2ao = 0;
		Real dis(0.0);
		numeric::deriv::distance_f1_f2_deriv( pAO, pH, dis, f1ao, f2ao );

		res1_atom_derivs[ orbital_surrogate_atom_index ].f1() += weight * orbital_min_dist_tau->distance_derivative * f1ao;
		res1_atom_derivs[ orbital_surrogate_atom_index ].f2() += weight * orbital_min_dist_tau->distance_derivative * f2ao;

		res2_atom_derivs[ orbital_min_dist_tau->h_index ].f1() -= weight * orbital_min_dist_tau->distance_derivative * f1ao;
		res2_atom_derivs[ orbital_min_dist_tau->h_index ].f2() -= weight * orbital_min_dist_tau->distance_derivative * f2ao;


	}



}

core::Size
OrbitalsScore::version() const
{
	return 2; // Initial versioning
}



}
}
}
