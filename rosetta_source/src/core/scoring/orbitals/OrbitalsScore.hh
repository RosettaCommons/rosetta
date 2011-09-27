// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_scoring_orbitals_OrbitalsScore_hh
#define INCLUDED_core_scoring_orbitals_OrbitalsScore_hh

#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/pose/Pose.fwd.hh>
#include <numeric/xyzVector.hh>
#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/scoring/orbitals/OrbitalsScore.fwd.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/orbitals/OrbitalTypeMapper.fwd.hh>


#include <core/scoring/ScoreFunction.fwd.hh>
#include <map>
#include <list>

namespace core{
namespace scoring{
namespace orbitals{

struct orbital_info_for_score{ // including default values
	OrbitalsLookup::h_type htype;
	core::chemical::orbitals::orbital_type_enum orbital_type;
	core::Real distance;
	core::Real tau;
	core::Size atom_index;
	core::Size h_index;
	core::Size orbital_index;
	bool lookup_backbone_energy;
	core::Real energy;
	core::Real distance_derivative;
	core::Real angle_derivative;
};

class OrbitalsScore : public methods::ContextDependentTwoBodyEnergy {
public:
	typedef methods::ContextDependentTwoBodyEnergy parent;


//virtual functions from score functions
public:
	OrbitalsScore(methods::EnergyMethodOptions const &);

	//clone
	virtual methods::EnergyMethodOP clone() const;

	virtual
	void setup_for_scoring(pose::Pose & pose, ScoreFunction const &) const;

	virtual
	void setup_for_derivatives( pose::Pose &pose, ScoreFunction const &  ) const;

	virtual
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const &,
		pose::Pose const &, // provides context
		EnergyMap const &,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	virtual
	void eval_intrares_energy(
		core::conformation::Residue const &,
		core::pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const;

	 virtual
	 bool defines_intrares_energy(const core::scoring::EnergyMap&) const; //false

	 virtual
	 void residue_pair_energy(
		 core::conformation::Residue const & res1,
		 core::conformation::Residue const & res2,
		 core::pose::Pose const &,
		 core::scoring::ScoreFunction const &,
		 EnergyMap & emap
	 ) const;

	virtual
	core::Real atomic_interaction_cutoff() const; //set to default

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > &  ) const;

	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const {
		return false;
	}

public:
	std::list< orbital_info_for_score >
	get_E_haro_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & haro_E,
		bool check_derivative,
		bool & was_scored
	) const;

	std::list< orbital_info_for_score >
	get_E_hpol_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & hpol_E,
		core::Real & sc_bb_E,
		bool check_derivative,
		bool & was_scored
	) const;


	std::list< orbital_info_for_score >
	get_sc_bb_E_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & sc_bb_E,
		bool check_derivative,
		bool & was_scored
	) const;


	void add_energies_to_emap(EnergyMap & emap, std::list< orbital_info_for_score > & orbital_hydrogen_energies) const;

	///@brief places index of an atom with orbital into variable:  atom_index_w_orbital. Places index of Hpol atom into: H_atom_index.
	/// Based on shortest distance between atom and hydrogen
	/// will also trigger check orbital distance if distance between atom and H are within max_dist_squared_
	void get_bonded_orb_atom_and_Hpol_atom(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		std::list<std::pair<core::Size, utility::vector1<core::Size > >  > & atoms_with_orbs_and_hydrogen_pairs,
		bool & check_orbital_distance
	) const;

	void get_bonded_orb_atom_and_Haro_atom(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		std::list< std::pair< core::Size, utility::vector1< core::Size > > > & atoms_with_orbs_and_hydrogen_pairs,
		bool & check_orbital_distance
	) const;

	void get_bonded_bb_orb_atom_and_Hpol_atom(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		std::list< std::pair< core::Size, utility::vector1< core::Size > > > & atoms_with_orbs_and_hydrogen_pairs,
		bool & check_orbital_distance
	) const;


	///@brief given an atom with an orbital and a hydrogen, get the shortest distance between the orbitals on the atom and the hydrogen. Basically,
	/// iterate through all orbitals on a given atom with a H on other residue. Hydrogen and atom index determined from previous functions.
	std::list< orbital_info_for_score >
	get_orb_H_distance_and_angle(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		const std::list<std::pair<core::Size, utility::vector1<core::Size > >  > & atoms_with_orbs_and_hydrogen_pairs,
		core::Real & min_dist_squared,
		core::chemical::orbitals::orbital_type_enum & orbital_type,
		bool & get_score
	) const;


	///@brief the htype and orbital_type are dertimened by whether we are looking for hpol or haro atoms
	/// to a given orbital_type. orbital_type is currently limited to 7 different types and in this context is an enum
	/// min_dist_squared are determined from get_orb_H_distance_and_angle() function.
	void get_orbital_H_score(
		OrbitalsLookup::h_type htype,
		core::chemical::orbitals::orbital_type_enum orbital_type,
		core::Real const & min_dist_squared,
		core::Real const & tau,
		core::Real & energy,
		core::Real & distance_derivative,
		core::Real & angle_derivative,
		bool check_derivative
	) const;


	void assign_atom_derivatves(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap const & weights,
		std::list< orbital_info_for_score > min_orbital_dist_tau,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

//virtual private functions
private:
	virtual core::Size version() const;

private:
	OrbitalsLookup const & lookup_table_;
	//core::scoring::etable::EtableEnergy lookup_Etable_;
	core::Real max_dist_squared_; //the maximum distance squared which orbitals are scored
	//core::scoring::etable::EtableCAP etable_;
	core::Real max_orbital_dist_squared_;//defaault 4A or 16A squared. based on statistics
	core::Real nbr_distance_squared_;//default 10A or 100A squared
	std::map<core::chemical::orbitals::orbital_type_enum, core::Real> min_orb_dist_enum_map_;
	core::Real Cpisp2_min2_;
	core::Real Npisp2_min2_;
	core::Real Npsp2_min2_;
	core::Real Opisp2_min2_;
	core::Real Opsp2_min2_;
	core::Real Opsp3_min2_;
	core::Real Spsp3_min2_;

};





}//namespace orbitals
}//namespace scoring
}//namespace core



#endif /* ORBITALSSCORE_HH_ */
