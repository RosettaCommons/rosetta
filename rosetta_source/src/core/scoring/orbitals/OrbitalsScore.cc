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

#include <core/conformation/RotamerSetBase.hh>

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

// Auto-header: duplicate removed #include <core/pose/Pose.hh>
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
	sts.push_back(orbitals_hpol_bb);
	return sts;
}

static basic::Tracer TR("core.scoring.orbitals_hpol");

OrbitalsScore::OrbitalsScore(methods::EnergyMethodOptions const &) :
					parent( new OrbitalsScoreCreator ),
					lookup_table_(core::scoring::ScoringManager::get_instance()->get_OrbitalsLookupTable()),
					max_dist_squared_(36),
					max_orbital_dist_squared_(11.56)

{
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
OrbitalsScore::prepare_rotamers_for_packing(
	pose::Pose const & /*pose*/,
	conformation::RotamerSetBase & set
) const
{
	for ( Size ii = 1; ii <= set.num_rotamers(); ++ii ) {
		set.nonconst_rotamer( ii )->update_orbital_coords(); // the Rotamer set does not take responsibility for this; why not? -- maybe Residue could?
	}
}

void
OrbitalsScore::update_residue_for_packing(
	pose::Pose & pose,
	Size resid
) const
{
	pose.update_orbital_coords( resid );
}




core::Size
OrbitalsScore::version() const
{
	return 2; // Initial versioning
}

/////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
/////////////////////////////////

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
		if ( res1.is_aromatic() && res2.is_aromatic() ) {
			core::Real HARO_sc_H_sc_orb_E(0.0);
			core::Real HARO_DHO_angle_E(0.0);
			get_E_haro_one_way(res1, res2, HARO_sc_H_sc_orb_E, HARO_DHO_angle_E);
			emap[orbitals_haro] += HARO_sc_H_sc_orb_E/2;
			HARO_sc_H_sc_orb_E = 0.0;
			HARO_DHO_angle_E = 0.0;
			get_E_haro_one_way(res2, res1, HARO_sc_H_sc_orb_E, HARO_DHO_angle_E);
			emap[orbitals_haro] += HARO_sc_H_sc_orb_E/2;
		} else {
			//PROF_START( basic::ORBITALS_E_1_WAY );
			{
				core::Real HPOL_sc_H_sc_orb_E(0.0);
				core::Real HPOL_bb_H_sc_orb_energy(0.0);
				core::Real HPOL_sc_H_bb_orb_energy(0.0);
				core::Real HPOL_DHO_angle_E(0.0);
				get_E_hpol_one_way(res1, res2, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, HPOL_sc_H_bb_orb_energy, HPOL_DHO_angle_E);
				emap[orbitals_hpol] += HPOL_sc_H_sc_orb_E/2;
				emap[orbitals_hpol_bb] += (HPOL_bb_H_sc_orb_energy+HPOL_sc_H_bb_orb_energy)/4;


			}
			{
				core::Real HPOL_sc_H_sc_orb_E(0.0);
				core::Real HPOL_bb_H_sc_orb_energy(0.0);
				core::Real HPOL_sc_H_bb_orb_energy(0.0);
				core::Real HPOL_DHO_angle_E(0.0);
				get_E_hpol_one_way(res2, res1, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, HPOL_sc_H_bb_orb_energy, HPOL_DHO_angle_E);
				emap[orbitals_hpol] += HPOL_sc_H_sc_orb_E/2;//divide by two because DHO and AOH angle energies included
				emap[orbitals_hpol_bb] += (HPOL_bb_H_sc_orb_energy+HPOL_sc_H_bb_orb_energy)/4;//divide by 4 because DHO and AOH angle energies included
			}
		}
	}
}

void OrbitalsScore::get_E_haro_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & HARO_sc_H_sc_orb_E,
		core::Real & HARO_DHO_angle_E
) const
{
	core::Real dummy_E1(0.0);
	core::Real dummy_E2(0.0);
	core::Real max_dist_squared(max_dist_squared_);
	for (
			chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.atoms_with_orb_index().begin(),
			atoms_with_orb_index_end = res1.atoms_with_orb_index().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index
	)
	{
		if ( !res1.atom_is_backbone(*atoms_with_orb_index) ) {
			for (
					chemical::AtomIndices::const_iterator
					haro_index = res2.Haro_index().begin(),
					haro_end = res2.Haro_index().end();
					haro_index != haro_end; ++haro_index
			)
			{
				numeric::xyzVector<core::Real> const & atom_xyz = res1.atom(*atoms_with_orb_index).xyz();
				numeric::xyzVector<core::Real> const & H_xyz = res2.atom(*haro_index).xyz();
				core::Size donor_id(res2.bonded_neighbor(*haro_index)[1]);
				numeric::xyzVector<core::Real> const & donor_xyz = res2.xyz(donor_id);

				core::Real temp_dist = atom_xyz.distance_squared(H_xyz);
				if ( temp_dist < max_dist_squared ) {
					core::Size atom_index(*atoms_with_orb_index);
					get_orb_H_distance_and_energy(res1, atom_index, atom_xyz, H_xyz, donor_xyz, HARO_sc_H_sc_orb_E, dummy_E1, dummy_E2, HARO_DHO_angle_E, lookup_table_.haro, false, false);
				}

			}
		}
	}
}


void OrbitalsScore::get_E_hpol_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & HPOL_sc_H_sc_orb_E,
		core::Real & HPOL_bb_H_sc_orb_energy,
		core::Real & HPOL_sc_H_bb_orb_energy,
		core::Real & HPOL_DHO_angle_E
) const
{
	core::Real max_dist_squared(max_dist_squared_);
	for (
			chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.atoms_with_orb_index().begin(),
			atoms_with_orb_index_end = res1.atoms_with_orb_index().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index
	)
	{
		numeric::xyzVector<core::Real> const & atom_xyz = res1.atom(*atoms_with_orb_index).xyz();
			for (
					chemical::AtomIndices::const_iterator
					hpol_index = res2.Hpol_index().begin(),
					hpol_end = res2.Hpol_index().end();
					hpol_index != hpol_end; ++hpol_index
			)
			{
				//this check is to look at bb orbital bb hydrogen. This potential does not calculate it.
				//The hbond_lr_bb and hbond_sr_bb scoring terms look into this.
				if(res1.atom_is_backbone(*atoms_with_orb_index) && res2.atom_is_backbone(*hpol_index)){
					continue;
				}

				numeric::xyzVector<core::Real> const & H_xyz = res2.atom(*hpol_index).xyz();
				core::Size donor_id(res2.bonded_neighbor(*hpol_index)[1]);
				numeric::xyzVector<core::Real> const & donor_xyz = res2.xyz(donor_id);
				core::Real temp_dist = atom_xyz.distance_squared(H_xyz);
				if ( temp_dist < max_dist_squared ) {
					core::Size atom_index(*atoms_with_orb_index);
					core::Size H_index(*hpol_index);
					if(res1.atom_is_backbone(atom_index)){
						get_orb_H_distance_and_energy(res1, atom_index, atom_xyz, H_xyz, donor_xyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, HPOL_sc_H_bb_orb_energy, HPOL_DHO_angle_E, lookup_table_.bb, false, true);
					}
					if(res2.atom_is_backbone(*hpol_index)){
						get_orb_H_distance_and_energy(res1, atom_index, atom_xyz, H_xyz, donor_xyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, HPOL_sc_H_bb_orb_energy, HPOL_DHO_angle_E, lookup_table_.sc_orb_bb_H, true, false);
					}else{
						get_orb_H_distance_and_energy(res1, atom_index, atom_xyz, H_xyz, donor_xyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, HPOL_sc_H_bb_orb_energy, HPOL_DHO_angle_E, lookup_table_.hpol, false, false);
					}
				}
			}
	}
}

void OrbitalsScore::get_orb_H_distance_and_energy(
		core::conformation::Residue const & res1,
		core::Size const & atom_index,
		numeric::xyzVector<core::Real> const & atom_xyz,
		numeric::xyzVector<core::Real> const & H_xyz,
		numeric::xyzVector<core::Real> const & donor_xyz,
		core::Real & sc_energy,
		core::Real & bb_h_energy,
		core::Real & bb_orb_energy,
		core::Real & DHO_angle_E,
		OrbitalsLookup::h_type htype,
		bool bb_h_flag,
		bool bb_orb_flag
) const
{
	core::Real d_deriv(0.0);
	core::Real a_deriv(0.0);
	utility::vector1< core::Size > const & orbital_indices(res1.bonded_orbitals(atom_index));
	core::Real added_sc_energy(sc_energy);
	core::Real added_bb_h_energy(bb_h_energy);
	core::Real added_bb_orb_energy(bb_orb_energy);
	for(
			utility::vector1< core::Size >::const_iterator
			orbital_index = orbital_indices.begin(),
			orbital_index_end = orbital_indices.end();
			orbital_index != orbital_index_end; ++orbital_index
	)
	{
		numeric::xyzVector< core::Real > const & orbital_xyz(res1.orbital_xyz(*orbital_index) );
		core::Real temp_dist_squared = orbital_xyz.distance_squared( H_xyz );
		if(temp_dist_squared < max_orbital_dist_squared_){
			core::chemical::orbitals::orbital_type_enum orbital_type = res1.orbital_type(*orbital_index).orbital_enum();
			if(lookup_table_.check_distance(temp_dist_squared, htype, orbital_type)){
				core::Real AOH_angle(cos_of(atom_xyz, orbital_xyz, H_xyz));//atom-orbital-hydrogen angle
				if(lookup_table_.check_AOH_angle(AOH_angle, htype, orbital_type)){
					core::Real dist_squared(std::sqrt(temp_dist_squared));
					core::Real DHO_angle(cos_of(donor_xyz, H_xyz, orbital_xyz));
					if(bb_h_flag){
						lookup_table_.get_dist_AOH_angle_energy(htype, orbital_type, dist_squared, AOH_angle, bb_h_energy, d_deriv, a_deriv, false);
						lookup_table_.get_DHO_angle_energy(htype, orbital_type, DHO_angle, DHO_angle_E, a_deriv, false);
						added_bb_h_energy += bb_h_energy+DHO_angle_E;
					}else if(bb_orb_flag){
						lookup_table_.get_dist_AOH_angle_energy(htype, orbital_type, dist_squared, AOH_angle, bb_orb_energy, d_deriv, a_deriv, false);
						lookup_table_.get_DHO_angle_energy(htype, orbital_type, DHO_angle, DHO_angle_E, a_deriv, false);
						added_bb_orb_energy += bb_orb_energy+DHO_angle_E;
					}else{
						lookup_table_.get_dist_AOH_angle_energy(htype, orbital_type, dist_squared, AOH_angle, sc_energy, d_deriv, a_deriv, false);
						lookup_table_.get_DHO_angle_energy(htype, orbital_type, DHO_angle, DHO_angle_E, a_deriv, false);
						added_sc_energy += sc_energy+DHO_angle_E;
					}
				}
			}
		}
	}
	bb_h_energy = added_bb_h_energy;
	sc_energy = added_sc_energy;
	bb_orb_energy = added_bb_orb_energy;
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
		if ( res1.is_aromatic() && res2.is_aromatic() ) {
			assign_haro_derivs_one_way(res1, res2, weights, r1_atom_derivs, r2_atom_derivs);
			assign_haro_derivs_one_way(res2, res1, weights, r2_atom_derivs, r1_atom_derivs);

		} else {
			{
				assign_hpol_derivs_one_way(res1, res2, weights, r1_atom_derivs, r2_atom_derivs);
			}
			{
				assign_hpol_derivs_one_way(res2, res1, weights, r2_atom_derivs, r1_atom_derivs);
			}
		}
	}
}




void OrbitalsScore::assign_haro_derivs_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	core::Real max_dist_squared(max_dist_squared_);
	for (
			chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.atoms_with_orb_index().begin(),
			atoms_with_orb_index_end = res1.atoms_with_orb_index().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index
	)
	{
		if ( !res1.atom_is_backbone(*atoms_with_orb_index) ) {
			for (
					chemical::AtomIndices::const_iterator
					haro_index = res2.Haro_index().begin(),
					haro_end = res2.Haro_index().end();
					haro_index != haro_end; ++haro_index
			)
			{
				numeric::xyzVector<core::Real> const & atom_xyz = res1.atom(*atoms_with_orb_index).xyz();
				numeric::xyzVector<core::Real> const & H_xyz = res2.atom(*haro_index).xyz();
				core::Real temp_dist = atom_xyz.distance_squared(H_xyz);
				if ( temp_dist < max_dist_squared ) {
					core::Size atom_index(*atoms_with_orb_index);
					core::Size H_index(*haro_index);
					assign_orb_H_derivs(
							res1, res2, atom_index, atom_xyz, H_index, H_xyz,
							lookup_table_.haro, weights, r1_atom_derivs, r2_atom_derivs
					);
				}
			}
		}
	}
}

void OrbitalsScore::assign_hpol_derivs_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs

) const
{
	core::Real max_dist_squared(max_dist_squared_);
	for (
			chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.atoms_with_orb_index().begin(),
			atoms_with_orb_index_end = res1.atoms_with_orb_index().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index
	)
	{

			for (
					chemical::AtomIndices::const_iterator
					hpol_index = res2.Hpol_index().begin(),
					hpol_end = res2.Hpol_index().end();
					hpol_index != hpol_end; ++hpol_index
			)
			{
				if(res1.atom_is_backbone(*atoms_with_orb_index) && res2.atom_is_backbone(*hpol_index)){
					continue;
				}
				numeric::xyzVector<core::Real> const & atom_xyz = res1.atom(*atoms_with_orb_index).xyz();
				numeric::xyzVector<core::Real> const & H_xyz = res2.atom(*hpol_index).xyz();
				core::Real temp_dist = atom_xyz.distance_squared(H_xyz);
				if ( temp_dist < max_dist_squared ) {
					core::Size atom_index(*atoms_with_orb_index);
					core::Size H_index(*hpol_index);
					assign_orb_H_derivs(
							res1, res2, atom_index, atom_xyz, H_index, H_xyz,
							lookup_table_.hpol, weights, r1_atom_derivs, r2_atom_derivs
					);
				}
			}

	}
}


void OrbitalsScore::assign_orb_H_derivs(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Size const & atom_index,
		numeric::xyzVector<core::Real> const & atom_xyz,
		core::Size const & H_index,
		numeric::xyzVector<core::Real> const & H_xyz,
		OrbitalsLookup::h_type htype,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
)const
{
	core::Real energy(0.0);
	core::Real d_deriv(0.0); //distance derivative
	core::Real a_deriv(0.0); //angle derivative

	core::Real max_orbital_dist_squared=max_orbital_dist_squared_;
	utility::vector1< core::Size > const & orbital_indices(res1.bonded_orbitals(atom_index));
	core::Real added_energy(energy);
	for(
			utility::vector1< core::Size >::const_iterator
			orbital_index = orbital_indices.begin(),
			orbital_index_end = orbital_indices.end();
			orbital_index != orbital_index_end; ++orbital_index
	)
	{
		numeric::xyzVector< core::Real > const orbital_xyz(res1.orbital_xyz(*orbital_index) );
		core::Real temp_dist_squared = orbital_xyz.distance_squared( H_xyz );
		if(temp_dist_squared < max_orbital_dist_squared){
			core::chemical::orbitals::orbital_type_enum orbital_type = res1.orbital_type(*orbital_index).orbital_enum();
			core::Real angle(cos_of(atom_xyz, orbital_xyz, H_xyz));
			core::Real dist_squared(std::sqrt(temp_dist_squared));

			if(res2.atom_is_backbone(H_index)){
				htype=lookup_table_.sc_orb_bb_H;
			}
			if(res1.atom_is_backbone(atom_index)){
				htype=lookup_table_.bb;
			}


			lookup_table_.get_dist_AOH_angle_energy(htype, orbital_type, dist_squared, angle, energy, d_deriv, a_deriv, true);



			Vector pAB(atom_xyz);
			Vector pAO(orbital_xyz);
			Vector pH(H_xyz);

			//core::Size orbital_surrogate_atom_index = res1.bonded_neighbor(orbital_min_dist_tau->atom_index)[1];
			core::Size orbital_surrogate_atom_index(atom_index);



			Vector f1ab,f2ab;
			Vector f1ao,f2ao;
			Vector f1h,f2h;
			core::Real tau(0.0);

			numeric::deriv::angle_p1_deriv( pAB, pAO, pH, tau, f1ab, f2ab  );
			numeric::deriv::angle_p2_deriv( pAB, pAO, pH, tau, f1ao, f2ao  );
			numeric::deriv::angle_p1_deriv( pH, pAO, pAB, tau, f1h, f2h  );

			Real weight(0.0);
			if (htype == lookup_table_.haro) {
				weight = weights[orbitals_haro];
			}else if(htype== lookup_table_.bb || htype== lookup_table_.sc_orb_bb_H){
				weight = weights[orbitals_hpol_bb];
			}else if (htype == lookup_table_.hpol) {
				weight = weights[orbitals_hpol];
			}

			Real neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;

			r1_atom_derivs[ atom_index ].f1() += weight * neg_sine_tau * a_deriv * f1ab;
			r1_atom_derivs[ atom_index ].f2() += weight * neg_sine_tau * a_deriv * f2ab;

			r1_atom_derivs[ orbital_surrogate_atom_index ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
			r1_atom_derivs[ orbital_surrogate_atom_index ].f2() += weight * neg_sine_tau * a_deriv * f2ao;

			r2_atom_derivs[ H_index ].f1() += weight * neg_sine_tau * a_deriv * f1h;
			r2_atom_derivs[ H_index ].f2() += weight * neg_sine_tau * a_deriv * f2h;


			f1ao = f2ao = 0;
			Real dis(0.0);
			numeric::deriv::distance_f1_f2_deriv( pAO, pH, dis, f1ao, f2ao );

			r1_atom_derivs[ orbital_surrogate_atom_index ].f1() += weight * d_deriv * f1ao;
			r1_atom_derivs[ orbital_surrogate_atom_index ].f2() += weight * d_deriv * f2ao;

			r2_atom_derivs[ H_index ].f1() -= weight * d_deriv * f1ao;
			r2_atom_derivs[ H_index ].f2() -= weight * d_deriv * f2ao;





		}
	}


}


}
}
}
