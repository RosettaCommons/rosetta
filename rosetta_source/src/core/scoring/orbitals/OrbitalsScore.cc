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

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <core/id/types.hh>
#include <core/kinematics/Jump.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

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
					max_orbital_dist_squared_(9)

{
	if(basic::options::option[ basic::options::OptionKeys::in::add_orbitals] != 1){
		utility_exit_with_message( "Trying to run features test without orbitals! Pass the flag -add_orbitals!" );
	}
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

void OrbitalsScore::setup_for_derivatives(pose::Pose & pose , ScoreFunction const &) const
{
	for(core::Size resid=1; resid <= pose.n_residue(); ++resid){
		pose.update_orbital_coords(resid);
	}

}
void
OrbitalsScore::finalize_after_derivatives( pose::Pose & /*pose*/, ScoreFunction const &  ) const{
/*	for(core::Size resid=1; resid <= pose.n_residue(); ++resid){
		pose.update_orbital_coords(resid);
	}*/
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

bool OrbitalsScore::cation_pi_rules(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Size const & Aindex,
		core::Size const & Dindex
)const
{
	if(res1.atom_type(Aindex).atom_type_name() == "COO" && res2.atom_type(Dindex).atom_type_name() == "OH") return false;
	if(res1.atom_type(Aindex).atom_type_name() == "aroC" && res2.atom_type(Dindex).atom_type_name() == "OH") return false;
	return true;
}


void
OrbitalsScore::residue_pair_energy(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::pose::Pose const &,
		core::scoring::ScoreFunction const &,
		EnergyMap & emap
) const
{
	if(res1.is_aromatic()) {
		core::Real HARO_scHscOrb_E(0.0);
		get_E_haro_one_way(res2, res1, HARO_scHscOrb_E);
		emap[orbitals_haro] += HARO_scHscOrb_E;
	}
	if(res2.is_aromatic()){
		core::Real HARO_scHscOrb_E(0.0);
		get_E_haro_one_way(res1, res2, HARO_scHscOrb_E);
		emap[orbitals_haro] += HARO_scHscOrb_E;
	}

	{
		core::Real HPOL_scHscOrb_E(0.0);
		core::Real HPOL_bb_E(0.0);
		get_E_hpol_one_way(res1, res2, HPOL_scHscOrb_E, HPOL_bb_E);
		emap[orbitals_hpol] += HPOL_scHscOrb_E;
		emap[orbitals_hpol_bb] += HPOL_bb_E;
	}
	{
		core::Real HPOL_scHscOrb_E(0.0);
		core::Real HPOL_bb_E(0.0);
		get_E_hpol_one_way(res2, res1, HPOL_scHscOrb_E, HPOL_bb_E);
		emap[orbitals_hpol] += HPOL_scHscOrb_E;
		emap[orbitals_hpol_bb] += HPOL_bb_E;
	}

}

void OrbitalsScore::get_E_haro_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & HARO_scHscOrb_E
) const
{
	core::Real dummy_E1(0.0);//needed for generalized function get_orb_H_distance_and_energy
	core::Real max_dist_squared(max_dist_squared_);
	for (
			chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.atoms_with_orb_index().begin(),
			atoms_with_orb_index_end = res1.atoms_with_orb_index().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index
	)
	{
		numeric::xyzVector<core::Real> const & Axyz = res1.atom(*atoms_with_orb_index).xyz();//acceptor xyz
		if ( !res1.atom_is_backbone(*atoms_with_orb_index) ) {
			for (
					chemical::AtomIndices::const_iterator
					haro_index = res2.Haro_index().begin(),
					haro_end = res2.Haro_index().end();
					haro_index != haro_end; ++haro_index
			)
			{
				numeric::xyzVector<core::Real> const & Hxyz = res2.atom(*haro_index).xyz(); //hydrogen xyz
				core::Real temp_dist = Axyz.distance_squared(Hxyz);
				if ( temp_dist < max_dist_squared ) {
					core::Size Aindex(*atoms_with_orb_index); //acceptor index
					core::Size donor_index(res2.bonded_neighbor(*haro_index)[1]);
					numeric::xyzVector<core::Real> const & Dxyz(res2.xyz(donor_index)); //donor xyz
					get_orb_H_distance_and_energy(res1, Aindex, Axyz, Hxyz, Dxyz, HARO_scHscOrb_E, dummy_E1, lookup_table_.Haro_scOrbH, false);
				}
			}
		}
	}
}


void OrbitalsScore::get_E_hpol_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & HPOL_sc_H_sc_orb_E,
		core::Real & HPOL_bb_H_sc_orb_energy
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
		numeric::xyzVector<core::Real> const & Axyz = res1.atom(*atoms_with_orb_index).xyz(); //acceptor xyz
		for (
				chemical::AtomIndices::const_iterator
				hpol_index = res2.Hpol_index().begin(),
				hpol_end = res2.Hpol_index().end();
				hpol_index != hpol_end; ++hpol_index
		)
		{
			//this check is to look at bb orbital bb hydrogen. This potential does not calculate it.
			//The hbond_lr_bb and hbond_sr_bb scoring terms look into this.
			core::Size donor_index(res2.bonded_neighbor(*hpol_index)[1]);
			if(res1.atom_is_backbone(*atoms_with_orb_index) && res2.atom_is_backbone(donor_index)){
				continue;
			}
			if(cation_pi_rules(res1, res2, *atoms_with_orb_index, donor_index) == false){
				continue;
			}

			numeric::xyzVector<core::Real> const & Hxyz = res2.atom(*hpol_index).xyz(); //hydrogen xyz
			core::Size donor_id(res2.bonded_neighbor(*hpol_index)[1]);
			core::Real temp_dist = Axyz.distance_squared(Hxyz);
			if ( temp_dist < max_dist_squared ) {
				core::Size Aindex(*atoms_with_orb_index);
				numeric::xyzVector<core::Real> const & Dxyz(res2.xyz(donor_index)); //donor xyz
				if(res2.atom_is_backbone(donor_index) || res1.atom_is_backbone(Aindex)){
					get_orb_H_distance_and_energy(res1, Aindex, Axyz, Hxyz, Dxyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, lookup_table_.Hpol_bbOrbH, true);
				}else{
					get_orb_H_distance_and_energy(res1, Aindex, Axyz, Hxyz, Dxyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, lookup_table_.Hpol_scOrbH, false);
				}
			}
		}
	}
}

void OrbitalsScore::get_orb_H_distance_and_energy(
		core::conformation::Residue const & res1,
		core::Size const & Aindex,
		numeric::xyzVector<core::Real> const & Axyz, //acceptor xyz
		numeric::xyzVector<core::Real> const & Hxyz,//hydrogen xyz
		numeric::xyzVector<core::Real> const & Dxyz, //donor xyz
		core::Real & sc_energy,
		core::Real & bb_h_energy,
		OrbitalsLookup::h_type htype,
		bool bb_h_flag
) const
{
	core::Real d_deriv(0.0);
	core::Real a_deriv(0.0);
	utility::vector1< core::Size > const & orbital_indices(res1.bonded_orbitals(Aindex));
	core::Real added_sc_energy(sc_energy);
	core::Real added_bb_h_energy(bb_h_energy);
	for(
			utility::vector1< core::Size >::const_iterator
			orbital_index = orbital_indices.begin(),
			orbital_index_end = orbital_indices.end();
			orbital_index != orbital_index_end; ++orbital_index
	)
	{
		numeric::xyzVector< core::Real > const & Orbxyz(res1.orbital_xyz(*orbital_index) );
		core::Real temp_dist_squared = Orbxyz.distance_squared( Hxyz );
		if(temp_dist_squared < max_orbital_dist_squared_){
			core::chemical::orbitals::orbital_type_enum orbital_type = res1.orbital_type(*orbital_index).orbital_enum();
			if(orbital_type==core::chemical::orbitals::O_pi_sp2_bb){
				orbital_type=core::chemical::orbitals::O_pi_sp2;
			}
			if(orbital_type==core::chemical::orbitals::O_p_sp2_bb){
				orbital_type=core::chemical::orbitals::O_p_sp2;
			}
			core::Real cosDHO(cos_of(Dxyz, Hxyz, Orbxyz));//Donor - Hydrogen - Orbital angle
			core::Real cosAOH(cos_of(Axyz, Orbxyz, Hxyz));
			core::Real dist(std::sqrt(temp_dist_squared));
			if(bb_h_flag){
				lookup_table_.OrbHdist_cosDHO_energy(htype, orbital_type, dist, cosDHO, bb_h_energy, d_deriv, a_deriv, false);
				added_bb_h_energy += bb_h_energy;
				bb_h_energy=0;
				lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, dist, cosAOH, bb_h_energy, d_deriv, a_deriv, false);
				added_bb_h_energy += bb_h_energy;
			}else{
				lookup_table_.OrbHdist_cosDHO_energy(htype, orbital_type, dist, cosDHO, sc_energy, d_deriv, a_deriv, false);
				added_sc_energy += sc_energy;
				sc_energy=0;
				lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, dist, cosAOH, sc_energy, d_deriv, a_deriv, false);
				added_sc_energy += sc_energy;
			}
		}
	}
	bb_h_energy = added_bb_h_energy;
	sc_energy = added_sc_energy;
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

	if(res1.is_aromatic()) {
		assign_haro_derivs_one_way(res2, res1, weights, r2_atom_derivs, r1_atom_derivs);
	}
	if(res2.is_aromatic()){
		assign_haro_derivs_one_way(res1, res2, weights, r1_atom_derivs, r2_atom_derivs);
	}
	assign_hpol_derivs_one_way(res1, res2, weights, r1_atom_derivs, r2_atom_derivs);
	assign_hpol_derivs_one_way(res2, res1, weights, r2_atom_derivs, r1_atom_derivs);

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
				numeric::xyzVector<core::Real> const & Axyz = res1.atom(*atoms_with_orb_index).xyz();
				numeric::xyzVector<core::Real> const & Hxyz = res2.atom(*haro_index).xyz();
				core::Real temp_dist = Axyz.distance_squared(Hxyz);
				if ( temp_dist < max_dist_squared ) {
					core::Size atom_index(*atoms_with_orb_index);
					core::Size H_index(*haro_index);
					assign_orb_H_derivs(
							res1, res2, atom_index, Axyz, H_index, Hxyz,
							lookup_table_.Haro_scOrbH, weights, r1_atom_derivs, r2_atom_derivs
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
			core::Size donor_index(res2.bonded_neighbor(*hpol_index)[1]);
			if(res1.atom_is_backbone(*atoms_with_orb_index) && res2.atom_is_backbone(donor_index)){
				continue;
			}
			if(cation_pi_rules(res1, res2, *atoms_with_orb_index, donor_index) == false){
				continue;
			}
			numeric::xyzVector<core::Real> const & Axyz = res1.atom(*atoms_with_orb_index).xyz();
			numeric::xyzVector<core::Real> const & Hxyz = res2.atom(*hpol_index).xyz();
			core::Real temp_dist = Axyz.distance_squared(Hxyz);
			if ( temp_dist < max_dist_squared ) {
				core::Size atom_index(*atoms_with_orb_index);
				core::Size H_index(*hpol_index);
				if(res2.atom_is_backbone(donor_index) || res1.atom_is_backbone(atom_index)){
					assign_orb_H_derivs(
							res1, res2, atom_index, Axyz, H_index, Hxyz,
							lookup_table_.Hpol_bbOrbH, weights, r1_atom_derivs, r2_atom_derivs
					);
				}else{
					assign_orb_H_derivs(
							res1, res2, atom_index, Axyz, H_index, Hxyz,
							lookup_table_.Hpol_scOrbH, weights, r1_atom_derivs, r2_atom_derivs
					);
				}
			}
		}

	}
}


void OrbitalsScore::assign_orb_H_derivs(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Size  & atom_index,
		numeric::xyzVector<core::Real> const & Axyz,
		core::Size const & H_index,
		numeric::xyzVector<core::Real> const & Hxyz,
		OrbitalsLookup::h_type htype,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
)const
{
	core::Real energy(0.0);
	core::Real d_deriv(0.0); //distance derivative
	core::Real a_deriv(0.0); //angle derivative

	Real weight(0.0);
	if (htype == lookup_table_.Haro_scOrbH) {
		weight = weights[orbitals_haro];
	}else if(htype== lookup_table_.Hpol_bbOrbH ){
		weight = weights[orbitals_hpol_bb];
	}else if (htype == lookup_table_.Hpol_scOrbH) {
		weight = weights[orbitals_hpol];
	}

	utility::vector1< core::Size > const & orbital_indices(res1.bonded_orbitals(atom_index));
	for(
			utility::vector1< core::Size >::const_iterator
			orbital_index = orbital_indices.begin(),
			orbital_index_end = orbital_indices.end();
			orbital_index != orbital_index_end; ++orbital_index
	)
	{
		numeric::xyzVector< core::Real > const Orbxyz(res1.orbital_xyz(*orbital_index) );
		core::Real temp_dist_squared = Orbxyz.distance_squared( Hxyz );
		if(temp_dist_squared < max_orbital_dist_squared_){
			core::chemical::orbitals::orbital_type_enum orbital_type = res1.orbital_type(*orbital_index).orbital_enum();

			if(orbital_type==core::chemical::orbitals::O_pi_sp2_bb){
				orbital_type=core::chemical::orbitals::O_pi_sp2;
			}
			if(orbital_type==core::chemical::orbitals::O_p_sp2_bb){
				orbital_type=core::chemical::orbitals::O_p_sp2;
			}



			//This starts the DHO derivative calculation.
			core::Size donor_index(res2.bonded_neighbor(H_index)[1]);
			numeric::xyzVector<core::Real> const & Dxyz(res2.xyz(donor_index));
			core::Real cosDHO(cos_of(Dxyz, Hxyz, Orbxyz ));
			core::Real OrbHdist(std::sqrt(temp_dist_squared));

			lookup_table_.OrbHdist_cosDHO_energy(htype, orbital_type, OrbHdist, cosDHO, energy, d_deriv, a_deriv, true);


			Vector pD(Dxyz);
			Vector pDH(Hxyz);
			Vector pO(Orbxyz);


			Vector f1ab,f2ab;
			Vector f1ao,f2ao;
			Vector f1h,f2h;
			core::Real tau(0.0);





			core::Size orbital_surrogate_atom_index = atom_index;





			numeric::deriv::angle_p1_deriv( pD, pDH, pO, tau, f1ab, f2ab  );
			numeric::deriv::angle_p2_deriv( pD, pDH, pO, tau, f1ao, f2ao  );
			numeric::deriv::angle_p1_deriv( pO, pDH, pD, tau, f1h, f2h  );


			Real neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;


			r2_atom_derivs[ H_index ].f1() += weight * neg_sine_tau * a_deriv * f1ab;
			r2_atom_derivs[ H_index ].f2() += weight * neg_sine_tau * a_deriv * f2ab;

			r2_atom_derivs[ donor_index ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
			r2_atom_derivs[ donor_index ].f2() += weight * neg_sine_tau * a_deriv * f2ao;

			r1_atom_derivs[ atom_index ].f1() += weight * neg_sine_tau * a_deriv * f1h;
			r1_atom_derivs[ atom_index ].f2() += weight * neg_sine_tau * a_deriv * f2h;


/*			if(res1.atom_type(atom_index).name() == "OH"){
				r1_atom_derivs[atom_index].f1() = weight * neg_sine_tau * a_deriv * f1ab;
				r1_atom_derivs[atom_index].f2() = weight * neg_sine_tau * a_deriv * f2ab;
				r1_atom_derivs[atom_index].f1() = weight * neg_sine_tau * a_deriv * f1h;
				r1_atom_derivs[atom_index].f2() = weight * neg_sine_tau * a_deriv * f2h;
			}*/


			f1ao = f2ao = 0;
			Real dis(0.0);
			numeric::deriv::distance_f1_f2_deriv( pDH, pO, dis, f1ao, f2ao );

			r2_atom_derivs[ donor_index ].f1() += weight * d_deriv * f1ao;
			r2_atom_derivs[ donor_index ].f2() += weight * d_deriv * f2ao;

			r1_atom_derivs[ atom_index ].f1() -= weight * d_deriv * f1ao;
			r1_atom_derivs[ atom_index ].f2() -= weight * d_deriv * f2ao;


			/////////
			//this starts the AOH derivative calculation

			core::Real cosAOH(cos_of(Axyz, Orbxyz, Hxyz));
			lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, OrbHdist, cosAOH, energy, d_deriv, a_deriv, true);



			Vector pAB(Axyz);
			Vector pAO(Orbxyz);
			Vector pH(Hxyz);

			f1ab=f2ab=f1ao=f2ao=f1h=f2h=0;
			tau=0;







			numeric::deriv::angle_p1_deriv( pAB, pAO, pH, tau, f1ab, f2ab  );
			numeric::deriv::angle_p2_deriv( pAB, pAO, pH, tau, f1ao, f2ao  );
			numeric::deriv::angle_p1_deriv( pH, pAO, pAB, tau, f1h, f2h  );

			neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;core::Size orbital_surrogate_atom_index = atom_index;


			r1_atom_derivs[ atom_index ].f1() += weight * neg_sine_tau * a_deriv * f1ab;
			r1_atom_derivs[ atom_index ].f2() += weight * neg_sine_tau * a_deriv * f2ab;

			r1_atom_derivs[ orbital_surrogate_atom_index ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
			r1_atom_derivs[ orbital_surrogate_atom_index ].f2() += weight * neg_sine_tau * a_deriv * f2ao;

			r2_atom_derivs[ H_index ].f1() += weight * neg_sine_tau * a_deriv * f1h;
			r2_atom_derivs[ H_index ].f2() += weight * neg_sine_tau * a_deriv * f2h;

/*			if(res1.atom_type(atom_index).name() == "OH"){
				r1_atom_derivs[res1.atom_base(atom_index)].f1() = weight * neg_sine_tau * a_deriv * f1ab;
				r1_atom_derivs[res1.atom_base(atom_index)].f2() = weight * neg_sine_tau * a_deriv * f2ab;
				r1_atom_derivs[res1.atom_base(atom_index)].f1() = weight * neg_sine_tau * a_deriv * f1h;
				r1_atom_derivs[res1.atom_base(atom_index)].f2() = weight * neg_sine_tau * a_deriv * f2h;

			}*/

			f1ao = f2ao = 0;
			dis= 0.0;
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
