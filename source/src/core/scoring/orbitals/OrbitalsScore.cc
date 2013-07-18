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
#include <core/conformation/Conformation.hh>

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
#include <core/id/AtomID.hh>
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
	sts.push_back( pci_cation_pi );
	sts.push_back( pci_pi_pi );
	sts.push_back(pci_hbond);
	sts.push_back(orbitals_hpol_bb);
	sts.push_back(pci_salt_bridge);
	return sts;
}

static basic::Tracer TR("core.scoring.orbitals_hpol");

//Because we don't really use the energy method options anyways

OrbitalsScore::OrbitalsScore() :
							parent( new OrbitalsScoreCreator ),
							lookup_table_(core::scoring::ScoringManager::get_instance()->get_OrbitalsLookupTable()),
							max_orbital_dist_squared_(9),
							max_dist_squared_(36)
{
	if(basic::options::option[ basic::options::OptionKeys::in::add_orbitals] != 1){
		utility_exit_with_message( "Trying to run features test without orbitals! Pass the flag -add_orbitals!" );
	}
}

OrbitalsScore::OrbitalsScore(methods::EnergyMethodOptions const &) :
							parent( new OrbitalsScoreCreator ),
							lookup_table_(core::scoring::ScoringManager::get_instance()->get_OrbitalsLookupTable()),
							max_orbital_dist_squared_(9),
							max_dist_squared_(36)

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

void OrbitalsScore::setup_for_scoring(pose::Pose & pose, ScoreFunction const & weights) const
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
OrbitalsScore::setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & minmap,
		ResSingleMinimizationData const & res1_data_cache,
		ResSingleMinimizationData const & res2_data_cache,
		ResPairMinimizationData & data_cache
) const{
	//std::cout << "we got to setup for minimizeing" << std::endl;
	conformation::Residue *res1_ptr = const_cast<conformation::Residue *>(&rsd1);
	res1_ptr->update_orbital_coords();
	conformation::Residue *res2_ptr = const_cast<conformation::Residue *>(&rsd2);
	res2_ptr->update_orbital_coords();
	//rsd1.update_orbital_coords();
	//rsd2.update_orbital_coords();
}

void
OrbitalsScore::finalize_after_derivatives( pose::Pose & pose, ScoreFunction const &  ) const{
	for(core::Size resid=1; resid <= pose.n_residue(); ++resid){
		pose.update_orbital_coords(resid);
	}
}


void
OrbitalsScore::finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &
) const{
	for(core::Size resid=1; resid <= pose.n_residue(); ++resid){
		pose.update_orbital_coords(resid);
	}
}


void
OrbitalsScore::setup_for_minimizing(
		pose::Pose & pose,
		ScoreFunction const & ,
		kinematics::MinimizerMapBase const &
) const{
	for(core::Size resid=1; resid <= pose.n_residue(); ++resid){
		pose.update_orbital_coords(resid);
	}
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
	return  7.0;
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



//orbital types found in
void OrbitalsScore::scfxn_rules_for_energy(
		bool hydrogen_interaction,
		core::Size orbtype1,
		OrbitalsLookup::h_type htype,
		core::Size orbtype2,
		core::Real energy,
		EnergyMap & emap
) const{
	if ( hydrogen_interaction ) {
		if ( orbtype1 == 1 ) {//c.pi.sp2 orbital type
			if ( htype == lookup_table_.Hpol_scOrbH ) emap[pci_cation_pi] += energy; //polar hydrogen
			if ( htype == lookup_table_.Haro_scOrbH ) emap[pci_pi_pi] += energy; //aromatic hydrogen
			if ( htype == lookup_table_.Hpol_bbOrbH ) emap[orbitals_hpol_bb] += energy; //bb hydrogen
		}
		if ( orbtype1 == 2 ) {//n.pi.sp2 orbtal type. Do nothing with hydrogen...no energy associated with this
		}
		if ( orbtype1 == 3 ) {//n.p.sp2 orbital type
			if ( htype == lookup_table_.Hpol_scOrbH ) emap[pci_hbond] += energy; //polar hydrogen
			//if ( htype == 2 ) emap[pci_pi_pi] += energy; //aromatic hydrogen
			//if ( htype == 3 ) emap[pci_pi_pi] += energy; //bb hydrogen
		}
    if ( orbtype1 == 4 ) {
      if ( htype == lookup_table_.Hpol_scOrbH ) emap[pci_salt_bridge] += energy;
    }
		if ( orbtype1 == 5 ) {
			if ( htype == lookup_table_.Hpol_scOrbH ) emap[pci_salt_bridge] += energy;
		}
		if ( orbtype1 == 6 ) {
			if ( htype == lookup_table_.Hpol_scOrbH ) emap[pci_hbond] += energy;
		}
		if ( orbtype1 == 9 ) {
			if ( htype == lookup_table_.Hpol_bbOrbH ) emap[orbitals_hpol_bb] += energy;
		}
	} else {
		if ( orbtype1 == static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			if(orbtype2 == static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2) ){
				emap[pci_pi_pi] += (energy);
			} //c.pi.sp2 to c.pi.sp2
		}
		if ( orbtype1 == 1 && orbtype2 == 2 ){ emap[pci_cation_pi] += energy;}//c.pi.sp2 to n.pi.sp2
		if ( orbtype1 == 2 && orbtype2 == 1 ){ emap[pci_cation_pi] += energy;}//c.pi.sp2 to n.pi.sp2
		if ( orbtype1 == 1 && orbtype2 == 3 ){ emap[pci_cation_pi] += energy;}//c.pi.sp2 to n.p.sp2
		if ( orbtype1 == 3 && orbtype2 == 1 ){ emap[pci_cation_pi] += energy;}//c.pi.sp2 to n.p.sp2

	}
}
core::Real OrbitalsScore::scfxn_rules_for_weight(
		bool hydrogen_interaction,
		core::Size orbtype1,
		OrbitalsLookup::h_type htype,
		core::Size orbtype2,
		EnergyMap const & emap
) const{
	if ( hydrogen_interaction ) {
		if ( orbtype1 == 1 ) {//c.pi.sp2 orbital type
			if ( htype == lookup_table_.Hpol_scOrbH ) return emap[pci_cation_pi]; //polar hydrogen
			if ( htype == lookup_table_.Haro_scOrbH ) return emap[pci_pi_pi] ; //aromatic hydrogen
			if ( htype == lookup_table_.Hpol_bbOrbH ) return emap[orbitals_hpol_bb] ; //bb hydrogen
		}
		if ( orbtype1 == 2 ) {//n.pi.sp2 orbtal type. Do nothing with hydrogen...no energy associated with this
		}
		if ( orbtype1 == 3 ) {//n.p.sp2 orbital type
			if ( htype == lookup_table_.Hpol_scOrbH ) return emap[pci_hbond] ; //polar hydrogen
			//if ( htype == 2 ) emap[pci_pi_pi] += energy; //aromatic hydrogen
			//if ( htype == 3 ) emap[pci_pi_pi] += energy; //bb hydrogen
		}
    if ( orbtype1 == 4 ) {
      if ( htype == lookup_table_.Hpol_scOrbH ) return emap[pci_salt_bridge] ;
    }

		if ( orbtype1 == 5 ) {
			if ( htype == lookup_table_.Hpol_scOrbH ) return emap[pci_salt_bridge] ;
		}
		if ( orbtype1 == 6 ) {
			if ( htype == lookup_table_.Hpol_scOrbH ) return emap[pci_hbond];
		}
		if ( orbtype1 == 9 ) {
			if ( htype == lookup_table_.Hpol_bbOrbH ) return emap[orbitals_hpol_bb];
		}
	} else {
		if ( orbtype1 == 1 && orbtype2 == 1 ) return emap[pci_pi_pi]; //c.pi.sp2 to c.pi.sp2
		if ( orbtype1 == 1 && orbtype2 == 2 ) return emap[pci_cation_pi]; //c.pi.sp2 to n.pi.sp2
		if ( orbtype1 == 2 && orbtype2 == 1 ) return emap[pci_cation_pi]; //c.pi.sp2 to n.pi.sp2
		if ( orbtype1 == 1 && orbtype2 == 3 ) return emap[pci_cation_pi]; //c.pi.sp2 to n.p.sp2
		if ( orbtype1 == 3 && orbtype2 == 1 ) return emap[pci_cation_pi]; //c.pi.sp2 to n.p.sp2
	}
	return 0.0;
}



bool OrbitalsScore::orb_orb_rules(
		const core::Size atype1,
		const core::Size atype2
)const
{

	if(atype1==6){//atype 6 is aroC, which would be for C.pi -> C.pi interactions
		if(atype2==11 || atype2==6)     {//atype 11 is Narg, which would be for N.pi -> C.pi interactions
			return true;
		}else {return false;}
	}
	if(atype1==11){
		if(atype2==11 || atype2==6)     {//atype 11 is Narg, which would be for N.pi -> C.pi interactions
			return true;
		}else {return false;}   }
	if(atype2==6){//atype 6 is aroC, which would be for C.pi -> C.pi interactions
		if(atype1==11 || atype1==6)     {//atype 11 is Narg, which would be for N.pi -> C.pi interactions
			return true;
		}else {return false;}
	}
	if(atype2==11){
		if(atype1==11 || atype1==6)     {//atype 11 is Narg, which would be for N.pi -> C.pi interactions
			return true;
		}else {return false;}   }
	return false;
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
	get_E_haro_one_way(res2, res1, emap);
	get_E_haro_one_way(res1, res2, emap);

	get_E_hpol_one_way(res1, res2, emap);
	get_E_hpol_one_way(res2, res1, emap);

	get_orb_orb_E(res1, res2, emap);
	get_orb_orb_E(res2, res1, emap);
}

//TODO This function duplicates code from the one directly below it.  The function below could
// Just call this function, provided that the extra function call, AtomID construction and
// if statement evaluation don't cause measurable slow down.  It would improve the maintainability of the code.
void OrbitalsScore::get_orb_orb_E(
		core::pose::Pose const & pose,
		core::id::AtomID const & atom1,
		core::id::AtomID const & atom2,
		EnergyMap & emap
)const
{
	core::Real orb_orb_E(0.0);
	core::conformation::Residue const & res1(pose.conformation().residue(atom1.rsd()));
	core::conformation::Residue const & res2(pose.conformation().residue(atom2.rsd()));
	if(
		res1.type().atom_is_backbone(atom1.atomno()) &&
		res2.type().atom_is_backbone(atom2.atomno()) &&
		res1.type().atom_type(atom1.atomno()).atom_has_orbital() &&
		res2.type().atom_type(atom2.atomno()).atom_has_orbital() &&
		orb_orb_rules(res1.atom_type_index(atom1.atomno()), res2.atom_type_index(atom2.atomno()))
	)
	{
		utility::vector1< core::Size > const & res1_orbs(res1.bonded_orbitals(atom1.atomno()));
		for(
				utility::vector1< core::Size >::const_iterator
				res1_orb = res1_orbs.begin(),
				res1_orb_end = res1_orbs.end();
				res1_orb != res1_orb_end; ++res1_orb
		){

			utility::vector1< core::Size > const & res2_orbs(res2.bonded_orbitals(atom2.atomno()));
			for(
					utility::vector1< core::Size >::const_iterator
					res2_orb = res2_orbs.begin(),
					res2_orb_end = res2_orbs.end();
					res2_orb != res2_orb_end; ++res2_orb
			){
				numeric::xyzVector< core::Real > const & res1_Orbxyz(res1.orbital_xyz(*res1_orb) );
				numeric::xyzVector< core::Real > const & res2_Orbxyz(res2.orbital_xyz(*res2_orb) );
				core::Real const orb1_orb2_dist= res1_Orbxyz.distance_squared(res2_Orbxyz);
				if(orb1_orb2_dist < 16){
					core::Size const & orbital_type1(res1.orbital_type_index(*res1_orb));
					core::Size const & orbital_type2(res2.orbital_type_index(*res2_orb));
					core::Real const dist(std::sqrt(orb1_orb2_dist));
					numeric::xyzVector< core::Real > const & Axyz(res1.xyz(atom1.atomno()));
					numeric::xyzVector< core::Real > const & Dxyz(res2.xyz(atom2.atomno()));
					core::Real const cosAOD(cos_of(Axyz, res1_Orbxyz, Dxyz));
					core::Real const cosDOA(cos_of(Dxyz, res2_Orbxyz, Axyz));
					core::Real d_deriv(0.0);
					core::Real a_deriv(0.0);
					lookup_table_.OrbOrbDist_cosAOD_energy(orbital_type1, orbital_type2, dist, cosAOD, orb_orb_E, d_deriv, a_deriv, false);
					scfxn_rules_for_energy(false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, orb_orb_E, emap ); //dummy value for htype given
					orb_orb_E=0.0;
					lookup_table_.OrbOrbDist_cosDOA_energy(orbital_type1, orbital_type2, dist, cosDOA, orb_orb_E, d_deriv, a_deriv, false);
					scfxn_rules_for_energy(false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, orb_orb_E,  emap ); //dummy value for htype given
					orb_orb_E=0.0;
				}
			}
		}
	}
}

void OrbitalsScore::get_orb_orb_E(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap & emap
)const{
	core::Real orb_orb_E(0.0);
	for (
			chemical::AtomIndices::const_iterator
			Aindex = res1.atoms_with_orb_index().begin(),
			Aindex_end = res1.atoms_with_orb_index().end();
			Aindex != Aindex_end; ++Aindex
	){
		if ( !res1.atom_is_backbone(*Aindex) ) {
			for (
					chemical::AtomIndices::const_iterator
					Dindex = res2.atoms_with_orb_index().begin(),
					Dindex_end = res2.atoms_with_orb_index().end();
					Dindex != Dindex_end; ++Dindex
			)
			{
				if(!res2.atom_is_backbone(*Dindex)){
					if(orb_orb_rules(res1.atom_type_index(*Aindex),res2.atom_type_index(*Dindex) )){
						utility::vector1< core::Size > const & res1_orbs(res1.bonded_orbitals(*Aindex));
						for(
								utility::vector1< core::Size >::const_iterator
								res1_orb = res1_orbs.begin(),
								res1_orb_end = res1_orbs.end();
								res1_orb != res1_orb_end; ++res1_orb
						){

							utility::vector1< core::Size > const & res2_orbs(res2.bonded_orbitals(*Dindex));
							for(
									utility::vector1< core::Size >::const_iterator
									res2_orb = res2_orbs.begin(),
									res2_orb_end = res2_orbs.end();
									res2_orb != res2_orb_end; ++res2_orb
							){
								numeric::xyzVector< core::Real > const & res1_Orbxyz(res1.orbital_xyz(*res1_orb) );
								numeric::xyzVector< core::Real > const & res2_Orbxyz(res2.orbital_xyz(*res2_orb) );
								core::Real const orb1_orb2_dist= res1_Orbxyz.distance_squared(res2_Orbxyz);
								if(orb1_orb2_dist < 16){
									core::Size const & orbital_type1(res1.orbital_type_index(*res1_orb));
									core::Size const & orbital_type2(res2.orbital_type_index(*res2_orb));
									core::Real const dist(std::sqrt(orb1_orb2_dist));
									numeric::xyzVector< core::Real > const & Axyz(res1.xyz(*Aindex));
									numeric::xyzVector< core::Real > const & Dxyz(res2.xyz(*Dindex));
									core::Real const cosAOD(cos_of(Axyz, res1_Orbxyz, Dxyz));
									core::Real const cosDOA(cos_of(Dxyz, res2_Orbxyz, Axyz));
									core::Real d_deriv(0.0);
									core::Real a_deriv(0.0);
									lookup_table_.OrbOrbDist_cosAOD_energy(orbital_type1, orbital_type2, dist, cosAOD, orb_orb_E, d_deriv, a_deriv, false);
									scfxn_rules_for_energy(false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, orb_orb_E, emap ); //dummy value for htype given
									orb_orb_E=0.0;
									lookup_table_.OrbOrbDist_cosDOA_energy(orbital_type1, orbital_type2, dist, cosDOA, orb_orb_E, d_deriv, a_deriv, false);
									scfxn_rules_for_energy(false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, orb_orb_E,  emap ); //dummy value for htype given
									orb_orb_E=0.0;
								}
							}
						}
					}
				}
			}
		}
	}
}


//TODO This function duplicates code from the one directly below it.  The function below could
// Just call this function, provided that the extra function call, AtomID construction and
// if statement evalution don't cause measurable slow down.  It would improve the maintainability of the code.
void OrbitalsScore::get_E_haro_one_way(
	core::pose::Pose const & pose,
	core::id::AtomID const & atom1,
	core::id::AtomID const & atom2,
	EnergyMap & emap
	) const
{
	core::conformation::Residue const & res1(pose.conformation().residue(atom1.rsd()));
	core::conformation::Residue const & res2(pose.conformation().residue(atom2.rsd()));
	
	core::Real dummy_E1(0.0);//needed for generalized function get_orb_H_distance_and_energy
	core::Real energy(0.0);

	if(
		res1.type().atom_type(atom1.atomno()).atom_has_orbital() && // Acceptor has orbital
		res2.type().atom_type(atom2.atomno()).name() == "Haro" &&  // Donor atom is Haro
		res1.atom_is_backbone(atom1.atomno()) //Acceptor is backbone
	)
	{
		numeric::xyzVector<core::Real> const & Axyz = res1.atom(atom1.atomno()).xyz();
		numeric::xyzVector<core::Real> const & Hxyz = res2.atom(atom2.atomno()).xyz();
		core::Real const temp_dist = Axyz.distance_squared(Hxyz);
		if (temp_dist < max_dist_squared_) {
			core::Size const donor_index(res2.bonded_neighbor(atom2.atomno())[1]);
			numeric::xyzVector<core::Real> const & Dxyz(res2.xyz(donor_index));
			get_orb_H_distance_and_energy(res1, atom1.atomno(), Axyz, Hxyz, Dxyz, energy, dummy_E1, lookup_table_.Haro_scOrbH, false, emap);
		}
	}
}

void OrbitalsScore::get_E_haro_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap & emap
) const
{

	core::Real dummy_E1(0.0);//needed for generalized function get_orb_H_distance_and_energy
	core::Real energy(0.0);
	for(
			chemical::AtomIndices::const_iterator
			atoms_with_orb_index = res1.atoms_with_orb_index().begin(),
			atoms_with_orb_index_end = res1.atoms_with_orb_index().end();
			atoms_with_orb_index != atoms_with_orb_index_end; ++atoms_with_orb_index
	)
	{
		core::Size const Aindex(*atoms_with_orb_index); //acceptor index
		if ( !res1.atom_is_backbone(Aindex) ) {
			numeric::xyzVector<core::Real> const & Axyz = res1.atom(*atoms_with_orb_index).xyz();//acceptor xyz
			for (
					chemical::AtomIndices::const_iterator
					haro_index = res2.Haro_index().begin(),
					haro_end = res2.Haro_index().end();
					haro_index != haro_end; ++haro_index
			)
			{
				numeric::xyzVector<core::Real> const & Hxyz = res2.atom(*haro_index).xyz(); //hydrogen xyz
				core::Real const temp_dist = Axyz.distance_squared(Hxyz);
				if ( temp_dist < max_dist_squared_ ) {
					core::Size const donor_index(res2.bonded_neighbor(*haro_index)[1]);
					numeric::xyzVector<core::Real> const & Dxyz(res2.xyz(donor_index)); //donor xyz
					get_orb_H_distance_and_energy(res1, Aindex, Axyz, Hxyz, Dxyz, energy, dummy_E1, lookup_table_.Haro_scOrbH, false, emap);
				}
			}
		}
	}
}

//TODO This function duplicates code from the one directly below it.  The function below could
// Just call this function, provided that the extra function call, AtomID construction and
// if statement evalution don't cause measurable slow down.  It would improve the maintainability of the code.
void OrbitalsScore::get_E_hpol_one_way(
		core::pose::Pose const & pose,
		core::id::AtomID const & atom1,
		core::id::AtomID const & atom2,
		EnergyMap & emap
) const
{
	core::Real HPOL_sc_H_sc_orb_E(0.0);
	core::Real HPOL_bb_H_sc_orb_energy(0.0);

	core::conformation::Residue const & res1(pose.conformation().residue(atom1.rsd()));
	core::conformation::Residue const & res2(pose.conformation().residue(atom2.rsd()));

	core::Size const donor_index(res2.bonded_neighbor(atom2.atomno())[1]);
	if(
		res1.type().atom_type(atom1.atomno()).atom_has_orbital() && // Acceptor has orbital
		res2.type().atom_type(atom2.atomno()).name() == "Hpol" &&  // Donor atom is Haro
		!(
			res1.atom_is_backbone(atom1.atomno()) &&
			res2.atom_is_backbone(donor_index)
		)
	)
	{
		numeric::xyzVector<core::Real> const & Axyz(res1.atom(atom1.atomno()).xyz());
		numeric::xyzVector<core::Real> const & Hxyz(res2.atom(atom2.atomno()).xyz());

		if( Axyz.distance_squared(Hxyz) < max_dist_squared_)
		{
			numeric::xyzVector<core::Real> const & Dxyz(res2.xyz(donor_index));
			if(res2.atom_is_backbone(donor_index) || res1.atom_is_backbone(atom1.atomno())){
				get_orb_H_distance_and_energy(res1, atom1.atomno(), Axyz, Hxyz, Dxyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, lookup_table_.Hpol_bbOrbH, true, emap);
			}else{
				get_orb_H_distance_and_energy(res1, atom1.atomno(), Axyz, Hxyz, Dxyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, lookup_table_.Hpol_scOrbH, false, emap);
			}
		}

	}
}


void OrbitalsScore::get_E_hpol_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap & emap
) const
{
	core::Real HPOL_sc_H_sc_orb_E(0.0);
	core::Real HPOL_bb_H_sc_orb_energy(0.0);
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
			core::Size const donor_index(res2.bonded_neighbor(*hpol_index)[1]);
			if(res1.atom_is_backbone(*atoms_with_orb_index) && res2.atom_is_backbone(donor_index)){
				continue;
			}
			numeric::xyzVector<core::Real> const & Hxyz = res2.atom(*hpol_index).xyz(); //hydrogen xyz
			if ( Axyz.distance_squared(Hxyz) < max_dist_squared_ ) {
				core::Size const Aindex(*atoms_with_orb_index);
				numeric::xyzVector<core::Real> const & Dxyz(res2.xyz(donor_index)); //donor xyz
				if(res2.atom_is_backbone(donor_index) || res1.atom_is_backbone(Aindex)){
					get_orb_H_distance_and_energy(res1, Aindex, Axyz, Hxyz, Dxyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, lookup_table_.Hpol_bbOrbH, true, emap);
				}else{
					get_orb_H_distance_and_energy(res1, Aindex, Axyz, Hxyz, Dxyz, HPOL_sc_H_sc_orb_E, HPOL_bb_H_sc_orb_energy, lookup_table_.Hpol_scOrbH, false, emap);
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
		bool bb_h_flag,
		EnergyMap & emap
) const
{
	core::Real d_deriv(0.0);
	core::Real a_deriv(0.0);
	utility::vector1< core::Size > const & orbital_indices(res1.bonded_orbitals(Aindex));
	for(
			utility::vector1< core::Size >::const_iterator
			orbital_index = orbital_indices.begin(),
			orbital_index_end = orbital_indices.end();
			orbital_index != orbital_index_end; ++orbital_index
	)
	{
		numeric::xyzVector< core::Real > const & Orbxyz(res1.orbital_xyz(*orbital_index) );
		core::Real const temp_dist_squared = Orbxyz.distance_squared( Hxyz );
		if(temp_dist_squared < max_orbital_dist_squared_){
			core::Size orbital_type= res1.orbital_type_index(*orbital_index);
			core::Real const cosDHO(cos_of(Dxyz, Hxyz, Orbxyz));//Donor - Hydrogen - Orbital angle
			core::Real const cosAOH(cos_of(Axyz, Orbxyz, Hxyz));
			core::Real const dist(std::sqrt(temp_dist_squared));
			//change the orbital type to regular orbs for e calculation
			if(orbital_type==static_cast<core::Size>(core::chemical::orbitals::O_pi_sp2_bb)){
				orbital_type=static_cast<core::Size>(core::chemical::orbitals::O_pi_sp2);
			}
			if(orbital_type==static_cast<core::Size>(core::chemical::orbitals::O_p_sp2_bb)){
				orbital_type=static_cast<core::Size>(core::chemical::orbitals::O_p_sp2);
			}
			if(bb_h_flag){
				lookup_table_.OrbHdist_cosDHO_energy(htype, orbital_type, dist, cosDHO, bb_h_energy, d_deriv, a_deriv, false);
				emap[orbitals_hpol_bb] += bb_h_energy;
				bb_h_energy=0;
				lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, dist, cosAOH, bb_h_energy, d_deriv, a_deriv, false, false);
				emap[orbitals_hpol_bb] += bb_h_energy;
				bb_h_energy=0;
			}else{
				lookup_table_.OrbHdist_cosDHO_energy(htype, orbital_type, dist, cosDHO, sc_energy, d_deriv, a_deriv, false);
				scfxn_rules_for_energy(true, orbital_type, htype, 0, sc_energy, emap);
				sc_energy=0;
				//a little confusing without context. This checks to see if the residue is an aromatic residue. If it is an aromatic residue
				//then we need to check if the orbital we are looking at is the action center orbital. If it is, then we use a separate
				//energy than if it were. This is done because the action center orbital has a different angle associated with the acceptor
				//orbital hydrogen angle. why? because there is no index for action centers, thefore the Acceptor is the first action center atom
				if(res1.aa() == chemical::aa_tyr || res1.aa() == chemical::aa_phe || res1.aa() == chemical::aa_tyr){
					if(*orbital_index<=2){
						lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, dist, cosAOH, sc_energy, d_deriv, a_deriv, false, true);
						scfxn_rules_for_energy(true, orbital_type, htype, 0, sc_energy, emap);
						sc_energy =0;
					}else{
						lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, dist, cosAOH, sc_energy, d_deriv, a_deriv, false, false);
						scfxn_rules_for_energy(true, orbital_type, htype, 0, sc_energy, emap);
						sc_energy =0;
					}
				}else{
					lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, dist, cosAOH, sc_energy, d_deriv, a_deriv, false, false);
					scfxn_rules_for_energy(true, orbital_type, htype, 0, sc_energy, emap);
					sc_energy =0;
				}
			}
		}
	}

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

	assign_haro_derivs_one_way(res2, res1, weights, r2_atom_derivs, r1_atom_derivs);
	assign_haro_derivs_one_way(res1, res2, weights, r1_atom_derivs, r2_atom_derivs);
	assign_hpol_derivs_one_way(res1, res2, weights, r1_atom_derivs, r2_atom_derivs);
	assign_hpol_derivs_one_way(res2, res1, weights, r2_atom_derivs, r1_atom_derivs);
	assign_orb_orb_derivs(res1, res2, weights, r1_atom_derivs, r2_atom_derivs);
	assign_orb_orb_derivs(res2, res1, weights, r2_atom_derivs, r1_atom_derivs);

}




void OrbitalsScore::assign_haro_derivs_one_way(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
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
				if ( Axyz.distance_squared(Hxyz) < max_dist_squared_ ) {
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
			core::Size const donor_index(res2.bonded_neighbor(*hpol_index)[1]);
			if(res1.atom_is_backbone(*atoms_with_orb_index) && res2.atom_is_backbone(donor_index)){
				continue;
			}

			numeric::xyzVector<core::Real> const & Axyz = res1.atom(*atoms_with_orb_index).xyz();
			numeric::xyzVector<core::Real> const & Hxyz = res2.atom(*hpol_index).xyz();
			if ( Axyz.distance_squared(Hxyz) < max_dist_squared_ ) {
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


void
OrbitalsScore::assign_orb_orb_derivs(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
)const {
	core::Real d_deriv(0.0);
	core::Real a_deriv(0.0);
	core::Real orb_orb_E(0.0);
	for (
			chemical::AtomIndices::const_iterator
			Aindex = res1.atoms_with_orb_index().begin(),
			Aindex_end = res1.atoms_with_orb_index().end();
			Aindex != Aindex_end; ++Aindex
	){
		if ( !res1.atom_is_backbone(*Aindex) ) {
			for (
					chemical::AtomIndices::const_iterator
					Dindex = res2.atoms_with_orb_index().begin(),
					Dindex_end = res2.atoms_with_orb_index().end();
					Dindex != Dindex_end; ++Dindex
			)
			{
				if(!res2.atom_is_backbone(*Dindex)){

					utility::vector1< core::Size > const & res1_orbs(res1.bonded_orbitals(*Aindex));
					for(
							utility::vector1< core::Size >::const_iterator
							res1_orb = res1_orbs.begin(),
							res1_orb_end = res1_orbs.end();
							res1_orb != res1_orb_end; ++res1_orb
					){
						numeric::xyzVector< core::Real > const res1_Orbxyz(res1.orbital_xyz(*res1_orb) );
						utility::vector1< core::Size > const & res2_orbs(res2.bonded_orbitals(*Dindex));
						for(
								utility::vector1< core::Size >::const_iterator
								res2_orb = res2_orbs.begin(),
								res2_orb_end = res2_orbs.end();
								res2_orb != res2_orb_end; ++res2_orb
						){

							if(res1.atom_type_index(*Aindex),res2.atom_type_index(*Dindex)  ){
								numeric::xyzVector< core::Real > const res2_Orbxyz(res2.orbital_xyz(*res2_orb) );
								core::Real const orb1_orb2_dist= res1_Orbxyz.distance_squared(res2_Orbxyz);
								if(orb1_orb2_dist < 16){
									core::Size orbital_type1 = res1.orbital_type_index(*res1_orb);
									core::Size orbital_type2 = res2.orbital_type_index(*res2_orb);
									core::Real dist(std::sqrt(orb1_orb2_dist));
									numeric::xyzVector< core::Real > const Axyz(res1.xyz(*Aindex));
									numeric::xyzVector< core::Real > const Dxyz(res2.xyz(*Dindex));
									core::Real cosAOD(cos_of(Axyz, res1_Orbxyz, Dxyz));
									core::Real cosDOA(cos_of(Dxyz, res2_Orbxyz, Axyz));
									lookup_table_.OrbOrbDist_cosDOA_energy(orbital_type1, orbital_type2, dist, cosDOA, orb_orb_E, d_deriv, a_deriv, true);
									core::Real weight = scfxn_rules_for_weight(false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, weights); //dummy for htype



									Vector pD(Dxyz);
									Vector pDH(res2_Orbxyz);
									Vector pO(Axyz);




									Vector f1ab,f2ab;
									Vector f1ao,f2ao;
									Vector f1h,f2h;
									core::Real tau(0.0);

									numeric::deriv::angle_p1_deriv( pD, pDH, pO, tau, f1ab, f2ab  );
									numeric::deriv::angle_p2_deriv( pD, pDH, pO, tau, f1ao, f2ao  );
									numeric::deriv::angle_p1_deriv( pO, pDH, pD, tau, f1h, f2h  );


									Real neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;


									r2_atom_derivs[ *Dindex ].f1() += weight * neg_sine_tau * a_deriv * f1ab;
									r2_atom_derivs[ *Dindex ].f2() += weight * neg_sine_tau * a_deriv * f2ab;

									r2_atom_derivs[ *Dindex ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
									r2_atom_derivs[ *Dindex ].f2() += weight * neg_sine_tau * a_deriv * f2ao;

									r1_atom_derivs[ *Aindex ].f1() += weight * neg_sine_tau * a_deriv * f1h;
									r1_atom_derivs[ *Aindex ].f2() += weight * neg_sine_tau * a_deriv * f2h;


									f1ao = f2ao = 0;
									Real dis(0.0);
									numeric::deriv::distance_f1_f2_deriv( pDH, res1_Orbxyz, dis, f1ao, f2ao );

									r2_atom_derivs[ *Dindex ].f1() += weight * d_deriv * f1ao;
									r2_atom_derivs[ *Dindex ].f2() += weight * d_deriv * f2ao;

									r1_atom_derivs[ *Aindex ].f1() -= weight * d_deriv * f1ao;
									r1_atom_derivs[ *Aindex ].f2() -= weight * d_deriv * f2ao;


									/////////
									//this starts the AOH derivative calculation

									d_deriv=0;
									a_deriv=0;
									lookup_table_.OrbOrbDist_cosAOD_energy(orbital_type1, orbital_type2, dist, cosAOD, orb_orb_E, d_deriv, a_deriv, true);



									Vector pAB(Axyz);
									Vector pAO(res1_Orbxyz);
									Vector pH(Dxyz);

									f1ab=f2ab=f1ao=f2ao=f1h=f2h=0;
									tau=0;



									numeric::deriv::angle_p1_deriv( pAB, pAO, pH, tau, f1ab, f2ab  );
									numeric::deriv::angle_p2_deriv( pAB, pAO, pH, tau, f1ao, f2ao  );
									numeric::deriv::angle_p1_deriv( pH, pAO, pAB, tau, f1h, f2h  );

									neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;core::Size orbital_surrogate_atom_index = atom_index;


									r1_atom_derivs[ *Aindex ].f1() += weight * neg_sine_tau * a_deriv * f1ab;
									r1_atom_derivs[ *Aindex ].f2() += weight * neg_sine_tau * a_deriv * f2ab;

									r1_atom_derivs[ *Aindex ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
									r1_atom_derivs[ *Aindex ].f2() += weight * neg_sine_tau * a_deriv * f2ao;

									r2_atom_derivs[ *Dindex ].f1() += weight * neg_sine_tau * a_deriv * f1h;
									r2_atom_derivs[ *Dindex ].f2() += weight * neg_sine_tau * a_deriv * f2h;

									f1ao = f2ao = 0;
									dis= 0.0;
									numeric::deriv::distance_f1_f2_deriv( pAO, res2_Orbxyz, dis, f1ao, f2ao );

									r1_atom_derivs[ *Aindex ].f1() += weight * d_deriv * f1ao;
									r1_atom_derivs[ *Aindex ].f2() += weight * d_deriv * f2ao;

									r2_atom_derivs[ *Dindex ].f1() -= weight * d_deriv * f1ao;
									r2_atom_derivs[ *Dindex ].f2() -= weight * d_deriv * f2ao;


								}
							}

						}
					}


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

	utility::vector1< core::Size > const & orbital_indices(res1.bonded_orbitals(atom_index));
	for(
			utility::vector1< core::Size >::const_iterator
			orbital_index = orbital_indices.begin(),
			orbital_index_end = orbital_indices.end();
			orbital_index != orbital_index_end; ++orbital_index
	)
	{
		numeric::xyzVector< core::Real > const Orbxyz(res1.orbital_xyz(*orbital_index) );
		core::Real const temp_dist_squared = Orbxyz.distance_squared( Hxyz );
		if(temp_dist_squared < max_orbital_dist_squared_){
			core::Size orbital_type = res1.orbital_type_index(*orbital_index);
			if(orbital_type==static_cast<core::Size>(core::chemical::orbitals::O_pi_sp2_bb)){
				orbital_type=static_cast<core::Size>(core::chemical::orbitals::O_pi_sp2);
			}
			if(orbital_type==static_cast<core::Size>(core::chemical::orbitals::O_p_sp2_bb)){
				orbital_type=static_cast<core::Size>(core::chemical::orbitals::O_p_sp2);
			}

			core::Real weight =  scfxn_rules_for_weight(true, orbital_type, htype, 0, weights); //check weights before we reassign orbitals type


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


			f1ao = f2ao = 0;
			Real dis(0.0);
			numeric::deriv::distance_f1_f2_deriv( pDH, pO, dis, f1ao, f2ao );

			r2_atom_derivs[ donor_index ].f1() += weight * d_deriv * f1ao;
			r2_atom_derivs[ donor_index ].f2() += weight * d_deriv * f2ao;

			r1_atom_derivs[ atom_index ].f1() -= weight * d_deriv * f1ao;
			r1_atom_derivs[ atom_index ].f2() -= weight * d_deriv * f2ao;


			/////////
			//this starts the AOH derivative calculation


			d_deriv =0;
			a_deriv =0;

			core::Real cosAOH(cos_of(Axyz, Orbxyz, Hxyz));
			//a little confusing without context. This checks to see if the residue is an aromatic residue. If it is an aromatic residue
			//then we need to check if the orbital we are looking at is the action center orbital. If it is, then we use a separate
			//energy than if it were. This is done because the action center orbital has a different angle associated with the acceptor
			//orbital hydrogen angle. why? because there is no index for action centers, thefore the Acceptor is the first action center atom

			if(res1.aa() == chemical::aa_tyr || res1.aa() == chemical::aa_phe || res1.aa() == chemical::aa_tyr){
				if(*orbital_index<=2){
					lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, OrbHdist, cosAOH, energy, d_deriv, a_deriv, true, true);
				}else{
					lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, OrbHdist, cosAOH, energy, d_deriv, a_deriv, true, false);
				}
			}else{
				lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, OrbHdist, cosAOH, energy, d_deriv, a_deriv, true, false);
			}





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
