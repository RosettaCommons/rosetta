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
	sts.push_back( pci_hbond );
	sts.push_back( orbitals_hpol_bb );
	sts.push_back( pci_salt_bridge );
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

OrbitalsScore::OrbitalsScore( methods::EnergyMethodOptions const & ) :
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

void OrbitalsScore::setup_for_scoring(pose::Pose & pose, ScoreFunction const & ) const
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
		pose::Pose const &,
		ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData &
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

ScoreType
score_type_for_orb_params(
	OrbitalsLookup const & lookup_table,
	bool hydrogen_interaction,
	bool backbone, // if the is a hydrogen interaction involving any backbone atoms, either the hydrogen or the acceptor orbital, then use the orbitals_hpol_bb ScoreType
	core::Size orbtype1,
	OrbitalsLookup::h_type htype,
	core::Size orbtype2
)
{
	if ( hydrogen_interaction ) {
		if ( backbone ) {
			return orbitals_hpol_bb;
		}

		if ( orbtype1 == 1 ) {//c.pi.sp2 orbital type
			if ( htype == lookup_table.Hpol_scOrbH ) return pci_cation_pi;  //polar hydrogen
			if ( htype == lookup_table.Haro_scOrbH ) return pci_pi_pi; //aromatic hydrogen
			if ( htype == lookup_table.Hpol_bbOrbH ) return orbitals_hpol_bb; //bb hydrogen
		}
		if ( orbtype1 == 2 ) {//n.pi.sp2 orbtal type. Do nothing with hydrogen...no energy associated with this
		}
		if ( orbtype1 == 3 ) {//n.p.sp2 orbital type
			if ( htype == lookup_table.Hpol_scOrbH ) return pci_hbond; //polar hydrogen
			//if ( htype == 2 ) return pci_pi_pi; //aromatic hydrogen
			//if ( htype == 3 ) return pci_pi_pi; //bb hydrogen
		}
    if ( orbtype1 == 4 ) {
      if ( htype == lookup_table.Hpol_scOrbH ) return pci_salt_bridge;
    }
		if ( orbtype1 == 5 ) {
			if ( htype == lookup_table.Hpol_scOrbH ) return pci_salt_bridge;
		}
		if ( orbtype1 == 6 ) {
			if ( htype == lookup_table.Hpol_scOrbH ) return pci_hbond;
		}
		if ( orbtype1 == 9 ) {
			if ( htype == lookup_table.Hpol_bbOrbH ) return orbitals_hpol_bb;
		}
	} else {
		if ( orbtype1 == static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2)){
			if(orbtype2 == static_cast <core::Size>(core::chemical::orbitals::C_pi_sp2) ){
				return pci_pi_pi;
			} //c.pi.sp2 to c.pi.sp2
		}
		if ( orbtype1 == 1 && orbtype2 == 2 ){ return pci_cation_pi;}//c.pi.sp2 to n.pi.sp2
		if ( orbtype1 == 2 && orbtype2 == 1 ){ return pci_cation_pi;}//c.pi.sp2 to n.pi.sp2
		if ( orbtype1 == 1 && orbtype2 == 3 ){ return pci_cation_pi;}//c.pi.sp2 to n.p.sp2
		if ( orbtype1 == 3 && orbtype2 == 1 ){ return pci_cation_pi;}//c.pi.sp2 to n.p.sp2

	}
	return fa_atr;
}


//orbital types found in
void OrbitalsScore::scfxn_rules_for_energy(
	bool hydrogen_interaction,
	bool backbone,
	core::Size orbtype1,
	OrbitalsLookup::h_type htype,
	core::Size orbtype2,
	core::Real energy,
	EnergyMap & emap
) const
{
	ScoreType which_st = score_type_for_orb_params( lookup_table_, hydrogen_interaction, backbone, orbtype1, htype, orbtype2 );
	if ( which_st != fa_atr ) {
		//		if ( energy != 0 ) {
		//	std::cout << "scfxn_rules_for_energy: " << hydrogen_interaction << " " << backbone << " " << orbtype1 << " " << orbtype2 << " " << energy << " " << which_st << std::endl;
		//}
		emap[ which_st ] += energy;
	}
}

core::Real OrbitalsScore::scfxn_rules_for_weight(
	bool hydrogen_interaction,
	bool backbone,
	core::Size orbtype1,
	OrbitalsLookup::h_type htype,
	core::Size orbtype2,
	EnergyMap const & emap
) const{
	ScoreType which_st = score_type_for_orb_params( lookup_table_, hydrogen_interaction, backbone, orbtype1, htype, orbtype2 );
	if ( which_st != fa_atr ) {
		//std::cout << "scfxn_rules_for_weight: " << hydrogen_interaction << " " << backbone << " " << orbtype1 << " " << orbtype2 << " " << which_st << std::endl;
		return emap[ which_st ];
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

bool
OrbitalsScore::requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & ) const { return true; }

/// @details OH, WOE IS ME!  Const casting an input Residue object? This is an awful thing for a score term to do.
void
OrbitalsScore::setup_for_scoring_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const &,
	ScoreFunction const &,
	ResSingleMinimizationData &
) const
{
	( const_cast< conformation::Residue * > (&rsd) )->update_orbital_coords();
}

bool
OrbitalsScore::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const & ) const { return true; }

/// @details OH, WOE IS ME!  Const casting an input Residue object? This is an awful thing for a score term to do.
void
OrbitalsScore::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const &,
	ScoreFunction const &,
	ResSingleMinimizationData &
) const
{
	( const_cast< conformation::Residue * > (&rsd) )->update_orbital_coords();
}


void OrbitalsScore::compute_orb_orb_E(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	core::Size res1_orb,
	core::Size res2_orb,
	core::Size res1_atomno,
	core::Size res2_atomno,
	EnergyMap & emap
) const
{
	core::Real orb_orb_E(0.0);
	numeric::xyzVector< core::Real > const & res1_Orbxyz(res1.orbital_xyz(res1_orb) );
	numeric::xyzVector< core::Real > const & res2_Orbxyz(res2.orbital_xyz(res2_orb) );
	core::Real const orb1_orb2_dist= res1_Orbxyz.distance_squared(res2_Orbxyz);
	if(orb1_orb2_dist < 9){
		core::Size const & orbital_type1(res1.orbital_type_index(res1_orb));
		core::Size const & orbital_type2(res2.orbital_type_index(res2_orb));
		core::Real const dist(std::sqrt(orb1_orb2_dist));
		numeric::xyzVector< core::Real > const & Axyz(res1.xyz(res1_atomno));
		numeric::xyzVector< core::Real > const & Dxyz(res2.xyz(res2_atomno));
		core::Real const cosAOD(cos_of(Axyz, res1_Orbxyz, Dxyz));
		core::Real const cosDOA(cos_of(Dxyz, res2_Orbxyz, Axyz));
		core::Real d_deriv(0.0);
		core::Real a_deriv(0.0);
		lookup_table_.OrbOrbDist_cosAOD_energy(orbital_type1, orbital_type2, dist, cosAOD, orb_orb_E, d_deriv, a_deriv, false);
		scfxn_rules_for_energy(false, false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, orb_orb_E, emap ); //dummy value for htype given
		orb_orb_E=0.0;
		lookup_table_.OrbOrbDist_cosDOA_energy(orbital_type1, orbital_type2, dist, cosDOA, orb_orb_E, d_deriv, a_deriv, false);
		scfxn_rules_for_energy(false, false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, orb_orb_E,  emap ); //dummy value for htype given
		orb_orb_E=0.0;
	}
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
		!res1.type().atom_is_backbone(atom1.atomno()) &&
		!res2.type().atom_is_backbone(atom2.atomno()) &&
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

				compute_orb_orb_E(res1,res2,*res1_orb,*res2_orb,atom1.atomno(),atom2.atomno(),emap);
			}
		}
	}
}


/*
void OrbitalsScore::cycle_through_orb_orb_interactions(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		EnergyMap & emap
){

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
									calculate_orb_orb_info_for_E(res1, res2, *res1_orb, *res2_orb, *Aindex, *Dindex, orb1_orb2_dist, emap);

								}
							}
						}
					}
				}
			}
		}
	}
}


void OrbitalsScore::calculate_orb_orb_info_for_E(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Size const & res1_orb,
		core::Size const & res2_orb,
		core::Size const & Aindex,
		core::Size const & Dindex,
		core::Real const & dist_squared,
		EnergyMap & emap
){
	numeric::xyzVector< core::Real > const & res1_Orbxyz(res1.orbital_xyz(res1_orb) );
	numeric::xyzVector< core::Real > const & res2_Orbxyz(res2.orbital_xyz(res2_orb) );
	core::Size const & orbital_type1(res1.orbital_type_index(res1_orb));
	core::Size const & orbital_type2(res2.orbital_type_index(res2_orb));
	core::Real const dist(std::sqrt(dist_squared));
	numeric::xyzVector< core::Real > const & Axyz(res1.xyz(Aindex));
	numeric::xyzVector< core::Real > const & Dxyz(res2.xyz(Dindex));
	core::Real const cosAOD(cos_of(Axyz, res1_Orbxyz, Dxyz));
	core::Real const cosDOA(cos_of(Dxyz, res2_Orbxyz, Axyz));
	core::Real d_deriv(0.0);
	core::Real a_deriv(0.0);
	lookup_table_.OrbOrbDist_cosAOD_energy(orbital_type1, orbital_type2, dist, cosAOD, orb_orb_E, d_deriv, a_deriv, false);
	scfxn_rules_for_energy(false, false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, orb_orb_E, emap ); //dummy value for htype given
	orb_orb_E=0.0;
	lookup_table_.OrbOrbDist_cosDOA_energy(orbital_type1, orbital_type2, dist, cosDOA, orb_orb_E, d_deriv, a_deriv, false);
	scfxn_rules_for_energy(false, false, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, orb_orb_E,  emap ); //dummy value for htype given
	orb_orb_E=0.0;

}
*/


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
								compute_orb_orb_E(res1,res2,*res1_orb,*res2_orb,*Aindex,*Dindex,emap);
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

core::chemical::orbitals::orbital_type_enum
lookup_type_for_orbital_type( core::chemical::orbitals::orbital_type_enum orbital_type )
{
	switch ( orbital_type ) {
		case core::chemical::orbitals::O_pi_sp2_bb : return core::chemical::orbitals::O_pi_sp2;
		case core::chemical::orbitals::O_p_sp2_bb : return core::chemical::orbitals::O_p_sp2;
		default : return orbital_type;
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
		if ( temp_dist_squared < max_orbital_dist_squared_ ) {
			core::Size orbital_type = res1.orbital_type_index(*orbital_index);
			core::Size orbital_lookup_type = lookup_type_for_orbital_type( chemical::orbitals::orbital_type_enum( orbital_type ) );

			core::Real const cosDHO(cos_of(Dxyz, Hxyz, Orbxyz)); // cos of Donor - Hydrogen - Orbital angle
			core::Real const cosAOH(cos_of(Axyz, Orbxyz, Hxyz)); // cos of Acceptor - Orbital - Hydrogen angle
			core::Real const dist(std::sqrt(temp_dist_squared));

			if ( bb_h_flag ) {
				lookup_table_.OrbHdist_cosDHO_energy(htype, orbital_lookup_type, dist, cosDHO, bb_h_energy, d_deriv, a_deriv, false);
				scfxn_rules_for_energy(true, bb_h_flag, orbital_type, htype, 0, bb_h_energy, emap);
				bb_h_energy=0;
				lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_lookup_type, dist, cosAOH, bb_h_energy, d_deriv, a_deriv, false, false);
				scfxn_rules_for_energy(true, bb_h_flag, orbital_type, htype, 0, bb_h_energy, emap);
				bb_h_energy=0;
			} else {
				lookup_table_.OrbHdist_cosDHO_energy(htype, orbital_type, dist, cosDHO, sc_energy, d_deriv, a_deriv, false);
				scfxn_rules_for_energy(true, bb_h_flag, orbital_type, htype, 0, sc_energy, emap);
				sc_energy=0;
				//a little confusing without context. This checks to see if the residue is an aromatic residue. If it is an aromatic residue
				//then we need to check if the orbital we are looking at is the action center orbital. If it is, then we use a separate
				//energy than if it were not. This is done because the action center orbital has a different angle associated with the acceptor
				//orbital hydrogen angle. why? because there is no index for action centers, thefore the Acceptor is the first action center atom
				bool const working_with_action_center_orbital = false;(
					(( res1.aa() == chemical::aa_tyr || res1.aa() == chemical::aa_phe || res1.aa() == chemical::aa_trp ) && *orbital_index <= 2 ));
				lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_type, dist, cosAOH, sc_energy, d_deriv, a_deriv, false, working_with_action_center_orbital );
				scfxn_rules_for_energy(true, bb_h_flag, orbital_type, htype, 0, sc_energy, emap);
				sc_energy=0;
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


core::Size
surrogate_atom_for_orbital(
	core::conformation::Residue const & res,
	core::Size acceptor_index,
	core::chemical::orbitals::orbital_type_enum orbital_type
)
{
	using core::chemical::orbitals::O_p_sp3;
	using core::chemical::orbitals::S_p_sp3;

	if ( orbital_type == O_p_sp3 || orbital_type == S_p_sp3 ) {
		// This makes an atom-tree assumption that the tree always goes out
		// from the acceptor atom towards atoms controlled by a higher "chi"
		// (or, if it's on the backbone, towards a higher-indexed backbone
		// dihedral, i.e. assuming an N-to-C fold tree). This assumption will
		// be violated frequently!  There doesn't seem to be a good way to fix
		// this issue except by 1) ignoring it, because it really is not a
		// problem for proteins, unless the atom tree is set so that a jump
		// lands on the HG and travels down towards the backbone along
		// OG --> CB --> CA, which I don't believe is practiced by anyone, or
		// 2) creating explicit orbitals psuedo atoms and assigning the f1/f2
		// vectors to those atoms, letting the atom tree sort out which DOFs
		// control the orbitals locations and not requiring the OrbitalsScore
		// term to either make assumptions about the kinematics of the system.
		if ( res.type().atom_is_backbone( acceptor_index ) ) {
			// uh oh!
		} else {
			// cases:
			// 1) there are no bonded atoms to this atom besides a single heavy atom
			// 2) there are no bonded atoms to this atom besides one or more hydrogen atoms (SER, THR, TYR, H2O)
			// 3) there are two or more bonded heavy atoms
			// 4) no bound atoms -- don't know what to do here!
			Size const n_bound_hydrogens = res.type().number_bonded_hydrogens( acceptor_index );
			Size const n_bound_hvyatoms = res.type().number_bonded_heavyatoms( acceptor_index );
			Size const n_bound_vrtatoms = res.type().nbrs( acceptor_index ).size() - n_bound_hydrogens - n_bound_hvyatoms;
			if ( n_bound_hydrogens == 0 ) {
				if ( n_bound_hvyatoms == 0 ) {
					if ( n_bound_vrtatoms == 0 ) {
						// case 4.  This is a weird case!  Give up and return the acceptor_index
						return acceptor_index;
					} else {
						// no bonded hydrogens or heavy atoms, but yet still bound to a virtual atom (there really ought to be at least two )
						// so let's assign the orbital derivatives to the first virtual atom.
						return res.type().bonded_neighbor( acceptor_index )[1];
					}
				} else if ( n_bound_hvyatoms == 1 ) {
					if ( n_bound_vrtatoms == 0 ) {
						// case 3a: in this case, we're talking about maybe a charged oxygen atom.  In the absence of
						// explicit orbital / virtual atom representation in the atom tree, there's no way to allow conformational
						// flexibility for the orbitals -- they can't spin! -- and so their positions must be calculated
						// strictly from the coordinate frame at the oxygen.  So in this case, we'll assign the
						// f1/f2 derivatives to the acceptor atom.
						return acceptor_index;
					} else {
						// case 3b: again, maybe a charged oxygen atom but now we have virtual atom(s)
						// giving us a coordinate frame with which we can reasonably allow the orbitals to spin.
						// Assign the orbital derivatives to the first virtual atom, since this atom is almost surely
						// being used to calculate the orbital locations
						return res.type().bonded_neighbor( acceptor_index )[1];
					}
				}
			} else {
				// there are some bound hydrogens
				if ( n_bound_hydrogens == 1 ) {
					if ( n_bound_hvyatoms == 0 ) {
						if ( n_bound_vrtatoms == 0 ) {
							// case 5a: OH- hydroxide floating around in space? There aren't enough atoms to define a coordinate
							// frame on the oxygen, so just assign the derivatives to the acceptor atom
							return acceptor_index;
						} else {
							// case 5b: OH- Hydroxide floating around in space with enough atoms to define a coordinate frame;
							// Assign the orbital derivatives to the first bonded atom which is either the H or the virt atom.
							return res.type().bonded_neighbor( acceptor_index )[1];
						}
					} else {
						// case 2: most common case.  SER/THR/TYR.  Assign the orbital derivatives to the hydrogen
						// also covers a case with R-OH-R where the hydrogen can be assumed a) to not be the
						// upstream atom (which is surely going to be one of the two heavy atoms) and b) to be
						// controlled by the same chi which is rotating the downstream heavy atom.
						return res.type().attached_H_begin( acceptor_index );
					}
				} else {
					// multiple bound hydrogens? Like NH2:
					// just return the first hydrogen atom
					return res.type().attached_H_begin( acceptor_index );
				}
			}
		}
	}

	// the default case is to assign the derivatives to the atom to which the orbital belongs
	return acceptor_index;
}

void
OrbitalsScore::assign_orb_orb_derivs(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
)const {
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
					if(orb_orb_rules(res1.atom_type_index(*Aindex),res2.atom_type_index(*Dindex))  ){
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
								numeric::xyzVector< core::Real > const pD(res2.xyz(*Dindex)); // the xyz of the donor
								numeric::xyzVector< core::Real > const pA(res1.xyz(*Aindex)); // the xyz of the acceptor
								numeric::xyzVector< core::Real > const pAO(res1.orbital_xyz(*res1_orb) ); // the xyz of the acceptor orbital
								core::Size orbital_type1 = res1.orbital_type_index(*res1_orb);
								core::Size const orbital1_surrogate_atom = surrogate_atom_for_orbital( res1, *Aindex, chemical::orbitals::orbital_type_enum( orbital_type1 ) );


								numeric::xyzVector< core::Real > const pDO(res2.orbital_xyz(*res2_orb) ); // the xyz of the donor orbital
								core::Size orbital_type2 = res2.orbital_type_index(*res2_orb);
								core::Size const orbital2_surrogate_atom = surrogate_atom_for_orbital( res2, *Dindex, chemical::orbitals::orbital_type_enum( orbital_type2 ) );
								core::Real const orb1_orb2_dist= pAO.distance_squared(pDO);
								if(orb1_orb2_dist < 9){
									core::Real const dist(std::sqrt(orb1_orb2_dist));

									core::Real const cosAOD(cos_of(pA, pAO, pD));
									core::Real const cosDOA(cos_of(pD, pDO, pA));

									core::Real const weight = scfxn_rules_for_weight(false /*is_hydrogen*/, false /*is_backbone*/, orbital_type1, lookup_table_.Hpol_scOrbH, orbital_type2, weights); //dummy for htype

									{ // scope

										/// first evaluate the derivative for the pAO-pDO distance and the "D--DO--A" angle
										Real d_deriv=0;
										Real a_deriv=0;
										lookup_table_.OrbOrbDist_cosDOA_energy(orbital_type1, orbital_type2, dist, cosDOA, orb_orb_E, d_deriv, a_deriv, true);

										Vector f1ad(0),f2ad(0);
										Vector f1ao(0),f2ao(0);
										Vector f1aa(0),f2aa(0);
										core::Real tau(0.0);
										numeric::deriv::angle_p1_deriv( pD, pDO, pA, tau, f1ad, f2ad  );
										numeric::deriv::angle_p2_deriv( pD, pDO, pA, tau, f1ao, f2ao  );
										numeric::deriv::angle_p1_deriv( pA, pDO, pD, tau, f1aa, f2aa  );

										Real const neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;

										r2_atom_derivs[ *Dindex ].f1() += weight * neg_sine_tau * a_deriv * f1ad;
										r2_atom_derivs[ *Dindex ].f2() += weight * neg_sine_tau * a_deriv * f2ad;

										r2_atom_derivs[ orbital2_surrogate_atom ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
										r2_atom_derivs[ orbital2_surrogate_atom ].f2() += weight * neg_sine_tau * a_deriv * f2ao;

										r1_atom_derivs[ *Aindex ].f1() += weight * neg_sine_tau * a_deriv * f1aa;
										r1_atom_derivs[ *Aindex ].f2() += weight * neg_sine_tau * a_deriv * f2aa;

										Vector f1ddo(0),f2ddo(0);
										Real dis(0.0);
										numeric::deriv::distance_f1_f2_deriv( pDO, pAO, dis, f1ddo, f2ddo );

										r2_atom_derivs[ orbital2_surrogate_atom ].f1() += weight * d_deriv * f1ddo;
										r2_atom_derivs[ orbital2_surrogate_atom ].f2() += weight * d_deriv * f2ddo;

										r1_atom_derivs[ orbital1_surrogate_atom ].f1() -= weight * d_deriv * f1ddo;
										r1_atom_derivs[ orbital1_surrogate_atom ].f2() -= weight * d_deriv * f2ddo;
									}

									{ // scope

										/////////
										//this starts the AOH derivative calculation

										Real d_deriv=0;
										Real a_deriv=0;
										lookup_table_.OrbOrbDist_cosAOD_energy(orbital_type1, orbital_type2, dist, cosAOD, orb_orb_E, d_deriv, a_deriv, true);

										Vector f1ad(0),f2ad(0);
										Vector f1ao(0),f2ao(0);
										Vector f1aa(0),f2aa(0);
										Real tau=0;
										numeric::deriv::angle_p1_deriv( pA, pAO, pD, tau, f1aa, f2aa  );
										numeric::deriv::angle_p2_deriv( pA, pAO, pD, tau, f1ao, f2ao  );
										numeric::deriv::angle_p1_deriv( pD, pAO, pA, tau, f1ad, f2ad  );

										Real const neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;core::Size orbital_surrogate_atom_index = atom_index;

										r1_atom_derivs[ *Aindex ].f1() += weight * neg_sine_tau * a_deriv * f1aa;
										r1_atom_derivs[ *Aindex ].f2() += weight * neg_sine_tau * a_deriv * f2aa;

										r1_atom_derivs[ orbital1_surrogate_atom ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
										r1_atom_derivs[ orbital1_surrogate_atom ].f2() += weight * neg_sine_tau * a_deriv * f2ao;

										r2_atom_derivs[ *Dindex ].f1() += weight * neg_sine_tau * a_deriv * f1ad;
										r2_atom_derivs[ *Dindex ].f2() += weight * neg_sine_tau * a_deriv * f2ad;

										Vector f1dao(0), f2dao(0);
										Real dis= 0.0;
										numeric::deriv::distance_f1_f2_deriv( pAO, pDO, dis, f1dao, f2dao );

										r1_atom_derivs[ orbital1_surrogate_atom ].f1() += weight * d_deriv * f1dao;
										r1_atom_derivs[ orbital1_surrogate_atom ].f2() += weight * d_deriv * f2dao;

										r2_atom_derivs[ orbital2_surrogate_atom ].f1() -= weight * d_deriv * f1dao;
										r2_atom_derivs[ orbital2_surrogate_atom ].f2() -= weight * d_deriv * f2dao;
									}

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
		core::Size  & atom_index, // the atom on residue 1 with the set of orbitals we are considering
		numeric::xyzVector<core::Real> const & Axyz,
		core::Size const & H_index, // the hydrogen atom on residue 2
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

	bool const backbone_interaction = res1.atom_is_backbone( atom_index ) || res2.atom_is_backbone( H_index );

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
		if ( temp_dist_squared < max_orbital_dist_squared_ ) {
			core::Size const orbital_type = res1.orbital_type_index(*orbital_index);
			core::Size const orbital_lookup_type = lookup_type_for_orbital_type( chemical::orbitals::orbital_type_enum( orbital_type ) );
			// which atom do we assign the f1/f2 vectors to?  For most cases, it will be the acceptor atom, but
			// for sp3 hybridized atoms, it will be the attached atom of the acceptor which is not the atoms' base;
			// that is, for SER/THR, it will be HG, for an ether, it will be the carbon downstream of the acceptor
			core::Size const orbital_surrogate_atom = surrogate_atom_for_orbital( res1, atom_index, chemical::orbitals::orbital_type_enum( orbital_type ) );

			core::Real weight = scfxn_rules_for_weight(true, backbone_interaction, orbital_type, htype, 0, weights);
			if ( weight == 0 ) return;

			//std::cout << "Weight " << backbone_interaction << " " << orbital_type << " " << htype << " " << score_type_for_orb_params( lookup_table_, true, backbone_interaction, orbital_type, htype, chemical::orbitals::C_pi_sp2 ) << " : " << weight << std::endl;

			//This starts the DHO derivative calculation.
			core::Size donor_index(res2.bonded_neighbor(H_index)[1]);
			numeric::xyzVector<core::Real> const & Dxyz(res2.xyz(donor_index));
			core::Real cosDHO(cos_of(Dxyz, Hxyz, Orbxyz ));
			core::Real OrbHdist(std::sqrt(temp_dist_squared));

			lookup_table_.OrbHdist_cosDHO_energy(htype, orbital_lookup_type, OrbHdist, cosDHO, energy, d_deriv, a_deriv, true);
			//std::cout << "DHO: " << res1.seqpos() << " " << res2.seqpos() << " " << weight << " " << energy << " " << d_deriv << " " << a_deriv << std::endl;

			Vector pD(Dxyz);
			Vector pH(Hxyz);
			Vector pO(Orbxyz);

			// a=angle, and then d=donor,h=hydrogen,o=orbital
			Vector f1ad,f2ad;
			Vector f1ah,f2ah;
			Vector f1ao,f2ao;
			core::Real tau(0.0);

			//core::Size orbital_surrogate_atom_index = atom_index;

			numeric::deriv::angle_p1_deriv( pD, pH, pO, tau, f1ad, f2ad );
			numeric::deriv::angle_p2_deriv( pD, pH, pO, tau, f1ah, f2ah );
			numeric::deriv::angle_p1_deriv( pO, pH, pD, tau, f1ao, f2ao );


			Real neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;

			r2_atom_derivs[ H_index ].f1() += weight * neg_sine_tau * a_deriv * f1ah;
			r2_atom_derivs[ H_index ].f2() += weight * neg_sine_tau * a_deriv * f2ah;

			r2_atom_derivs[ donor_index ].f1() += weight * neg_sine_tau * a_deriv * f1ad;
			r2_atom_derivs[ donor_index ].f2() += weight * neg_sine_tau * a_deriv * f2ad;

			r1_atom_derivs[ orbital_surrogate_atom ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
			r1_atom_derivs[ orbital_surrogate_atom ].f2() += weight * neg_sine_tau * a_deriv * f2ao;


			Vector f1dh( 0 ), f2dh( 0 );
			Real dis(0.0);
			numeric::deriv::distance_f1_f2_deriv( pH, pO, dis, f1dh, f2dh );

			r2_atom_derivs[ H_index ].f1() += weight * d_deriv * f1dh;
			r2_atom_derivs[ H_index ].f2() += weight * d_deriv * f2dh;

			r1_atom_derivs[ orbital_surrogate_atom ].f1() -= weight * d_deriv * f1dh;
			r1_atom_derivs[ orbital_surrogate_atom ].f2() -= weight * d_deriv * f2dh;

			/////////
			//this starts the AOH derivative calculation

			d_deriv =0;
			a_deriv =0;

			core::Real cosAOH(cos_of(Axyz, Orbxyz, Hxyz));
			//a little confusing without context. This checks to see if the residue is an aromatic residue. If it is an aromatic residue
			//then we need to check if the orbital we are looking at is the action center orbital. If it is, then we use a separate
			//energy than if it were not. This is done because the action center orbital has a different angle associated with the acceptor
			//orbital hydrogen angle. why? because there is no index for action centers, thefore the Acceptor is the first action center atom

			bool const working_with_action_center_orbital = (
				(( res1.aa() == chemical::aa_tyr || res1.aa() == chemical::aa_phe || res1.aa() == chemical::aa_trp ) && *orbital_index <= 2 ));
			lookup_table_.OrbHdist_cosAOH_energy(htype, orbital_lookup_type, OrbHdist, cosAOH, energy, d_deriv, a_deriv, true, working_with_action_center_orbital );
			//std::cout << "AOH: " << res1.seqpos() << " " << res2.seqpos() << " " << weight << " " << energy << " " << d_deriv << " " << a_deriv << std::endl;

			Vector pA(Axyz);
			//Vector pO(Orbxyz);
			//Vector pH(Hxyz);

			Vector f1aa(0),f2aa(0);
			//f1ad=f2ab
			//=f1ao=f2ao=f1ah=f2ah=0;
			//tau=0;

			numeric::deriv::angle_p1_deriv( pA, pO, pH, tau, f1aa, f2aa  );
			numeric::deriv::angle_p2_deriv( pA, pO, pH, tau, f1ao, f2ao  );
			numeric::deriv::angle_p1_deriv( pH, pO, pA, tau, f1ah, f2ah  );

			neg_sine_tau = -sin( tau ); // d cos(theta) / d theta;core::Size orbital_surrogate_atom_index = atom_index;

			r1_atom_derivs[ atom_index ].f1() += weight * neg_sine_tau * a_deriv * f1aa;
			r1_atom_derivs[ atom_index ].f2() += weight * neg_sine_tau * a_deriv * f2aa;

			r1_atom_derivs[ orbital_surrogate_atom ].f1() += weight * neg_sine_tau * a_deriv * f1ao;
			r1_atom_derivs[ orbital_surrogate_atom ].f2() += weight * neg_sine_tau * a_deriv * f2ao;

			r2_atom_derivs[ H_index ].f1() += weight * neg_sine_tau * a_deriv * f1ah;
			r2_atom_derivs[ H_index ].f2() += weight * neg_sine_tau * a_deriv * f2ah;

			f1dh = f2dh = 0;
			dis= 0.0;
			numeric::deriv::distance_f1_f2_deriv( pH, pO, dis, f1dh, f2dh );

			r1_atom_derivs[ orbital_surrogate_atom ].f1() -= weight * d_deriv * f1dh;
			r1_atom_derivs[ orbital_surrogate_atom ].f2() -= weight * d_deriv * f2dh;

			r2_atom_derivs[ H_index ].f1() += weight * d_deriv * f1dh;
			r2_atom_derivs[ H_index ].f2() += weight * d_deriv * f2dh;
		}
	}
}


}
}
}
