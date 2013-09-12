// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_FullAtomVDW_BasePhosphate
/// @brief  RNA_FullAtomVDW_BasePhosphate energy method class implementation
/// @author Parin Sripakdeevong

// Unit Headers
#include <core/scoring/rna/RNA_FullAtomVDW_BasePhosphate.hh>
#include <core/scoring/rna/RNA_FullAtomVDW_BasePhosphateCreator.hh>
#include <core/scoring/etable/Etable.hh>

// Package Headers
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/types.hh>


// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/conversions.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/NeighborList.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
// C++


using namespace core::chemical::rna;


namespace core {
namespace scoring {
namespace rna {



/// @details This must return a fresh instance of the RNA_FullAtomVDW_BasePhosphateCreator class
/// never an instance already in use
methods::EnergyMethodOP
RNA_FullAtomVDW_BasePhosphateCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {


	etable::Etable const & etable_in= *( ScoringManager::get_instance()->etable( options.etable_type() ) );

	etable::TableLookupEtableEnergy etable_energy_in( etable_in, options );

	// [Note, previously created pointer to a regular object instead of OP -- energy leak.]
	//	etable::EtableEnergy & etable_energy_in= (*(new etable::EtableEnergy( etable_in, options ) ) );

	//RNA_FullAtomVDW_BasePhosphateCreator is friend of BaseEtableEnergy.hh and hence can access its private variables.
	// No.  Don't go around declaring other classes your friends.  That's just pisses on everyone else's efforts to create
	// an object oriented program.
	etable_energy_in.intrares_evaluator().set_scoretypes(
		fa_intra_RNA_base_phos_atr,
		fa_intra_RNA_base_phos_rep,
		fa_intra_RNA_base_phos_sol );

	return new RNA_FullAtomVDW_BasePhosphate( etable_energy_in, etable_in);
}


ScoreTypes
RNA_FullAtomVDW_BasePhosphateCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_intra_RNA_base_phos_atr );
	sts.push_back( fa_intra_RNA_base_phos_rep );
	sts.push_back( fa_intra_RNA_base_phos_sol );
	return sts;
}


/// constructor
RNA_FullAtomVDW_BasePhosphate::RNA_FullAtomVDW_BasePhosphate(
	etable::TableLookupEtableEnergy const & etable_energy_in,
	etable::Etable const & etable_in
):
	parent( new RNA_FullAtomVDW_BasePhosphateCreator ),
	etable_energy_( etable_energy_in), //Hacky thing, created the etable_energy_ energy_method to get access to its function.
	//etable_energy_( *(new etable::EtableEnergy(etable_energy_in)) ),
	//etable_energy_( *( new etable::EtableEnergy( etable_in, opts ) ) ),
	etable_(etable_in)
{

	//etable_energy_.set_scoretypes( fa_intra_RNA_base_phos_atr, fa_intra_RNA_base_phos_rep, unfolded);
	//	etable_energy_( *( new etable::BaseEtableEnergy< etable::EtableEnergy >(new etable::EtableEnergyCreator, etable_in, options, fa_atr, fa_rep, fa_sol ) ) ),
	//	etable_energy_( etable::BaseEtableEnergy< etable::EtableEnergy > (new etable::EtableEnergyCreator, etable_in, options, unfolded, unfolded, unfolded ) )
	//etable_energy_( etable::EtableEnergyCreator::create_energy_method(options) )

	if ( basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ] ){
		utility_exit_with_message(  "RNA_LJ_BaseEnergy not compatible with analytic_etable_evaluation yet -- rerun with flag -analytic_etable_evaluation false." );
	}

}


RNA_FullAtomVDW_BasePhosphate::~RNA_FullAtomVDW_BasePhosphate() {}


/// clone
methods::EnergyMethodOP
RNA_FullAtomVDW_BasePhosphate::clone() const
{

	return new RNA_FullAtomVDW_BasePhosphate( *this );
}




void
RNA_FullAtomVDW_BasePhosphate::residue_fast_pair_energy_attached_H(
	conformation::Residue const & res1,
	int const atomno1,
	conformation::Residue const & res2,
	Size const atomno2,
	Size const at1hbegin, //at1hbegin and at1hend define a range of hydrogen atom indices -- those h's bound to at1
	Size const at1hend,
	Size const at2hbegin,
	Size const at2hend,
	EnergyMap & emap
) const
{
	using conformation::Atom;

	Weight weight( 1.0 );

	Atom const & atom1( res1.atom( atomno1 ) );
	Atom const & atom2( res2.atom( atomno2 ) );


	// Heavy Atom in res1 to Hs in res2
	for ( Size i = at2hbegin; i<= at2hend; ++i )
	{
		Atom const & H2( res2.atom( i ) );
		weight = 1.0;
		etable_energy_.pair_energy_H( atom1, H2, weight, emap );
	}


	// Hs in res1 to heavy Atom and Hs in res2
	for ( Size i = at1hbegin; i<= at1hend; ++i )
	{
		Atom const & H1( res1.atom(i) );
		weight = 1.0;
		// H in res1 to heavy Atom in res2
		etable_energy_.pair_energy_H( H1, atom2, weight, emap );

		// H in res1 to Hs in res2
		for ( Size j = at2hbegin; j<= at2hend; ++j ) {
			Atom const & H2( res2.atom(j) );
			weight = 1.0f;
			etable_energy_.pair_energy_H( H1, H2, weight, emap );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void
RNA_FullAtomVDW_BasePhosphate::residue_energy(
		conformation::Residue const & rsd,
		EnergyMap & emap 	) const
{
	using conformation::Atom;

	//runtime_assert(false);
	//utility_exit_with_message("Why am I here?");

	if(rsd.is_RNA()==false) return;


	DistanceSquared dsq;

	//Weight const weight=1.0;
	Real const weight=1.0;

	// get hydrogen interaction cutoff
	Real const Hydrogen_interaction_cutoff2=( etable_energy_.hydrogen_interaction_cutoff2() );

	EnergyMap mock_emap;

	typedef utility::vector1< Size > const & vect;

	vect rhbegin( rsd.attached_H_begin() );
	vect rhend(   rsd.attached_H_end()   );

	Size const rsdnheavyatoms = rsd.nheavyatoms();

	// Atom pairs
	for ( Size i=1; i <= rsdnheavyatoms; ++i ) {
		Atom const & atom1( rsd.atom(i) );
		for ( Size j=i+1; j <= rsdnheavyatoms; ++j ) {
			Atom const & atom2( rsd.atom(j) );

			bool base_phosphate_atom_pair=false;

			if( rsd.RNA_type().atom_is_phosphate( i ) && rsd.RNA_type().is_RNA_base_atom( j ) ) base_phosphate_atom_pair=true;
			if( rsd.RNA_type().atom_is_phosphate( j ) && rsd.RNA_type().is_RNA_base_atom( i ) ) base_phosphate_atom_pair=true;

			if(base_phosphate_atom_pair==false) continue;

			if(rsd.path_distance( i, j ) < 4){
				 std::cout << "rsd.name3()=" << rsd.name3() << " rsd.seqpos() = " << rsd.seqpos() << std::endl;
				 std::cout << "i= " << i << " j=" << j << std::endl;
				 std::cout << "rsd.atom_name(i)= " << rsd.atom_name(i) << " rsd.atom_name(j)=" << rsd.atom_name(j) << std::endl;
				 std::cout << "rsd.path_distance( i, j )=" << rsd.path_distance( i, j ) << std::endl;

				 std::cout << "rsdnheavyatoms= " << rsdnheavyatoms << std::endl;
				 std::cout << "rsd.first_sidechain_atom()= " << rsd.first_sidechain_atom() << std::endl;
				 for ( Size atomno=1; atomno <= rsdnheavyatoms; ++atomno ) {
				 	 std::cout << "rsd.atom_name(" << atomno << " )= " << rsd.atom_name(atomno) << std::endl;
				 }

				 std::cout << "rsd.nchi()=" << rsd.nchi() << std::endl;

				 utility::vector1< utility::vector1< Size > > const & chi_atoms_list=rsd.type().chi_atoms();
				 utility::vector1< utility::vector1< Size > > const & all_chi_atoms_list=rsd.type().atoms_last_controlled_by_chi(); /*chi # 1 must be nucleic acid "chi"*/


				 for(Size chi_id=1; chi_id<=chi_atoms_list.size(); chi_id++){
					 std::cout << std::endl;
					 std::cout << "chi_atoms_list[" << chi_id<< "].size()=" << chi_atoms_list[chi_id].size() << std::endl;
					 for(Size n=1; n<=chi_atoms_list[chi_id].size(); n++){
						 Size const atomno=chi_atoms_list[chi_id][n];
					 	 std::cout << "rsd.atom_name(" << atomno << " )= " << rsd.atom_name(atomno) << std::endl;
					 }

				 }
				 for(Size chi_id=1; chi_id<=chi_atoms_list.size(); chi_id++){
					 std::cout << std::endl;
					 std::cout << "all_chi_atoms_list[" << chi_id<< "].size()=" << all_chi_atoms_list[chi_id].size() << std::endl;
					 for(Size n=1; n<=all_chi_atoms_list[chi_id].size(); n++){
						 Size const atomno=all_chi_atoms_list[chi_id][n];
					 	 std::cout << "rsd.atom_name(" << atomno << " )= " << rsd.atom_name(atomno) << std::endl;
					 }

				 }

				 utility_exit_with_message("rsd.path_distance( i, j ) < 4");
			}

			Real const START_fa_atr=mock_emap[fa_intra_RNA_base_phos_atr];
			Real const START_fa_rep=mock_emap[fa_intra_RNA_base_phos_rep];

			etable_energy_.atom_pair_energy( atom1, atom2, weight, mock_emap, dsq );

			if(rsd.is_virtual(i) || rsd.is_virtual(j)){
				Real const diff_fa_atr_E=mock_emap[fa_intra_RNA_base_phos_atr]-START_fa_atr;
				Real const diff_fa_rep_E=mock_emap[fa_intra_RNA_base_phos_rep]-START_fa_rep;

				if(diff_fa_atr_E>0.001 || diff_fa_atr_E<-0.001){
					std::cout << "At least one of the atoms in the pair " << i << "," << j << " | res= " << rsd.seqpos() << " is virtual";
					std::cout << " but diff_fa_atr_E= " << diff_fa_atr_E << " is non-zero!" <<std::endl;
					utility_exit_with_message("at least one of the atom in the pair is virtual BUT diff_fa_atr_E>0.001 || diff_fa_atr_E<-0.001");
				}



				if(diff_fa_rep_E>0.001 || diff_fa_rep_E<-0.001){
					std::cout << "At least one of the atoms in the pair " << i << "," << j << " | res= " << rsd.seqpos() << " is virtual";
					std::cout << " but diff_fa_rep_E= " << diff_fa_rep_E << " is non-zero!" <<std::endl;
					utility_exit_with_message("at least one of the atom in the pair is virtual BUT diff_fa_rep_E>0.001 || diff_fa_rep_E<-0.001");
				}
			}

			if ( dsq < Hydrogen_interaction_cutoff2 ) {
				residue_fast_pair_energy_attached_H(rsd, i, rsd, j, rhbegin[ i ], rhend[ i ], rhbegin[ j ], rhend[ j ], mock_emap);
			}
		}
	}

	emap[ fa_intra_RNA_base_phos_atr ] += mock_emap[ fa_intra_RNA_base_phos_atr ]; //Fix from = to += on Jan 10, 2012, since now emap is the total_energy, accumulated over all residues)
	emap[ fa_intra_RNA_base_phos_rep ] += mock_emap[ fa_intra_RNA_base_phos_rep ]; //Fix from = to += on Jan 10, 2012, since now emap is the total_energy, accumulated over all residues)

	//Note that right now fa_intra_RNA_base_phos_sol is not include since it is not currently used in the RNA force-field. The term have also not yet been tested!//
	//If you want to implement this term, please ensure your implementation works properly (i.e. perform numerical_derivative_check() and etc),
	//before committing the code to TRUNK!)
	//Parin S. (sripakpa@stanford.edu). Jan 11, 2012

}


///////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomVDW_BasePhosphate::residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		EnergyMap & emap 	) const {

	return residue_energy( rsd, emap );
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomVDW_BasePhosphate::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & /*sfxn*/, // needed for non-nblist minimization
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	if(weights[ fa_intra_RNA_base_phos_sol] != 0.0){
		//Please refer to paragraph at the end of RNA_FullAtomVDW_BasePhosphate::residue_energy for explanation.
		//Again, if you want to implement this term, please ensure your implementation works properly (i.e. perform numerical_derivative_check() and etc),
		//before committing the code to TRUNK!)
		// Parin S. (sripakpa@stanford.edu). Jan 11, 2012
		utility_exit_with_message("weights[ fa_intra_RNA_base_phos_sol ] != 0.0, but this term is not yet implemented!");
	}

	Size const seq_num = id.rsd();

	if(pose.residue( seq_num ).is_RNA()==false) return;

	conformation::Residue const & rsd= pose.residue( seq_num );

	conformation::Atom const & atom1( rsd.atom( id.atomno() ) );

	Real const cp_weight=1.0;

	//std::cout << "eval_atom_deriv: seq_num=" << id.rsd() << " id.atomno()=" << id.atomno() << "[" << rsd.atom_name(id.atomno()) <<"]" << std::endl;


	Vector f1,f2;
	for ( Size nbr_atomno=1; nbr_atomno<=rsd.natoms(); nbr_atomno++ ) {

 		if( rsd.is_virtual( id.atomno() ) ) continue; //Is this necessary?
 		if( rsd.is_virtual( nbr_atomno  ) ) continue; //Is this necessary?

		bool base_phosphate_atom_pair=false;

		if( rsd.RNA_type().atom_is_phosphate( id.atomno() ) && rsd.RNA_type().is_RNA_base_atom( nbr_atomno  ) ) base_phosphate_atom_pair=true;
		if( rsd.RNA_type().atom_is_phosphate( nbr_atomno  ) && rsd.RNA_type().is_RNA_base_atom( id.atomno() ) ) base_phosphate_atom_pair=true;

		if(base_phosphate_atom_pair==false) continue;

		if(rsd.path_distance( nbr_atomno, id.atomno() ) < 4) utility_exit_with_message("rsd.path_distance( nbr_atomno, id.atomno() ) < 4");

		conformation::Atom const & atom2( rsd.atom( nbr_atomno ) );

		Real const dE_dR_over_r= etable_energy_.intrares_evaluator().eval_dE_dR_over_r( atom1, atom2, weights, f1, f2 );

		if ( dE_dR_over_r != 0.0 ) {
			F1 += dE_dR_over_r * cp_weight * f1;
			F2 += dE_dR_over_r * cp_weight * f2;
		}

	}

}

// check compatibility with atomtypeset
///////////////////////////////////////////////////////////////////////////////

void
RNA_FullAtomVDW_BasePhosphate::setup_for_scoring( pose::Pose & pose, scoring::ScoreFunction const & /*scfxn*/ ) const
{


	if (pose.total_residue()>0) {
		if ( pose.residue(1).type().atom_type_set_ptr() != etable_.atom_set() ) {
			utility_exit_with_message( "Illegal attempt to score with non-identical atom set between pose and etable " );
		}
	}

}


///@details Make sure that the neighborlist is up-to-date before evaluating derivatives
void
RNA_FullAtomVDW_BasePhosphate::setup_for_derivatives( pose::Pose & pose, scoring::ScoreFunction const & /*scfxn*/ ) const
{

	if (pose.total_residue()>0) {
		if ( pose.residue(1).type().atom_type_set_ptr() != etable_.atom_set() ) {
			utility_exit_with_message( "Illegal attempt to score with non-identical atom set between pose and etable " );
		}
	}

}

///@details Make sure that the neighborlist is up-to-date before evaluating derivatives
/*
void
RNA_FullAtomVDW_BasePhosphate::setup_for_packing( pose::Pose & pose, pack::task::PackerTask const & ) const
{

	//NOT SURE IF THIS WORKS!


	if (pose.total_residue()>0) {
		if ( pose.residue(1).type().atom_type_set_ptr() != etable_.atom_set() ) {
			utility_exit_with_message( "Illegal attempt to score with non-identical atom set between pose and etable " );
		}
	}

}
*/

void RNA_FullAtomVDW_BasePhosphate::indicate_required_context_graphs( utility::vector1< bool > & ) const{}


core::Size
RNA_FullAtomVDW_BasePhosphate::version() const
{
	return 1; // First version, created by Parin Sripakdeevong (sripakpa@stanford.edu), Jan 2012.
}




} // rna
} // scoring
} // core

