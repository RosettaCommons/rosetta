// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/elec/RNA_FA_ElecEnergy.cc
/// @brief  Electrostatics energy method for RNA class implementation
/// @author Rhiju Das
/// Edited by Joseph Yesselm (9.6.13)


// Unit headers
#include <core/scoring/elec/RNA_FA_ElecEnergy.hh>
#include <core/scoring/elec/RNA_FA_ElecEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

#include <core/kinematics/MinimizerMapBase.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/DerivVectorPair.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/RotamerSetBase.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


// ObjexxFCL headers


// C++

// Open questions:
// Why not just use the existing FA_ElecEnergy plus add terms to it?

namespace core {
namespace scoring {
namespace elec {


/// @details This must return a fresh instance of the RNA_FA_ElecEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_FA_ElecEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new RNA_FA_ElecEnergy( options ) );
}

ScoreTypes
RNA_FA_ElecEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_elec_rna_phos_phos );
	sts.push_back( fa_elec_rna_phos_sugr );
	sts.push_back( fa_elec_rna_phos_base );
	sts.push_back( fa_elec_rna_sugr_sugr );
	sts.push_back( fa_elec_rna_sugr_base );
	sts.push_back( fa_elec_rna_base_base );
	return sts;
}

////////////////////////////////////////////////////////////////////////////
RNA_FA_ElecEnergy::RNA_FA_ElecEnergy(
	methods::EnergyMethodOptions const & options
):
	parent( options )
{
	set_score_types( methods::EnergyMethodCreatorOP( new RNA_FA_ElecEnergyCreator ) );
}


////////////////////////////////////////////////////////////////////////////
RNA_FA_ElecEnergy::RNA_FA_ElecEnergy( RNA_FA_ElecEnergy const & src ):
	parent( src )
{
	set_score_types( methods::EnergyMethodCreatorOP( new RNA_FA_ElecEnergyCreator ) );
}

/// clone
methods::EnergyMethodOP
RNA_FA_ElecEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNA_FA_ElecEnergy( *this ) );
}

//////////////////////////////////////
bool is_phosphate_2( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_phosphate_2( rsd, rsd.atom_base( i ) );
	} else {
		//MAGIC NUMBERS! BAD!
		return ( i==1 || i==2 || i==3 || i==4 || i==9 ); //P, OP2, OP1, O5', O3'
	}
}

//////////////////////////////////////
bool is_sugar_2( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_sugar_2( rsd, rsd.atom_base( i ) );
	} else {
		return ( ( i < rsd.first_sidechain_atom() && !is_phosphate_2(rsd, i) ) ||
			i == 12 /*hey this is really bad!*/ );
	}
}

//////////////////////////////////////
bool is_base_2( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_base_2( rsd, rsd.atom_base( i ) );
	} else {
		return (!is_sugar_2(rsd, i) && !is_phosphate_2(rsd, i) );
	}
}

RNAAtomType
assign_rna_atom_type(
	conformation::Residue const & rsd,
	Size const i
) {
	if ( rsd.atom_is_hydrogen(i) ) {
		return assign_rna_atom_type( rsd, rsd.atom_base( i ) );
	}

	// AMW TODO: this cannot stand.
	if ( i == 1 || i == 2 || i == 3 || i == 4 || i == 9 ) { //P, OP2, OP1, O5', O3'
		return PHOSPHATE;
	} else if ( ( i < rsd.first_sidechain_atom() && !is_phosphate_2(rsd, i) ) || i == 12 ) {
		return SUGAR;
	} else {
		return BASE;
	}
}

void
RNA_FA_ElecEnergy::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( ! pose.energies().use_nblist() ) return;

	// stash our nblist inside the pose's energies object
	Energies & energies( pose.energies() );

	// setup the atom-atom nblist
	NeighborListOP nblist;
	Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? option[ run::nblist_autoupdate_narrow ] : 1.5;
	Real const XX = coulomb().max_dis() + 2 * tolerated_motion;
	nblist = NeighborListOP( new NeighborList( min_map.domain_map(), XX*XX, XX*XX, XX*XX) );
	if ( pose.energies().use_nblist_auto_update() ) {
		nblist->set_auto_update( tolerated_motion );
	}
	// this partially becomes the EtableEnergy classes's responsibility
	nblist->setup( pose, sfxn, *this);
	energies.set_nblist( EnergiesCacheableDataType::RNA_ELEC_NBLIST, nblist );
}

void
RNA_FA_ElecEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sfxn ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	pose.update_residue_neighbors();
	setup_for_scoring( pose, sfxn );
}

void
RNA_FA_ElecEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const
{
	//std::cout << "amw In setup_for_scoring, hoping nblist is already there." << std::endl;
	pose.update_residue_neighbors();

	if ( ! pose.energies().use_nblist() ) return;

	NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::RNA_ELEC_NBLIST ) );
	nblist.prepare_for_scoring( pose, scfxn, *this );
}


// The FA_ElectEnergy method stores a vector of rotamer trie objects in the Energies
// object for use in rapid rotamer/background energy calculations.  Overrides default
// do-nothing behavior.
void
RNA_FA_ElecEnergy::setup_for_packing(
	pose::Pose &,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}

// @brief Creates a rotamer trie for the input set of rotamers and stores the trie
// in the rotamer set.
void
RNA_FA_ElecEnergy::prepare_rotamers_for_packing(
	pose::Pose const &,
	conformation::RotamerSetBase &
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}


// @brief Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
RNA_FA_ElecEnergy::update_residue_for_packing(
	pose::Pose &,
	Size
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


///
void
RNA_FA_ElecEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scfxn,
	EnergyMap & emap
) const {
	if ( pose.energies().use_nblist() ) return;
	if ( ! rsd1.is_RNA() || ! rsd2.is_RNA() ) { return; }

	if ( scfxn.has_nonzero_weight(fa_elec_rna_phos_phos) ) {
		emap[ fa_elec_rna_phos_phos ] += rna_fa_elec_one_way(rsd1,rsd2,PHOSPHATE,PHOSPHATE);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_sugr_sugr) ) {
		emap[ fa_elec_rna_sugr_sugr ] += rna_fa_elec_one_way(rsd1,rsd2,SUGAR,SUGAR);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_base_base) ) {
		emap[ fa_elec_rna_base_base ] += rna_fa_elec_one_way(rsd1,rsd2,BASE,BASE);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_phos_sugr) ) {
		emap[ fa_elec_rna_phos_sugr ] += rna_fa_elec_one_way(rsd1,rsd2,PHOSPHATE,SUGAR) +
			rna_fa_elec_one_way(rsd1,rsd2,SUGAR,PHOSPHATE);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_phos_base) ) {
		emap[ fa_elec_rna_phos_base ] += rna_fa_elec_one_way(rsd1,rsd2,PHOSPHATE,BASE) +
			rna_fa_elec_one_way(rsd1,rsd2,BASE,PHOSPHATE);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_sugr_base) ) {
		emap[ fa_elec_rna_sugr_base ] += rna_fa_elec_one_way(rsd1,rsd2,SUGAR,BASE) +
			rna_fa_elec_one_way(rsd1,rsd2,BASE,SUGAR);
	}
}

bool
RNA_FA_ElecEnergy::minimize_in_whole_structure_context( pose::Pose const & pose ) const
{
	return pose.energies().use_nblist_auto_update();
}

void
RNA_FA_ElecEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const {
	if ( ! rsd1.is_RNA() || ! rsd2.is_RNA() ) { return; }
	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( rna_elec_pair_nblist ) ));
	auto const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( rna_elec_pair_nblist ) ) );
	Real dsq, score( 0.0 );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	for ( auto const & neighb : neighbs ) {
		score += score_atom_pair( rsd1, rsd2, neighb.atomno1(), neighb.atomno2(), emap, sfxn, neighb.weight(), dsq );
	}
	// Nothing done with score, only emap
}

void
RNA_FA_ElecEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData & pair_data
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( pose.energies().use_nblist_auto_update() ) return;

	etable::count_pair::CountPairFunctionCOP count_pair = get_count_pair_function( rsd1, rsd2 );
	debug_assert( rsd1.seqpos() < rsd2.seqpos() );

	// update the existing nblist if it's already present in the min_data object
	ResiduePairNeighborListOP nblist( utility::pointer::static_pointer_cast< core::scoring::ResiduePairNeighborList > ( pair_data.get_data( rna_elec_pair_nblist ) ));
	if ( ! nblist ) nblist = ResiduePairNeighborListOP( new ResiduePairNeighborList );

	/// STOLEN CODE!
	Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
	Real const XX2 = std::pow( coulomb().max_dis() + 2*tolerated_narrow_nblist_motion, 2 );

	nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair );

	pair_data.set_data( rna_elec_pair_nblist, nblist );
}

Real
RNA_FA_ElecEnergy::score_atom_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Size const at1,
	Size const at2,
	EnergyMap & emap,
	ScoreFunction const & sfxn,
	Real const cpweight,
	Real & d2
) const {
	Real energy = cpweight * coulomb().eval_atom_atom_fa_elecE(
		rsd1.xyz(at1), rsd1.atomic_charge(at1),
		rsd2.xyz(at2), rsd2.atomic_charge(at2), d2);

	//std::cout << "energy from atoms " << rsd1.seqpos() << " " << at1 << "  " << rsd2.seqpos() << " " << at2 << " is " << energy << std::endl;

	RNAAtomType const t1 = assign_rna_atom_type( rsd1, at1 );
	RNAAtomType const t2 = assign_rna_atom_type( rsd2, at2 );

	if ( t1 == PHOSPHATE && t2 == PHOSPHATE ) {
		if ( sfxn.has_nonzero_weight( fa_elec_rna_phos_phos ) ) emap[ fa_elec_rna_phos_phos ] += energy;
	}

	if ( t1 == SUGAR && t2 == SUGAR ) {
		if ( sfxn.has_nonzero_weight( fa_elec_rna_sugr_sugr ) ) emap[ fa_elec_rna_sugr_sugr ] += energy;
	}

	if ( t1 == BASE && t2 == BASE ) {
		if ( sfxn.has_nonzero_weight( fa_elec_rna_base_base ) ) emap[ fa_elec_rna_base_base ] += energy;
	}

	if ( ( t1 == PHOSPHATE && t2 == SUGAR ) || ( ( t2 == PHOSPHATE && t1 == SUGAR ) ) ) {
		if ( sfxn.has_nonzero_weight( fa_elec_rna_phos_sugr ) ) emap[ fa_elec_rna_phos_sugr ] += energy;
	}

	if ( ( t1 == PHOSPHATE && t2 == BASE ) || ( ( t2 == PHOSPHATE && t1 == BASE ) ) ) {
		if ( sfxn.has_nonzero_weight( fa_elec_rna_phos_base ) ) emap[ fa_elec_rna_phos_base ] += energy;
	}

	if ( ( t1 == SUGAR && t2 == BASE ) || ( ( t2 == SUGAR && t1 == BASE ) ) ) {
		if ( sfxn.has_nonzero_weight( fa_elec_rna_sugr_base ) ) emap[ fa_elec_rna_sugr_base ] += energy;
	}

	return energy;
}


//////////////////////////////////////////////////////////////////////////////////
// Sept 5th 2013, implemetnation of an optmized design function only phos-phos interactions are evaluated!
//////////////////////////////////////////////////////////////////////////////////

void
RNA_FA_ElecEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const {
	if ( set1.num_rotamers() >= 1 && set2.num_rotamers() >= 1 &&
			set1.rotamer(1)->is_RNA() && set2.rotamer(1)->is_RNA() ) {
		grandparent::evaluate_rotamer_pair_energies( set1, set2, pose, sfxn, weights, energy_table );
	} // else, non rna rna interaction; early return
}

void
RNA_FA_ElecEnergy::backbone_backbone_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	emap[ fa_elec_rna_phos_phos ] +=  rna_fa_elec_one_way(rsd1, rsd2, PHOSPHATE, PHOSPHATE);
}


Real
RNA_FA_ElecEnergy::rna_fa_elec_one_way(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	RNAAtomType const & type1 ,
	RNAAtomType const & type2
) const {
	Real energy( 0.0 );

	using namespace etable::count_pair;

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	for ( Size i=1, i_end = rsd1.natoms(); i<= i_end; ++i ) {

		RNAAtomType const atom_type_1 = assign_rna_atom_type(rsd1, i);
		if ( atom_type_1 != type1 ) continue;

		Vector const & i_xyz( rsd1.xyz(i) );
		Real const i_charge( rsd1.atomic_charge(i) );

		for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {

			RNAAtomType const atom_type_2 = assign_rna_atom_type(rsd2, j);
			if ( atom_type_2 != type2 ) continue;

			Real weight(1.0);
			Size path_dist( 0 );
			if ( ! cpfxn->count( i, j, weight, path_dist ) ) continue;

			energy += weight * coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), rsd2.atomic_charge(j));
		}
	}
	return energy;
}

void
RNA_FA_ElecEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const {
	if ( set.num_rotamers() >= 1 && set.rotamer(1)->is_RNA() && residue.is_RNA() ) {
		grandparent::evaluate_rotamer_background_energies( set, residue, pose, sfxn, weights, energy_vector );
	} // else, non rna rna interaction; early return
}

void
RNA_FA_ElecEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	if ( ! pose.energies().use_nblist() || ! pose.energies().use_nblist_auto_update() ) return;

	NeighborList const & nblist
		( pose.energies().nblist( EnergiesCacheableDataType::RNA_ELEC_NBLIST ) );

	nblist.check_domain_map( pose.energies().domain_map() );
	//utility::vector1< conformation::Residue const * > resvect;
	//resvect.reserve( pose.total_residue() );
	//for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
	// resvect.push_back( & pose.residue( ii ) );
	//}

	//Real total_score( 0.0 );
	for ( Size i = 1, i_end = pose.size(); i<= i_end; ++i ) {
		conformation::Residue const & ires = pose.residue( i );
		for ( Size ii=1, ii_end=ires.natoms(); ii<= ii_end; ++ii ) {
			AtomNeighbors const & nbrs( nblist.upper_atom_neighbors(i,ii) );

			bool const atom1_is_base = is_base_2(ires, ii);
			bool const atom1_is_sugar = is_sugar_2(ires, ii);
			bool const atom1_is_phosphate = is_phosphate_2(ires, ii);
			debug_assert( atom1_is_base || atom1_is_sugar || atom1_is_phosphate );

			for ( AtomNeighbor const & nbr : nbrs ) {
				Size const  j( nbr.rsd() );
				Size const jj( nbr.atomno() );

				conformation::Residue const & jres = pose.residue( j );

				bool const atom2_is_base = is_base_2(jres, jj);
				bool const atom2_is_sugar = is_sugar_2(jres, jj);
				bool const atom2_is_phosphate = is_phosphate_2(jres, jj);
				debug_assert( atom2_is_base || atom2_is_sugar || atom2_is_phosphate );

				core::Real wt_envdep = 1.0;
				Real score = nbr.weight() *
					coulomb().eval_atom_atom_fa_elecE( ires.xyz(ii), ires.atomic_charge(ii), jres.xyz(jj), jres.atomic_charge(jj) );

				if ( atom1_is_base && atom2_is_base ) {
					totals[ fa_elec_rna_base_base ] += wt_envdep*score;
				} else if (  (atom1_is_base && atom2_is_sugar)  || (atom1_is_sugar && atom2_is_base) ) {
					totals[ fa_elec_rna_sugr_base ] += wt_envdep*score;
				} else if (  (atom1_is_base && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_base) ) {
					totals[ fa_elec_rna_phos_base ] += wt_envdep*score;
				} else if (  (atom1_is_sugar && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_sugar) ) {
					totals[ fa_elec_rna_phos_sugr ] += wt_envdep*score;
				} else if (  (atom1_is_sugar && atom2_is_sugar)  ) {
					totals[ fa_elec_rna_sugr_sugr ] += wt_envdep*score;
				} else if (  (atom1_is_phosphate && atom2_is_phosphate)  ) {
					totals[ fa_elec_rna_phos_phos ] += wt_envdep*score;
				}
			}
		}
	}
}

Real rna_elec_weight(
	EnergyMap const & weights,
	bool const b1,
	bool const s1,
	bool const p1,
	bool const b2,
	bool const s2,
	bool const p2
) {
	if ( b1 && b2 ) {
		return weights[ fa_elec_rna_base_base ];
	}
	if ( s1 && s2 ) {
		return weights[ fa_elec_rna_sugr_sugr ];
	}
	if ( p1 && p2 ) {
		return weights[ fa_elec_rna_phos_phos ];
	}
	if ( ( b1 && s2 ) || ( s1 && b2 ) ) {
		return weights[ fa_elec_rna_sugr_base ];
	}
	if ( ( b1 && p2 ) || ( p1 && b2 ) ) {
		return weights[ fa_elec_rna_phos_base ];
	}
	if ( ( s1 && p2 ) || ( p1 && s2 ) ) {
		return weights[ fa_elec_rna_phos_sugr ];
	}

	utility_exit_with_message( "Atom is neither base nor sugar nor phosphate?" );
}

void
RNA_FA_ElecEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( rna_elec_pair_nblist ) ));

	auto const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( rna_elec_pair_nblist ) ) );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	for ( auto const & neighb : neighbs ) {
		Vector const & atom1xyz( rsd1.xyz( neighb.atomno1() ) );
		Vector const & atom2xyz( rsd2.xyz( neighb.atomno2() ) );

		Real const at1_charge( rsd1.atomic_charge( neighb.atomno1() ) );
		Real const at2_charge( rsd2.atomic_charge( neighb.atomno2() ) );

		Vector f2 = ( atom1xyz - atom2xyz );
		Real const dis2( f2.length_squared() );
		Real const dE_dr_over_r = neighb.weight() * coulomb().eval_dfa_elecE_dr_over_r( dis2, at1_charge, at2_charge );

		if ( dE_dr_over_r == 0.0 ) continue;

		bool const atom1_is_base = is_base_2(rsd1, neighb.atomno1());
		bool const atom1_is_sugar = is_sugar_2(rsd1, neighb.atomno1());
		bool const atom1_is_phosphate = is_phosphate_2(rsd1, neighb.atomno1());
		debug_assert( atom1_is_base || atom1_is_sugar || atom1_is_phosphate );

		bool const atom2_is_base = is_base_2(rsd2, neighb.atomno2());
		bool const atom2_is_sugar = is_sugar_2(rsd2, neighb.atomno2());
		bool const atom2_is_phosphate = is_phosphate_2(rsd2, neighb.atomno2());
		debug_assert( atom2_is_base || atom2_is_sugar || atom2_is_phosphate );

		Real sfxn_weight = rna_elec_weight( weights, atom1_is_base, atom1_is_sugar, atom1_is_phosphate, atom2_is_base, atom2_is_sugar, atom2_is_phosphate );
		if ( sfxn_weight == 0 ) continue;

		Vector f1 = atom1xyz.cross( atom2xyz );
		f1 *= dE_dr_over_r * sfxn_weight;
		f2 *= dE_dr_over_r * sfxn_weight;
		r1_atom_derivs[ neighb.atomno1() ].f1() += f1;
		r1_atom_derivs[ neighb.atomno1() ].f2() += f2;
		r2_atom_derivs[ neighb.atomno2() ].f1() -= f1;
		r2_atom_derivs[ neighb.atomno2() ].f2() -= f2;
	}
}

void
RNA_FA_ElecEnergy::eval_atom_derivative_RNA(
	conformation::Residue const & rsd1,
	Size const i,
	conformation::Residue const & rsd2,
	Size const j,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	using namespace etable::count_pair;

	CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	Real const i_charge( rsd1.atomic_charge( i ) );
	Vector const & i_xyz( rsd1.xyz(i) );

	// NOTE THAT is_base, is_sugar, and is_phosphate contains MAGIC NUMBERS (for speed!)
	// Might be better to make it precomputed as part of the residue_type definition?
	//This repeats some unnecessary stuff, but hey, derivatives don't have to be that fast.
	bool const atom1_is_base = is_base_2(rsd1, i);
	bool const atom1_is_sugar = is_sugar_2(rsd1, i);
	bool const atom1_is_phosphate = is_phosphate_2(rsd1, i);
	debug_assert( atom1_is_base || atom1_is_sugar || atom1_is_phosphate );

	Real const j_charge( rsd2.atomic_charge(j) );
	if ( j_charge == 0.0 ) return;
	Real weight(1.0);
	Size path_dist( 0 );
	if ( ! cpfxn->count( i, j, weight, path_dist ) ) return;

	Vector const & j_xyz( rsd2.xyz(j) );
	Vector const f2( i_xyz - j_xyz );
	Real const dis2( f2.length_squared() );
	Real const dE_dr_over_r = weight *
		coulomb().eval_dfa_elecE_dr_over_r( dis2, i_charge, j_charge );
	if ( dE_dr_over_r == 0.0 ) return;

	Vector const f1( i_xyz.cross( j_xyz ) );

	bool const atom2_is_base = is_base_2(rsd2, j);
	bool const atom2_is_sugar = is_sugar_2(rsd2, j);
	bool const atom2_is_phosphate = is_phosphate_2(rsd2, j);
	debug_assert( atom2_is_base || atom2_is_sugar || atom2_is_phosphate );

	if ( atom1_is_base && atom2_is_base ) {
		F1 += weights[ fa_elec_rna_base_base ] * dE_dr_over_r * f1;
		F2 += weights[ fa_elec_rna_base_base ] * dE_dr_over_r * f2;
	} else if (  (atom1_is_base && atom2_is_sugar)  || (atom1_is_sugar && atom2_is_base) ) {
		F1 += weights[ fa_elec_rna_sugr_base ] * dE_dr_over_r * f1;
		F2 += weights[ fa_elec_rna_sugr_base ] * dE_dr_over_r * f2;
	} else if (  (atom1_is_base && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_base) ) {
		F1 += weights[ fa_elec_rna_phos_base ] * dE_dr_over_r * f1;
		F2 += weights[ fa_elec_rna_phos_base ] * dE_dr_over_r * f2;
	} else if (  (atom1_is_sugar && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_sugar) ) {
		F1 += weights[ fa_elec_rna_phos_sugr ] * dE_dr_over_r * f1;
		F2 += weights[ fa_elec_rna_phos_sugr ] * dE_dr_over_r * f2;
	} else if (  (atom1_is_sugar && atom2_is_sugar)  ) {
		F1 += weights[ fa_elec_rna_sugr_sugr ] * dE_dr_over_r * f1;
		F2 += weights[ fa_elec_rna_sugr_sugr ] * dE_dr_over_r * f2;
	} else if (  (atom1_is_phosphate && atom2_is_phosphate)  ) {
		F1 += weights[ fa_elec_rna_phos_phos ] * dE_dr_over_r * f1;
		F2 += weights[ fa_elec_rna_phos_phos ] * dE_dr_over_r * f2;
	}
}


void
RNA_FA_ElecEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	// Only use when autoupdating
	if ( ! pose.energies().use_nblist_auto_update() ) return;

	using namespace etable::count_pair;

	// what is my charge?
	Size const pos1( atom_id.rsd() );
	Size const i   ( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( pos1 ) );

	if ( ! rsd1.is_RNA() ) return;
	Real const i_charge( rsd1.atomic_charge( i ) );
	if ( i_charge == 0.0 ) return;

	int const pos1_map( domain_map( pos1 ) );
	bool const pos1_fixed( pos1_map != 0 );

	debug_assert( pose.energies().use_nblist() );

	// cached energies object
	NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::RNA_ELEC_NBLIST ) );
	AtomNeighbors const & nbrs( nblist.atom_neighbors(pos1,i) );

	for ( scoring::AtomNeighbor const & nbr : nbrs ) {
		Size const pos2( nbr.rsd() );
		Size const j( nbr.atomno() );

		if ( pos1_fixed && pos1_map == domain_map( pos2 ) ) continue; // fixed wrt one another

		conformation::Residue const & rsd2( pose.residue( pos2 ) );
		debug_assert( pos2 != pos1 );

		if ( rsd2.is_RNA() ) {
			eval_atom_derivative_RNA( rsd1, i, rsd2, j, weights, F1, F2 );
		}
	}
}


/// @brief RNA_FA_ElecEnergy2 is context independent; no context graphs required
void
RNA_FA_ElecEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}
core::Size
RNA_FA_ElecEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
