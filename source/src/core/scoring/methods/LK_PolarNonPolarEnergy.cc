// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/LK_PolarNonPolarEnergy.hh
/// @author Rhiju Das


// Unit headers
#include <core/scoring/methods/LK_PolarNonPolarEnergy.hh>
#include <core/scoring/methods/LK_PolarNonPolarEnergyCreator.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>

#include <ObjexxFCL/format.hh>

#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>

using namespace core::chemical::rna;

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.LK_PolarNonPolarEnergy", basic::t_info );

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is largely deprecated now. The 'lk_costheta' functionality was not found to be any better
//  than geom_sol for RNA applications, and other solvation options like lk_ball and occ_sol are
//  more likely to continue development in Rosetta.
//
// Was mainly in use for lk_nonpolar, but that is more efficiently computed through normal fa_sol/etable machinery
//  using the NO_LK_POLAR_DESOLVATION flag. I checked that we get the same energies
//  when using TableLookupEvaluator etables. The numbers are different by ~0.1% for AnalyticEtableEvaluator
//  because of the way the spline fits to the separate solvation components (atom1 on atom2; atom2 to atom1) do
//  not perfectly add up to a spline fit on their sum. The former is the way we calculate below,
//  and the latter is the way used in the normal fa_sol/etable machinery.
//
// Was also in use in parin's SWA runs for lk_nonpolar_intra_RNA, but I am replacing that, again, with standard
//  fa_sol calculations with PUT_INTRA_INTO_TOTAL score terms.
//
// Still, do NOT remove for a while -- lk_nonpolar has been in wide use in RNA code, and this is still
//  the way to play with rebalancing fa_sol & lk_nonpolar separately. Revisit & consider removal in early 2016.
//
//   -- rhiju, 2014
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace methods {

using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the LK_PolarNonPolarEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
LK_PolarNonPolarEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return EnergyMethodOP( new LK_PolarNonPolarEnergy( *(ScoringManager::get_instance()->etable( options ).lock()),
		options.analytic_etable_evaluation() ) );
}

ScoreTypes
LK_PolarNonPolarEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( lk_costheta );
	sts.push_back( lk_polar );
	sts.push_back( lk_nonpolar );
	sts.push_back( lk_polar_intra_RNA );
	sts.push_back( lk_nonpolar_intra_RNA );
	return sts;
}


LK_PolarNonPolarEnergy::LK_PolarNonPolarEnergy( etable::Etable const & etable_in, bool const analytic_etable_evaluation ):
	parent( methods::EnergyMethodCreatorOP( new LK_PolarNonPolarEnergyCreator ) ),
	safe_max_dis2_( etable_in.get_safe_max_dis2() ),
	max_dis_( etable_in.max_dis() ),
	verbose_( false )
{
	if ( analytic_etable_evaluation ) {
		etable_evaluator_ = etable::EtableEvaluatorOP( new etable::AnalyticEtableEvaluator( etable_in ) );
	} else {
		etable_evaluator_ = etable::EtableEvaluatorOP( new etable::TableLookupEvaluator( etable_in ) );
	}
}

////////////////////////////////////////////////
LK_PolarNonPolarEnergy::LK_PolarNonPolarEnergy( LK_PolarNonPolarEnergy const & src ):
	parent( src ),
	etable_evaluator_( src.etable_evaluator_ ),
	safe_max_dis2_( src.safe_max_dis2_ ),
	max_dis_( src.max_dis_ ),
	verbose_( src.verbose_ )
{}


void
LK_PolarNonPolarEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const{

	if ( rsd.is_RNA()==false ) return;

	Real lk_polar_intra_RNA_score, lk_nonpolar_intra_RNA_score, lk_costheta_intra_RNA_score;
	bool const compute_polar    =  scorefxn.has_nonzero_weight( lk_polar ) || scorefxn.has_nonzero_weight( lk_costheta );
	bool const compute_nonpolar =  scorefxn.has_nonzero_weight( lk_nonpolar );

	get_residue_energy_RNA_intra( rsd, pose, lk_polar_intra_RNA_score, lk_nonpolar_intra_RNA_score, lk_costheta_intra_RNA_score, compute_polar, compute_nonpolar );
	emap[ lk_polar_intra_RNA ]    += lk_polar_intra_RNA_score;
	emap[ lk_nonpolar_intra_RNA ] += lk_nonpolar_intra_RNA_score;

}

bool
LK_PolarNonPolarEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	bool method_1 = ( weights[lk_polar_intra_RNA] > 0.0 || weights[lk_nonpolar_intra_RNA] > 0.0 ) ? true : false;

	return method_1;
}


Distance
LK_PolarNonPolarEnergy::atomic_interaction_cutoff() const
{
	return max_dis_;
}

/// clone
EnergyMethodOP
LK_PolarNonPolarEnergy::clone() const
{
	return EnergyMethodOP( new LK_PolarNonPolarEnergy( *this ) );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const
{
	//TR << "residue_pair_energy() was called..." << std::endl;
	//if ( pose.energies().use_nblist() ) return;
	Real lk_polar_score, lk_nonpolar_score, lk_costheta_score;

	bool const compute_polar    =  scorefxn.has_nonzero_weight( lk_polar ) || scorefxn.has_nonzero_weight( lk_costheta );
	bool const compute_nonpolar =  scorefxn.has_nonzero_weight( lk_nonpolar );

	get_residue_pair_energy_one_way( rsd1, rsd2, pose, lk_polar_score, lk_nonpolar_score, lk_costheta_score,
		compute_polar, compute_nonpolar );
	emap[ lk_polar ]    += lk_polar_score;
	emap[ lk_nonpolar ] += lk_nonpolar_score;
	emap[ lk_costheta ] += lk_costheta_score;

	get_residue_pair_energy_one_way( rsd2, rsd1, pose, lk_polar_score, lk_nonpolar_score, lk_costheta_score,
		compute_polar, compute_nonpolar );
	emap[ lk_polar ]    += lk_polar_score;
	emap[ lk_nonpolar ] += lk_nonpolar_score;
	emap[ lk_costheta ] += lk_costheta_score;

}

////////////////////////////////////////////////
Vector
LK_PolarNonPolarEnergy::get_base_vector( conformation::Residue const & rsd1, Size const i, pose::Pose const & pose ) const
{
	// Use DB/APL prescription of looking through bonded neighbors.
	Size non_H_neighbors = 0;
	Vector  base_pseudo_atom(0);
	for ( Size ii = 1; ii <=rsd1.bonded_neighbor(i).size(); ++ii ) {
		Size neighbor_id = rsd1.bonded_neighbor(i)[ii];
		if ( !  rsd1.atom_is_hydrogen(neighbor_id) ) {
			base_pseudo_atom += rsd1.xyz(neighbor_id);
			non_H_neighbors++;
		}
	}

	if ( rsd1.type().n_residue_connections_for_atom(i) > 0  ) {
		/// CONTEXT DEPENDENCY HERE -- e.g. if c_prev moves, rsd1 needs to be rescored.  Fortunately,
		/// if c_prev moves, the internal "psi" for rsd1 will be updated, and this residue will be rescored.

		for ( Size ii = 1; ii <= rsd1.type().residue_connections_for_atom(i).size(); ++ii ) {
			chemical::ResConnID const ii_conn = rsd1.connect_map( rsd1.type().residue_connections_for_atom(i)[ ii ] );
			Size const neighbor_res_id(  ii_conn.resid() );
			if ( neighbor_res_id < 1 ) continue;
			Size const nieghbor_atom_id( pose.residue( ii_conn.resid() ).residue_connection( ii_conn.connid() ).atomno() );
			if ( ! pose.residue( neighbor_res_id ).atom_is_hydrogen( nieghbor_atom_id ) ) {
				base_pseudo_atom += pose.residue( neighbor_res_id ).xyz( nieghbor_atom_id );
				non_H_neighbors++;
			}
		}
	}

	if ( non_H_neighbors > 0 )  base_pseudo_atom /= non_H_neighbors;

	// Note -- probably should have a big WARNING show up
	// if we're trying to compute this scorefunction for water or something
	// where there is no base atom.  That's where non_H_neighbors is 0.

	Vector res1_base_vector =  rsd1.xyz(i) - base_pseudo_atom;
	Vector res1_base_vector_norm = (res1_base_vector).normalized();
	return res1_base_vector_norm;
}

////////////////////////////////////////////////
// Why is all this code copied? Argh. -- rhiju
void
LK_PolarNonPolarEnergy::get_residue_energy_RNA_intra(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	Real & lk_polar_intra_RNA_score,
	Real & lk_nonpolar_intra_RNA_score,
	Real & lk_costheta_intra_RNA_score,
	bool const compute_polar,
	bool const compute_nonpolar
) const
{

	using namespace etable::count_pair;

	lk_polar_intra_RNA_score =  0.0;
	lk_nonpolar_intra_RNA_score =  0.0;
	lk_costheta_intra_RNA_score =  0.0;

	conformation::Residue const & rsd1=rsd; //Assume rsd1 contributes polar atoms.
	conformation::Residue const & rsd2=rsd; //Assume rsd2 contributes occluding atoms.

	CountPairFunctionOP cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd, CP_CROSSOVER_4); //intra_res version!

	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {

		Vector const heavy_atom_i( rsd1.xyz( i ) );

		//Just compute cos(theta) for polars.
		bool is_polar( false );
		if ( rsd1.atom_type(i).is_acceptor() || rsd1.atom_type(i).is_donor() ) is_polar = true;

		if ( is_polar  && !compute_polar    ) continue;
		if ( !is_polar && !compute_nonpolar ) continue;

		//Need to figure out "direction"...
		Vector const res1_base_vector_norm = is_polar ? get_base_vector( rsd1, i, pose ) : Vector( 0.0 );

		for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( !cpfxn->count( i, j, cp_weight, path_dist ) ) continue;

			if ( !is_base_phosphate_atom_pair(rsd1, rsd2, i, j) ) continue;

			Vector const heavy_atom_j( rsd2.xyz( j ) );

			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			Real dotprod( 1.0 );
			Real dummy_deriv( 0.0 );
			Real temp_score = cp_weight * eval_lk( rsd1.atom( i ), rsd2.atom( j ), dummy_deriv, false);

			if ( is_polar ) {
				lk_polar_intra_RNA_score += temp_score;

				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Vector const d_ij_norm = d_ij.normalized();
				dotprod = dot( res1_base_vector_norm, d_ij_norm );
				temp_score *= dotprod;
				lk_costheta_intra_RNA_score += temp_score;

				if ( verbose_ && std::abs( temp_score ) > 0.1 ) {
					TR << "Occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( i ) << " covered by " <<
						rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(j) << " ==> " << F(8,5,temp_score) << ' ' << F(8,3,dotprod) << std::endl;
				}

			} else {
				lk_nonpolar_intra_RNA_score += temp_score;
				if ( verbose_ && std::abs( temp_score ) > 0.1 ) {
					TR << "Nonpolar occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( i ) << " covered by " <<
						rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(j) << " ==> " << F(8,5,temp_score) << ' ' << F(8,3,dotprod) << std::endl;
				}
			}

		} // j
	} // i
}

////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::get_residue_pair_energy_one_way(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	Real & lk_polar_score,
	Real & lk_nonpolar_score,
	Real & lk_costheta_score,
	bool const compute_polar,
	bool const compute_nonpolar
) const
{
	using namespace etable::count_pair;
	using basic::options::option;
	using basic::options::OptionKeys::score::lk_polar_without_proline_N;

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	lk_polar_score =  0.0;
	lk_nonpolar_score =  0.0;
	lk_costheta_score =  0.0;

	//Assume rsd1 contributes polar atoms.
	//Assume rsd2 contributes occluding atoms.
	for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {

		Vector const & heavy_atom_i( rsd1.xyz( i ) );

		//Just compute cos(theta) for polars.
		bool const is_acceptor = rsd1.atom_type(i).is_acceptor();
		bool const is_donor = option[lk_polar_without_proline_N] ? // bazzoli
			rsd1.heavyatom_has_polar_hydrogens(i) : rsd1.atom_type(i).is_donor();
		bool const is_polar = is_acceptor || is_donor;

		if ( is_polar  && !compute_polar    ) continue;
		if ( !is_polar && !compute_nonpolar ) continue;

		//Need to figure out "direction"...
		Vector const res1_base_vector_norm = is_polar ? get_base_vector( rsd1, i, pose ) : Vector( 0.0 );

		for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {

			Real cp_weight = 1.0;
			Size path_dist( 0 );
			if ( !cpfxn->count( i, j, cp_weight, path_dist ) ) continue;

			Vector const heavy_atom_j( rsd2.xyz( j ) );

			Vector const d_ij = heavy_atom_j - heavy_atom_i;
			Real const d2 = d_ij.length_squared();

			if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

			Real dotprod( 1.0 );
			Real dummy_deriv( 0.0 );
			Real temp_score = cp_weight * eval_lk( rsd1.atom( i ), rsd2.atom( j ), dummy_deriv, false);

			if ( is_polar ) {
				lk_polar_score += temp_score;

				Vector const d_ij = heavy_atom_j - heavy_atom_i;
				Vector const d_ij_norm = d_ij.normalized();
				dotprod = dot( res1_base_vector_norm, d_ij_norm );
				temp_score *= dotprod;
				lk_costheta_score += temp_score;

				if ( verbose_ && std::abs( temp_score ) > 0.1 ) {
					TR << "Occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( i ) << " covered by " <<
						rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(j) << " ==> " << F(8,5,temp_score) << ' ' << F(8,3,dotprod) << std::endl;
				}

			} else {
				lk_nonpolar_score += temp_score;
				if ( verbose_ && std::abs( temp_score ) > 0.1 ) {
					TR << "Nonpolar occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( i ) << " covered by " <<
						rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(j) << " ==> " << temp_score << ' ' << F(8,3,dotprod) << std::endl;
				}
			}

		} // j
	} // i
}

/////////////////////////////////////////////////////////////////////////////
// derivatives
/////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const &
) const
{
	//TR << "setup_for_derivatives() was called..." << std::endl;
	pose.update_residue_neighbors();
}

///////
void
LK_PolarNonPolarEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn) const
{
	// We need the H-bond set -- well, at least the backbone/backbone h-bonds
	// when computing geometric solvation scores.
	// Since this is probably being computed elsewhere, might make sense
	// to have a "calculated" flag.
	// But, anyway, the geometric sol calcs take way longer than this.
	//std::cout<<"test";
	//TR << "setup_for_scoring() was called..." << std::endl;
	pose.update_residue_neighbors();

	if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::LK_POLARNONPOLAR_NBLIST ) );
		nblist.prepare_for_scoring( pose, scfxn, *this );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//set_nres_mono(pose);
	//TR << "setup_for_minimizing() was called..." << std::endl;

	//    if ( pose.energies().use_nblist() ) {
	//        TR << "Using neighborlist..." << std::endl;
	//    }
	if ( !pose.energies().use_nblist() ) return;

	// stash our nblist inside the pose's energies object
	Energies & energies( pose.energies() );

	// setup the atom-atom nblist
	NeighborListOP nblist;
	Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? option[ run::nblist_autoupdate_narrow ] : 1.5;
	Real const XX = max_dis_ + 2 * tolerated_motion;
	nblist = NeighborListOP( new NeighborList( min_map.domain_map(), XX*XX, XX*XX, XX*XX) );
	if ( pose.energies().use_nblist_auto_update() ) {
		//TR << "Using neighborlist auto-update..." << std::endl;
		nblist->set_auto_update( tolerated_motion );
	}
	// this partially becomes the EtableEnergy classes's responsibility
	nblist->setup( pose, sfxn, *this);
	energies.set_nblist( EnergiesCacheableDataType::LK_POLARNONPOLAR_NBLIST, nblist );
}

///////////////////////////////////////////////////////////////////////////////
bool
LK_PolarNonPolarEnergy::defines_score_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool res_moving_wrt_eachother
) const
{
	if ( rsd1.seqpos() == rsd2.seqpos() ) {
		return false;
	}
	return res_moving_wrt_eachother;
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
LK_PolarNonPolarEnergy::get_count_pair_function(
	Size const res1,
	Size const res2,
	pose::Pose const & pose,
	ScoreFunction const &
) const
{
	using namespace etable::count_pair;
	if ( res1 == res2 ) {
		return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );
	}

	conformation::Residue const & rsd1( pose.residue( res1 ) );
	conformation::Residue const & rsd2( pose.residue( res2 ) );
	return get_count_pair_function( rsd1, rsd2 );
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
LK_PolarNonPolarEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	using namespace etable::count_pair;

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );

	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		return CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	}
	return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairAll ) );
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
LK_PolarNonPolarEnergy::get_intrares_countpair(
	conformation::Residue const & rsd1,
	pose::Pose const &,
	ScoreFunction const &
) const
{
	using namespace etable::count_pair;
	return CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_4);
}

///////////////////////////////////////////////////////////////////////////////
bool
LK_PolarNonPolarEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

////////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData & pair_data
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//TR << "setup_for_minimizing_for_residue_pair() was called..." << std::endl;

	etable::count_pair::CountPairFunctionCOP count_pair = get_count_pair_function( rsd1, rsd2 );

	// update the existing nblist if it's already present in the min_data object
	ResiduePairNeighborListOP nblist( utility::pointer::static_pointer_cast< core::scoring::ResiduePairNeighborList > ( pair_data.get_data( lk_PolarNonPolar_pair_nblist ) ));
	if ( ! nblist ) nblist = ResiduePairNeighborListOP( new ResiduePairNeighborList );

	/// STOLEN CODE!
	Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
	Real const XX2 = std::pow( max_dis_ + 2*tolerated_narrow_nblist_motion, 2 );

	nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair );

	pair_data.set_data( lk_PolarNonPolar_pair_nblist, nblist );
}

////////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const
{
	//TR << "residue_pair_energy_ext() was called..." << std::endl;

	using namespace etable::count_pair;
	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	Real lk_polar_score( 0.0 );
	Real lk_nonpolar_score( 0.0 );
	Real lk_costheta_score ( 0.0 );
	bool const compute_polar    =  scorefxn.has_nonzero_weight( lk_polar ) || scorefxn.has_nonzero_weight( lk_costheta );
	bool const compute_nonpolar =  scorefxn.has_nonzero_weight( lk_nonpolar );


	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( lk_PolarNonPolar_pair_nblist ) ) );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	Size m;
	Size n;
	Real cp_weight;
	Size path_dist;
	Real dotprod;
	Real dummy_deriv;
	Real temp_score;
	std::pair<Real,Real> scores;

	for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
		m = neighbs[ ii ].atomno1();
		if ( m > rsd1.nheavyatoms() ) continue;
		n = neighbs[ ii ].atomno2();
		if ( n > rsd2.nheavyatoms() ) continue;
		cp_weight = 1.0;
		path_dist = 0;
		if ( !cpfxn->count( m, n, cp_weight, path_dist ) ) continue;
		bool const is_polar_m = ( rsd1.atom_type(m).is_acceptor() || rsd1.atom_type(m).is_donor());
		bool const is_polar_n = ( rsd2.atom_type(n).is_acceptor() || rsd2.atom_type(n).is_donor());
		if ( !compute_polar && is_polar_m && is_polar_n ) continue;
		if ( !compute_nonpolar && !is_polar_m && !is_polar_n ) continue;
		//Real const d2 = d_ij.length_squared();
		//if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

		dummy_deriv = 0.0;
		scores = eval_lk_efficient( rsd1.atom( m ), rsd2.atom( n ), dummy_deriv, false );
		//temp_score = cp_weight * eval_lk( rsd1.atom( m ), rsd2.atom( n ), dummy_deriv, false);
		temp_score = cp_weight * scores.first;

		if ( compute_polar && is_polar_m ) {
			Vector const & heavy_atom_m( rsd1.xyz( m ) );
			Vector const & heavy_atom_n( rsd2.xyz( n ) );
			Vector const d_ij = heavy_atom_n - heavy_atom_m;
			Vector const d_ij_norm = d_ij.normalized();
			lk_polar_score += temp_score;
			Vector const res1_base_vector_norm = get_base_vector( rsd1, m, pose );
			dotprod = dot( res1_base_vector_norm, d_ij_norm );
			temp_score *= dotprod;
			lk_costheta_score += temp_score;
		} else if ( !is_polar_m && compute_nonpolar ) {
			lk_nonpolar_score += temp_score;
		}

		temp_score = cp_weight * scores.second;

		if ( is_polar_n && compute_polar ) {
			Vector const & heavy_atom_m( rsd1.xyz( m ) );
			Vector const & heavy_atom_n( rsd2.xyz( n ) );
			Vector const d_ij = heavy_atom_n - heavy_atom_m;
			Vector const d_ij_norm = d_ij.normalized();
			lk_polar_score += temp_score;
			Vector const res2_base_vector_norm = get_base_vector( rsd2, n, pose );
			dotprod = -1.0*dot( res2_base_vector_norm, d_ij_norm );
			temp_score *= dotprod;
			lk_costheta_score += temp_score;
		} else if ( !is_polar_n && compute_nonpolar ) {
			lk_nonpolar_score += temp_score;
		}
	}
	emap[ lk_polar ]    += lk_polar_score;
	emap[ lk_nonpolar ] += lk_nonpolar_score;
	emap[ lk_costheta ] += lk_costheta_score;
}

////////////////////////////////////////////////
Real
LK_PolarNonPolarEnergy::eval_lk(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real & deriv,
	bool const & eval_deriv
) const
{

	Real temp_score( 0.0 ); //, temp_score2( 0.0 );
	deriv = 0.0;

	etable_evaluator_->atom_pair_lk_energy_and_deriv_v( atom1, atom2, temp_score, deriv, eval_deriv );

	return temp_score; // original -- works with nonanalytic.
}

////////////////////////////////////////////////
std::pair<Real,Real>
LK_PolarNonPolarEnergy::eval_lk_efficient(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real & deriv,
	bool const & eval_deriv
) const
{

	Real score_1( 0.0 );
	Real score_2( 0.0 );
	deriv = 0.0;

	etable_evaluator_->atom_pair_lk_energy_and_deriv_v_efficient( atom1, atom2, score_1, score_2, deriv, eval_deriv );

	return std::make_pair(score_1,score_2);
}

//////////////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::eval_atom_derivative_intra_RNA(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & /*domain_map*/,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	bool do_eval_intra_RNA= (weights[lk_polar_intra_RNA]>0.0 || weights[lk_nonpolar_intra_RNA]>0.0) ? true : false;

	if ( do_eval_intra_RNA==false ) return; //early return.

	Size const i( atom_id.rsd() );
	Size const j( atom_id.rsd() );

	conformation::Residue const & rsd1( pose.residue( i ) );
	conformation::Residue const & rsd2( pose.residue( j ) );

	if ( rsd1.is_RNA()==false ) return;
	if ( rsd2.is_RNA()==false ) return; //no effect!

	Size const m( atom_id.atomno() );

	if ( m > rsd1.nheavyatoms() ) return;

	if ( verbose_ ) {
		std::cout << "Start LK_PolarNonPolarEnergy::eval_atom_derivative, intra_res case, res= " << i << " atomno= " << m << "[" << rsd1.atom_name(m) <<  "]" << std::endl;
	}

	//std::cout << "LK_PolarNonPolarEnergy::eval_atom_derivative early return since pos1_fixed && domain_map(i=" << i << ") == domain_map(j=" << j << ")" << std::endl;
	//  return; //Fixed w.r.t. one another.
	//}

	Vector const heavy_atom_i( rsd1.xyz( m ) );

	Real deriv( 0.0 );
	Vector const res1_base_vector_norm = get_base_vector( rsd1, m, pose );

	using namespace etable::count_pair;

	CountPairFunctionOP cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_4); //intra_res version!

	bool atom1_is_polar( false );
	if ( rsd1.atom_type(m).is_acceptor() || rsd1.atom_type(m).is_donor() ) atom1_is_polar = true;

	for ( Size n = 1; n <= rsd2.nheavyatoms(); ++n ) {

		Real cp_weight = 1.0;
		Size path_dist( 0 );
		if ( !cpfxn->count(m, n, cp_weight, path_dist ) ) continue;

		if ( is_base_phosphate_atom_pair(rsd1, rsd2, m, n)==false ) continue;

		Vector const heavy_atom_j( rsd2.xyz( n ) );
		Vector const d_ij = heavy_atom_j - heavy_atom_i;
		Real const d2 = d_ij.length_squared();
		Real const d = std::sqrt( d2 );
		Vector const d_ij_norm = d_ij.normalized();

		if ( ( d2 >= safe_max_dis2_) || ( d2 == Real(0.0) ) ) continue;

		bool atom2_is_polar( false );
		if ( rsd2.atom_type(n).is_acceptor() || rsd2.atom_type(n).is_donor() ) atom2_is_polar = true;

		Vector const res2_base_vector_norm = get_base_vector( rsd2, n, pose );

		Vector f1_fwd( 0.0 ), f2_fwd( 0.0 ), f1_bkd( 0.0 ), f2_bkd( 0.0 );
		Real lk_score1( 0.0 ), lk_score2( 0.0 );
		Real dotprod_fwd( 1.0 ), dotprod_bkd( 1.0 );

		//Forward direction first.
		lk_score1 = cp_weight * eval_lk( rsd1.atom(m), rsd2.atom(n), deriv, true );

		f2_fwd =   -1.0 * cp_weight * deriv * d_ij_norm;
		f1_fwd =   1.0 * cross( f2_fwd, heavy_atom_j );

		if ( atom1_is_polar ) {

			F1 += weights[ lk_polar_intra_RNA ] * f1_fwd;
			F2 += weights[ lk_polar_intra_RNA ] * f2_fwd;

			dotprod_fwd = dot( res1_base_vector_norm, d_ij_norm );

			f2_fwd *= dotprod_fwd;
			f2_fwd -= lk_score1 * ( 1/d ) *  (res1_base_vector_norm  - dotprod_fwd * d_ij_norm );

			f1_fwd =   1.0 * cross( f2_fwd, heavy_atom_j );

			lk_score1 *= dotprod_fwd; //to check later (verbose)
		} else {
			F1 += weights[ lk_nonpolar_intra_RNA  ] * f1_fwd;
			F2 += weights[ lk_nonpolar_intra_RNA  ] * f2_fwd;
		}

		/////////////////////////////////
		// Backwards
		Vector d_ji_norm = -d_ij_norm;

		lk_score2 = cp_weight * eval_lk( rsd2.atom(n), rsd1.atom(m),  deriv, true );

		f2_bkd =   -1.0 * deriv * cp_weight * d_ji_norm;
		f1_bkd =   1.0 * cross( f2_bkd, heavy_atom_i );


		if ( atom2_is_polar ) {

			F1 -= weights[ lk_polar_intra_RNA ] * f1_bkd;
			F2 -= weights[ lk_polar_intra_RNA ] * f2_bkd;

			dotprod_bkd = dot( res2_base_vector_norm, d_ji_norm );

			f2_bkd *= dotprod_bkd;
			f2_bkd -= lk_score2 * ( 1/d ) *  (res2_base_vector_norm  - dotprod_bkd * d_ji_norm );

			f1_bkd =   1.0 * cross( f2_bkd, heavy_atom_i );

			lk_score2 *= dotprod_bkd; //to check later (verbose)
		} else {
			F1 -= weights[ lk_nonpolar_intra_RNA ] * f1_bkd;
			F2 -= weights[ lk_nonpolar_intra_RNA  ] * f2_bkd;
		}

		if ( verbose_  &&  (std::abs( lk_score1 ) > 0.1 || std::abs( lk_score2 ) > 0.1 ) ) {
			std::cout << "Occlusion penalty: " << rsd1.name1() << rsd1.seqpos() << " " << rsd1.atom_name( m ) << " covered by " <<
				rsd2.name1() << rsd2.seqpos() << " " << rsd2.atom_name(n) <<
				" " << F(8,3,lk_score1) << " " << F(8,3,lk_score2) <<
				" ==> " << " DERIV " <<
				F(8,6,f2_fwd( 1 ) ) <<  ' ' << F(8,6,f2_bkd(1) ) <<
				' ' << std::sqrt( d2 ) << " " << cp_weight << " " <<
				F(8,3,dotprod_fwd) << " " << F(8,3,dotprod_bkd) << std::endl;
		}
	}

	if ( verbose_ ) {
		std::cout << "LK_PolarNonPolarEnergy::eval_atom_derivative, intra_res :";
		std::cout << " F1= " << F1[0] << " " << F1[1] << " " << F1[2];
		std::cout << " F2= " << F2[0] << " " << F2[1] << " " << F2[2] << std::endl;
		std::cout << "Finish LK_PolarNonPolarEnergy::eval_atom_derivative, intra_res case, res= " << i << " atomno= " << m << "[" << rsd1.atom_name(m) <<  "]" << std::endl;
	}
}

//////////////////////////////////////////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & scorefxn,// sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( ! pose.energies().use_nblist_auto_update() ) return;
	//TR << "eval_atom_derivative() was called..." << std::endl;
	if ( defines_intrares_energy( weights ) ) {
		eval_atom_derivative_intra_RNA(atom_id, pose, domain_map, weights, F1, F2);
	}

	Size const i( atom_id.rsd() );
	bool const pos1_fixed( domain_map( i ) != 0 );
	Size const m( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( i ) );
	if ( m > rsd1.nheavyatoms() ) return;
	if ( rsd1.is_virtual( m ) )  return;
	Vector const & heavy_atom_m( rsd1.xyz( m ) );
	bool const is_polar_m = ( rsd1.atom_type(m).is_acceptor() || rsd1.atom_type(m).is_donor());

	// Size const nres = pose.total_residue();
	NeighborList const & nblist
		( pose.energies().nblist( EnergiesCacheableDataType::LK_POLARNONPOLAR_NBLIST ) );
	//TR << "checkpoint..." << std::endl;
	AtomNeighbors const & nbrs( nblist.atom_neighbors(i,m) );

	using namespace etable::count_pair;
	CountPairFunctionOP cpfxn;

	bool const compute_polar    =  scorefxn.has_nonzero_weight( lk_polar ) || scorefxn.has_nonzero_weight( lk_costheta );
	bool const compute_nonpolar =  scorefxn.has_nonzero_weight( lk_nonpolar );

	Real cp_weight;
	Size path_dist;
	Real dotprod;
	Real deriv ( 0.0 );
	Real lk_score;
	Vector f1;
	Vector f2;

	for ( scoring::AtomNeighbors::const_iterator it2=nbrs.begin(),
			it2e=nbrs.end(); it2 != it2e; ++it2 ) {
		scoring::AtomNeighbor const & nbr( *it2 );
		Size const j( nbr.rsd() );
		if ( i == j ) continue;
		if ( pos1_fixed && domain_map(i) == domain_map(j) ) continue;
		Size const n( nbr.atomno() );
		conformation::Residue const & rsd2( pose.residue( j ) );

		if ( n > rsd2.nheavyatoms() ) continue;
		if ( rsd2.is_virtual( n ) ) continue;

		cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
		cp_weight = 1.0;
		path_dist = 0;
		if ( !cpfxn->count( m, n, cp_weight, path_dist ) ) continue;

		bool const is_polar_n = ( rsd2.atom_type(n).is_acceptor() || rsd2.atom_type(n).is_donor());
		if ( is_polar_m && is_polar_n  && !compute_polar ) continue;
		if ( !compute_nonpolar && !is_polar_m && !is_polar_n ) continue;

		Vector const & heavy_atom_n( rsd2.xyz( n ) );
		Vector const d_ij = heavy_atom_n - heavy_atom_m;
		Vector const d_ij_norm = d_ij.normalized();

		lk_score = cp_weight * eval_lk( rsd1.atom( m ), rsd2.atom( n ), deriv, true);
		f2 = -1.0 * cp_weight * deriv * d_ij_norm;
		f1 = cross( f2, heavy_atom_n );
		if ( compute_polar && is_polar_m ) {
			Real const d = d_ij.length();
			F1 += weights[ lk_polar ] * f1;
			F2 += weights[ lk_polar ] * f2;
			Vector const res1_base_vector_norm = get_base_vector( rsd1, m, pose );
			dotprod = dot( res1_base_vector_norm, d_ij_norm );
			f2 *= dotprod;
			f2 -= lk_score * (1/d) * (res1_base_vector_norm  - dotprod * d_ij_norm );
			f1 = cross( f2, heavy_atom_n );
			F1 += weights[ lk_costheta ] * f1;
			F2 += weights[ lk_costheta ] * f2;
		} else if ( !is_polar_m && compute_nonpolar ) {
			F1 += weights[ lk_nonpolar ] * f1;
			F2 += weights[ lk_nonpolar ] * f2;
		}

		lk_score = cp_weight * eval_lk( rsd2.atom( n ), rsd1.atom( m ), deriv, true);
		f2 = cp_weight * deriv * d_ij_norm;
		f1 = cross( f2, heavy_atom_m );

		if ( compute_polar && is_polar_n ) {
			Real const d = d_ij.length();
			F1 -= weights[ lk_polar ] * f1;
			F2 -= weights[ lk_polar ] * f2;
			Vector const res2_base_vector_norm = get_base_vector( rsd2, n, pose );
			dotprod = dot( res2_base_vector_norm, -d_ij_norm );
			f2 *= dotprod;
			f2 -= lk_score * ( 1/d ) *  (res2_base_vector_norm  + dotprod * d_ij_norm );
			f1 = cross( f2, heavy_atom_m );
			F1 -= weights[ lk_costheta ] * f1;
			F2 -= weights[ lk_costheta ] * f2;
		} else if ( !is_polar_n && compute_nonpolar ) {
			F1 -= weights[ lk_nonpolar ] * f1;
			F2 -= weights[ lk_nonpolar ] * f2;
		}
	}
}

////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /* context_graphs_required */ ) const
{}

////////////////////////////////////////////////
void
LK_PolarNonPolarEnergy::finalize_total_energy(
	pose::Pose & /*pose*/,
	ScoreFunction const & /*scorefxn*/,
	EnergyMap & /*totals*/
) const
{
	//TR << "finalize_total_energy() was called..." << std::endl;
	//TR << pose.energies().use_nblist() << std::endl;
	return;
}

/////////////////////////////////////////
core::Size
LK_PolarNonPolarEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}


