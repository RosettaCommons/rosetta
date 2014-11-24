// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilot/chrisk/rotamer_analysis.cc
/// @brief

//core library
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/chemical/util.hh>

#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/types.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/electron_density/util.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/tools/make_map.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rtmin.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <basic/Tracer.hh>

//protocols library (Movers)
#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/electron_density/util.hh>

//utilities

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rot_anl.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

//local options
namespace basic{ namespace options{ namespace OptionKeys{
}}}//basic::options::OptionKeys

////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace basic;
using namespace id;
using namespace pose;
using namespace pack;
using namespace conformation;
using namespace scoring;
using namespace options;
using namespace basic::options::OptionKeys;
using namespace optimization;
namespace OK = OptionKeys;
using utility::vector1;
using std::string;
using import_pose::pose_from_pdb;
using io::pdb::dump_pdb; // deprecated though
using namespace ObjexxFCL;
using basic::T;
using basic::Warning;
using basic::Error;

static thread_local basic::Tracer TRnat( "rot_anl.nat" );
static thread_local basic::Tracer TRmin( "rot_anl.min" );
static thread_local basic::Tracer TRrtmin( "rot_anl.rtmin" );
static thread_local basic::Tracer TRscmove( "rot_anl.scmove" );

//Map of the atoms that are used in each motif
std::map < std::string, utility::vector1< std::string > > sc_rmsd_AtomIDs(
	utility::tools::make_map(
	std::string(" DA"), utility::tools::make_vector1( std::string("N6"), std::string("C5"), std::string("N7") ),
	std::string(" DC"), utility::tools::make_vector1( std::string("N4"), std::string("C4"), std::string("C5") ),
	std::string(" DG"), utility::tools::make_vector1( std::string("O6"), std::string("C5"), std::string("N7") ),
	std::string(" DT"), utility::tools::make_vector1( std::string("O4"), std::string("C5"), std::string("C7") ),
	std::string("GLY"), utility::tools::make_vector1( std::string("N"), std::string("CA"), std::string("C") ),
	std::string("ALA"), utility::tools::make_vector1( std::string("CB"), std::string("CA"), std::string("N") ),
	std::string("CYS"), utility::tools::make_vector1( std::string("SG"), std::string("CB"), std::string("CA") ),
	std::string("ASP"), utility::tools::make_vector1( std::string("OD1"), std::string("CG"), std::string("OD2") ),
	std::string("GLU"), utility::tools::make_vector1( std::string("OE1"), std::string("CD"), std::string("OE2") ),
	std::string("PHE"), utility::tools::make_vector1( std::string("CE1"), std::string("CZ"), std::string("CE2") ),
	std::string("HIS"), utility::tools::make_vector1( std::string("ND1"), std::string("CE1"), std::string("NE2") ),
	std::string("ILE"), utility::tools::make_vector1( std::string("CD1"), std::string("CG1"), std::string("CG2") ),
	std::string("LYS"), utility::tools::make_vector1( std::string("NZ"), std::string("CE"), std::string("CD") ),
	std::string("LEU"), utility::tools::make_vector1( std::string("CG"), std::string("CD1"), std::string("CD2") ),
	std::string("MET"), utility::tools::make_vector1( std::string("CG"), std::string("SD"), std::string("CE") ),
	std::string("ASN"), utility::tools::make_vector1( std::string("OD1"), std::string("CG"), std::string("ND2") ),
	std::string("GLN"), utility::tools::make_vector1( std::string("OE1"), std::string("CD"), std::string("NE2") ),
	std::string("PRO"), utility::tools::make_vector1( std::string("CG"), std::string("CD"), std::string("NV") ),
	std::string("ARG"), utility::tools::make_vector1( std::string("NH1"), std::string("CZ"), std::string("NH2") ),
	std::string("SER"), utility::tools::make_vector1( std::string("OG"), std::string("CB"), std::string("CA") ),
	std::string("THR"), utility::tools::make_vector1( std::string("OG1"), std::string("CB"), std::string("CG2") ),
	std::string("VAL"), utility::tools::make_vector1( std::string("CG1"), std::string("CB"), std::string("CG2") ),
	std::string("TRP"), utility::tools::make_vector1( std::string("NE1"), std::string("CZ2"), std::string("CE3") ),
	std::string("TYR"), utility::tools::make_vector1( std::string("OH"), std::string("CE1"), std::string("CE2") )
	) );

//set occupancy and bfactor data from native pose into pose
void
set_pose_occ_and_bfac(
	Pose & pose,
	Pose const native_pose
)
{
	for( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ){
		if( pose.residue( seqpos ).name3().compare( native_pose.residue( seqpos ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
		Residue rsd( pose.residue( seqpos ) );
		for( Size ii = 1; ii <= rsd.natoms(); ++ii ){
			if( rsd.atom_is_hydrogen( ii ) ) continue;
			pose.pdb_info()->occupancy( seqpos, ii, native_pose.pdb_info()->occupancy( seqpos, ii ) );
			pose.pdb_info()->temperature( seqpos, ii, native_pose.pdb_info()->temperature( seqpos, ii ) );
		}
	}
}

//get lk energy of one atom by a given residue
//graciously stolen from pbradley
Real
get_atom_lk_energy_by_residue_no_count_pair(
		scoring::ScoreFunctionOP scorefxn,
		conformation::Residue const & rsd1,
		Size const atom1,
		conformation::Residue const & rsd2
		)
{
//  solv1             (  ScoringManager::get_instance()->etable( options.etable_type() )->solv1()),
//  safe_max_dis2     (  ScoringManager::get_instance()->etable( options.etable_type() )->get_safe_max_dis2() ),
//  etable_bins_per_A2(  ScoringManager::get_instance()->etable( options.etable_type() )->get_bins_per_A2() ),

	ObjexxFCL::FArray3D< Real > const & solv1( ScoringManager::get_instance()->etable(
		scorefxn->energy_method_options().etable_type() )->solv1() );
	Real const safe_max_dis2( ScoringManager::get_instance()->etable(
		scorefxn->energy_method_options().etable_type() )->get_safe_max_dis2() );
	Real const etable_bins_per_A2( ScoringManager::get_instance()->etable(
		scorefxn->energy_method_options().etable_type() )->get_bins_per_A2() );

	// setup residue information
	Vector const & atom1_xyz( rsd1.xyz( atom1 ) );
	Size const atom1_type_index( rsd1.atom( atom1 ).type() );

	Real total_lk_energy( 0.0 );
	for ( Size atom2=1; atom2<= rsd2.nheavyatoms(); ++atom2 ) {
		Vector const & atom2_xyz( rsd2.xyz( atom2 ) );

		Real const d2( atom1_xyz.distance_squared( atom2_xyz ) );

		if ( ( d2 >= safe_max_dis2) || ( d2 < 1e-3 ) ) continue; // exclude self...

		// setup for solvation Etable lookups
		Size const atom2_type_index( rsd2.atom( atom2 ).type() );
		Real const d2_bin = d2 * etable_bins_per_A2;
		int disbin = static_cast< int >( d2_bin ) + 1;
		Real  frac = d2_bin - ( disbin - 1 );
		int const l1 = solv1.index( disbin, atom2_type_index, atom1_type_index );

		Real const lk_energy_of_atom1_by_atom2
			( ( ( 1. - frac ) * solv1[ l1 ] + frac * solv1[ l1+1 ] ) );

		total_lk_energy += lk_energy_of_atom1_by_atom2;
//		std::cout << rsd1.name3() + " " + rsd1.atom_name( atom1 ) + " <- " + rsd2.name3() + " " + rsd2.atom_name( atom2 ) + "\t" + string_of( lk_energy_of_atom1_by_atom2 ) + "\n";

	}
	return total_lk_energy;
}

/*
//total LK burial for single atom, sums over self residue atoms and energy graph nbr atoms
Real
calc_lk_burial_for_single_atom(
		Size const atom1,
		conformation::Residue const & rsd1,
		pose::Pose const & pose
		)
{


	PROF_START( util::CALC_LK_BURIAL_FOR_SINGLE_ATOM );

	/// this could be bad if atom type sets don't match up
	static chemical::ResidueTypeSet const * fa_standard_rsd_set
		( & ( *chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) ) );

	static methods::EnergyMethodOptions energy_method_options;

	static methods::LK_BallEnergy lk_ball_energy( energy_method_options );


	if ( & ( rsd1.residue_type_set() ) != fa_standard_rsd_set ) {
		TR.Trace << "skipping lk desolvation calculation for non-fa-standard residue" << std::endl;
		return 0.0;
	}

	/// what is this atom's desolvation parameter?

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	Size const pos1( rsd1.seqpos() );

	Real weighted_desolvation_no_count_pair( 0.0 );
	Real const lk_dgfree( -1.0 * rsd1.atom_type( atom1 ).lk_dgfree() );
	for ( graph::Graph::EdgeListConstIter
			ir  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			ire = energy_graph.get_node( pos1 )->const_edge_list_end();
			ir != ire; ++ir ) {
		EnergyEdge const * edge( static_cast< EnergyEdge const *> (*ir) );
		Size const pos2( edge->get_other_ind( pos1 ) );
		conformation::Residue const & rsd2( pose.residue( pos2 ) );
		if ( rsd2.aa() == chemical::aa_h2o || rsd2.aa() == chemical::aa_vrt ) continue;
		assert( pos2 != pos1 );
		weighted_desolvation_no_count_pair +=
			lk_ball_energy.calculate_lk_desolvation_of_single_atom_by_residue_no_count_pair( atom1, rsd1, rsd2 );
	}

	/// add something here to check for water rotamers built where there was a virtual residue -- will have no neighbors!
	if ( energy_graph.get_node( pos1 )->const_edge_list_begin() ==
			energy_graph.get_node( pos1 )->const_edge_list_end() ) {
		TR.Trace << "calc_lk_desolvation_for_single_atom: no nbrs!" << std::endl;
	}

	weighted_desolvation_no_count_pair +=
		lk_ball_energy.calculate_lk_desolvation_of_single_atom_by_residue_no_count_pair( atom1, rsd1, rsd1 );

	PROF_STOP( util::CALC_LK_BURIAL_FOR_SINGLE_ATOM );

	return weighted_desolvation_no_count_pair / lk_dgfree;

}
*/

Real
get_atom_lk_energy(
	pose::Pose const pose,
	scoring::ScoreFunctionOP scorefxn,
	Size const seqpos,
	Size const iatom
)
{
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	Real fa_sol_atom = 0.;
	core::conformation::Residue const & resl( pose.residue( seqpos ) );
	core::conformation::Atom const & atoml( resl.atom(iatom) );
	scoring::etable::TableLookupEtableEnergy const etable_energy(
		*ScoringManager::get_instance()->etable( scorefxn->energy_method_options().etable_type() ),
		scorefxn->energy_method_options()
	);
	//over all other residue nbrs in energy graph
	for ( graph::Graph::EdgeListConstIter
			ir  = energy_graph.get_node( seqpos )->const_edge_list_begin(),
			ire = energy_graph.get_node( seqpos )->const_edge_list_end();
			ir != ire; ++ir ) {
		EnergyEdge const * edge( static_cast< EnergyEdge const *> (*ir) );
		Size const j( edge->get_other_ind( seqpos ) );
		core::conformation::Residue const & resu( pose.residue( j ) );

		/*
		//this just gets symmetrized energy!!
		//over all atoms in nbr residue
		for ( Size jatom=1, jatom_end = resu.nheavyatoms(); jatom <= jatom_end; ++jatom ) {
			core::conformation::Atom const & atomu( resu.atom(jatom) );
			core::DistanceSquared dsq;
			core::scoring::Weight weight = 1.0;
			core::Size path_dist(0);
			core::Energy atr, rep, solv, bb;
			dsq = atoml.xyz().distance_squared( atomu.xyz() );

			core::scoring::etable::count_pair::CountPairFunctionOP cpfxn = core::scoring::etable::count_pair::CountPairFactory::create_count_pair_function( resl, resu, scoring::etable::count_pair::CP_CROSSOVER_4 );
			if ( cpfxn->count( iatom, jatom, weight, path_dist ) ) {
				etable_energy.atom_pair_energy(atoml,atomu,weight,atr,rep,solv,bb,dsq);
				fa_sol_atom += solv;
			}
		}
		*/
		fa_sol_atom += get_atom_lk_energy_by_residue_no_count_pair( scorefxn, resl, iatom, resu );
	}

	/*
	//this just gets symmetrized energy!!
	//over all other atoms in self residue
	for ( Size jatom=1, jatom_end = resl.nheavyatoms(); jatom <= jatom_end; ++jatom ) {
		//skip self
		if( jatom == iatom ) continue;
		core::conformation::Atom const & atomu( resl.atom(jatom) );
		core::DistanceSquared dsq;
		core::scoring::Weight weight = 1.0;
		core::Size path_dist(0);
		core::Energy atr, rep, solv, bb;
		dsq = atoml.xyz().distance_squared( atomu.xyz() );

		core::scoring::etable::count_pair::CountPairFunctionOP cpfxn = core::scoring::etable::count_pair::CountPairFactory::create_count_pair_function( resl, resl, scoring::etable::count_pair::CP_CROSSOVER_4 );
		if ( cpfxn->count( iatom, jatom, weight, path_dist ) ) {
			etable_energy.atom_pair_energy(atoml,atomu,weight,atr,rep,solv,bb,dsq);
			fa_sol_atom += solv;
		}
	}
	*/
	fa_sol_atom += get_atom_lk_energy_by_residue_no_count_pair( scorefxn, resl, iatom, resl );

	return fa_sol_atom;
}

Real
get_atom_lk_burial(
	pose::Pose const pose,
	scoring::ScoreFunctionOP scorefxn,
	Size const seqpos,
	Size const iatom
)
{
	Real const atom_lk_energy( get_atom_lk_energy( pose, scorefxn, seqpos, iatom ) );
	Real const lk_dgfree( pose.residue( seqpos ).atom_type( iatom ).lk_dgfree() );
	//HACKATTACK! carbonyl carbons have 0 dgfree! returns bogus -1 if dgfree == 0
	if( lk_dgfree == 0 ) return -1;
	Real atom_lk_burial( atom_lk_energy / ( -1 * lk_dgfree ) );
//	std::cout << pose.residue( seqpos ).atom_name( iatom ) + "\t" << atom_lk_burial << "\t" + string_of( lk_dgfree ) + "\n";
	return atom_lk_burial;
}

//add up atomic lk burials, normalize by nheavyatoms
Real
get_res_avg_lk_burial(
	pose::Pose const pose,
	scoring::ScoreFunctionOP scorefxn,
	Size const seqpos,
	bool incl_bb,
	bool incl_sc
)
{
	Real res_lk_burial( 0. );
	Size n_atoms( 0 );
	//HACKATTACK! carbonyl carbons have 0 dgfree! skip val if is NAN
	for( Size iatom = 1; iatom <= pose.residue( seqpos ).nheavyatoms(); ++iatom ){
		//exclude bb or sc atoms?
		if( !incl_bb && iatom < pose.residue( seqpos ).first_sidechain_atom() ) continue;
		if( !incl_sc && iatom >= pose.residue( seqpos ).first_sidechain_atom() ) continue;
		Real atom_lk_burial( get_atom_lk_burial( pose, scorefxn, seqpos, iatom ) );
		//check for bogus neg value
		if( atom_lk_burial < 0 ) continue;
		res_lk_burial += get_atom_lk_burial( pose, scorefxn, seqpos, iatom );
		++n_atoms;
	}
	//return 0 if skipped all atoms
	if( n_atoms == 0 ) return 0;
	return res_lk_burial / n_atoms;
}

//get sum of sr and lr hbond_bb by res
Real
get_res_hbond_bb_score_raw(
	pose::Pose const pose,
	Size const seqpos
)
{
	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	HBondSet hbond_bb_set;
	//exclude bbsc, scbb, and scsc
	fill_hbond_set( pose, false, hbond_bb_set, false, true, true, true );
	Real rsd_hbond_bb_energy( 0.0 );
	for( Size i_hb = 1; i_hb <= Size ( hbond_bb_set.nhbonds() ); ++i_hb ){
			HBond hb( hbond_bb_set.hbond( i_hb ) );
			if ( hb.don_res() == seqpos || hb.acc_res() == seqpos ) {
					//divide by 2 for half of interacting pair
					rsd_hbond_bb_energy += ( hb.energy() / 2.0 );
			}
	}
	return rsd_hbond_bb_energy;
}

//get # sidechain hbonds
Size
get_n_bb_hbonds(
		Pose const pose,
		Size const seqpos
		)
{
	using namespace core::scoring::hbonds;
	Size n_hbonds( 0 );
	HBondSet hbset;
	hbset.setup_for_residue_pair_energies( pose, false, false );
	for( Size i = 1; i <= hbset.nhbonds(); ++i ){
		HBond hb( hbset.hbond( i ) );
		if( hb.don_res() == seqpos && hb.don_hatm_is_protein_backbone() ) ++n_hbonds;
		else if( hb.acc_res() == seqpos && hb.acc_atm_is_protein_backbone() ) ++n_hbonds;
	}
	return n_hbonds;
}

//get sidechain avg bfactor of seqpos, ignore hydrogens
//get # sidechain hbonds
Size
get_n_sc_hbonds(
		Pose const pose,
		Size const seqpos
		)
{
	using namespace core::scoring::hbonds;
	Size n_hbonds( 0 );
	HBondSet hbset;
	hbset.setup_for_residue_pair_energies( pose, false, false );
	for( Size i = 1; i <= hbset.nhbonds(); ++i ){
		HBond hb( hbset.hbond( i ) );
		if( hb.don_res() == seqpos && !hb.don_hatm_is_protein_backbone() ) ++n_hbonds;
		else if( hb.acc_res() == seqpos && !hb.acc_atm_is_protein_backbone() ) ++n_hbonds;
	}
	return n_hbonds;
}

//get sidechain avg bfactor of seqpos, ignore hydrogens
Real
get_sc_bfactor(
		Pose const native_pose,
		Size const seqpos
		)
{
	Residue native_rsd( native_pose.residue( seqpos ) );
	Real sc_bfactor( 0.0 );
	Size count( 0 );
	//skip if no sidechain
	if( native_rsd.first_sidechain_atom() > native_rsd.natoms() ) return 0;
	//get bfacs
	for( Size ii = native_rsd.first_sidechain_atom(); ii <= native_rsd.natoms(); ++ii ){
		if( native_rsd.atom_is_hydrogen( ii ) ) continue;
		if( native_pose.pdb_info()->occupancy( seqpos, ii ) >= 0.5 ){
			sc_bfactor += native_pose.pdb_info()->temperature( seqpos, ii );
			++count;
		}
	}
	sc_bfactor /= count;
	return sc_bfactor;
}

//get sidechain rmsd of seqpos, ignore hydrogens and low occ
Real
get_sc_rmsd(
		Pose const pose,
		Pose const native_pose,
		Size const seqpos
		)
{
	if( pose.residue( seqpos ).name3().compare( native_pose.residue( seqpos ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
	Residue rsd( pose.residue( seqpos ) );
	Residue native_rsd( native_pose.residue( seqpos ) );
	Real sc_rmsd( 0.0 );
	Size count( 0 );
	for( Size ii = rsd.first_sidechain_atom(); ii <= rsd.natoms(); ++ii ){
		if( rsd.atom_is_hydrogen( ii ) ) continue;
		if( native_pose.pdb_info()->occupancy( seqpos, ii ) >= 0.5 ){
			sc_rmsd += rsd.xyz( ii ).distance_squared( native_rsd.xyz( ii ) );
			++count;
		}
	}
	sc_rmsd = std::sqrt( sc_rmsd / count );
	return sc_rmsd;
}

//get automorphic rmsd
Real
get_sc_automorphic_rmsd(
		Pose pose,
		Pose ref_pose,
		Size const seqpos
		)
{
	if( pose.residue( seqpos ).name3().compare( ref_pose.residue( seqpos ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
	//add variant type to virt bb atoms
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_BB", seqpos );
	pose::add_variant_type_to_pose_residue( ref_pose, "VIRTUAL_BB", seqpos );

	Residue rsd( pose.residue( seqpos ) );
	Residue native_rsd( ref_pose.residue( seqpos ) );

	Real sc_rmsd( scoring::automorphic_rmsd( pose.residue( seqpos ), ref_pose.residue( seqpos ), false /*superpose*/ ) );
	return sc_rmsd;
}

//get sidechain rmsd of seqpos, ignore hydrogens and low occ
Real
get_three_atom_sc_rmsd(
		Pose const pose,
		Pose const native_pose,
		Size const seqpos
		)
{
	if( pose.residue( seqpos ).name3().compare( native_pose.residue( seqpos ).name3() ) != 0 ) utility_exit_with_message( "Native residue type mismatch at " + string_of( seqpos ) + "\n" );
	Residue rsd( pose.residue( seqpos ) );
	Residue native_rsd( native_pose.residue( seqpos ) );
	vector1< std::string > atom_ids( sc_rmsd_AtomIDs[ rsd.name3() ] );
	Real sc_rmsd( 0.0 );
	Size count( 0 );
	for( Size ii = 1; ii <= atom_ids.size(); ++ii ){
		sc_rmsd += rsd.xyz( rsd.atom_index( atom_ids[ ii ] ) ).distance_squared( native_rsd.xyz( native_rsd.atom_index( atom_ids[ ii ] ) ) );
		++count;
	}
	sc_rmsd = std::sqrt( sc_rmsd / count );
	return sc_rmsd;
}

void
split_fa_dun(
		Pose const pose,
		Size const seqpos,
		Real & fa_dun_rot,
		Real & fa_dun_dev,
		pack::dunbrack::RotVector & rotvec
		)
{

	using namespace scoring;
	using namespace pack::dunbrack;

	conformation::Residue const & rsd( pose.residue( seqpos ) );
	bool is_dun02( true );
        if( option[ corrections::correct ].user() || option[ corrections::score::dun10 ].user() ) is_dun02 = false;
	RotamerLibrary const & rlib( * RotamerLibrary::get_instance() );
	RotamerLibraryScratchSpace scratch;

	Real fa_dun_tot = rlib.rotamer_energy( rsd, scratch );
	fa_dun_rot = scratch.fa_dun_rot();
	fa_dun_dev = scratch.fa_dun_dev();
	//semirotameric returns 0 for fa_dun_rot!
	if( std::abs( fa_dun_tot - fa_dun_rot - fa_dun_dev ) > 0.01 ) fa_dun_rot = ( fa_dun_tot - fa_dun_dev );

	rotamer_from_chi( rsd, rotvec );
}

std::string
get_aa_torsion_string(
		Pose const pose,
		Size const seqpos
		)
{
	std::string tor_anl;
	tor_anl += ( "phi: " + string_of( pose.phi( seqpos ) ) + " " );
	tor_anl += ( "psi: " + string_of( pose.psi( seqpos ) ) + " " );
	if( pose.residue( seqpos ).nchi() < 1 ){
		tor_anl += ( "chi1: 0 chi2: 0 chi3: 0 chi4: 0 " );
		return tor_anl;
	}
	tor_anl += ( "chi1: " + string_of( pose.chi( 1, seqpos ) ) + " " );
	if( pose.residue( seqpos ).nchi() < 2 ){
		tor_anl += ( "chi2: 0 chi3: 0 chi4: 0 " );
		return tor_anl;
	}
	tor_anl += ( "chi2: " + string_of( pose.chi( 2, seqpos ) ) + " " );
	if( pose.residue( seqpos ).nchi() < 3 ){
		tor_anl += ( "chi3: 0 chi4: 0 " );
		return tor_anl;
	}
	tor_anl += ( "chi3: " + string_of( pose.chi( 3, seqpos ) ) + " " );
	if( pose.residue( seqpos ).nchi() < 4 ){
		tor_anl += ( "chi4: 0 " );
		return tor_anl;
	}
	tor_anl += ( "chi4: " + string_of( pose.chi( 4, seqpos ) ) + " " );
	return tor_anl;
}

std::string
res_lvl_analysis(
		Pose pose,
		Pose const native_pose,
		Size const seqpos,
		ScoreFunctionOP scorefxn,
		ScoreFunctionOP scorefxn_edens,
		bool const do_sc_rmsd
		)
{

	scorefxn->score( pose );

	//extra data
	Real total_score( pose.energies().residue_total_energies( seqpos ).dot( scorefxn->weights() ) );
	Real sc_auto_rmsd_nat( 0 );
	Real sc_auto_rmsd_min( 0 );
	if( do_sc_rmsd ){
		sc_auto_rmsd_nat = get_sc_automorphic_rmsd( pose, native_pose, seqpos );
	}
	Real sc_bfactor( get_sc_bfactor( native_pose, seqpos ) );
	Size n_sc_hbonds( get_n_sc_hbonds( pose, seqpos ) );
	Real res_lk_burial( get_res_avg_lk_burial( pose, scorefxn, seqpos, true, true ) );
	Real edens_score( 0.0 );
	if( option[ edensity::mapfile ].user() ){
		scorefxn_edens->score( pose );
		edens_score = pose.energies().residue_total_energies( seqpos )[ elec_dens_window ];
	}

	//dunbrack analysis
	Real fa_dun_rot, fa_dun_dev;
	pack::dunbrack::RotVector rotvec;
	split_fa_dun( pose, seqpos, fa_dun_rot, fa_dun_dev, rotvec );
	Real fa_dun_wt( scorefxn->get_weight( fa_dun ) );
	fa_dun_rot *= fa_dun_wt;
	fa_dun_dev *= fa_dun_wt;
	std::string rotbin_str( "" );
	for( Size i_rotvec = 1; i_rotvec <= 4; ++i_rotvec ){
		if( i_rotvec <= rotvec.size() ) rotbin_str += ( string_of( rotvec[ i_rotvec ] ) );
		else rotbin_str += "0";
	}

	//output
	std::string score_data( pose.energies().residue_total_energies( seqpos ).weighted_string_of( scorefxn->weights() ) + " total_score: " + string_of( total_score ) );
	std::string extra_data( " n_sc_hbonds: " + string_of( n_sc_hbonds ) + " scauto_rmsd_nat: " + string_of( sc_auto_rmsd_nat ) );
	if( option[ edensity::mapfile ].user() ) extra_data += " elec_dens_window: " + string_of( edens_score );
	extra_data += " sc_bfactor: " + string_of( sc_bfactor ) + " fa_dun_rot: " + string_of( fa_dun_rot ) + " fa_dun_dev: " + string_of( fa_dun_dev ) + " rotbin: " + rotbin_str + " res_lk_burial: " + string_of( res_lk_burial );
	std::string torsion_data( get_aa_torsion_string( pose, seqpos ) );
	return score_data + " " + extra_data + " " + torsion_data;
}

void
get_res_data_ss(
		io::silent::SilentStructOP & ss,
		Pose pose,
		Pose const native_pose,
		Size const seqpos,
		ScoreFunctionOP scorefxn,
		ScoreFunctionOP scorefxn_edens,
		bool const do_sc_rmsd
		)
{
	using namespace core;
	using namespace core::scoring;

	std::map< std::string, Real > res_data_map;

	scorefxn->score( pose );
    EnergyMap weights( pose.energies().weights() );
	EnergyMap rsd_energies( pose.energies().residue_total_energies( seqpos )  );

	Real total_score( pose.energies().residue_total_energies( seqpos ).dot( scorefxn->weights() ) );
	ss->add_energy( "score", total_score );

	for ( int ii = 1; ii <= scoring::n_score_types; ++ii ) {
		if ( weights[ ScoreType(ii) ] != 0.0 ) {
			Real const value( rsd_energies[ ScoreType(ii) ] );
			std::string const scorename( name_from_score_type( ScoreType(ii) ) );
			ss->add_energy( scorename, value );
		}
	}

	//dunbrack split
	Real fa_dun_rot, fa_dun_dev;
	pack::dunbrack::RotVector rotvec;
	split_fa_dun( pose, seqpos, fa_dun_rot, fa_dun_dev, rotvec );
	Real fa_dun_wt( scorefxn->get_weight( fa_dun ) );
	fa_dun_rot *= fa_dun_wt;
	fa_dun_dev *= fa_dun_wt;
	ss->add_energy( "fa_dun_rot", fa_dun_rot );
	ss->add_energy( "fa_dun_dev", fa_dun_dev );

	//res hbond_bb energy
	Real hbond_bb_score_raw( get_res_hbond_bb_score_raw( pose, seqpos ) );
	ss->add_energy( "hbond_bb_raw", hbond_bb_score_raw );

	//avg atomic lk burial in residue
	Real res_lk_burial( get_res_avg_lk_burial( pose, scorefxn, seqpos, true, true ) );
	Real bb_lk_burial( get_res_avg_lk_burial( pose, scorefxn, seqpos, true, false ) );
	Real sc_lk_burial( get_res_avg_lk_burial( pose, scorefxn, seqpos, false, true ) );
	ss->add_energy( "res_lk_burial", res_lk_burial );
	ss->add_energy( "bb_lk_burial", bb_lk_burial );
	ss->add_energy( "sc_lk_burial", sc_lk_burial );

	//electron density res score
	Real edens_score( 0.0 );
	if( option[ edensity::mapfile ].user() ){
		scorefxn_edens->score( pose );
		edens_score = pose.energies().residue_total_energies( seqpos )[ elec_dens_window ];
	}
	ss->add_energy( "elec_dens_window", edens_score );

	//n hbonds
	Size n_sc_hbonds( get_n_sc_hbonds( pose, seqpos ) );
	Size n_bb_hbonds( get_n_bb_hbonds( pose, seqpos ) );
	ss->add_energy( "n_sc_hbonds", static_cast< Real >( n_sc_hbonds ) );
	ss->add_energy( "n_bb_hbonds", static_cast< Real >( n_bb_hbonds ) );

	//things after this point should not be delta'ed

	//this will xform 4 rotbin indices into a 4 digit number
	Real rotbin_val( 0 );
	Size n_rotbins( 4 );
	for( Size i_rotvec = 1; i_rotvec <= n_rotbins; ++i_rotvec ){
		if( i_rotvec <= rotvec.size() ) rotbin_val += ( ( rotvec[ i_rotvec ] + 1 ) * std::pow( 10.0, static_cast< Real >( 4 - i_rotvec ) ) );
		else rotbin_val += std::pow( 10.0, static_cast< Real >( 4 - i_rotvec ) );
	}
	ss->add_energy( "rotbin", static_cast< Real >( rotbin_val ) );

	//torsion angles
	ss->add_energy( "phi", pose.phi( seqpos ) );
	ss->add_energy( "psi", pose.psi( seqpos ) );
	ss->add_energy( "chi1", pose.chi( 1, seqpos ) );
	ss->add_energy( "chi2", pose.chi( 2, seqpos ) );
	ss->add_energy( "chi3", pose.chi( 3, seqpos ) );
	ss->add_energy( "chi4", pose.chi( 4, seqpos ) );

	//sidechain bfactor from native
	Real sc_bfactor( get_sc_bfactor( native_pose, seqpos ) );
	ss->add_energy( "sc_bfactor", sc_bfactor );

	//sidechain rmsd
	Real sc_auto_rmsd_nat( 0 );
	if( do_sc_rmsd ){
		sc_auto_rmsd_nat = get_sc_automorphic_rmsd( pose, native_pose, seqpos );
	}
	ss->add_energy( "sc_rmsd", sc_auto_rmsd_nat );

	Size residx( pose.residue( seqpos ).aa() );
	ss->add_energy( "residx", static_cast< Real >( residx ) );

}

//minimizes sidechains
void
minimize_all_sidechains(
		Pose & pose,
		scoring::ScoreFunctionOP const scorefxn,
		bool const min_dna
		)
{
	using namespace protocols;
	using namespace protocols::moves;

	( *scorefxn )( pose );
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	for( Size i = 1; i <= pose.total_residue(); ++i ){
		if( pose.residue( i ).is_protein() ) mm->set_chi( i, true );
		//min dna?
		else if( min_dna && pose.residue( i ).is_DNA() ){
			mm->set_chi( i, true );
			//and sugar
			mm->set( id::TorsionID( i, id::CHI, 2 ), true );
			mm->set( id::TorsionID( i, id::CHI, 3 ), true );
			mm->set( id::TorsionID( i, id::CHI, 4 ), true );
		}
	}

	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( mm, scorefxn, "dfpmin", 0.001, true ) );
	min_mover->apply( pose );
}

void
minimize_sidechain(
		Pose & pose,
		scoring::ScoreFunctionOP const scorefxn,
		Size const seqpos
		)
{
	using namespace protocols;
	using namespace protocols::moves;

	( *scorefxn )( pose );

	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_chi( false );
	mm->set_chi( seqpos, true );
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( mm, scorefxn, "dfpmin", 0.001, true ) );
	min_mover->apply( pose );

}

//gen packer task for repacking one res
core::pack::task::PackerTaskOP
single_res_task(
		Pose pose,
		Size const seqpos,
		bool const no_incl_curr,
		bool design
		)
{
	//define task
	core::pack::task::PackerTaskOP packer_task( pack::task::TaskFactory::create_packer_task( pose ) );
	packer_task->initialize_from_command_line();
	if( !no_incl_curr ) packer_task->or_include_current( true );
	if( !design ) packer_task->restrict_to_repacking();
	vector1< bool > repack_this( pose.total_residue(), false );
	repack_this[ seqpos ]  = true;
	packer_task->restrict_to_residues( repack_this );
	return packer_task;
}

void
scmove_residue(
		Pose & pose,
		scoring::ScoreFunctionOP const scorefxn,
		Size const seqpos,
		bool const no_incl_curr
		)
{
	using namespace protocols;
	using namespace protocols::moves;

	( *scorefxn )( pose );

	//scmover doesnt handle design (??)
	core::pack::task::PackerTaskOP packer_task( single_res_task( pose, seqpos, no_incl_curr, false ) );

	//create movers and apply
	protocols::simple_moves::sidechain_moves::SidechainMoverOP scmover( new protocols::simple_moves::sidechain_moves::SidechainMover() );
	scmover->set_task( packer_task );
	scmover->init_task( pose );
	scmover->set_prob_uniform( 0.1 );
	scmover->set_prob_withinrot( 0.0 );
	scmover->set_prob_random_pert_current( 0.0 );
	scmover->set_preserve_detailed_balance( false );
	MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, 1.0 ) );

	//and apply in loop
	Size n = option[ rot_anl::nloop_scmove ];
	Size nloop = static_cast< Size >( std::pow( static_cast< Real >( n ), static_cast< Real >( pose.residue( seqpos ).nchi() ) ) );
	for( Size i = 1; i <= nloop; ++i ){
		scmover->apply( pose );
		minimize_sidechain( pose, scorefxn, seqpos );
		mc->boltzmann( pose );
	}
	mc->recover_low( pose );
}

void
rottrial_residue(
		Pose & pose,
		scoring::ScoreFunctionOP const scorefxn,
		Size const seqpos,
		bool const no_incl_curr,
		bool const design
		)
{
	using namespace protocols;
	using namespace protocols::moves;

	( *scorefxn )( pose );

	core::pack::task::PackerTaskOP packer_task( single_res_task( pose, seqpos, no_incl_curr, design ) );

	//create movers
	protocols::simple_moves::RotamerTrialsMoverOP rottrial ( new protocols::simple_moves::RotamerTrialsMover( scorefxn, *packer_task ) );

	//and apply
	rottrial->apply( pose );
}

void
rottrialmin_residue(
		Pose & pose,
		scoring::ScoreFunctionOP const scorefxn,
		Size const seqpos,
		bool const design
		)
{
	using namespace protocols;
	using namespace protocols::moves;

	( *scorefxn )( pose );
	Pose pose_in( pose );
	Real total_score_in( scorefxn->score( pose ) );

	//no_incl_curr is ignored by rtmin, resets to include current
	core::pack::task::PackerTaskOP packer_task( single_res_task( pose, seqpos, false, design ) );

	//create movers
	protocols::simple_moves::RotamerTrialsMinMoverOP rottrialmin ( new protocols::simple_moves::RotamerTrialsMinMover( scorefxn, *packer_task ) );

	//and apply
	rottrialmin->apply( pose );
	//double check score did not incr.
	if( total_score_in < scorefxn->score( pose ) ) pose = pose_in;
}

///////////////////////////////////////////////////////////////////////////////
	void
RotamerAnalysis()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//tag is "." or "tag."
	std::string tag( option[ rot_anl::tag ] );
	if( option[ rot_anl::tag ].user() ) tag = tag + ".";

	//allow design in packer tasks?
	bool const design( option[ rot_anl::design ] );

	std::string pdbname( start_file() );
	Size pdbnamestart( 0 );
	if( pdbname.find_last_of( "/" ) < ( pdbname.size() - 1 ) ) pdbnamestart = pdbname.find_last_of( "/" ) + 1;
	std::string pdbnametag( pdbname, pdbnamestart, pdbname.size() - pdbnamestart - 4 );
	Pose pose;
	pose_from_pdb( pose, pdbname );

	//create a ScoreFunction from commandline options
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::scoring::ScoreFunctionOP scorefxn_edens( ( *scorefxn ).clone() );

	//add edens score?
	if( option[ edensity::mapfile ].user() ){
		protocols::electron_density::SetupForDensityScoringMoverOP edens_mover( new protocols::electron_density::SetupForDensityScoringMover );
		edens_mover->apply( pose );
		std::string const pdbname_out( pdbnametag + ".edock." + "pdb" );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn_edens );
		pose.dump_scored_pdb( pdbname_out, *scorefxn_edens );
	}

	//set non-edens scorefxn edens wt to 0.0 so can rtmin and stuff
	scorefxn->set_weight( elec_dens_window, 0.0 );

	std::string native_pdbname( pdbname );
	Pose native_pose( pose );
	if( option[ in::file::native ].user() ){
		native_pdbname = option[ in::file::native ]();
		pose_from_pdb( native_pose, native_pdbname );
		if( option[ edensity::mapfile ].user() ) core::scoring::calpha_superimpose_pose( native_pose, pose );
	}
	set_pose_occ_and_bfac( pose, native_pose );

	scorefxn->score( pose );
	scorefxn->score( native_pose );

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
	//if premin flag, min and dump
	if( option[ rot_anl::premin ] ){
		minimize_all_sidechains( pose, scorefxn, false );
		std::string const pdbname_out( pdbnametag + "." + tag + "pdb" );
		pose.dump_scored_pdb( pdbname_out, *scorefxn );
	}
*/

	//silent-type output
	io::silent::SilentFileData sfd;

	//create a task factory to get resfile data
	core::pack::task::TaskFactoryOP task_factory = new core::pack::task::TaskFactory;
	//use resfile to note residues to analyze...
	Pose start_pose( pose );
	if( option[ packing::resfile ].user() ) task_factory->push_back( new core::pack::task::operation::ReadResfile );
	else task_factory->push_back( new core::pack::task::operation::RestrictToRepacking() );
	core::pack::task::PackerTaskOP packer_task( task_factory->create_task_and_apply_taskoperations( pose ) );
	//each residue
	for( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ){
		if( !pose.residue( seqpos ).is_protein() ) continue;
		if( packer_task->being_packed( seqpos ) ){
			pose = start_pose;
			std::string reschain( pose.pdb_info()->pose2pdb( seqpos ) );
			std::string pdbname_reschain( reschain );
			if( option[ edensity::mapfile ].user() ) reschain = pdbname_reschain = string_of( seqpos );
			else pdbname_reschain = reschain.replace( reschain.find( " " ), 1, "_" );

			{
			//res analysis
			//			TRnat << pdbname + " " + reschain + " " + pose.residue( seqpos ).name3() + res_lvl_analysis( pose, native_pose, seqpos, scorefxn, scorefxn_edens, !design ) << std::endl;
			//silent file stylee
					io::silent::SilentStructOP ss( new io::silent::ScoreFileSilentStruct );
					std::string const decoytag_out( pdbnametag + "." + tag + string_of( seqpos ) );
					ss->decoy_tag( decoytag_out );
					get_res_data_ss( ss, pose, native_pose, seqpos, scorefxn, scorefxn_edens, !design );
					std::string const scname_out( pdbnametag + "." + tag + "sc" );
					sfd.write_silent_struct( *ss, scname_out );
			}

			//begin perturbations
			if( option[ rot_anl::min ] ){
				//if edens, min+repack "minimized" rotamer
				if( option[ edensity::mapfile ].user() ){
					minimize_sidechain( pose, scorefxn_edens, seqpos );
					rottrial_residue( pose, scorefxn_edens, seqpos, false, design );
					minimize_sidechain( pose, scorefxn_edens, seqpos );
				}
				else	minimize_sidechain( pose, scorefxn, seqpos );

				std::string const pdbname_out( pdbnametag + "." + tag + "min." + pdbname_reschain + ".pdb" );
				//dump_pdb
				if( option[ rot_anl::dump_pdb ] ) pose.dump_scored_pdb( pdbname_out, *scorefxn );
				//res analysis
//				TRmin << pdbname_out + " " + reschain + " " + pose.residue( seqpos ).name3() + res_lvl_analysis( pose, native_pose, seqpos, scorefxn, scorefxn_edens, !design ) << std::endl;
			}
			Pose const min_pose( pose );

			if( option[ rot_anl::repack ] ){
				//if did edens, reset pose so can incl native rotamer!
				if( option[ edensity::mapfile ].user() ) pose = native_pose;
				minimize_sidechain( pose, scorefxn, seqpos );
				rottrial_residue( pose, scorefxn, seqpos, false, design );
				minimize_sidechain( pose, scorefxn, seqpos );
				std::string const pdbname_out( pdbnametag + "." + tag + "rt." + pdbname_reschain + ".pdb" );
				//dump_pdb
				pose.dump_scored_pdb( pdbname_out, *scorefxn );
				//res analysis
//				TRrtmin << pdbname_out + " " + reschain + " " + pose.residue( seqpos ).name3() + res_lvl_analysis( pose, native_pose, seqpos, scorefxn, scorefxn_edens, !design ) << std::endl;
			}

			if( option[ rot_anl::rtmin ] ){
				//if did edens, reset pose so can incl native rotamer!
				if( option[ edensity::mapfile ].user() ) pose = native_pose;
				rottrialmin_residue( pose, scorefxn, seqpos, design );
				std::string const pdbname_out( pdbnametag + "." + tag + "rtmin." + pdbname_reschain + ".pdb" );
				//dump_pdb
				if( option[ rot_anl::dump_pdb ] ) pose.dump_scored_pdb( pdbname_out, *scorefxn );
				//res analysis
//				TRrtmin << pdbname_out + " " + reschain + " " + pose.residue( seqpos ).name3() + res_lvl_analysis( pose, native_pose, seqpos, scorefxn, scorefxn_edens, !design ) << std::endl;
			}

			if( option[ rot_anl::scmove ] ){
				scmove_residue( pose, scorefxn, seqpos, false );
				std::string const pdbname_out( pdbnametag + "." + tag + "scmove." + pdbname_reschain + ".pdb" );
				//dump_pdb
				if( option[ rot_anl::dump_pdb ] ) pose.dump_scored_pdb( pdbname_out, *scorefxn );
				//res analysis
//				TRscmove << pdbname_out + " " + reschain + " " + pose.residue( seqpos ).name3() + res_lvl_analysis( pose, native_pose, seqpos, scorefxn, scorefxn_edens, !design ) << std::endl;
			}


			/*
			//dump struct if passes score, rmsd diff filter
			//TODO: don't calc autormsd twice
			if( min_pose.energies().residue_total_energies( seqpos ).dot( scorefxn->weights() ) - pose.energies().residue_total_energies( seqpos ).dot( scorefxn->weights() ) > option[ rot_anl::score_tol ]
			&& get_sc_automorphic_rmsd( pose, native_pose, seqpos ) > option[ rot_anl::rmsd_tol ] ){
			//dump pdb if given flag
			if( option[ rot_anl::dump_pdb ] ) pose.dump_scored_pdb( pdbname_out, *scorefxn );
			}
			 */



		}
	}
}


	void*
my_main( void*)
{

	RotamerAnalysis();
	exit(0);

}

	int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


