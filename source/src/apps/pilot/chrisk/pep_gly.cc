// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief


// libRosetta headers
 #include <protocols/frags/VallData.hh>
 #include <protocols/frags/TorsionFragment.hh>

 #include <core/scoring/dna/setup.hh>
 #include <core/scoring/dna/base_geometry.hh>
 #include <core/scoring/dna/BasePartner.hh>
 #include <core/scoring/GenBornPotential.hh>
 #include <core/scoring/LREnergyContainer.hh>
 #include <core/scoring/methods/Methods.hh>

 #include <protocols/simple_moves/BackboneMover.hh>
 #include <protocols/simple_moves/MinMover.hh>
 #include <protocols/moves/MonteCarlo.hh>
 #include <protocols/moves/Mover.hh>
 #include <protocols/moves/MoverContainer.hh>
 #include <protocols/moves/OutputMovers.hh>
 #include <protocols/rigid/RigidBodyMover.hh>
 // #include <protocols/moves/rigid_body_moves.hh>
 #include <protocols/moves/TrialMover.hh>
 #include <protocols/simple_moves/PackRotamersMover.hh>
 #include <protocols/simple_moves/RotamerTrialsMover.hh>
 #include <protocols/moves/RepeatMover.hh>

 #include <protocols/loops/ccd_closure.hh>
 #include <protocols/loops/loops_main.hh>

 #include <protocols/viewer/viewers.hh>

 #include <core/types.hh>

 #include <core/scoring/sasa.hh>

// #include <basic/prof.hh> // profiling
// #include <basic/CacheableData.hh> // profiling

 #include <core/id/SequenceMapping.hh>

 #include <core/chemical/AtomTypeSet.hh>
 #include <core/chemical/MMAtomTypeSet.hh>

 #include <core/chemical/AA.hh>
 #include <core/conformation/Residue.hh>
 #include <core/conformation/ResidueMatcher.hh>
 #include <core/pack/rotamer_set/RotamerCouplings.hh>
 #include <core/chemical/ResidueTypeSet.hh>
 #include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
 #include <core/chemical/VariantType.hh>

 #include <core/chemical/ChemicalManager.hh>

 #include <core/scoring/etable/Etable.hh>
 #include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
 #include <core/scoring/Ramachandran.hh>
 #include <core/pack/dunbrack/RotamerLibrary.hh>
 #include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
 #include <core/scoring/hbonds/HBondSet.hh>
 #include <core/scoring/hbonds/hbonds.hh>
 #include <core/scoring/etable/count_pair/CountPairFunction.hh>

 #include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/TaskOperation.hh>

 #include <core/kinematics/FoldTree.hh>
 #include <protocols/viewer/visualize.hh>
#include <core/kinematics/MoveMap.hh>
 #include <core/kinematics/util.hh>
 #include <core/id/AtomID_Map.hh>
// #include <core/id/AtomID_Map.Pose.hh>

 #include <core/mm/MMTorsionLibrary.hh>
 #include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/util.hh>//option.hh>
// #include <basic/options/after_opts.hh>

 #include <basic/basic.hh>

 #include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

 #include <numeric/xyzVector.hh>
 #include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// //REMOVE LATER!
// #include <utility/io/izstream.hh>

#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/PeriodicFunc.hh>
#include <core/scoring/constraints/DOF_Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/fragment/picking/FragmentLibraryManager.hh>
#include <core/fragment/picking/vall/VallLibrarian.hh>
#include <core/fragment/picking/vall/gen/LengthGen.hh>
#include <core/fragment/picking/vall/scores/VallFragmentScore.hh>

// // C++ headers
 #include <cstdlib>
 #include <fstream>
 #include <iostream>
 #include <string>

//silly using/typedef
#include <basic/Tracer.hh>


#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pep_spec.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>



 using basic::T;
 using basic::Error;
 using basic::Warning;

 static numeric::random::RandomGenerator RG(16621);

using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
//using namespace protocols;

using utility::vector1;
using std::string;
using io::pdb::dump_pdb;

using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace pose;
using namespace chemical;
using namespace conformation;
using namespace scoring;
using namespace optimization;
using namespace kinematics;
using namespace protocols::moves;
using namespace id;
using namespace protocols::frags;

Size
get_n_pep_nbrs(
	pose::Pose const & pose,
	vector1< bool > const is_pep,
	Real const cutoff_cg
)
{
	Size n_pep_nbrs( 0 );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !is_pep[i] ) continue;
		bool cg_res_has_nbr( false );
		Residue const & rsd1( pose.residue(i) );
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			if( cg_res_has_nbr ) break;
			if ( is_pep[j] ) continue;
			Residue const & rsd2( pose.residue(j) );
			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				if( cg_res_has_nbr ) break;
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
					if ( rsd1.xyz(ii).distance( rsd2.xyz(jj) ) < cutoff_cg ) {
						cg_res_has_nbr = true;
						break;
					}
				}
			}
		}
		if( cg_res_has_nbr ) ++n_pep_nbrs;
	}
	return n_pep_nbrs;
}



void
dump_efactor_pdb(
	pose::Pose & pose,
	core::scoring::ScoreFunctionOP const scorefxn,
	std::string const & tag
)
{
	Size const nres( pose.total_residue() );
//	id::AtomID_Mask const & mask;
//	id::initialize( mask, pose );

	( *scorefxn )( pose );

	char const *filename = tag.c_str();
	std::fstream out( filename, std::ios::out );


	Size number(0);

	static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );

	out << "MODEL     " << tag << "\n";
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		Real residue_total_energy( pose.energies().residue_total_energies( i ).dot( scorefxn->weights() ) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );

 //			if ( ! mask[ id::AtomID( j,i ) ] ) continue;

			//skip outputing virtual atom unless specified
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::output_virtual ]() &&
				(int)atom.type() == (int)rsd.atom_type_set().n_atomtypes() ) continue;

			++number;
			assert( rsd.chain() < int(chains.size()) ); // silly restriction
			char const chain( chains[ rsd.chain() ] );
			out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
				rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
				F(8,3,atom.xyz()(1)) <<
				F(8,3,atom.xyz()(2)) <<
				F(8,3,atom.xyz()(3)) <<
				F(6,2,1.0) << F(6,2, residue_total_energy ) << '\n';
		}
	}
	out << "ENDMDL\n";
}

void
gen_fold_tree_for_nbr_segments(
	pose::Pose & pose,
	FoldTree & ftree,
	vector1< bool > const & is_ligand,
	vector1< bool > const & is_skipped,
	Real const & nbr_cutoff,
	vector1< bool > & is_nbr
)
{
	Size nres( pose.total_residue() );
	//define neighbors, all protein residues within cutoff from ligand, excluding is_skipped
	set_ss_from_phipsi( pose );
	for( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( !is_ligand[ i ] ) continue;
		Residue const & rsd1( pose.residue(i) );
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			Residue const & rsd2( pose.residue(j) );
			if( is_ligand[ j ] || is_skipped[ j ] || !rsd2.is_protein() || !( pose.secstruct( j ) == 'L' ) ) continue;
			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				if( is_nbr[ j ] ) break;
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
					if ( rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) ) < nbr_cutoff*nbr_cutoff ) {
						is_nbr[ j ] = true;
						break;
					}
				}
			}
		}
	}

	for( Size i = 2, j = 1; i <= pose.total_residue(); ++i, ++j ){
		if( is_nbr[ i ] && !is_nbr[ i-1 ] ) j = 1;
		if( is_nbr[ i ] && !is_nbr[ i+1 ] ){
			ftree.new_jump( i - j, i - static_cast< int >( j / 2 ), i - j );
			ftree.new_jump( i - j, i + 1, i );
		}
	}

}

bool
has_clash(
		pose::Pose & pose,
		vector1< bool > is_checked,
		scoring::ScoreFunctionOP const & scorefxn,
		Real const clash_threshold
	 )
{
	using namespace scoring;
	using namespace chemical;

	using namespace ObjexxFCL::format; // I and F

	bool is_clash( false );
	for( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ){
		if( !is_checked[ seqpos ] ) continue;

		( *scorefxn )( pose );


		// cached energies object
		Energies & energies( pose.energies() );

		// the neighbor/energy links
		EnergyGraph & energy_graph( energies.energy_graph() );

		// search upstream
		for ( graph::Graph::EdgeListIter
				iru  = energy_graph.get_node( seqpos )->lower_edge_list_begin(),
				irue = energy_graph.get_node( seqpos )->lower_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
//			Size const j( edge->get_first_node_ind() );

			// the pair energies cached in the link
			EnergyMap const & emap( edge->fill_energy_map());
			Real const clash( emap[ fa_rep ] );
			if ( clash > clash_threshold ){
				is_clash = true;
				break;
			}
		}
		if( is_clash == true ) break;

		// and downstream
		for ( graph::Graph::EdgeListIter
				iru  = energy_graph.get_node( seqpos )->upper_edge_list_begin(),
				irue = energy_graph.get_node( seqpos )->upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
//			Size const j( edge->get_second_node_ind() );

			// the pair energies cached in the link
			EnergyMap const & emap( edge->fill_energy_map());
			Real const clash( emap[ fa_rep ] );
			if ( clash > clash_threshold ){
				is_clash = true;
				break;
			}
		}
		if( is_clash == true ) break;
	}
	return is_clash;
}
/*
std::string
pep_rmsd_analysis(
	pose::Pose pose,
	Size prot_begin,
	Size prot_end,
	Size pep_begin,
	Size pep_anchor,
	Size pep_end
)
{

	pose::Pose ref_pose;
	std::string ref_name( option[ pep_spec::ref_pose ] );
	core::import_pose::pose_from_pdb( ref_pose, ref_name );

	Size ref_pep_anchor_in( option[ pep_spec::pep_anchor ] );
	if( option[ pep_spec::ref_pep_anchor ].user() ) ref_pep_anchor_in = option[ pep_spec::ref_pep_anchor ];
	std::string ref_pep_chain_in( option[ pep_spec::pep_chain ] );
	if( option[ pep_spec::ref_pep_chain ].user() ) ref_pep_chain_in = option[ pep_spec::ref_pep_chain ];
	Size ref_pep_anchor( ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_in ) );
	Size ref_pep_chain( ref_pose.chain( ref_pep_anchor ) );
	Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pep_chain ) );
	Size ref_pep_end( ref_pose.conformation().chain_end( ref_pep_chain ) );

	Size ref_prot_chain;
	for( Size i = 1; i <= ref_pose.conformation().num_chains(); ++i ){
		if( i == ref_pep_chain ) continue;
		if( !( ref_pose.residue( ref_pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			ref_prot_chain = i;
			break;
		}
	}
	Size ref_prot_begin( ref_pose.conformation().chain_begin( ref_prot_chain ) );
	Size ref_prot_end( ref_pose.conformation().chain_end( ref_prot_chain ) );
	Size ref_prot_anchor( ref_prot_begin );

	Size pep_nterm( pep_anchor - pep_begin );
	Size ref_pep_nterm( ref_pep_anchor - ref_pep_begin );
	Size pep_cterm( pep_end - pep_anchor );
	Size ref_pep_cterm( ref_pep_end - ref_pep_anchor );

	Size nterm( pep_nterm );
	Size cterm( pep_cterm );
	if( pep_nterm < ref_pep_nterm ) ref_pep_begin += ref_pep_nterm - pep_nterm;
	else if( pep_nterm > ref_pep_nterm ){
		pep_begin += pep_nterm - ref_pep_nterm;
		nterm = ref_pep_nterm;
	}
	if( pep_cterm < ref_pep_cterm ) ref_pep_end -= ref_pep_cterm - pep_cterm;
	else if( pep_cterm > ref_pep_cterm ){
		pep_end -= pep_cterm - ref_pep_cterm;
		cterm = ref_pep_cterm;
	}
	//superpose if needed
	if( option[ pep_spec::ref_align ] ){
		id::AtomID_Map< id::AtomID > atom_map;
		id::initialize( atom_map, pose, id::BOGUS_ATOM_ID );

		for ( Size i = prot_begin; i <= prot_end; ++i ) {
			id::AtomID const id1( pose.residue( i ).atom_index( "CA" ), i );
			id::AtomID const id2( ref_pose.residue( i + static_cast< int >( ref_prot_begin ) - static_cast< int >( prot_begin ) ).atom_index( "CA" ), i + static_cast< int >( ref_prot_begin ) - static_cast< int >( prot_begin ) );
			atom_map[ id1 ] = id2;
		}
		core::scoring::ScoreFunctionOP full_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::wts ] ) );
		superimpose_pose( pose, ref_pose, atom_map );

	}
	Real total_rmsd( 0 );
	std::string rmsd_analysis( "" );
	for( int i_seq = 0; i_seq <= nterm + cterm; ++i_seq ){
		Size ref_pep_seqpos = ref_pep_begin + i_seq;
		Size pep_seqpos = pep_begin + i_seq;
		Real sd( ref_pose.residue( ref_pep_seqpos ).xyz( "CA" ).distance_squared( pose.residue( pep_seqpos ).xyz( "CA" ) ) );
		total_rmsd += sd;
		rmsd_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_rmsd:\t" + string_of( std::sqrt( sd ) ) + "\t";
	}
	total_rmsd = std::sqrt( total_rmsd / ( nterm+cterm+1 ) );
	rmsd_analysis += "total_rmsd:\t" + string_of( total_rmsd ) + "\t";
	return rmsd_analysis;

}

std::string
pep_phipsi_analysis(
	pose::Pose pose,
	Size prot_begin,
	Size prot_end,
	Size pep_begin,
	Size pep_anchor,
	Size pep_end
)
{

	pose::Pose ref_pose;
	std::string ref_name( option[ pep_spec::ref_pose ] );
	core::import_pose::pose_from_pdb( ref_pose, ref_name );

	Size ref_pep_anchor_in( option[ pep_spec::pep_anchor ] );
	if( option[ pep_spec::ref_pep_anchor ].user() ) ref_pep_anchor_in = option[ pep_spec::ref_pep_anchor ];
	std::string ref_pep_chain_in( option[ pep_spec::pep_chain ] );
	if( option[ pep_spec::ref_pep_chain ].user() ) ref_pep_chain_in = option[ pep_spec::ref_pep_chain ];
	Size ref_pep_anchor( ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_in ) );
	Size ref_pep_chain( ref_pose.chain( ref_pep_anchor ) );
	Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pep_chain ) );
	Size ref_pep_end( ref_pose.conformation().chain_end( ref_pep_chain ) );

	Size ref_prot_chain;
	for( Size i = 1; i <= ref_pose.conformation().num_chains(); ++i ){
		if( i == ref_pep_chain ) continue;
		if( !( ref_pose.residue( ref_pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			ref_prot_chain = i;
			break;
		}
	}
	Size ref_prot_begin( ref_pose.conformation().chain_begin( ref_prot_chain ) );
	Size ref_prot_end( ref_pose.conformation().chain_end( ref_prot_chain ) );
	Size ref_prot_anchor( ref_prot_begin );

	Size pep_nterm( pep_anchor - pep_begin );
	Size ref_pep_nterm( ref_pep_anchor - ref_pep_begin );
	Size pep_cterm( pep_end - pep_anchor );
	Size ref_pep_cterm( ref_pep_end - ref_pep_anchor );

	Size nterm( pep_nterm );
	Size cterm( pep_cterm );
	if( pep_nterm < ref_pep_nterm ) ref_pep_begin += ref_pep_nterm - pep_nterm;
	else if( pep_nterm > ref_pep_nterm ){
		pep_begin += pep_nterm - ref_pep_nterm;
		nterm = ref_pep_nterm;
	}
	if( pep_cterm < ref_pep_cterm ) ref_pep_end -= ref_pep_cterm - pep_cterm;
	else if( pep_cterm > ref_pep_cterm ){
		pep_end -= pep_cterm - ref_pep_cterm;
		cterm = ref_pep_cterm;
	}

	Real total_phi( 0 );
	Real total_psi( 0 );
	Real total_omega( 0 );
	std::string phipsi_analysis( "" );
	for( int i_seq = 0; i_seq <= nterm + cterm; ++i_seq ){
		Size ref_pep_seqpos = ref_pep_begin + i_seq;
		Size pep_seqpos = pep_begin + i_seq;
		Real ramadev( 0 );

		//delta phi, psi, omega
		Real phi( std::abs( ref_pose.phi( ref_pep_seqpos ) - pose.phi( pep_seqpos ) ) );
		if( phi > 180 ) phi = std::abs( 360 -phi );
		ramadev += ( phi * phi );
		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_phi:\t" + string_of( phi ) + "\t";

		Real psi( std::abs( ref_pose.psi( ref_pep_seqpos ) - pose.psi( pep_seqpos ) ) );
		if( psi > 180 ) psi = std::abs( 360 - psi );
		ramadev += ( psi * psi );
		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_psi:\t" + string_of( psi ) + "\t";

		Real omega( std::abs( ref_pose.omega( ref_pep_seqpos ) - pose.omega( pep_seqpos ) ) );
		if( omega > 180 ) omega = std::abs( 360 - omega );
		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_omega:\t" + string_of( omega ) + "\t";

		phipsi_analysis += "p" + string_of( i_seq - static_cast< int >( nterm ) ) + "_ramadev:\t" + string_of( std::sqrt( ramadev / 2 ) ) + "\t";
	}
	return phipsi_analysis;

}
*/
void
pep_scan_analysis(
	pose::Pose pose,
	core::scoring::ScoreFunctionOP soft_scorefxn,
	core::scoring::ScoreFunctionOP full_scorefxn,
	std::string filename,
	Size pep_begin,
	Size pep_anchor,
	Size pep_end
)
{
	( *full_scorefxn )( pose );
	Pose start_pose( pose );

	vector1< bool > is_pep( pose.total_residue(), false );
	for( Size i = 1; i <= pose.total_residue(); ++i ){
		if( i >= pep_begin && i <= pep_end ) is_pep[ i ] = true;
	}

	for( Size pep_seqpos = pep_begin; pep_seqpos <= pep_end; ++pep_seqpos ){
		if( pep_seqpos == pep_anchor ) continue;
		pose = start_pose;

		char pep_aa( pose.residue( pep_seqpos ).name1() );

		float cutoff( option[ pep_spec::interface_cutoff ] );
		vector1< bool > is_mut_nbr( pose.total_residue(), false );
		Residue const & rsd1( pose.residue( pep_seqpos ) );
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			Residue const & rsd2( pose.residue(j) );
			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
					if ( rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) ) < cutoff*cutoff ) {
						is_mut_nbr[j] = true;
						break;
					}
				}
			}
		}

		{
			// the movable dof's
			kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
			mm_min->set_chi( is_pep );
			mm_min->set_chi( is_mut_nbr );
			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "dfpmin", 0.001, true );

			//define design task and repack task
			pack::task::RestrictResidueToRepackingOperationOP restrict_to_repack_taskop( new pack::task::RestrictResidueToRepackingOperation() );
			pack::task::PreventRepackingOperationOP prevent_repack_taskop( new pack::task::PreventRepackingOperation() );
			pack::task::PackerTaskOP rp_task( pack::task::TaskFactory::create_packer_task( pose ));
			rp_task->initialize_from_command_line().or_include_current( true );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( is_mut_nbr[i] ) {
					restrict_to_repack_taskop->include_residue( i );
					rp_task->nonconst_residue_task( i ).restrict_to_repacking();
				} else {
					rp_task->nonconst_residue_task( i ).prevent_repacking();
					prevent_repack_taskop->include_residue( i );
				}
			}
			pack::task::TaskFactoryOP rottrial_task_factory( new pack::task::TaskFactory );
			rottrial_task_factory->push_back( new pack::task::InitializeFromCommandlineOperation() );
			rottrial_task_factory->push_back( new pack::task::IncludeCurrentOperation() );
			rottrial_task_factory->push_back( restrict_to_repack_taskop );
			rottrial_task_factory->push_back( prevent_repack_taskop );

			protocols::simple_moves::PackRotamersMoverOP pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, rp_task, 1 ) );
			protocols::simple_moves::RotamerTrialsMoverOP rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, rottrial_task_factory ) );
			SequenceMoverOP design_seq = new SequenceMover;
			if( !option[ pep_spec::test_no_pack ] ){
				design_seq->add_mover( pack );
				design_seq->add_mover( rottrial );
			}
			if( !option[ pep_spec::test_no_min ] ){
				design_seq->add_mover( min_mover );
			}

			design_seq->apply( pose );
			( *full_scorefxn )( pose );

		}

		Pose mut_start_pose( pose );
		EMapVector start_emap( pose.energies().total_energies() );
		Real start_total( start_emap.dot( full_scorefxn->weights() ) );
		EMapVector res_avg_emap;
		Real res_avg_total( 0 );
		for( Size resindex = 1; resindex <=20; ++resindex ){
			pose = mut_start_pose;
			chemical::make_sequence_change( pep_seqpos, chemical::AA(resindex), pose );
			if( pose.residue( pep_seqpos ).name1() == 'P' ) continue;
			char mut_aa( pose.residue( pep_seqpos ).name1() );

			// the movable dof's
			kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
			mm_min->set_chi( is_pep );
			mm_min->set_chi( is_mut_nbr );
			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "dfpmin", 0.001, true );

			//define design task and repack task
			pack::task::RestrictResidueToRepackingOperationOP restrict_to_repack_taskop( new pack::task::RestrictResidueToRepackingOperation() );
			pack::task::PreventRepackingOperationOP prevent_repack_taskop( new pack::task::PreventRepackingOperation() );
			pack::task::PackerTaskOP rp_task( pack::task::TaskFactory::create_packer_task( pose ));
			rp_task->initialize_from_command_line().or_include_current( true );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( is_mut_nbr[i] ) {
					restrict_to_repack_taskop->include_residue( i );
					rp_task->nonconst_residue_task( i ).restrict_to_repacking();
				} else {
					rp_task->nonconst_residue_task( i ).prevent_repacking();
					prevent_repack_taskop->include_residue( i );
				}
			}
			pack::task::TaskFactoryOP rottrial_task_factory( new pack::task::TaskFactory );
			rottrial_task_factory->push_back( new pack::task::InitializeFromCommandlineOperation() );
			rottrial_task_factory->push_back( new pack::task::IncludeCurrentOperation() );
			rottrial_task_factory->push_back( restrict_to_repack_taskop );
			rottrial_task_factory->push_back( prevent_repack_taskop );

			protocols::simple_moves::PackRotamersMoverOP pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, rp_task, 1 ) );
			protocols::simple_moves::RotamerTrialsMoverOP rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, rottrial_task_factory ) );
			SequenceMoverOP design_seq = new SequenceMover;
			if( !option[ pep_spec::test_no_pack ] ){
				design_seq->add_mover( pack );
				design_seq->add_mover( rottrial );
			}
			if( !option[ pep_spec::test_no_min ] ){
				design_seq->add_mover( min_mover );
			}

			design_seq->apply( pose );
			( *full_scorefxn )( pose );

			EMapVector diff_emap( pose.energies().total_energies() );
			diff_emap -= start_emap;
			Real diff_total( pose.energies().total_energies().dot( full_scorefxn->weights() ) - start_total );

			res_avg_emap += diff_emap;
			res_avg_total += diff_total;

			std::cout << pep_aa << "\t" << string_of( pep_seqpos ) << "\t" << mut_aa << "\tspec_res:\t" << diff_emap.weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< diff_total << "\n";

		}

		for( EMapVector::iterator itr = res_avg_emap.begin(); itr != res_avg_emap.end(); ++itr ){
			*itr *= ( 1.0 / 19 );
		}
		res_avg_total = res_avg_total / 19;

		std::cout << pep_aa << "\t" << string_of( pep_seqpos ) << "\t" << "X" << "\tspec_res_avg:\t" << res_avg_emap.weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< res_avg_total << "\n";

	}

}

void
pep_energies_analysis(
	pose::Pose pose,
	core::scoring::ScoreFunctionOP full_scorefxn,
	std::string filename,
	core::scoring::EMapVector & emap_total_avg,
	vector1< core::scoring::EMapVector > & emap_res_avg,
	vector1< core::scoring::EMapVector > & emap_res_sd
)
{
	( *full_scorefxn )( pose );
	EMapVector emap_total( pose.energies().total_energies() );
	emap_total_avg += emap_total;

	for( Size i_seq = 1; i_seq <= pose.total_residue(); ++i_seq ){
		emap_res_avg[ i_seq ] += pose.energies().residue_total_energies( i_seq );
		char pep_aa( pose.residue( i_seq ).name1() );
		std::cout << filename << "\t" << pep_aa << "\t" << string_of( i_seq ) << "\ttotal_res\t" << pose.energies().residue_total_energies( i_seq ).weighted_string_of( full_scorefxn->weights() ) <<"\ttotal_score:\t"<< pose.energies().residue_total_energies( i_seq ).dot( full_scorefxn->weights() ) << "\n";
	}
	std::cout << std::endl;
}

void add_pep_buffer_res(
	pose::Pose & pose,
	Size & pep_begin,
	Size & pep_anchor,
	Size & pep_end,
	bool add_nterm,
	bool add_cterm
)
{
		ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
		ResidueOP ala( ResidueFactory::create_residue( rsd_set.name_map( "ALA" ) ) );

		if( add_cterm ){
			pose.conformation().safely_append_polymer_residue_after_seqpos( *ala, pep_end, true );
			pep_end = pep_end + 1;
			Size nres_pep = nres_pep + 1;
			pose.set_omega( pep_end - 1, 180.0 );
		}

		if( add_nterm ){
			pose.conformation().safely_prepend_polymer_residue_before_seqpos( *ala, pep_begin, true );
			pep_end = pep_end + 1;
			pep_anchor = pep_anchor + 1;
			Size nres_pep = nres_pep + 1;
			pose.set_omega( pep_begin, 180.0 );
		}
}

void remove_pep_buffer_res(
	pose::Pose & pose,
	Size & pep_begin,
	Size & pep_anchor,
	Size & pep_end,
	bool add_nterm,
	bool add_cterm
)
{
		if( add_cterm ){
			pose.conformation().delete_residue_slow( pep_end );
			pep_end = pep_end - 1;
		}

		if( add_nterm ){
			pose.conformation().delete_residue_slow( pep_begin );
			pep_anchor = pep_anchor - 1;
			pep_end = pep_end - 1;
		}
}

//make backbone from random frag phi/psi insertions//
pose::Pose
gen_pep_bb_frag(
	pose::Pose pose,
	Size pep_begin,
	Size pep_end,
	scoring::ScoreFunctionOP cen_scorefxn,
	protocols::frags::TorsionFragmentLibrary lib
)
{
	Size nres_pep( pep_end - pep_begin + 1 );

	( *cen_scorefxn )( pose );
	Real best_score( pose.energies().total_energies().dot( cen_scorefxn->weights() ) );
	Size n_build_loop( option[ pep_spec::n_build_loop ] );
	MonteCarloOP mc_frag ( new MonteCarlo( pose, *cen_scorefxn, 2.0 ) );
	for ( Size build_loop_inner = 1; build_loop_inner <= n_build_loop; ++build_loop_inner ) {
		//choose an insertion position
		Size pos(  static_cast< int > ( nres_pep * RG.uniform() ) + 1 );
		Size lib_pos(  static_cast< int > ( 20 * RG.uniform() ) + 1 );

		Size const nfrags( lib[ lib_pos ].size() );
		int const frag_index( static_cast< int >( nfrags * RG.uniform() + 1 ) );
		lib[ lib_pos ][ frag_index ].insert( pose, pep_begin - 1 + pos );
		if( mc_frag->boltzmann( pose ) ){
			Real test_score( pose.energies().total_energies().dot( cen_scorefxn->weights() ) );
			//need to do my own eval for <= because MC mover is < only
			if( test_score <= best_score ){
				best_score = test_score;
				mc_frag->reset( pose );
			}
		}
	}
	mc_frag->recover_low( pose );
	return pose;
}

//make backbone from random rama phi_psi insertions//
pose::Pose
gen_pep_bb_rama(
	pose::Pose pose,
	Size pep_begin,
	Size pep_end,
	scoring::ScoreFunctionOP cen_scorefxn
)
{
	Size nres_pep( pep_end - pep_begin + 1 );

	vector1< core::scoring::Ramachandran > rama_movers;
	core::scoring::Ramachandran loop_rama_mover( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
	loop_rama_mover.init_rama_sampling_table( 3 );
	rama_movers.push_back( loop_rama_mover );
	if( option[ pep_spec::random_rama_ss_type ] ){
		core::scoring::Ramachandran helix_rama_mover( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
		helix_rama_mover.init_rama_sampling_table( 1 );
		rama_movers.push_back( helix_rama_mover );
		core::scoring::Ramachandran sheet_rama_mover( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
		sheet_rama_mover.init_rama_sampling_table( 2 );
		rama_movers.push_back( sheet_rama_mover );
	}

	( *cen_scorefxn )( pose );
	Size n_build_loop( option[ pep_spec::n_build_loop ] );
	MonteCarloOP mc_rama ( new MonteCarlo( pose, *cen_scorefxn, 2.0 ) );
	for ( Size build_loop_inner = 1; build_loop_inner <= n_build_loop; ++build_loop_inner ) {
		// choose an insertion position
		Size pos(  static_cast< int > ( nres_pep * RG.uniform() ) + pep_begin );
		Size rama_mover_index( 1 );
		if( option[ pep_spec::random_rama_ss_type ] ) rama_mover_index = static_cast< int > ( RG.uniform() * rama_movers.size() + 1 );
		Real rama_phi, rama_psi;
		if( option[ pep_spec::use_input_seq ] ){
			chemical::AA actual_aa( aa_from_oneletter_code( pose.residue( pos ).name1() ) );
			rama_movers[ rama_mover_index ].random_phipsi_from_rama( actual_aa, rama_phi, rama_psi );
		}
		else{
			//use ala rama map 95% of time, use gly map 5%
//			int resindex( static_cast< int > ( 20 * RG.uniform() + 1 ) );
			chemical::AA aa( aa_from_oneletter_code( 'A' ) );
			if( RG.uniform() > 0.95 ) aa = aa_from_oneletter_code( 'G' );
			rama_movers[ rama_mover_index ].random_phipsi_from_rama( aa, rama_phi, rama_psi );
		}
		pose.set_phi( pos, rama_phi );
		pose.set_psi( pos, rama_psi );
		mc_rama->boltzmann( pose );
	}
	mc_rama->recover_low( pose );
	return pose;
}


pose::Pose
perturb_pep_bb(
	pose::Pose pose,
	kinematics::MoveMapOP mm_move,
	scoring::ScoreFunctionOP cen_scorefxn,
	Size n_iter
)
{
	( *cen_scorefxn )( pose );
	protocols::simple_moves::SmallMoverOP cg_small( new protocols::simple_moves::SmallMover( mm_move, 2.0, 1 ) );
	cg_small->angle_max( 'H', 2.0 );
	cg_small->angle_max( 'E', 2.0 );
	cg_small->angle_max( 'L', 2.0 );

	MonteCarloOP mc_cg ( new MonteCarlo( pose, *cen_scorefxn, 2.0  ) );
	TrialMoverOP cg_trial = new TrialMover( cg_small, mc_cg );
	RepeatMoverOP cg_cycle = new RepeatMover( cg_trial, 100 );
	for( Size cgmove_loop = 1; cgmove_loop <= n_iter; cgmove_loop++ ){
		mc_cg->reset( pose );
		cg_cycle->apply( pose );
	}
	mc_cg->recover_low( pose );
	return pose;
}

pose::Pose
packmin_unbound_pep(
	pose::Pose pose,
	scoring::ScoreFunctionOP full_scorefxn,
	scoring::ScoreFunctionOP soft_scorefxn
)
{
	pose.update_residue_neighbors();
	( *full_scorefxn )( pose );
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().or_include_current( true );
	task->restrict_to_repacking();
	protocols::simple_moves::PackRotamersMoverOP pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, task, 1 ) );
	pack->apply( pose );
	if( !option[ pep_spec::test_no_min ] ){
		kinematics::MoveMapOP mm ( new kinematics::MoveMap );
		mm->set_chi( true );
		protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm, full_scorefxn, "dfpmin", 0.001, true );
		min_mover->apply( pose );
	}
	return pose;
}

void
RunPepSpec()
{
	//load loop lengths//
	Size n_peptides( option[ pep_spec::n_peptides ] );

	std::string input_name ( basic::options::start_file() );
	Pose start_pose;
	core::import_pose::pose_from_pdb( start_pose, input_name );
	ResidueTypeSet const & rsd_set( start_pose.residue(1).residue_type_set() );

	//data out
	std::string out_nametag( "data" );
	if( option[ out::file::o ].user() ) out_nametag = option[ out::file::o ];

	//create output file
	std::string out_file_name_str( out_nametag + ".spec" );
	char const *out_file_name = out_file_name_str.c_str();
	std::fstream out_file( out_file_name, std::ios::out );

	std::string output_seq;

	//define scoring functions//
	core::scoring::ScoreFunctionOP soft_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::soft_wts ] ) );
	core::scoring::ScoreFunctionOP full_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::wts ] ) );
	core::scoring::ScoreFunctionOP cen_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::cen_wts ] ) );

	VallData vall( option[ pep_spec::vall ] );
	TorsionFragmentLibrary lib;

	lib.resize( 20 );
	//rsd type
	for( Size i = 1; i <= 20; ++i ){
		char aa( oneletter_code_from_aa( chemical::AA( i ) ) );
		std::string frag_seq;
		frag_seq.push_back( aa );
		std::string ss_seq( "S" );
		Real seq_wt( 1.0 );
		Real ss_wt( 0.0 );
		vall.get_frags( 10000, frag_seq, ss_seq, seq_wt, ss_wt, false, false, true, lib[ i ] );
	}

	std::string cg_res_type( option[ pep_spec::cg_res_type ] );

	for( Size peptide_loop = 1; peptide_loop <= n_peptides; ++peptide_loop ){	//LOOP how many final peps

		Pose & pose( start_pose );
		ResidueOP ala( ResidueFactory::create_residue( rsd_set.name_map( "ALA" ) ) );
//		Size n_append( static_cast< int >( 6 * RG.uniform() + 1 ) );
//		Size n_prepend( static_cast< int >( 6 * RG.uniform() + 1 ) );
/*
		Size n_append( 6 );
		Size n_prepend( 6 );
		for( Size i = 1; i <= n_append; ++i ){
			pose.conformation().safely_append_polymer_residue_after_seqpos( *ala, i, true );
		}
		for( Size i = 1; i <= n_prepend; ++i ){
			pose.conformation().safely_prepend_polymer_residue_before_seqpos( *ala, 1, true );
		}
		for( Size i_omega = 1; i_omega <= pose.total_residue() - 1; ++i_omega ){
			pose.set_omega( i_omega, 180.0 );
		}
*/
		Size pep_begin( 1 );
		Size pep_end( pose.total_residue() );
/*
		//rsd type
		for( Size mut_site = pep_begin; mut_site <= pep_end; mut_site++ ){
			if( option[ pep_spec::cg_res_type ].user() ) chemical::make_sequence_change( mut_site, chemical::aa_from_name( cg_res_type ), pose );
			else chemical::make_sequence_change( mut_site, chemical::AA( static_cast< int > ( 20 * RG.uniform() + 1 ) ), pose );
		}
*/
		//gen fold tree//
		FoldTree f( pose.total_residue() );
		pose.fold_tree( f );

		//convert to CG residues//
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
////////////////////////////////////////////
		Size pep_anchor( 2 );
		if( option[ pep_spec::add_buffer_res ] ) add_pep_buffer_res( pose, pep_begin, pep_anchor, pep_end, true, false );

		if( option[ pep_spec::cen_bb_frag ] ) pose = gen_pep_bb_frag( pose, pep_begin, pep_end, cen_scorefxn, lib );
		else if( option[ pep_spec::cen_bb_rama ] ) pose = gen_pep_bb_rama( pose, pep_begin, pep_end, cen_scorefxn );

		//min movemap//
		kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
		mm_min->set_chi( true );
		kinematics::MoveMapOP mm_move ( new kinematics::MoveMap );
		mm_move->set_bb( true );

		//small CG backbone moves//
		Size n_cgperturb_iter( 1 );
		if( !option[ pep_spec::test_no_cgrelax ] ) pose = perturb_pep_bb( pose, mm_move, cen_scorefxn, n_cgperturb_iter );

		//debug dump CG data and pdb//
		if( option[ pep_spec::test_frag_only ] ){
			if( option[ pep_spec::test_dump_all ] ) dump_efactor_pdb( pose, cen_scorefxn, "pdbs/" + out_nametag + "_" + string_of( peptide_loop ) + "_cen.pdb" );
			out_file << out_nametag + "_" + string_of( peptide_loop ) + "_cen.pdb" +
				"\t"<<pose.energies().total_energies().weighted_string_of( cen_scorefxn->weights() );
			out_file<<"\ttotal_score:\t"<<pose.energies().total_energies()[ total_score ]<<"\t";
			out_file<<std::endl;

			continue;
		}

		//switch back to fullatom
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );

/*
		//randomize pep sequence//
		for(Size mut_site = pep_begin; mut_site <= pep_end; mut_site++){ //over all pep positions
			int resindex;
			resindex = static_cast< int > ( 20 * RG.uniform() + 1 );
			chemical::make_sequence_change( mut_site, chemical::AA(resindex), pose );
		}
*/
		//define movers//
		protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "dfpmin", 0.001, true );

		//define design task and repack task
		pack::task::PackerTaskOP dz_task( pack::task::TaskFactory::create_packer_task( pose ));
		dz_task->initialize_from_command_line();
		dz_task->restrict_to_repacking();

		pack::task::TaskFactoryOP task_factory( new pack::task::TaskFactory );
		task_factory->push_back( new pack::task::InitializeFromCommandlineOperation() );
		task_factory->push_back( new pack::task::RestrictToRepackingOperation() );

		protocols::simple_moves::PackRotamersMoverOP dz_pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, dz_task, 1 ) );
		protocols::simple_moves::RotamerTrialsMoverOP dz_rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, task_factory ) );
		SequenceMoverOP design_seq = new SequenceMover;
		design_seq->add_mover( dz_pack );
		design_seq->add_mover( dz_rottrial );
//		design_seq->add_mover( min_mover );

		EMapVector emap_total_avg;
		vector1< EMapVector > emap_res_avg( pose.total_residue() );
		vector1< EMapVector > emap_res_sd( pose.total_residue() );

		//replace termini
		core::pose::add_lower_terminus_type_to_pose_residue( pose, pep_begin );
		core::pose::add_upper_terminus_type_to_pose_residue( pose, pep_end );
		pose.conformation().update_polymeric_connection( pep_begin );
		pose.conformation().update_polymeric_connection( pep_begin + 1 );
		pose.conformation().update_polymeric_connection( pep_end );
		pose.conformation().update_polymeric_connection( pep_end - 1 );

		//RUN MOVERS//
		( *soft_scorefxn )( pose );
		design_seq->apply( pose );
/////////////////////////////
		if( option[ pep_spec::energies_analysis ] ) pep_energies_analysis( pose, full_scorefxn, "WITH BUFFER", emap_total_avg, emap_res_avg, emap_res_sd );
		if( option[ pep_spec::scan_analysis ] ) pep_scan_analysis( pose, soft_scorefxn, full_scorefxn, "WITH BUFFER", pep_begin, pep_anchor, pep_end );

		if( option[ pep_spec::add_buffer_res ] ) remove_pep_buffer_res( pose, pep_begin, pep_anchor, pep_end, true, false );



		//replace termini
		core::pose::add_lower_terminus_type_to_pose_residue( pose, pep_begin );
		core::pose::add_upper_terminus_type_to_pose_residue( pose, pep_end );
		pose.conformation().update_polymeric_connection( pep_begin );
		pose.conformation().update_polymeric_connection( pep_begin + 1 );
		pose.conformation().update_polymeric_connection( pep_end );
		pose.conformation().update_polymeric_connection( pep_end - 1 );

		( *full_scorefxn )( pose );

		//Analysis//
		output_seq.clear();
		for(Size i = pep_begin; i <= pep_end; i++){
			output_seq.append( 1, pose.residue( i ).name1() );
		}
		Real total_score( pose.energies().total_energies().dot( full_scorefxn->weights() ) );

		//Output//
		out_file << out_nametag + "_" + string_of( peptide_loop ) + ".pdb" + "\t" << output_seq << "\t" << pose.energies().total_energies().weighted_string_of( full_scorefxn->weights() );
		out_file<<"\ttotal_score:\t"<<total_score<<"\t";
		out_file<<std::endl;

		if( option[ pep_spec::energies_analysis ] ) pep_energies_analysis( pose, full_scorefxn, "NO BUFFER", emap_total_avg, emap_res_avg, emap_res_sd );
		if( option[ pep_spec::scan_analysis ] ) pep_scan_analysis( pose, soft_scorefxn, full_scorefxn, "NO BUFFER", pep_begin, pep_anchor, pep_end );
		if( option[ pep_spec::test_dump_all ] ) dump_efactor_pdb( pose, full_scorefxn, "pdbs/" + out_nametag + "_" + string_of( peptide_loop ) + ".pdb" );

	}
}


int main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init(argc, argv);

		RunPepSpec();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
