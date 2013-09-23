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
#include <utility/file/file_sys_util.hh>

 #include <numeric/xyzVector.hh>
 #include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// //REMOVE LATER!
// #include <utility/io/izstream.hh>

#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/PeriodicFunc.hh>
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
aa2index(
	chemical::AA aa
)
{
	Size index( 0 );
	switch( aa ){
		case chemical::aa_ala :
			index = 1;
			break;
		case chemical::aa_cys :
			index = 2;
			break;
		case chemical::aa_asp :
			index = 3;
			break;
		case chemical::aa_glu :
			index = 4;
			break;
		case chemical::aa_phe :
			index = 5;
			break;
		case chemical::aa_gly :
			index = 6;
			break;
		case chemical::aa_his :
			index = 7;
			break;
		case chemical::aa_ile :
			index = 8;
			break;
		case chemical::aa_lys :
			index = 9;
			break;
		case chemical::aa_leu :
			index = 10;
			break;
		case chemical::aa_met :
			index = 11;
			break;
		case chemical::aa_asn :
			index = 12;
			break;
		case chemical::aa_pro :
			index = 13;
			break;
		case chemical::aa_gln :
			index = 14;
			break;
		case chemical::aa_arg :
			index = 15;
			break;
		case chemical::aa_ser :
			index = 16;
			break;
		case chemical::aa_thr :
			index = 17;
			break;
		case chemical::aa_val :
			index = 18;
			break;
		case chemical::aa_trp :
			index = 19;
			break;
		case chemical::aa_tyr :
			index = 20;
			break;
		default :
			break;
	}
	return index;
}

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

	std::cout << "sel flex_prot, ";
	for( Size i = 1; i <= pose.total_residue(); ++i ){
		if( is_nbr[ i ] ) std::cout << "resi " << string_of( i ) << " or ";
	}
	std::cout << std::endl;
}

bool
has_clash(
		pose::Pose pose,
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
			Size const j( edge->get_first_node_ind() );

			// the pair energies cached in the link
			EnergyMap const & emap( edge->fill_energy_map());
			Real const clash( emap[ fa_rep ] );
			if ( clash > clash_threshold ){
//				std::cout<< "fa_rep: " << string_of( clash ) << " at " << pose.residue( j ).name1() << string_of( j ) << "-" << pose.residue( seqpos ).name1() << string_of( seqpos ) << std::endl;
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
			Size const j( edge->get_second_node_ind() );

			// the pair energies cached in the link
			EnergyMap const & emap( edge->fill_energy_map());
			Real const clash( emap[ fa_rep ] );
			if ( clash > clash_threshold ){
//				std::cout<< "fa_rep: " << string_of( clash ) << " at " << pose.residue( j ).name1() << string_of( j ) << "-" << pose.residue( seqpos ).name1() << string_of( seqpos ) << std::endl;
				is_clash = true;
				break;
			}
		}
		if( is_clash == true ) break;
	}
	return is_clash;
}

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
		core::pose::initialize_atomid_map( atom_map, pose, id::BOGUS_ATOM_ID );

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

void add_termini(
	pose::Pose & pose,
	Size & pep_begin,
	Size & pep_end
)
{
	core::pose::add_lower_terminus_type_to_pose_residue( pose, pep_begin );
	core::pose::add_upper_terminus_type_to_pose_residue( pose, pep_end );
	pose.conformation().update_polymeric_connection( pep_begin );
	pose.conformation().update_polymeric_connection( pep_begin + 1 );
	pose.conformation().update_polymeric_connection( pep_end );
	pose.conformation().update_polymeric_connection( pep_end - 1 );

}

void add_pep_res(
	pose::Pose & pose,
	Size prot_anchor,
	Size & pep_begin,
	Size & pep_anchor,
	Size & pep_end,
	bool add_nterm,
	bool add_cterm
)
{
		ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
		ResidueOP gly( ResidueFactory::create_residue( rsd_set.name_map( "GLY" ) ) );

		if( add_cterm ){
			pose.conformation().safely_append_polymer_residue_after_seqpos( *gly, pep_end, true );
			pep_end = pep_end + 1;
			pose.set_omega( pep_end - 1, 180.0 );
			pose.conformation().update_polymeric_connection( pep_end );
			pose.conformation().update_polymeric_connection( pep_end - 1 );
		}

		if( add_nterm ){
			pose.conformation().safely_prepend_polymer_residue_before_seqpos( *gly, pep_begin, true );
			pep_end = pep_end + 1;
			pep_anchor = pep_anchor + 1;
			pose.set_omega( pep_begin, 180.0 );
			pose.conformation().update_polymeric_connection( pep_begin );
			pose.conformation().update_polymeric_connection( pep_begin + 1 );
		}

		//replace termini
		add_termini( pose, pep_begin, pep_end );
}

void remove_pep_res(
	pose::Pose & pose,
	Size prot_anchor,
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
			pose.conformation().update_polymeric_connection( pep_end );
			pose.conformation().update_polymeric_connection( pep_end - 1 );
		}

		if( add_nterm ){
			pose.conformation().delete_residue_slow( pep_begin );
			pep_anchor = pep_anchor - 1;
			pep_end = pep_end - 1;
			pose.conformation().update_polymeric_connection( pep_begin );
			pose.conformation().update_polymeric_connection( pep_begin + 1 );
		}

		//replace termini
		add_termini( pose, pep_begin, pep_end );
}

//make backbone from random frag phi/psi insertions//
void
gen_pep_bb_frag(
	pose::Pose & pose,
	Size pep_begin,
	Size pep_anchor,
	Size pep_end,
	scoring::ScoreFunctionOP cen_scorefxn,
	protocols::frags::TorsionFragmentLibrary lib
)
{
	Size nres_pep( pep_end - pep_begin + 1 );

	//randomize backbone
	for( Size seqpos = pep_begin; seqpos <= pep_end; ++seqpos ){
		pose.set_phi( seqpos, RG.uniform() * 360.0 );
		pose.set_psi( seqpos, RG.uniform() * 360.0 );
		pose.set_omega( seqpos, 180.0 );
	}

	( *cen_scorefxn )( pose );
	Real best_score( pose.energies().total_energies().dot( cen_scorefxn->weights() ) );
	Size n_build_loop( option[ pep_spec::n_build_loop ] );
	MonteCarloOP mc_frag ( new MonteCarlo( pose, *cen_scorefxn, 3.0 ) );
//	mc_frag->set_autotemp( true, 1.0 );
	for ( Size build_loop_inner = 1; build_loop_inner <= n_build_loop; ++build_loop_inner ) {

		//choose an insertion position, index from 1
		Size pos(  static_cast< int > ( nres_pep * RG.uniform() ) + 1 );
		Size lib_pos(  static_cast< int > ( 20 * RG.uniform() ) + 1 );
		//use actual aa if is anchor
		if( option[ pep_spec::test_no_design ] || option[ pep_spec::use_input_seq ] || pos == ( pep_anchor - pep_begin + 1 ) ){
			lib_pos = aa2index( pose.residue( pep_begin + pos - 1 ).aa() );
		}

		Size const nfrags( lib[ lib_pos ].size() );
		int const frag_index( static_cast< int >( nfrags * RG.uniform() + 1 ) );
		lib[ lib_pos ][ frag_index ].insert( pose, pep_begin + pos - 1 );
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
}

void
gen_pep_bb_frag_old(
	pose::Pose & pose,
	Size pep_begin,
	Size pep_anchor,
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
//	mc_frag->set_autotemp( true, 1.0 );
	for ( Size build_loop = 1; build_loop <= static_cast< int >( std::sqrt( n_build_loop ) ); ++build_loop ) {
		mc_frag->reset( pose );
		for ( Size build_loop_inner = 1; build_loop_inner <= static_cast< int >( std::sqrt( n_build_loop ) ); ++build_loop_inner ) {
			// choose an insertion position
			Size pos(  static_cast< int > ( nres_pep * RG.uniform() ) + 1 );

			Size const nfrags( lib[ pos ].size() );
			int const frag_index( static_cast< int >( nfrags * RG.uniform() + 1 ) );
			lib[ pos ][ frag_index ].insert( pose, pep_begin - 1 + pos );

			mc_frag->boltzmann( pose );
		}

	}
	mc_frag->recover_low( pose );

}

//make backbone from random rama phi_psi insertions//
void
gen_pep_bb_rama(
	pose::Pose & pose,
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
}


void
perturb_pep_bb(
	pose::Pose & pose,
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
}

void
refine_fa_pep_bb(
	Pose & pose,
	Size pep_anchor,
	vector1< bool > is_pep,
	vector1< bool > is_pep_nbr,
	scoring::ScoreFunctionOP full_scorefxn
)
{
	for( Size mut_site = 1; mut_site <= pose.total_residue(); mut_site++ ){
		if( !( pose.residue( mut_site ).is_protein() ) ) continue;
		else if( is_pep[ mut_site ] && mut_site != pep_anchor ){
			chemical::make_sequence_change( mut_site, chemical::aa_gly, pose );
		}
		else if( is_pep_nbr[ mut_site ] && pose.residue( mut_site ).name1() != 'G' ){
			chemical::make_sequence_change( mut_site, chemical::aa_ala, pose );
		}
	}
	core::scoring::ScoreFunctionOP rep_scorefxn( new ScoreFunction );
	rep_scorefxn->set_weight( fa_rep, full_scorefxn->get_weight( fa_rep ) );
	rep_scorefxn->set_weight( rama, 1.0 );
	( *rep_scorefxn )( pose );

	kinematics::MoveMapOP mm ( new kinematics::MoveMap );
	mm->set_bb( is_pep );

	if( has_clash( pose, is_pep, rep_scorefxn, option[ pep_spec::clash_fa_rep_tol ] ) ){
		protocols::simple_moves::SmallMoverOP rep_small_mover( new protocols::simple_moves::SmallMover( mm, 1.0, 1 ) );
		rep_small_mover->angle_max( 'H', 1.0 );
		rep_small_mover->angle_max( 'E', 1.0 );
		rep_small_mover->angle_max( 'L', 1.0 );
		MonteCarloOP mc_rep ( new MonteCarlo( pose, *rep_scorefxn, 1.0 ) );
		TrialMoverOP rep_trial = new TrialMover( rep_small_mover, mc_rep );
		RepeatMoverOP rep_cycle = new RepeatMover( rep_trial, 1000 );

		rep_cycle->apply( pose );
		mc_rep->recover_low( pose );
	}
	if( has_clash( pose, is_pep, rep_scorefxn, option[ pep_spec::clash_fa_rep_tol ] ) ){
		protocols::simple_moves::MinMoverOP rep_min_mover = new protocols::simple_moves::MinMover( mm, rep_scorefxn, "dfpmin", 0.001, true );
		rep_min_mover->apply( pose );
	}
}

void
packmin_unbound_pep(
	pose::Pose & pose,
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
}

void
RunPepSpec()
{
	//load loop lengths//
	Size n_peptides( option[ pep_spec::n_peptides ] );
	Size n_cgrelax_loop( option[ pep_spec::n_cgrelax_loop ]  );
	Size n_farelax_loop( option[ pep_spec::n_farelax_loop ]  );

	//data out
	std::string out_nametag( "data" );
	if( option[ out::file::o ].user() ) out_nametag = option[ out::file::o ];
	std::string pdb_dir( out_nametag + ".pdbs" );
	utility::file::create_directory( pdb_dir );
	std::string out_file_name_str( out_nametag + ".spec" );
	char const *out_file_name = out_file_name_str.c_str();
	std::fstream out_file( out_file_name, std::ios::out );

	//load input structure pdbs/csts//
	pose::Pose pose;
	vector1< std::string > pdb_filenames;
	if( option[ pep_spec::pdb_list ].user() ){
		std::string pdb_list_filename( option[ pep_spec::pdb_list ] );
		std::ifstream pdb_list_data( pdb_list_filename.c_str() );
		if ( !pdb_list_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + pdb_list_filename + '\n' );
		}
		std::string pdb_list_line;
		while( !getline( pdb_list_data, pdb_list_line, '\n' ).eof() ) {
			std::string this_filename( pdb_list_line );
			pdb_filenames.push_back( this_filename );
		}
	}
	else{
		pdb_filenames.push_back( basic::options::start_file() );

	}

	core::import_pose::pose_from_pdb( pose, pdb_filenames[ 1 ] );
	ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );

	//convert user input to internal values//
	Size const pep_anchor_in( option[ pep_spec::pep_anchor ] );
	std::string const pep_chain_in( option[ pep_spec::pep_chain ] );


	Size start_pep_anchor( pose.pdb_info()->pdb2pose( pep_chain_in[ 0 ], pep_anchor_in ) );
	if( start_pep_anchor == 0 ) utility_exit_with_message( "pep_chain / start_pep_anchor combo not found\n" );
	Size pep_chain( pose.chain( start_pep_anchor ) );
	Size start_pep_begin( pose.conformation().chain_begin( pep_chain ) );
	Size start_pep_end( pose.conformation().chain_end( pep_chain ) );
	Size prot_chain;
	for( Size i = 1; i <= pose.conformation().num_chains(); ++i ){
		if( i == pep_chain ) continue;
		if( !( pose.residue( pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			prot_chain = i;
			break;
		}
	}
	Size prot_begin( pose.conformation().chain_begin( prot_chain ) );
	Size prot_end( pose.conformation().chain_end( prot_chain ) );
	Size prot_anchor( prot_begin );

	Real cutoff( option[ pep_spec::interface_cutoff ] );
	std::string output_seq;

	//gen fold tree//
	FoldTree f( pose.total_residue() );

	std::string start_pep_anchor_root( "CB" );
	if( pose.residue( start_pep_anchor ).name1() == 'G' ) start_pep_anchor_root = "CA";
	Size pep_jump( 2 );
	if( prot_chain < pep_chain ){
		pep_jump = f.new_jump( prot_anchor, start_pep_anchor, start_pep_begin - 1 );
		f.set_jump_atoms( pep_jump, "CA", start_pep_anchor_root );
		if( start_pep_end != pose.total_residue() ) f.new_jump( prot_anchor, start_pep_end + 1, start_pep_end );
		if( prot_end + 1 != start_pep_begin ) f.new_jump( prot_anchor, prot_end + 1, prot_end );
	}
	else{
		pep_jump = f.new_jump( start_pep_anchor, prot_anchor, prot_begin - 1 );
		f.set_jump_atoms( pep_jump, start_pep_anchor_root, "CA" );
		if( prot_end != pose.total_residue() ) f.new_jump( prot_anchor, prot_end + 1, prot_end );
		if( start_pep_end + 1 != prot_begin ) f.new_jump( start_pep_anchor, start_pep_end + 1, start_pep_end );
	}
	f.reorder( prot_anchor );
	std::cout << f << "\n";

	//define scoring functions//
	core::scoring::ScoreFunctionOP soft_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::soft_wts ] ) );
	core::scoring::ScoreFunctionOP full_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::wts ] ) );
	core::scoring::ScoreFunctionOP cen_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::cen_wts ] ) );

	//////////////fragment library/////////////////
/*
	VallLibrarian librarian;
	librarian.add_fragment_gen( new LengthGen( 1 ) );
	librarian.catalog( FragmentLibraryManager::get_instance()->get_Vall() );
*/
	VallData vall( option[ pep_spec::vall ] );
	TorsionFragmentLibrary lib;

	lib.resize( 20 );
	//lib indexed same as AA enum
	for( Size i = 1; i <= 20; ++i ){
		char aa( oneletter_code_from_aa( chemical::AA( i ) ) );
		std::string frag_seq;
		frag_seq.push_back( aa );
		std::string ss_seq( "S" );
		Real seq_wt( 1.0 );
		Real ss_wt( 0.0 );
		vall.get_frags( 10000, frag_seq, ss_seq, seq_wt, ss_wt, false, false, true, lib[ i ] );
	}
/*
	TorsionFragmentLibrary old_lib;
	old_lib.resize( start_pep_end - start_pep_begin + 1 );
	for ( Size i=1; i<= start_pep_end - start_pep_begin + 1; ++i ) {
		Size seqpos( i + start_pep_begin - 1 );

		std::string frag_seq;
		std::string ss_seq;
		for( Size ii = 1; ii <= 1; ++ii ){
			frag_seq.append( 1, pose.residue( seqpos ).name1() );
			ss_seq.append( 1, 'S' );
		}

		Real seq_wt( 0.0 );
		vall.get_frags( 10000, frag_seq, ss_seq, seq_wt, 1.0, false, false, true, old_lib[i] );
	}
*/

	for( Size peptide_loop = 1; peptide_loop <= n_peptides; ++peptide_loop ){

		Size pep_begin( start_pep_begin );
		Size pep_anchor( start_pep_anchor );
		Size pep_end( start_pep_end );

		//load random start pdb//
		int pose_index( static_cast< int >( RG.uniform() * pdb_filenames.size() + 1 ) );
		std::string pdb_filename( pdb_filenames[ pose_index ] );
		std::cout<<"Initializing "<< out_nametag + "_" + string_of( peptide_loop ) + " with " + pdb_filename << std::endl;
		core::import_pose::pose_from_pdb( pose, pdb_filename );
		pose::Pose start_pose( pose );

		protocols::viewer::add_conformation_viewer( pose.conformation(), "pep_spec_pose" );

		//set foldtree in the pose//
		pose.fold_tree( f );

//		rigid::RigidBodyPerturbMoverOP rb_mover = new rigid::RigidBodyPerturbMover( pep_jump, 5.0, 0.0 );
//		rb_mover->apply( pose );

		//remove pep_anchor termini
		if( pep_begin == pep_anchor && pep_end == pep_anchor ){
			std::string pep_anchor_type( pose.residue( pep_anchor ).name3() );
			Residue rsd1( pose.residue( prot_anchor ) );
			Residue rsd2( pose.residue( pep_anchor ) );
			StubID upstubid(  AtomID( rsd1.atom_index( "N" ), prot_anchor ) ,
					AtomID( rsd1.atom_index( "CA" ), prot_anchor ) ,
					AtomID( rsd1.atom_index( "C" ), prot_anchor ) ) ;

			StubID downstubid(  AtomID( rsd2.atom_index( "N" ), pep_anchor ) ,
					AtomID( rsd2.atom_index( "CA" ), pep_anchor ) ,
					AtomID( rsd2.atom_index( "C" ), pep_anchor ) ) ;
			RT const rt( pose.conformation().get_stub_transform( upstubid, downstubid ) );

			pose = pose.split_by_chain( prot_chain );
			ResidueOP pep_anchor_res_ptr( ResidueFactory::create_residue( pose.residue( 1 ).residue_type_set().name_map( pep_anchor_type ) ) );
			pose.append_residue_by_jump( *pep_anchor_res_ptr, prot_anchor, "", "", true );
			pose.conformation().set_stub_transform( upstubid, downstubid, rt );
			for( Size i_chi = 1; i_chi <= pose.residue( pep_anchor ).nchi(); ++i_chi ){
				pose.set_chi( i_chi, pep_anchor, start_pose.residue( pep_anchor ).chi( i_chi ) );
			}
		}
		Residue pep_anchor_res( pose.residue( pep_anchor ) );

		ResidueOP gly( ResidueFactory::create_residue( rsd_set.name_map( "GLY" ) ) );
		for( Size ii = 1; ii <= option[ pep_spec::n_append ]; ++ii ){
			pose.append_polymer_residue_after_seqpos( *gly, pep_end, true );
			pep_end = pep_end + 1;
			pose.set_omega( pep_end - 1, 180.0 );
			pose.conformation().update_polymeric_connection( pep_end );
			pose.conformation().update_polymeric_connection( pep_end - 1 );
		}
		for( Size ii = 1; ii <= option[ pep_spec::n_prepend ]; ++ii ){
			pose.prepend_polymer_residue_before_seqpos( *gly, pep_begin, true );
			pep_end = pep_end + 1;
			pep_anchor = pep_anchor + 1;
			pose.set_omega( pep_begin, 180.0 );
			pose.conformation().update_polymeric_connection( pep_begin );
			pose.conformation().update_polymeric_connection( pep_begin + 1 );
		}

		//use input sequence//
		if( option[ pep_spec::use_input_seq ] ){
			std::string const input_seq( option[ pep_spec::input_seq ] );
			for( Size mut_site = pep_begin; mut_site <= pep_end; mut_site++ ){
				if(mut_site==pep_anchor) continue;
				chemical::AA input_aa( aa_from_oneletter_code( input_seq[ mut_site - pep_begin ] ) );
				chemical::make_sequence_change( mut_site, input_aa, pose );
			}
		}
		//or erase identity//
		else if( !( option[ pep_spec::test_no_design ] ) ){
			for( Size mut_site = pep_begin; mut_site <= pep_end; mut_site++ ){
				if(mut_site==pep_anchor) continue;
				chemical::make_sequence_change( mut_site, chemical::aa_ala, pose );
			}
		}
		//convert to CG residues//
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

		//add buffer res
		bool add_nterm( true );
		if( pep_anchor == pep_begin ) add_nterm = false;
		bool add_cterm( true );
		if( pep_anchor == pep_end ) add_cterm = false;
		if( option[ pep_spec::add_buffer_res ] ) add_pep_res( pose, prot_anchor, pep_begin, pep_anchor, pep_end, add_nterm, add_cterm );

		if( option[ pep_spec::cen_bb_frag ] ) gen_pep_bb_frag( pose, pep_begin, pep_anchor, pep_end, cen_scorefxn, lib );
		else if( option[ pep_spec::cen_bb_rama ] ) gen_pep_bb_rama( pose, pep_begin, pep_end, cen_scorefxn );

		//define prot, pep, and anchor residues//
		vector1< bool > is_pep( pose.total_residue(), false ), is_pep_not_anchor( pose.total_residue(), false ), is_prot( pose.total_residue(), false ), is_anchor( pose.total_residue(), false );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			is_pep[i] = ( i >= pep_begin && i <= pep_end );
			is_prot[i] = ( i >= prot_begin && i <= prot_end );
			is_anchor[i] = ( ( i == pep_anchor || i == prot_anchor ) );
			is_pep_not_anchor[i] = ( i >= pep_begin && i <= pep_end && !(i==pep_anchor) );
		}

		vector1< bool > is_flex_prot( pose.total_residue(), false );
		if( option[ pep_spec::flex_prot_bb ] ) gen_fold_tree_for_nbr_segments( pose, f, is_pep, is_anchor, cutoff, is_flex_prot );

		//make peptide from CG backbone//
		Pose restart_cgrelax_pose( pose );
		Size restart_pep_begin( pep_begin );
		Size restart_pep_anchor( pep_anchor );
		Size restart_pep_end( pep_end );
		for(Size cgrelax_loop = 1; cgrelax_loop <= n_cgrelax_loop; cgrelax_loop++){
			pose = restart_cgrelax_pose;
			pep_begin = restart_pep_begin;
			pep_anchor = restart_pep_anchor;
			pep_end = restart_pep_end;

			//min movemap//
			kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
			mm_min->set_chi( is_pep );
			if( !option[ pep_spec::test_no_min_bb ] ){
				mm_min->set_bb( is_pep );
			}
			if( option[ pep_spec::flex_prot_bb ] ) mm_min->set_chi( is_flex_prot );
			if( option[ pep_spec::flex_prot_bb ] ) mm_min->set_bb( is_flex_prot );

			//small-shear movemap//
			kinematics::MoveMapOP mm_move ( new kinematics::MoveMap );
			if( !option[ pep_spec::test_no_min_bb ] ){
				mm_move->set_bb( is_pep );
				if( option[ pep_spec::flex_prot_bb ] ) mm_move->set_bb( is_flex_prot );
			}

			//small CG backbone moves//
			Size n_cgperturb_iter( 1 );
			if( cgrelax_loop > 1 && cgrelax_loop <= 10 ) n_cgperturb_iter = 10;
			else n_cgperturb_iter = 100;
			if( !option[ pep_spec::test_no_cgrelax ] ) perturb_pep_bb( pose, mm_move, cen_scorefxn, n_cgperturb_iter );

			//debug dump CG data and pdb//
			if( option[ pep_spec::dump_cg_bb ] ){
				//remove buffer res
				if( option[ pep_spec::add_buffer_res ] ) remove_pep_res( pose, prot_anchor, pep_begin, pep_anchor, pep_end, add_nterm, add_cterm );

				if( option[ pep_spec::test_dump_all ] ) dump_efactor_pdb( pose, cen_scorefxn, pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + "_cen.pdb" );
				out_file << out_nametag + "_" + string_of( peptide_loop ) + "_cen.pdb" +
					"\t"<<pose.energies().total_energies().weighted_string_of( cen_scorefxn->weights() );
				out_file<<"\ttotal_score:\t"<<pose.energies().total_energies()[ total_score ]<<"\t";
				if( option[ pep_spec::rmsd_analysis ] ) out_file << pep_rmsd_analysis( pose, prot_begin, prot_end, pep_begin, pep_anchor, pep_end );
				if( option[ pep_spec::phipsi_analysis ] ) out_file << pep_phipsi_analysis( pose, pep_begin, pep_anchor, pep_end );

				out_file<<std::endl;

				continue;
			}

			//switch back to fullatom
			core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );

			//replace prot residues w/ original rotamers
			for(Size resnum = prot_begin; resnum <= prot_end; ++resnum){
				pose.replace_residue( resnum, start_pose.residue( resnum ), false );
			}
			pose.replace_residue( pep_anchor, pep_anchor_res, true );

			//replace termini
			add_termini( pose, pep_begin, pep_end );

			//turn on constraints//
			using namespace core::scoring::constraints;
			if( option[ pep_spec::constrain_pep_anchor ] ){
				soft_scorefxn->set_weight( atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ] );
				full_scorefxn->set_weight( atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ] );

				ConstraintCOP cst1( new AtomPairConstraint( AtomID( pose.residue( prot_anchor ).atom_index( "N" ), prot_anchor ), AtomID( pose.residue( pep_anchor ).atom_index( "N" ), pep_anchor ), new HarmonicFunc( pose.residue( prot_anchor ).xyz( "N" ).distance( pose.residue( pep_anchor ).xyz( "N" ) ), 0.1 ) ) );
				ConstraintCOP cst2( new AtomPairConstraint( AtomID( pose.residue( prot_anchor ).atom_index( "CA" ), prot_anchor ), AtomID( pose.residue( pep_anchor ).atom_index( "CA" ), pep_anchor ), new HarmonicFunc( pose.residue( prot_anchor ).xyz( "CA" ).distance( pose.residue( pep_anchor ).xyz( "CA" ) ), 0.1 ) ) );
				ConstraintCOP cst3( new AtomPairConstraint( AtomID( pose.residue( prot_anchor ).atom_index( "C" ), prot_anchor ), AtomID( pose.residue( pep_anchor ).atom_index( "C" ), pep_anchor ), new HarmonicFunc( pose.residue( prot_anchor ).xyz( "C" ).distance( pose.residue( pep_anchor ).xyz( "C" ) ), 0.1 ) ) );
				pose.add_constraint( cst1 );
				pose.add_constraint( cst2 );
				pose.add_constraint( cst3 );
			}

			//define neighbors
			vector1< bool > is_pep_nbr( pose.total_residue(), false );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const & rsd1( pose.residue(i) );
				if ( is_pep[i] ) continue;
				for ( Size j=1; j<= pose.total_residue(); ++j ) {
					Residue const & rsd2( pose.residue(j) );
					if ( !is_pep[j] ) continue;
					if( is_pep_nbr[i] ) break;
					for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
						if( is_pep_nbr[i] ) break;
						for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
							if ( rsd1.xyz(ii).distance( rsd2.xyz(jj) ) < cutoff ) {
								is_pep_nbr[i] = true;
								break;
							}
						}
					}
				}
			}
			mm_min->set_chi( is_pep_nbr );

			//poly-G FA refinement//
			if( !( option[ pep_spec::test_no_design ] || option[ pep_spec::test_no_farelax ] || option[ pep_spec::use_input_seq ] ) ){
				refine_fa_pep_bb( pose,	pep_anchor, is_pep, is_pep_nbr, full_scorefxn );
			}

			//replace prot residues w/ original rotamers
			for(Size resnum = prot_begin; resnum <= prot_end; ++resnum){
				pose.replace_residue( resnum, start_pose.residue( resnum ), false );
			}

			//allow jump relaxation//
			if( option[ pep_spec::constrain_pep_anchor ] ) mm_min->set_jump( 1, true );

			//define movers//
			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "dfpmin", 0.001, true );

			protocols::simple_moves::ShearMoverOP shear_mover( new protocols::simple_moves::ShearMover( mm_move, 5.0, 1 ) );	//LOOP
			shear_mover->angle_max( 'H', 5.0 );
			shear_mover->angle_max( 'E', 5.0 );
			shear_mover->angle_max( 'L', 5.0 );

			//define design task and repack task
			pack::task::RestrictResidueToRepackingOperationOP restrict_to_repack_taskop( new pack::task::RestrictResidueToRepackingOperation() );
			pack::task::PreventRepackingOperationOP prevent_repack_taskop( new pack::task::PreventRepackingOperation() );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( is_pep[ i ] && i != pep_anchor ) {
					if( ( option[ pep_spec::test_no_design ] || option[ pep_spec::use_input_seq ] ) ){
						restrict_to_repack_taskop->include_residue( i );
					}
				} else if ( is_pep_nbr[i] && is_prot[i] ) {
					restrict_to_repack_taskop->include_residue( i );
				} else {
					prevent_repack_taskop->include_residue( i );
				}
			}
			pack::task::TaskFactoryOP rottrial_task_factory( new pack::task::TaskFactory );
			rottrial_task_factory->push_back( new pack::task::InitializeFromCommandlineOperation() );
			rottrial_task_factory->push_back( new pack::task::IncludeCurrentOperation() );
			rottrial_task_factory->push_back( restrict_to_repack_taskop );
			rottrial_task_factory->push_back( prevent_repack_taskop );
			pack::task::PackerTaskOP dz_task( rottrial_task_factory->create_task_and_apply_taskoperations( pose ));

			protocols::simple_moves::PackRotamersMoverOP dz_pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, dz_task, 1 ) );
			protocols::simple_moves::RotamerTrialsMoverOP dz_rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, rottrial_task_factory ) );
			SequenceMoverOP design_seq = new SequenceMover;
			if( !option[ pep_spec::test_no_pack ] ){
				design_seq->add_mover( dz_pack );
				design_seq->add_mover( dz_rottrial );
			}
			if( !option[ pep_spec::test_no_min ] ){
				design_seq->add_mover( min_mover );
			}

			SequenceMoverOP relax_seq = new SequenceMover;
			relax_seq->add_mover( shear_mover );
			relax_seq->add_mover( dz_rottrial );
			if( !option[ pep_spec::test_no_min ] ){
				relax_seq->add_mover( min_mover );
			}

			//RUN MOVERS//
			( *soft_scorefxn )( pose );
			design_seq->apply( pose );

			//remove buffer res
			if( option[ pep_spec::add_buffer_res ] ) remove_pep_res( pose, prot_anchor, pep_begin, pep_anchor, pep_end, add_nterm, add_cterm );

/*
			MonteCarloOP mc_relax ( new MonteCarlo( pose, *full_scorefxn, 1.0 ) );
			TrialMoverOP relax_trial = new TrialMover( relax_seq, mc_relax );
			RepeatMoverOP relax_cycle = new RepeatMover( relax_trial, n_farelax_loop );
			if( !option[ pep_spec::test_no_farelax ] ){
				relax_cycle->apply( pose );
			}
			mc_relax->recover_low( pose );
*/

			//must remove constraints!
			pose.constraint_set( 0 );

			( *full_scorefxn )( pose );

			//Analysis//
			output_seq.clear();
			for(Size i = pep_begin; i <= pep_end; i++){
				output_seq.append( 1, pose.residue( i ).name1() );
			}
			Real total_score( pose.energies().total_energies().dot( full_scorefxn->weights() ) );

			//score pep unbound state//
			Real pep_score( 0 );
			Real bind_score( 0 );
			Pose pep_pose( pose.split_by_chain( pep_chain ) );
			if( option[ pep_spec::score_binding ] ){
				packmin_unbound_pep( pep_pose, full_scorefxn, soft_scorefxn );
				( *full_scorefxn )( pep_pose );
				pep_score = pep_pose.energies().total_energies().dot( full_scorefxn->weights() );
				bind_score = total_score - pep_score;
			}

			//Output//
			std::string pdb_name( pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + ".pdb" );
			if( n_cgrelax_loop > 1 ) pdb_name = pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + "_" + string_of( cgrelax_loop ) + ".pdb";

			out_file << pdb_name + "\t" << output_seq << "\t" << pose.energies().total_energies().weighted_string_of( full_scorefxn->weights() );
			out_file<<"\ttotal_score:\t"<<total_score<<"\t";
			if( option[ pep_spec::score_binding ] ) out_file<<"binding_score:\t"<<bind_score<<"\t";
			if( option[ pep_spec::rmsd_analysis ] ) out_file << pep_rmsd_analysis( pose, prot_begin, prot_end, pep_begin, pep_anchor, pep_end );
			if( option[ pep_spec::phipsi_analysis ] ) out_file << pep_phipsi_analysis( pose, pep_begin, pep_anchor, pep_end );
			out_file<<std::endl;

			if( option[ pep_spec::test_dump_all ] ) dump_efactor_pdb( pose, full_scorefxn, pdb_name );

/*
			if( option[ pep_spec::score_binding ] ){
				std::cout<<"!!pep "+ out_nametag + "_" + string_of( peptide_loop ) + "_" + string_of( cgrelax_loop ) + "_pep.pdb" + "\t"<<output_seq<<"\t"<<pep_pose.energies().total_energies().weighted_string_of( full_scorefxn->weights() );
			std::cout<<"\ttotal_score:\t"<<pep_score<<std::endl;
				if( option[ pep_spec::test_dump_all ] ) dump_pdb( pep_pose, out_nametag + "_" + string_of( peptide_loop ) + "_" + string_of( cgrelax_loop ) + "_pep.pdb" );
			}
*/
		}
	}
}


void*
my_main( void*)
{

	RunPepSpec();
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
	}

	return 0;
}
