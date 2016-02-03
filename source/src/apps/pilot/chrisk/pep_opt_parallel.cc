// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
 #include <protocols/frags/VallData.hh>
 #include <protocols/frags/TorsionFragment.hh>
 #include <numeric/random/random_permutation.hh>

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
 #include <core/chemical/ResidueTypeSelector.hh>
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

#include <core/io/pdb/pdb_writer.hh>

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

#include <mpi.h>

//Auto Headers
#include <core/import_pose/import_pose.hh>
 using basic::T;
 using basic::Error;
 using basic::Warning;


using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
//using namespace protocols;

using utility::vector1;
using std::string;
using io::pdb::old_dump_pdb;

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

void
normalize_res_freqs(
		vector1< Real > & res_freqs,
		Size n
)
{
	for( Size i = 1; i <= res_freqs.size(); ++i ){
		res_freqs[ i ] = res_freqs[ i ] / n;
	}
}

void
minimize_ref_weights(
	vector1< Real > & ref_weights,
	vector1< Real > res_freqs,
	vector1< Real > ideal_res_freqs,
	Real stepsize
)
{
	Real tol( 0.005 );
	// if majority are low, only change wts for aa's that are high, vice versa
	Size n_under( 0 );
	for( Size i = 1; i <= ref_weights.size(); ++i ){
		if( res_freqs[ i ] - ideal_res_freqs[ i ] <= 0.000 ) ++n_under;
	}
	// move wts up/down
	Real stepsize_max( 0.3 );
	Real stepsize_min( 0.01 );
	for( Size i = 1; i <= ref_weights.size(); ++i ){
		Real d_freq( res_freqs[ i ] - ideal_res_freqs[ i ] );
		if( d_freq < tol && n_under > ( ref_weights.size() / 2.0 ) ) continue;
		else if( d_freq > tol && n_under < ( ref_weights.size() / 2.0 ) ) continue;

		stepsize = ( res_freqs[ i ] < 0.001 ? stepsize_max : stepsize * std::abs( d_freq ) );
		if( stepsize > stepsize_max ) stepsize = stepsize_max;
		if( stepsize < stepsize_min ) stepsize = stepsize_min;
		if( d_freq > tol ) ref_weights[ i ] += stepsize;
		else if( d_freq < tol ) ref_weights[ i ] -= stepsize;
	}
}

void
print_res_freqs(
		vector1< Real > res_freqs,
		vector1< Real > ideal_res_freqs,
		std::string res_freqs_file_name_str
)
{
	char const *res_freqs_file_name = res_freqs_file_name_str.c_str();
	std::fstream res_freqs_file( res_freqs_file_name, std::ios::out );

	Real rmsd( 0 );
	for( Size i = 1; i <= res_freqs.size(); ++i ){
		res_freqs_file << chemical::oneletter_code_from_aa( chemical::AA( i ) ) << "\t";
		res_freqs_file << string_of( res_freqs[ i ] ) << "\n";
		rmsd += ( std::pow( res_freqs[ i ] - ideal_res_freqs[ i ], 2 ) );
	}
	rmsd = std::sqrt( rmsd / res_freqs.size() );
	res_freqs_file << "RMSD\t" << string_of( rmsd ) << "\n";
}

vector1< Real >
operator +( vector1< Real > v1, vector1< Real > v2 )
{
	vector1< Real > vsum( v1.size() );
	for( Size i = 1; i <= v1.size(); ++i ){
		vsum[ i ] = v1[ i ] + v2[ i ];
	}
	return vsum;
}

void
print_scorefxn(
	ScoreFunctionOP scorefxn,
	vector1< Real > ref_weights,
	std::string wts_file_name_str
)
{
	char const *wts_file_name = wts_file_name_str.c_str();
	std::fstream wts_file( wts_file_name, std::ios::out );

	wts_file << "ETABLE " << scorefxn->energy_method_options().etable_type() << "\n";
	wts_file << "METHOD_WEIGHTS ref  ";
	for( Size ii = 1; ii <= ref_weights.size(); ++ii ){
		wts_file << ref_weights[ ii ] << " ";
	}
	wts_file << "\n";
	for( Size ii = 1; ii <= core::scoring::n_score_types; ++ii ) {
		if ( scorefxn->get_weight( ScoreType( ii ) ) != 0 ) {
			wts_file << name_from_score_type( ScoreType( ii ) ) << " " << scorefxn->get_weight( ScoreType( ii ) ) << "\n";
		}
	}
}

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

void
send_string_to_node( int destination, std::string const & string_to_send )
{
	int tag( 1 );
	int len( string_to_send.size() );
	MPI_Send( &len, 1, MPI_INT, destination, tag, MPI_COMM_WORLD );
	MPI_Send( const_cast< char * > (string_to_send.c_str()), len, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
}

std::string
receive_string_from_node( int source )
{
	int len( 0 );
	int tag( 1 );
	MPI_Status stat;
	MPI_Recv( &len, 1, MPI_INT, source, tag, MPI_COMM_WORLD, & stat );
	char * str = new char[ len + 1 ];
	str[ len ] = '\0'; // ? do I need null terminated strings?
	MPI_Recv( str, len, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
	std::string return_string( str, len );
	delete [] str;
	return return_string;

}

void
RunPepSpec()
{

	int my_rank( 0 ), nprocs( 1 ), tag_( 1 );
	MPI_Status stat_;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);/* get current process id */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);/* get number of processes */

/////////////////shuffle all_pdb_filenames @ each iter so none get too far ahead?///////////////
	//split up pdbs
	vector1< std::string > all_all_pdb_filenames;
	Size n_pdbs_per_block( 0 );
	if( my_rank == 0 ){

		std::string pdb_list_filename( option[ pep_spec::pdb_list ] );
		std::ifstream pdb_list_data( pdb_list_filename.c_str() );
		if ( !pdb_list_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + pdb_list_filename + '\n' );
		}
		std::string pdb_list_line;
		while( !getline( pdb_list_data, pdb_list_line, '\n' ).eof() ) {
			std::string this_filename( pdb_list_line );
			all_all_pdb_filenames.push_back( this_filename );
		}
		n_pdbs_per_block = all_all_pdb_filenames.size();

		if( option[ pep_spec::n_opt_pdbs ].user() ) n_pdbs_per_block = option[ pep_spec::n_opt_pdbs ];
		if( n_pdbs_per_block > all_all_pdb_filenames.size() ) utility_exit_with_message( "n_opt_pdbs exceeds size of pdb_list!\n" );
	}

	//data out
	std::string out_nametag( "opt" );
	if( option[ out::file::o ].user() ) out_nametag = option[ out::file::o ];

	//define scoring functions//

	vector1< Real > ideal_res_freqs( 20, ( 1.0 / 20.0 ) );
	Real stepsize( 2.0 );

	for( Size i_opt = 1; i_opt <= option[ pep_spec::n_opt_loop ] + 1; ++i_opt ){
		MPI_Barrier( MPI_COMM_WORLD );
		vector1< std::string > all_pdb_filenames;
		vector1< std::string > my_pdb_filenames;

		if( my_rank == 0 ){

			numeric::random::random_permutation( all_all_pdb_filenames, numeric::random::rg() );
			all_pdb_filenames = all_all_pdb_filenames;
			all_pdb_filenames.resize( n_pdbs_per_block );

			Size const num_pdbs_per_cpu = all_pdb_filenames.size() / nprocs;
			Size const nextra = all_pdb_filenames.size() - num_pdbs_per_cpu * nprocs;

			Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_pdbs_per_cpu;
			std::cout << "CPU " << my_rank << " loading:\t";
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				my_pdb_filenames.push_back( all_pdb_filenames[ ii ] );
				std::cout << my_pdb_filenames[ ii ] << " ";
			}
			std::cout << "\n";

			Size ii_offset = my_njobs;
			for ( Size ii = 1; ii < nprocs; ++ii ) {
				Size ii_njobs = ( nextra > ii ? 1 : 0 ) + num_pdbs_per_cpu;
				MPI_Send( & ii_njobs, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD );
				for ( Size jj = ii_offset + 1; jj <= ii_offset + ii_njobs; ++jj ) {
					send_string_to_node( ii, all_pdb_filenames[ jj ] );
				}
				ii_offset += ii_njobs;
			}
		}
		else{
			Size my_njobs;
			MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );
			std::cout << "CPU " << my_rank << " loading:\t";
			for ( Size ii = 1; ii <= my_njobs; ++ii ) {
				my_pdb_filenames.push_back( receive_string_from_node( 0 ) );
				std::cout << my_pdb_filenames[ ii ] << " ";
			}
			std::cout << "\n";
		}

		std::string scorefxn_filename;
		if( i_opt == 1 ) scorefxn_filename = std::string( option[ pep_spec::wts ] );
		else scorefxn_filename = out_nametag + "." + string_of( i_opt - 1 ) + ".wts";
		core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( scorefxn_filename ) );

		vector1< Real > res_freqs( chemical::num_canonical_aas, 0.0 );
		Size n_positions( 0 );

		for( Size peptide_loop = 1; peptide_loop <= my_pdb_filenames.size(); ++peptide_loop ){

			Pose pose;
			std::string pdb_filename( my_pdb_filenames[ peptide_loop ] );
			core::import_pose::pose_from_file( pose, pdb_filename , core::import_pose::PDB_file);

			if( option[ pep_spec::pep_chain ].user() ){
				Size pep_anchor_in( option[ pep_spec::pep_anchor ] );
				std::string const pep_chain_in( option[ pep_spec::pep_chain ] );
				Size pep_anchor( pose.pdb_info()->pdb2pose( pep_chain_in[0], pep_anchor_in ) );
				Size pep_chain( pose.chain( pep_anchor ) );
				pose = pose.split_by_chain( pep_chain );
			}

			Size pep_begin( 1 );
			Size pep_end( pose.total_residue() );
			//ignore termini
			n_positions += ( pose.total_residue() - 2 );

			//gen fold tree//
			FoldTree f( pose.total_residue() );
			pose.fold_tree( f );

			//randomize pep sequence//
			for(Size mut_site = pep_begin; mut_site <= pep_end; mut_site++){ //over all pep positions
				int resindex;
				resindex = static_cast< int > ( 20 * numeric::random::rg().uniform() + 1 );
				chemical::make_sequence_change( mut_site, chemical::AA(resindex), pose );
			}

			//define design task and repack task
			pack::task::PackerTaskOP dz_task( pack::task::TaskFactory::create_packer_task( pose ));
			dz_task->initialize_from_command_line();

			pack::task::TaskFactoryOP rottrial_task_factory( new pack::task::TaskFactory );
			rottrial_task_factory->push_back( new pack::task::InitializeFromCommandlineOperation() );

			protocols::simple_moves::PackRotamersMoverOP dz_pack( new protocols::simple_moves::PackRotamersMover( scorefxn, dz_task, 1 ) );
			protocols::simple_moves::RotamerTrialsMoverOP dz_rottrial ( new protocols::simple_moves::RotamerTrialsMover( scorefxn, rottrial_task_factory ) );
			SequenceMoverOP design_seq = new SequenceMover;
			design_seq->add_mover( dz_pack );
			design_seq->add_mover( dz_rottrial );

			//RUN MOVERS//
			( *scorefxn )( pose );
			design_seq->apply( pose );

			//Analysis//
			std::string output_seq;
			for(Size i = pep_begin + 1; i <= pep_end - 1; i++){
				output_seq.push_back( pose.residue( i ).name1() );
				res_freqs[ aa2index( pose.residue( i ).aa() ) ] = res_freqs[ aa2index( pose.residue( i ).aa() ) ] + 1;
			}

			//Output//
			if( option[ pep_spec::test_dump_all ] ){
				std::cout << out_nametag + "_" + string_of( peptide_loop ) + ".pdb" + "\t" << output_seq << "\n";
//				dump_efactor_pdb( pose, scorefxn, "pdbs/" + out_nametag + "_" + string_of( peptide_loop ) + ".pdb" );
			}
		}

		if( my_rank == 0 ){

			vector1< Real > all_res_freqs( res_freqs );
			Size all_n_positions( n_positions );

			for ( Size ii = 1; ii < nprocs; ++ii ) {
				vector1< Real > ii_res_freqs;
				for ( Size jj = 1; jj <= all_res_freqs.size(); ++jj ) {
					Real jj_res_freq( 0 );
					MPI_Recv( & jj_res_freq, 1, MPI_UNSIGNED_LONG, ii, tag_, MPI_COMM_WORLD, & stat_ );
					all_res_freqs[ jj ] += jj_res_freq;
				}

				Size ii_n_positions;
				MPI_Recv( & ii_n_positions, 1, MPI_INT, ii, tag_, MPI_COMM_WORLD, & stat_ );
				all_n_positions += ii_n_positions;
			}

			normalize_res_freqs( all_res_freqs, all_n_positions  );
			std::string all_res_freqs_file_name( out_nametag + "." + string_of( i_opt - 1 ) + ".freq" );
			print_res_freqs( all_res_freqs, ideal_res_freqs, all_res_freqs_file_name );

			if( i_opt < option[ pep_spec::n_opt_loop ] + 1 ){
				vector1< Real > ref_weights( scorefxn->energy_method_options().method_weights( ref ) );
				minimize_ref_weights( ref_weights, all_res_freqs, ideal_res_freqs, stepsize );
				std::string wts_file_name( out_nametag + "." + string_of( i_opt ) + ".wts" );
				print_scorefxn( scorefxn, ref_weights, wts_file_name );
			}

		}
		else{
			vector1< Real > my_res_freqs( res_freqs );
			for ( Size jj = 1; jj <= my_res_freqs.size(); ++jj ) {
				Real jj_res_freq( my_res_freqs[ jj ] );
				MPI_Send( & jj_res_freq, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD );
			}

			Size my_n_positions( n_positions );
			MPI_Send( & my_n_positions, 1, MPI_INT, 0, tag_, MPI_COMM_WORLD );
		}

	}

	MPI_Finalize();

}


int main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		MPI_Init(&argc, &argv);

		devel::init(argc, argv);

		RunPepSpec();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
