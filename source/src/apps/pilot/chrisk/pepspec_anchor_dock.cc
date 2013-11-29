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

 #include <core/scoring/constraints/CoordinateConstraint.hh>
 #include <core/scoring/constraints/ConstraintIO.hh>
 #include <core/scoring/func/FlatHarmonicFunc.hh>

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

// #include <core/util/prof.hh> // profiling
// #include <core/util/CacheableData.hh> // profiling


 #include <core/chemical/AA.hh>
 #include <core/chemical/AtomTypeSet.hh>

 #include <core/chemical/AA.hh>
 #include <core/conformation/Conformation.hh>
 #include <core/conformation/Residue.hh>
 #include <core/conformation/ResidueMatcher.hh>
 #include <core/conformation/ResidueFactory.hh>
 #include <core/pack/rotamer_set/RotamerCouplings.hh>
 #include <core/chemical/ResidueType.hh>
 #include <core/chemical/ResidueTypeSet.hh>
 #include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
 #include <core/chemical/VariantType.hh>
 #include <core/chemical/util.hh>
 #include <core/chemical/ChemicalManager.hh>

 #include <core/scoring/rms_util.hh>
 #include <core/scoring/EnergyMap.hh>
 #include <core/scoring/Energies.hh>
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
#include <core/pack/task/operation/TaskOperations.hh>

 #include <core/kinematics/RT.hh>
 #include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
 #include <core/kinematics/util.hh>
 #include <core/id/AtomID_Map.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <basic/options/util.hh>//option.hh>
// #include <basic/options/after_opts.hh>

 #include <basic/basic.hh>

 #include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/import_pose/import_pose.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

 #include <numeric/conversions.hh>
 #include <numeric/xyzVector.hh>
 #include <numeric/random/random.hh>
 #include <numeric/random/random_permutation.hh>
 #include <numeric/NumericTraits.hh>
 #include <ObjexxFCL/string.functions.hh>

// //REMOVE LATER!
// #include <utility/io/izstream.hh>

#include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/L1ScoringScheme.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

// // C++ headers
 #include <cstdlib>
 #include <fstream>
 #include <iostream>
 #include <string>
 #include <sstream>

//silly using/typedef

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/pepspec.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

 static numeric::random::RandomGenerator RG(16621);

using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using utility::vector1;
using std::string;
using io::pdb::dump_pdb;

basic::Tracer TR("apps.pilot.chrisk/pep_prep");

//parse cst line as: name1 res1 x y z x0 sd tol
struct pep_coord_cst {
	std::string atom_name;
	int pep_pos;
	Real x;
	Real y;
	Real z;
	Real x0;
	Real sd;
	Real tol;
};
vector1< pep_coord_cst > pep_coord_csts;

//print homol_csts data
void
print_pep_csts(
	vector1< Vector > pep_cst_all_vects,
	vector1< Real > pep_cst_all_tols,
	vector1< Real > pep_cst_all_sds,
	vector1< Size > n_pep_cs,
	vector1< Size > pep_cst_all_indices,
	vector1< std::string > pep_cst_all_atom_names
)
{
	Size min_cas_for_cst( 3 );
	Real min_tol_for_cst( 1.0 );
	Real max_pep_cs( *( std::max_element( n_pep_cs.begin(), n_pep_cs.end() ) ) );
	std::string cst_file_name_str( option[ out::file::o ]+".cst" );
	char const *cst_file_name = cst_file_name_str.c_str();
	std::fstream cst_file( cst_file_name, std::ios::out );
	for( Size i = 1; i <= pep_cst_all_vects.size(); ++i ){
		if( n_pep_cs[ i ] < min_cas_for_cst && option[ pepspec::p_homol_csts ].user() ) continue;
		Real this_tol( pep_cst_all_tols[ i ] + pep_cst_all_tols[ i ] * ( 1.0 - std::pow( ( n_pep_cs[ i ] / max_pep_cs ), 1.0 / 3.0 ) ) );
		Vector pep_cst_vect( pep_cst_all_vects[ i ] );
		int pep_pos( pep_cst_all_indices[ i ] );
		Real pep_cst_tol( pep_cst_all_tols[ i ] );
//		if( pep_cst_tol < min_tol_for_cst ) pep_cst_tol = min_tol_for_cst;
		cst_file << pep_cst_all_atom_names[ i ] + "\t" + string_of( pep_pos ) + "\t" + string_of( pep_cst_vect.x() ) + "\t" + string_of( pep_cst_vect.y() ) + "\t" + string_of( pep_cst_vect.z() ) + "\t0.0\t" + string_of( pep_cst_all_sds[ i ] ) + "\t" + string_of( this_tol ) + "\n";
	}
}

void
set_pep_cst(
	pose::Pose & pose,
	Size pep_anchor,
	pep_coord_cst pep_cst
)
{
	using namespace core::scoring::constraints;
	using namespace core::chemical;

	Size seqpos( pep_cst.pep_pos + pep_anchor );
	Size prot_anchor( 2 );
	Vector pep_cst_vector;
	pep_cst_vector.x( pep_cst.x );
	pep_cst_vector.y( pep_cst.y );
	pep_cst_vector.z( pep_cst.z );
	ConstraintCOP this_cst( new CoordinateConstraint( id::AtomID( pose.residue( seqpos ).atom_index( pep_cst.atom_name ), seqpos ), id::AtomID( pose.residue( prot_anchor ).atom_index( "CA" ), prot_anchor ), pep_cst_vector, new FlatHarmonicFunc( pep_cst.x0, pep_cst.sd, pep_cst.tol ) ) );
	pose.add_constraint( this_cst );
}

/// @details  This function will make a sequence mutation while trying to preserve the variants
void
make_sequence_change(
		Size const seqpos,
		chemical::AA const & new_aa,
		pose::Pose & pose
		)
{
	Size const which_his_variant = 1;
	conformation::Residue const & current_rsd( pose.residue( seqpos ) );
	if ( current_rsd.aa() == new_aa ) return; // already done

	chemical::ResidueTypeCOPs rsd_types
		( chemical::ResidueSelector().set_aa( new_aa ).match_variants( current_rsd.type() ).select( current_rsd.residue_type_set() ) );

	Size rsd_types_index( 1 );
	std::string const errmsg
		( "make_sequence_change failed: new_aa= "+chemical::name_from_aa(new_aa)+" rsd_types.size()= "+string_of( rsd_types.size() ) );

	if ( new_aa == chemical::aa_his ) {
		if ( rsd_types.size() != 2 || which_his_variant > 2 ) utility_exit_with_message( errmsg );
		rsd_types_index = which_his_variant;
	} else if ( rsd_types.size() != 1 ) {
		utility_exit_with_message( errmsg );
	}

	conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( *(rsd_types[ rsd_types_index ] ),
				current_rsd, pose.conformation() ) );
	pose.replace_residue( seqpos, *new_rsd, false );
}


bool
has_clash(
		pose::Pose pose,
		vector1< bool > check_seqpos,
		scoring::ScoreFunctionOP const & scorefxn,
		Real const clash_threshold,
		bool print_clash
	 )
{
	using namespace scoring;
	using namespace chemical;

	( *scorefxn )( pose );


	// cached energies object
	Energies & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph & energy_graph( energies.energy_graph() );
	bool is_clash( false );

	for( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ){
		if( !check_seqpos[ seqpos ] ) continue;
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
				if( print_clash ) TR<< "fa_rep: " << string_of( clash ) << " at " << pose.residue( j ).name1() << string_of( j ) << "-" << pose.residue( seqpos ).name1() << string_of( seqpos ) << std::endl;
				is_clash = true;
				return is_clash;
			}
		}

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
				if( print_clash ) TR<< "fa_rep: " << string_of( clash ) << " at " << pose.residue( j ).name1() << string_of( j ) << "-" << pose.residue( seqpos ).name1() << string_of( seqpos ) << std::endl;
				is_clash = true;
				return is_clash;
			}
		}
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
	std::string ref_name( option[ in::file::native ]() );
	import_pose::pose_from_pdb( ref_pose, ref_name );

	Size ref_pep_anchor( pep_anchor );
	if( option[ pepspec::native_pep_anchor ].user() ){
		Size ref_pep_anchor_in = option[ pepspec::native_pep_anchor ];
		std::string ref_pep_chain_in( option[ pepspec::pep_chain ] );
		if( option[ pepspec::native_pep_chain ].user() ) ref_pep_chain_in = option[ pepspec::native_pep_chain ];
		ref_pep_anchor = ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_in );
	}
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
	if( option[ pepspec::native_align ] ){
		id::AtomID_Map< id::AtomID > atom_map;
		pose::initialize_atomid_map( atom_map, pose, id::BOGUS_ATOM_ID );

		for ( Size i = prot_begin; i <= prot_end; ++i ) {
			id::AtomID const id1( pose.residue( i ).atom_index( "CA" ), i );
			id::AtomID const id2( ref_pose.residue( i + static_cast< int >( ref_prot_begin ) - static_cast< int >( prot_begin ) ).atom_index( "CA" ), i + static_cast< int >( ref_prot_begin ) - static_cast< int >( prot_begin ) );
			atom_map[ id1 ] = id2;
		}
		core::scoring::ScoreFunctionOP full_scorefxn(  core::scoring::getScoreFunction() );			
		core::scoring::superimpose_pose( pose, ref_pose, atom_map );

	}
	Real sd( 0.0 );
	std::string rmsd_analysis( "" );
	Size natoms( 0 );
	for( int i_seq = 0; i_seq <= nterm + cterm; ++i_seq ){
		Size ref_pep_seqpos = ref_pep_begin + i_seq;
		Size pep_seqpos = pep_begin + i_seq;

//		if( pose.residue( pep_seqpos ).natoms() != ref_pose.residue( ref_pep_seqpos ).natoms() ) utility_exit_with_message( "natoms mismatch for fa_rmsd\n" );
		for( Size atomno = 1; atomno <= pose.residue( pep_seqpos ).natoms(); ++atomno ){
				sd += ref_pose.residue( ref_pep_seqpos ).xyz( atomno ).distance_squared( pose.residue( pep_seqpos ).xyz( atomno ) );
				++natoms;
		}	

	}
	Real total_rmsd = std::sqrt( sd / natoms );
	rmsd_analysis += "fa_rmsd:\t" + string_of( total_rmsd ) + "\t";
	return rmsd_analysis;

}


Real
average(
	vector1< Real > real_vec
)
{
	Real avg( 0 );
	for( Size i = 1; i <= real_vec.size(); ++i ){
		avg += real_vec[ i ] / real_vec.size();
	}
	return avg;
}

//function for getting array of radian angles ready for averaging
vector1< Real >
shift_angles(
	vector1< Real > & angles
)
{
	const Real pi = numeric::NumericTraits<Real>::pi();
	Size n_angles( angles.size() );

	vector1< Real > pos_angles; 
	vector1< Real > neg_angles; 
	for( Size i = 1; i <= n_angles; ++i ){ 
		if( angles[ i ] >=0 ) pos_angles.push_back( angles[ i ] );
		else neg_angles.push_back( angles[ i ] );
	}
	Size n_pos_angles( pos_angles.size() );
	Size n_neg_angles( neg_angles.size() );
	if( n_pos_angles >= n_neg_angles ){
		Real avg_pos_angle( average( pos_angles ) );
		for( Size i = 1; i <= n_angles; ++i ){ 
			if( angles[ i ] < 0  && abs( angles[ i ] - avg_pos_angle ) > pi ) angles[ i ] += ( 2 * pi );
		}
	}
	else{
		Real avg_neg_angle( average( neg_angles ) );
		for( Size i = 1; i <= n_angles; ++i ){ 
			if( angles[ i ] >= 0  && abs( angles[ i ] - avg_neg_angle ) > pi ) angles[ i ] -= ( 2 * pi );
		}
	}

	return angles;
}


void
run_pep_prep()
{
	using namespace pose;
	using namespace chemical;
	using namespace conformation;
	using namespace scoring;
	using namespace optimization;
	using namespace kinematics;
	using namespace sequence;

        using namespace protocols::moves;
	
        using namespace id;
        using namespace protocols::frags;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	typedef  numeric::xyzVector< Real >  Vector; // DOUBLE!
	typedef  numeric::xyzMatrix< Real >  Matrix; // DOUBLE!

	std::string pdb_out_list_filename_str( option[ out::file::o ]+".pdblist" );
	char const *pdb_out_list_filename = pdb_out_list_filename_str.c_str();
	std::fstream pdb_out_list_file( pdb_out_list_filename, std::ios::out );
/*
	std::string cst_list_file_name_str( option[ out::file::o ]+".cstlist" );
	char const *cst_list_file_name = cst_list_file_name_str.c_str();
	std::fstream cst_list_file( cst_list_file_name, std::ios::out );
*/

	std::string out_nametag( "data" );
	if( option[ out::file::o ].user() ) out_nametag = option[ out::file::o ];
	std::string out_file_name_str( out_nametag + ".spec" );
	char const *out_file_name = out_file_name_str.c_str();
	std::fstream out_file( out_file_name, std::ios::out );

	vector1< string > ref_pdbs;
	vector1< string > ref_pep_chains;
	vector1< string > ref_pep_anchors_int;
	if( option[ pepspec::ref_pdb_list ].user() ){
		std::string ref_filename( option[ pepspec::ref_pdb_list ] );
		std::ifstream ref_file( ref_filename.c_str() );
		std::string line;
		std::string word;
		while( !getline( ref_file, line ).eof() ){
			std::istringstream line_parse( line );
			for( Size ii = 1; !getline( line_parse, word, '\t' ).eof(); ++ii ){
				if( ii == 1 ) ref_pdbs.push_back( word );
				else if( ii == 2 ) ref_pep_chains.push_back( word );
				else utility_exit_with_message( "ref_pdb_list has too many args @ line \"" + line + "\"\n" );
			}
			getline( line_parse, word );
			ref_pep_anchors_int.push_back( word );
		}

	}
	else utility_exit_with_message( "ref_pdb_list not defined!\n" );
	if( ref_pdbs.size() < 1 ) utility_exit_with_message( "no ref pdbs loaded!\n" );
	if( ref_pep_chains.size() < 1 ) utility_exit_with_message( "no ref chains loaded!\n" );

	Size n_ref_poses( ref_pdbs.size() );

	vector1< Pose > ref_poses( n_ref_poses );
	vector1< Size > ref_pep_anchors( n_ref_poses );
	vector1< Vector > jump_trans_vector( n_ref_poses );
	vector1< Matrix > jump_rot_matrix( n_ref_poses );
	vector1< Real > jump_bond_angle( n_ref_poses );
	vector1< Real > jump_tor1_angle( n_ref_poses );
	vector1< Real > jump_tor2_angle( n_ref_poses );

	std::string prot_anchor_stub1( "N" );
	std::string prot_anchor_stub2( "CA" );
	std::string prot_anchor_stub3( "C" );

	std::string pep_anchor_stub1( "CB" );
	std::string pep_anchor_stub2( "CA" );
	std::string pep_anchor_stub3( "N" );

	/////////////////////////////////////////////

	vector1< std::string > pdb_filenames;
	if( option[ pepspec::pdb_list ].user() ){
		std::string pdb_list_filename( option[ pepspec::pdb_list ] );
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
	Pose pose;
	import_pose::pose_from_pdb( pose, pdb_filenames[ 1 ] );

	ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );

	Size const nres( pose.total_residue() );

	Size prot_chain;
	for( Size i = 1; i <= pose.conformation().num_chains(); ++i ){
		if( !( pose.residue( pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			prot_chain = i;
			break;
		}
	}
	if( option[ pepspec::remove_input_bb ] ) pose = pose.split_by_chain( prot_chain );
	Size prot_begin( pose.conformation().chain_begin( prot_chain ) ); 
	Size prot_end( pose.conformation().chain_end( prot_chain ) ); 
	Size prot_anchor( prot_begin + 1 );


	// create a simple foldtree of the right size
	FoldTree f( pose.total_residue() );

	//gen fold tree
	if( prot_begin != 1 ) f.new_jump( prot_anchor, 1, 1 );
	if( prot_end != pose.total_residue() ) f.new_jump( prot_anchor, prot_end + 1, prot_end );

	core::scoring::ScoreFunctionOP scorefxn(  getScoreFunction() );
	core::scoring::ScoreFunctionOP soft_scorefxn(  ScoreFunctionFactory::create_score_function( option[ pepspec::soft_wts ] ) );

	//for each desired anchor-docked protein
	Size n_peptides( option[ pepspec::n_peptides ] );
	int pose_index( 0 );
	if( option[ pepspec::run_sequential ] ) n_peptides = pdb_filenames.size();
	for( Size peptide_loop = 1; peptide_loop <= n_peptides; ++peptide_loop ){

		if( option[ pepspec::run_sequential ] ) ++pose_index;
		else pose_index = static_cast< int >( RG.uniform() * pdb_filenames.size() + 1 );
		std::string pdb_filename( pdb_filenames[ pose_index ] );
		TR<<"Initializing "<< "prep_" + string_of( peptide_loop ) + ".pdb with " + pdb_filename << std::endl;
		import_pose::pose_from_pdb( pose, pdb_filename );

		// set the new foldtree in the pose
		pose.fold_tree( f );
		
		//prepack? 
		if( !option[ pepspec::no_prepack_prot ] ){
			( *soft_scorefxn )( pose );

			pack::task::PackerTaskOP prepack_task( pack::task::TaskFactory::create_packer_task( pose ));
			prepack_task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
			protocols::simple_moves::PackRotamersMoverOP prepack_mover( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, prepack_task, 1 ) );
			prepack_mover->apply( pose );

			pack::task::TaskFactoryOP prepack_task_factory( new pack::task::TaskFactory );
			prepack_task_factory->push_back( new pack::task::operation::InitializeFromCommandline() );
			prepack_task_factory->push_back( new pack::task::operation::IncludeCurrent() );
			prepack_task_factory->push_back( new pack::task::operation::RestrictToRepacking() );
			protocols::simple_moves::RotamerTrialsMoverOP prepack_rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, prepack_task_factory ) );
			prepack_rottrial->apply( pose );

			kinematics::MoveMapOP mm_prepack ( new kinematics::MoveMap );
			mm_prepack->set_chi( true );
			protocols::simple_moves::MinMoverOP prepack_min_mover = new protocols::simple_moves::MinMover( mm_prepack, scorefxn, "dfpmin", 0.001, true );
			prepack_min_mover->apply( pose );

			dump_pdb( pose, pdb_filename + ".prepack" );
		}
		( *scorefxn )( pose );
		Real prot_score( pose.energies().total_energies().dot( scorefxn->weights() ) );

		//align to other pre-aligned struct of same seq
		if( option[ pepspec::prep_align_prot_to ].user() ){
			std::string align_name( option[ pepspec::prep_align_prot_to ] );
			Pose align_pose;
			import_pose::pose_from_pdb( align_pose, align_name );
			id::AtomID_Map< id::AtomID > atom_map;
			pose::initialize_atomid_map( atom_map, pose, id::BOGUS_ATOM_ID );
			for ( Size i = prot_begin; i <= prot_end; ++i ) {
				id::AtomID const id1( pose.residue( i ).atom_index( "CA" ), i );
				id::AtomID const id2( align_pose.residue( i ).atom_index( "CA" ), i );
				atom_map[ id1 ] = id2;
			}
			superimpose_pose( pose, align_pose, atom_map );
		}

		// add anchor in random pos
		std::string const anchor_type( option[ pepspec::anchor_type ] );
		ResidueOP pep_anchor_res( ResidueFactory::create_residue( rsd_set.name_map( anchor_type ) ) );
		pose.append_residue_by_jump( *pep_anchor_res, prot_anchor, "CA", "CA", true );
		Size pep_anchor( pose.total_residue() );
		Size pep_jump( 1 );
		Pose start_pose( pose );
		/*
		   if( anchor_type == "PRO" ){
		   pep_anchor_stub2 = "C";
		   pep_anchor_stub3 = "O";
		   }
		 */
		FoldTree f_jump( pose.total_residue() );
		FoldTree f_orient( pose.total_residue() );

		f_jump.new_jump( prot_anchor, pep_anchor, prot_end );
		if( prot_end != pose.total_residue() - 1 ){
			f_jump.new_jump( prot_anchor, prot_end + 1, prot_end );
			f_orient.new_jump( prot_anchor, prot_end + 1, prot_end );
		}

		f_orient.new_chemical_bond( prot_anchor, pep_anchor, prot_anchor_stub2, pep_anchor_stub1 , prot_end );

		f_jump.reorder( prot_anchor );
		f_jump.set_jump_atoms( 1, prot_anchor_stub2, pep_anchor_stub1  );


		//homol_csts
		//calc vect sum, then normalize to avgs, then reiter over refs to get max distance
		//pep_ca_vects is peptides< seqpos< pep_pos, Vector > >
		vector1< vector1< std::pair < int, Vector > > > pep_ca_vects;
		vector1< vector1< std::pair < int, Vector > > > pep_cb_vects;
		int ref_pep_pos_min( 0 );
		int ref_pep_pos_max( 0 );

		/////////////////////////////////////////////////////////////////////////////////////////
		// copy anchor xyz coords from ref pdbs to pose and calc/store kinematic orientation
		bool pep_anchor_is_nterm( false );
		for( Size i = 1; i <= n_ref_poses; ++i ){

			//START CUT OUT, put this part of the ref prot loop outside the gen docked target loop

			//load ref pdb name and chain, anchor data
			std::string ref_input_name( ref_pdbs[ i ] );
			std::string const ref_pep_chain_in( ref_pep_chains[ i ] );
			std::istringstream ref_pep_anchor_in( ref_pep_anchors_int[ i ] );
			int ref_pep_anchor_int;
			ref_pep_anchor_in >> ref_pep_anchor_int;
			if( peptide_loop == 1 ) TR << "Loading " << ref_input_name << ": pep anchor @ chain " << ref_pep_chain_in << " seqpos " << string_of( ref_pep_anchor_int ) << std::endl;

			//load ref pdb
			Pose ref_pose;
			import_pose::pose_from_pdb( ref_pose, ref_input_name );

			Size const ref_pep_anchor( ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_int ) );
			//if gly, just change to ala, orient and calc data, then change back
			//can we just use one of the Halphas?
			//if( ref_pose.residue( ref_pep_anchor ).name3() == "GLY" ) utility_exit_with_message( "Cannot yet handle GLY anchors\n" );
			if( ref_pose.residue( ref_pep_anchor ).name3() == "GLY" ){ 
					make_sequence_change( ref_pep_anchor, chemical::aa_ala, ref_pose );
			}
			//store ref pose, anchor
			ref_poses[ i ] = ref_pose;
			ref_pep_anchors[ i ] = ref_pep_anchor;
			Size const ref_pep_chain( ref_pose.chain( ref_pep_anchor ) );
			Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pose.chain( ref_pep_anchor ) ) );
			Size ref_pep_end( ref_pose.conformation().chain_end( ref_pose.chain( ref_pep_anchor ) ) );
			//set ref prot chain  = first non-peptide protein chain in pose
			Size ref_prot_chain;
			for( Size ii = 1; ii <= ref_pose.conformation().num_chains(); ++ii ){
				if( ii == ref_pep_chain ) continue;
				if( !( ref_pose.residue( ref_pose.conformation().chain_begin( ii ) ).is_protein() ) ) continue;
				else{
					ref_prot_chain = ii;
					break;
				}
			}
			Size ref_prot_begin( ref_pose.conformation().chain_begin( ref_prot_chain ) ); 

			//sequence alignment of ref prot and target?
			if( option[ pepspec::seq_align ] ){
				utility::file::FileName blosum62( "/home/cking/blosum/BLOSUM62" );
				ScoringSchemeOP ss( new MatrixScoringScheme( -11, -1, blosum62 ) );
				NWAligner nw_aligner;
				SequenceOP seq1( new Sequence( ref_pose.chain_sequence( ref_prot_chain ), "ref", ref_prot_begin ) );
				SequenceOP seq2( new Sequence( pose.chain_sequence( prot_chain ), "target", prot_begin ) );
				SequenceAlignment global_align = nw_aligner.align( seq1, seq2, ss );
				SequenceMapping seq_map( global_align.sequence_mapping( 1, 2 ) );
				//seq_map.show( TR );
				vector1< Size > seq_align( seq_map.mapping() );
				id::AtomID_Map< id::AtomID > atom_map;
				pose::initialize_atomid_map( atom_map, pose, id::BOGUS_ATOM_ID );
				Real total_rmsd( 0 );
				Size npos( 0 );
				for( Size ii = 1; ii <= seq_align.size(); ++ii ){
					if( seq_align[ ii ] == 0 ) continue;
					id::AtomID const id1( ref_pose.residue( ii ).atom_index( "CA" ), ii );
					id::AtomID const id2( pose.residue( seq_align[ ii ] ).atom_index( "CA" ), seq_align[ ii ] );
					atom_map[ id1 ] = id2;
					total_rmsd += ref_pose.residue( ii ).xyz( "CA" ).distance_squared( pose.residue( seq_align[ ii ] ).xyz( "CA" ) );
					++npos;
				}
				total_rmsd = std::sqrt( total_rmsd / npos );
				TR << string_of( total_rmsd ) + " Ca RMSD over " + string_of( npos ) + " atoms\n";
				if( npos < pose.total_residue() / 2 ) utility_exit_with_message( "Alignment falied at " + ref_input_name + "\n" );
				//superimpose ref, target with atom map from seq alignment
				superimpose_pose( ref_pose, pose, atom_map );
			}

			//store calpha vectors for csts
			if( option[ pepspec::homol_csts ].user() && peptide_loop == 1 ){
				for( Size i_mut = ref_pep_begin; i_mut <= ref_pep_end; ++i_mut ){
					if( i_mut == ref_pep_anchor ) continue;
					make_sequence_change( i_mut, chemical::aa_ala, ref_pose );
				}
				vector1< std::pair< int, Vector > > pep_ca_vect;
				vector1< std::pair< int, Vector > > pep_cb_vect;
				for( Size ii = ref_pep_begin; ii <= ref_pep_end; ++ii ){
					int ref_pep_pos( ii - ref_pep_anchor );
					if( ref_pep_pos < ref_pep_pos_min ) ref_pep_pos_min = ref_pep_pos;
					if( ref_pep_pos > ref_pep_pos_max ) ref_pep_pos_max = ref_pep_pos;
					Vector pep_ca_coord( ref_pose.residue( ii ).xyz( "CA" ) );
					Vector pep_cb_coord( ref_pose.residue( ii ).xyz( "CB" ) );
					std::pair< int, Vector > pep_ca_pair( ref_pep_pos, pep_ca_coord );
					std::pair< int, Vector > pep_cb_pair( ref_pep_pos, pep_cb_coord );
					pep_ca_vect.push_back( pep_ca_pair );
					pep_cb_vect.push_back( pep_cb_pair );
				}
				//store ca, cb vecotrs for constraints
				pep_ca_vects.push_back( pep_ca_vect );
				pep_cb_vects.push_back( pep_cb_vect );
			}

			//END CUT OUT

			/////////////////////////////////////////////////////////////////////////////////////
			//set target anchor coords to those from ref anchor and store translations, rotations

			//check for noisy stub conformations due to no restraint at termini
			if( ref_pose.conformation().chain_begin( ref_pep_chain ) == ref_pep_anchor ){
				//check if previous refs were stored w/ nterm stub
				if( pep_anchor_is_nterm == false && i != 1 ) utility_exit_with_message( "nterm anchor must be first\n" );
				pep_anchor_is_nterm = true;
				//			pep_anchor_stub2 = "C";
				//			pep_anchor_stub3 = "O";
				pep_anchor_stub3 = "C";
			}
			if( ref_pose.conformation().chain_end( ref_pep_chain ) == ref_pep_anchor && pep_anchor_is_nterm == true ) utility_exit_with_message( "cannot have nterm and cterm ref pep anchors\n" ); 
			//TODO: WTF is this?
			if( pose.residue( pep_anchor ).name3() != "PRO" && ref_pose.residue( ref_pep_anchor ).name3() == "PRO" ) utility_exit_with_message( "cannot have PRO ref anchor for non-PRO anchor, need coords for H!\n" );

			if( pose.residue( pep_anchor ).name3() != "PRO" && !pep_anchor_is_nterm ) pose.set_xyz( AtomID( pose.residue( pep_anchor ).atom_index( "H" ), pep_anchor ), ref_pose.xyz( AtomID( ref_pose.residue( ref_pep_anchor ).atom_index( "H" ), ref_pep_anchor ) ) );
			pose.set_xyz( AtomID( pose.residue( pep_anchor ).atom_index( "N" ), pep_anchor ), ref_pose.xyz( AtomID( ref_pose.residue( ref_pep_anchor ).atom_index( "N" ), ref_pep_anchor ) ) );
			pose.set_xyz( AtomID( pose.residue( pep_anchor ).atom_index( "CA" ), pep_anchor ), ref_pose.xyz( AtomID( ref_pose.residue( ref_pep_anchor ).atom_index( "CA" ), ref_pep_anchor ) ) );
			pose.set_xyz( AtomID( pose.residue( pep_anchor ).atom_index( "CB" ), pep_anchor ), ref_pose.xyz( AtomID( ref_pose.residue( ref_pep_anchor ).atom_index( "CB" ), ref_pep_anchor ) ) );
			pose.set_xyz( AtomID( pose.residue( pep_anchor ).atom_index( "C" ), pep_anchor ), ref_pose.xyz( AtomID( ref_pose.residue( ref_pep_anchor ).atom_index( "C" ), ref_pep_anchor ) ) );
			pose.set_xyz( AtomID( pose.residue( pep_anchor ).atom_index( "HA" ), pep_anchor ), ref_pose.xyz( AtomID( ref_pose.residue( ref_pep_anchor ).atom_index( "HA" ), ref_pep_anchor ) ) );
			pose.set_xyz( AtomID( pose.residue( pep_anchor ).atom_index( "O" ), pep_anchor ), ref_pose.xyz( AtomID( ref_pose.residue( ref_pep_anchor ).atom_index( "O" ), ref_pep_anchor ) ) );


			pose.fold_tree( f_jump );

			Residue rsd1( pose.residue( prot_anchor ) );
			Residue rsd2( pose.residue( pep_anchor ) );
			StubID upstubid(  AtomID( rsd1.atom_index( prot_anchor_stub1 ), prot_anchor ) ,
					AtomID( rsd1.atom_index( prot_anchor_stub2 ), prot_anchor ) ,
					AtomID( rsd1.atom_index( prot_anchor_stub3 ), prot_anchor ) ) ;

			StubID downstubid(  AtomID( rsd2.atom_index( pep_anchor_stub1 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub2 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub3 ), pep_anchor ) ) ;

			RT const rt( pose.conformation().get_stub_transform( upstubid, downstubid ) );
			jump_trans_vector[i] = rt.get_translation();
			jump_rot_matrix[i] = rt.get_rotation();

			if( peptide_loop == 1 ) TR << "RB translation: x, y, z: " << jump_trans_vector[i].x() << ", "  << jump_trans_vector[i].y() << ", " << jump_trans_vector[i].z() << std::endl;
			///////////////////////

			pose.fold_tree( f_orient );

			jump_bond_angle[i] = pose.conformation().bond_angle( AtomID( rsd1.atom_index( prot_anchor_stub2 ), prot_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub1 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub2 ), pep_anchor ) );

			jump_tor1_angle[i] = pose.conformation().torsion_angle( AtomID( rsd1.atom_index( prot_anchor_stub2 ), prot_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub1 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub2 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub3 ), pep_anchor ) );

			jump_tor2_angle[i] = pose.conformation().torsion_angle( AtomID( rsd1.atom_index( prot_anchor_stub1 ), prot_anchor ) ,
					AtomID( rsd1.atom_index( prot_anchor_stub2 ), prot_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub1 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub2 ), pep_anchor ) );

			if( peptide_loop == 1 ) TR << "RB rotation: bond, tor1, tor2: " << degrees( jump_bond_angle[i] ) << ", "  << degrees( jump_tor1_angle[i] ) << ", " << degrees( jump_tor2_angle[i] ) << std::endl;

		}


		//homol_csts
		Size ref_pep_nres_max( ( -1 * ref_pep_pos_min ) + 1 + ref_pep_pos_max );
		vector1< Vector > pep_cst_ca_vects( ref_pep_nres_max, Vector( 0.0, 0.0, 0.0 ) );
		vector1< Vector > pep_cst_cb_vects( ref_pep_nres_max, Vector( 0.0, 0.0, 0.0 ) );
		vector1< Real > pep_cst_ca_tols( ref_pep_nres_max, 0.0 );
		vector1< Real > pep_cst_cb_tols( ref_pep_nres_max, 0.0 );
		vector1< Real > pep_cst_ca_sds( ref_pep_nres_max, 0.0 );
		vector1< Real > pep_cst_cb_sds( ref_pep_nres_max, 0.0 );
		vector1< Size > pep_cst_ca_indices( ref_pep_nres_max, 0 );
		vector1< Size > pep_cst_cb_indices( ref_pep_nres_max, 0 );
		if( option[ pepspec::homol_csts ].user() ){
			vector1< Size > n_pep_cas( ref_pep_nres_max, 0 );
			//sum pep_ca_vects -> pep_cst_ca_vects and n_pep_cas
			for( Size i = 1; i <= pep_ca_vects.size(); ++i ){
				vector1< std::pair< int, Vector > > pep_ca_vect( pep_ca_vects[ i ] );
				vector1< std::pair< int, Vector > > pep_cb_vect( pep_cb_vects[ i ] );
				for( Size ii = 1; ii <= pep_ca_vect.size(); ++ii ){
					int pep_pos( pep_ca_vect[ ii ].first );
					Vector pep_ca_coord( pep_ca_vect[ ii ].second );
					Vector pep_cb_coord( pep_cb_vect[ ii ].second );
					Size i_pep( pep_pos - ref_pep_pos_min + 1 );
					++n_pep_cas[ i_pep ];
					pep_cst_ca_vects[ i_pep ] += pep_ca_coord;
					pep_cst_cb_vects[ i_pep ] += pep_cb_coord;
				}
			}
			//normalize pep_cst_ca_vects -> avgs
			for( Size i = 1; i <= pep_cst_ca_vects.size(); ++i ){
				pep_cst_ca_vects[ i ] = pep_cst_ca_vects[ i ] / n_pep_cas[ i ];
				pep_cst_cb_vects[ i ] = pep_cst_cb_vects[ i ] / n_pep_cas[ i ];
			}
			//find max dists -> tols
			for( Size i = 1; i <= pep_ca_vects.size(); ++i ){
				vector1< std::pair< int, Vector > > pep_ca_vect( pep_ca_vects[ i ] );
				vector1< std::pair< int, Vector > > pep_cb_vect( pep_cb_vects[ i ] );
				for( Size ii = 1; ii <= pep_ca_vect.size(); ++ii ){
					int pep_pos( pep_ca_vect[ ii ].first );
					Vector pep_ca_coord( pep_ca_vect[ ii ].second );
					Vector pep_cb_coord( pep_cb_vect[ ii ].second );
					Size i_pep( pep_pos - ref_pep_pos_min + 1 );
					Real dist_to_ca_avg( pep_ca_coord.distance( pep_cst_ca_vects[ i_pep ] ) );
					Real dist_to_cb_avg( pep_cb_coord.distance( pep_cst_cb_vects[ i_pep ] ) );
					pep_cst_ca_sds[ i_pep ] += ( dist_to_ca_avg * dist_to_ca_avg );
					pep_cst_cb_sds[ i_pep ] += ( dist_to_cb_avg * dist_to_cb_avg );
					if( dist_to_ca_avg > pep_cst_ca_tols[ i_pep ] ) pep_cst_ca_tols[ i_pep ] = dist_to_ca_avg;
					if( dist_to_cb_avg > pep_cst_cb_tols[ i_pep ] ) pep_cst_cb_tols[ i_pep ] = dist_to_cb_avg;
				}
			}
			for( Size i = 1; i <= ref_pep_nres_max; ++i ){
				if( n_pep_cas[ i ] > 1 ) pep_cst_ca_sds[ i ] = std::sqrt( pep_cst_ca_sds[ i ] / ( n_pep_cas[ i ] - 1 ) );
				else pep_cst_ca_sds[ i ] = 100;
				if( n_pep_cas[ i ] > 1 ) pep_cst_cb_sds[ i ] = std::sqrt( pep_cst_cb_sds[ i ] / ( n_pep_cas[ i ] - 1 ) );
				else pep_cst_cb_sds[ i ] = 100;
				pep_cst_ca_indices[ i ] = ref_pep_pos_min + i - 1;
				pep_cst_cb_indices[ i ] = ref_pep_pos_min + i - 1;
			}
			//combine ca and cb data, build atom name vector
			vector1< Vector > pep_cst_all_vects( pep_cst_ca_vects );
			vector1< Real > pep_cst_all_tols( pep_cst_ca_tols );
			vector1< Real > pep_cst_all_sds( pep_cst_ca_sds );
			vector1< std::string > pep_cst_all_atom_names( ref_pep_nres_max, "CA" );
			vector1< Size > pep_cst_all_indices( pep_cst_ca_indices );
			vector1< Size > n_pep_cs( n_pep_cas );
			/*
			for( Size i = 1; i <= ref_pep_nres_max; ++i ){
				pep_cst_all_vects.push_back( pep_cst_cb_vects[ i ] );
				pep_cst_all_tols.push_back( pep_cst_cb_tols[ i ] );
				pep_cst_all_sds.push_back( pep_cst_cb_sds[ i ] );
				pep_cst_all_atom_names.push_back( "CB" );
				pep_cst_all_indices.push_back( pep_cst_cb_indices[ i ] );
				n_pep_cs.push_back( n_pep_cas[ i ] );
			}
			*/
			//now have pep_cst_ca_vects and pep_cst_ca_tols
			if( peptide_loop == 1 ) print_pep_csts( pep_cst_all_vects, pep_cst_all_tols, pep_cst_all_sds, n_pep_cs, pep_cst_all_indices, pep_cst_all_atom_names );
		}

		TR << "peptide stub atoms:\t" << pep_anchor_stub1 << "\t" << pep_anchor_stub2 << "\t" << pep_anchor_stub3 << "\n";    

		Size sd_denom( n_ref_poses );

		Real max_std_devs( option[ pepspec::n_anchor_dock_std_devs ] );

		// calc trans vector avg //
		Vector jump_trans_vector_avg( 0 );
		for( Size i = 1; i <= n_ref_poses; ++i ){ 
			jump_trans_vector_avg += ( jump_trans_vector[i] / n_ref_poses );
		}
		// ...and std dev //
		Real jump_trans_x_sd( 0 );
		Real jump_trans_y_sd( 0 );
		Real jump_trans_z_sd( 0 );
		if( option[ pepspec::prep_trans_std_dev ].user() ){
			jump_trans_x_sd = option[ pepspec::prep_trans_std_dev ];
			jump_trans_y_sd = option[ pepspec::prep_trans_std_dev ];
			jump_trans_z_sd = option[ pepspec::prep_trans_std_dev ];
		}
		else{
			for( Size i = 1; i<= n_ref_poses; ++i ){ 
				jump_trans_x_sd += numeric::square( jump_trans_vector[i].x() - jump_trans_vector_avg.x() ) / ( sd_denom );
				jump_trans_y_sd += numeric::square( jump_trans_vector[i].y() - jump_trans_vector_avg.y() ) / ( sd_denom );
				jump_trans_z_sd += numeric::square( jump_trans_vector[i].z() - jump_trans_vector_avg.z() ) / ( sd_denom );
			}
			jump_trans_x_sd = max_std_devs * std::sqrt( jump_trans_x_sd );
			jump_trans_y_sd = max_std_devs * std::sqrt( jump_trans_y_sd );
			jump_trans_z_sd = max_std_devs * std::sqrt( jump_trans_z_sd );
		}
		Vector jump_trans_vector_sd( jump_trans_x_sd, jump_trans_y_sd, jump_trans_z_sd );
		Real jump_trans_magn_sd( std::sqrt( jump_trans_x_sd * jump_trans_x_sd + jump_trans_y_sd * jump_trans_y_sd + jump_trans_z_sd * jump_trans_z_sd ) );

		if( peptide_loop == 1 ) TR << "RB translation avg: x, y, z: " << jump_trans_vector_avg.x() << ", "  << jump_trans_vector_avg.y() << ", " << jump_trans_vector_avg.z() << std::endl;
		if( peptide_loop == 1 ) TR << "Max RB translation deviation: x, y, z: " << jump_trans_x_sd << ", "  << jump_trans_y_sd << ", " << jump_trans_z_sd << std::endl;

		// calc orientation avgs //
		jump_bond_angle = shift_angles( jump_bond_angle );
		jump_tor1_angle = shift_angles( jump_tor1_angle );
		jump_tor2_angle = shift_angles( jump_tor2_angle );
		Real jump_bond_angle_avg( average( jump_bond_angle ) );
		Real jump_tor1_angle_avg( average( jump_tor1_angle ) );
		Real jump_tor2_angle_avg( average( jump_tor2_angle ) );
		// ...and std devs //
		Real jump_bond_angle_sd( 0 );
		Real jump_tor1_angle_sd( 0 );
		Real jump_tor2_angle_sd( 0 );
		if( option[ pepspec::prep_rot_std_dev ].user() ){
			jump_bond_angle_sd = radians( static_cast< Real >( option[ pepspec::prep_rot_std_dev ] ) );
			jump_tor1_angle_sd = radians( static_cast< Real >( option[ pepspec::prep_rot_std_dev ] ) );
			jump_tor2_angle_sd = radians( static_cast< Real >( option[ pepspec::prep_rot_std_dev ] ) );
		}
		else{
			for( Size i = 1; i <= n_ref_poses; ++i ){ 
				jump_bond_angle_sd += numeric::square( jump_bond_angle[i] - jump_bond_angle_avg ) / ( sd_denom );
				jump_tor1_angle_sd += numeric::square( jump_tor1_angle[i] - jump_tor1_angle_avg ) / ( sd_denom );
				jump_tor2_angle_sd += numeric::square( jump_tor2_angle[i] - jump_tor2_angle_avg ) / ( sd_denom );
			}
			jump_bond_angle_sd = max_std_devs * std::sqrt( jump_bond_angle_sd );
			jump_tor1_angle_sd = max_std_devs * std::sqrt( jump_tor1_angle_sd );
			jump_tor2_angle_sd = max_std_devs * std::sqrt( jump_tor2_angle_sd );
		}
		if( peptide_loop == 1 ) TR << "RB rotation avg: jump_bond_angle, jump_tor1, jump_tor2: " << degrees(jump_bond_angle_avg) << ", "  << degrees(jump_tor1_angle_avg) << ", " << degrees(jump_tor2_angle_avg) << std::endl;
		if( peptide_loop == 1 ) TR << "Max RB rotation deviation: jump_bond_angle, jump_tor1, jump_tor2: " << degrees(jump_bond_angle_sd) << ", "  << degrees(jump_tor1_angle_sd) << ", " << degrees(jump_tor2_angle_sd) << std::endl;

		////////////////////////////////////////////////////////////////////

		pose = start_pose;
		pose.fold_tree( f_jump );
		//set orientation to avg+-SD
		{
			Residue rsd1( pose.residue( prot_anchor ) );
			Residue rsd2( pose.residue( pep_anchor ) );
			StubID upstubid(  AtomID( rsd1.atom_index( prot_anchor_stub1 ), prot_anchor ) ,
					AtomID( rsd1.atom_index( prot_anchor_stub2 ), prot_anchor ) ,
					AtomID( rsd1.atom_index( prot_anchor_stub3 ), prot_anchor ) ) ;

			StubID downstubid(  AtomID( rsd2.atom_index( pep_anchor_stub1 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub2 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub3 ), pep_anchor ) ) ;

			Real jump_trans_x_div( jump_trans_vector_avg.x() + ( 2 * RG.uniform() - 1 ) * jump_trans_x_sd );
			Real jump_trans_y_div( jump_trans_vector_avg.y() + ( 2 * RG.uniform() - 1 ) * jump_trans_y_sd );
			Real jump_trans_z_div( jump_trans_vector_avg.z() + ( 2 * RG.uniform() - 1 ) * jump_trans_z_sd );

			Real jump_bond_angle_div( jump_bond_angle_avg + ( 2 * RG.uniform() - 1 ) * jump_bond_angle_sd );
			Real jump_tor1_angle_div( jump_tor1_angle_avg + ( 2 * RG.uniform() - 1 ) * jump_tor1_angle_sd );
			Real jump_tor2_angle_div( jump_tor2_angle_avg + ( 2 * RG.uniform() - 1 ) * jump_tor2_angle_sd );

			Vector jump_trans_vector_div( jump_trans_x_div, jump_trans_y_div, jump_trans_z_div );
			RT rt( jump_rot_matrix[1], jump_trans_vector_div );

			pose.conformation().set_stub_transform( upstubid, downstubid, rt );

			pose.fold_tree( f_orient );

			pose.conformation().set_bond_angle( AtomID( rsd1.atom_index(prot_anchor_stub2), prot_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub1 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub2 ), pep_anchor ) ,
					jump_bond_angle_div );
			pose.conformation().set_torsion_angle( AtomID( rsd1.atom_index(prot_anchor_stub2), prot_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub1 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub2 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub3 ), pep_anchor ),
					jump_tor1_angle_div );
			pose.conformation().set_torsion_angle( AtomID( rsd1.atom_index(prot_anchor_stub1), prot_anchor ) ,
					AtomID( rsd1.atom_index(prot_anchor_stub2), prot_anchor ),
					AtomID( rsd2.atom_index( pep_anchor_stub1 ), pep_anchor ) ,
					AtomID( rsd2.atom_index( pep_anchor_stub2 ), pep_anchor ) ,
					jump_tor2_angle_div );

			pose.fold_tree( f_jump );
		}
		/////////////////////////////////////////////////////////

		Size ref_index( static_cast< int >( RG.uniform() * n_ref_poses + 1 ) );
		if( option[ pepspec::prep_use_ref_rotamers ] ){
			pose.replace_residue( pep_anchor, ref_poses[ ref_index ].residue( ref_pep_anchors[ ref_index ] ), true );
		}

		//set p0 homol_csts
		else if( option[ pepspec::homol_csts ].user() ){
			scorefxn->set_weight( coordinate_constraint, 0.1 );
			Vector p0_vect( pep_cst_ca_vects[ -1 * ref_pep_pos_min + 1 ] );
			pep_coord_cst p0_cst;
			p0_cst.atom_name = "CA";
			p0_cst.pep_pos = 0;
			p0_cst.x = p0_vect.x();
			p0_cst.y = p0_vect.y();
			p0_cst.z = p0_vect.z();
			p0_cst.x0 = 0;
			p0_cst.sd = 1.0;
			p0_cst.tol = pep_cst_ca_tols[ -1 * ref_pep_pos_min + 1 ];
			set_pep_cst( pose, pep_anchor, p0_cst );
		} 

		vector1< bool > ignore_clash( pose.total_residue(), false );
		TR << "Ignoring clash residues:\t";
		for( Size i = 1; i <= pose.total_residue() - 1; ++i ){
			vector1< bool > this_clash( pose.total_residue(), false );
			this_clash[ i ] = true;
			if( has_clash( pose, this_clash, scorefxn, option[ pepspec::clash_cutoff ], false ) ){
				TR << string_of( i ) + ", ";
				ignore_clash[ i ] = true;
			}
		}
		TR << "\n";

		//Real Docking!!
		MonteCarloOP mc_dock2 ( new MonteCarlo( pose, *scorefxn, 0.8 ) );
		pose.update_residue_neighbors();
		Pose restart_pose( pose ); 
		for( Size i_dock = 1; i_dock <= option[ pepspec::n_dock_loop]; ++i_dock ){

			pose = restart_pose;
			( *scorefxn )( pose );
			MonteCarloOP mc_dock ( new MonteCarlo( pose, *scorefxn, 0.8 ) );
			vector1< bool > check_clash( pose.total_residue(), false );
			check_clash[ pep_anchor ] = true;
			rigid::RigidBodyPerturbNoCenterMoverOP rb_mover = new rigid::RigidPerturbNoCenterMover( pep_jump, 2.5, 0.25 );
			rb_mover->apply( pose );
			mc_dock->boltzmann( pose );

			//get seqpos nbrs from energy map and rottrials
			vector1< bool > is_pep_nbr( pose.total_residue(), false );
			EnergyGraph const & energy_graph( pose.energies().energy_graph() );
			for ( graph::Graph::EdgeListConstIter
					ir  = energy_graph.get_node( pep_anchor )->const_edge_list_begin(),
					ire = energy_graph.get_node( pep_anchor )->const_edge_list_end();
					ir != ire; ++ir ) {
				Size this_nbr( (*ir)->get_other_ind( pep_anchor ) );
				is_pep_nbr[ this_nbr ] = true;
				check_clash[ this_nbr ] = ( true && !ignore_clash[ this_nbr ] );
			}

			{
				kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
				//mm_min->set_jump( pep_jump, true );
				mm_min->set_chi( pep_anchor );
				mm_min->set_chi( is_pep_nbr );
				protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, scorefxn, "dfpmin", 0.001, true );
				min_mover->apply( pose );
			}
			mc_dock->boltzmann( pose );

			( *soft_scorefxn )( pose );
			pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
			task->initialize_from_command_line().or_include_current( true );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if( !option[ pepspec::anchor_type ].user() && i == pep_anchor ) continue;
				else if( i == pep_anchor && option[ pepspec::prep_use_ref_rotamers ] ) task->nonconst_residue_task( i ).prevent_repacking();
				else if( i == pep_anchor || is_pep_nbr[ i ] ) task->nonconst_residue_task( i ).restrict_to_repacking();
				else task->nonconst_residue_task( i ).prevent_repacking();
			}
			protocols::simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, task, 1 ) );
			packer->apply( pose );
			( *scorefxn )( pose );

			{
				kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
				//mm_min->set_jump( pep_jump, true );
				mm_min->set_chi( pep_anchor );
				mm_min->set_chi( is_pep_nbr );
				protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, scorefxn, "dfpmin", 0.001, true );
				min_mover->apply( pose );
			}
			mc_dock->boltzmann( pose );
			mc_dock->recover_low( pose );

			mc_dock2->boltzmann( pose );
		}
		mc_dock2->recover_low( pose );
		vector1< bool > check_clash( pose.total_residue(), false );
		EnergyGraph const & energy_graph( pose.energies().energy_graph() );
		for ( graph::Graph::EdgeListConstIter
				ir  = energy_graph.get_node( pep_anchor )->const_edge_list_begin(),
				ire = energy_graph.get_node( pep_anchor )->const_edge_list_end();
				ir != ire; ++ir ) {
			Size this_nbr( (*ir)->get_other_ind( pep_anchor ) );
			check_clash[ this_nbr ] = ( true && !ignore_clash[ this_nbr ] );
		}
		check_clash[ pep_anchor ] = true;
		if( has_clash( pose, check_clash, scorefxn, option[ pepspec::clash_cutoff ], true ) ){
			//dump_pdb( pose, "clash." + string_of( peptide_loop ) + ".pdb" );
			TR << "Failed to resolve clash in iter " + string_of( peptide_loop ) + "\n";
			if( peptide_loop > 0 ) --peptide_loop;
			continue;
		}


		std::string pdb_name ( option[ out::file::o ] );
		if( option[pepspec::n_peptides ] > 1 ) pdb_name += "_" + string_of( peptide_loop );
		pdb_name += ".pdb";

		( *scorefxn )( pose );
		out_file << pdb_name + "\t";
		out_file << pose.energies().total_energies().weighted_string_of( scorefxn->weights() );
		Real total_score( pose.energies().total_energies().dot( scorefxn->weights() ) );
		out_file << " total_score: " << total_score << "\t";
		Real total_prot_score( total_score - prot_score );
		out_file << " total-prot_score: " << total_prot_score << "\t";
		if( option[ pepspec::rmsd_analysis ] ) out_file << pep_rmsd_analysis( pose, prot_begin, prot_end, pep_anchor, pep_anchor, pep_anchor );
		out_file << std::endl; 

		////////////////////////////////////////////////////////

		Size n_append( option[ pepspec::n_append ] );
		Size n_prepend( option[ pepspec::n_prepend ] );

		ResidueOP ala( ResidueFactory::create_residue( rsd_set.name_map( "ALA" ) ) );

		Size pep_begin( pep_anchor );
		Size pep_end( pep_anchor );
		for( Size i = 1; i <= n_append; ++i ){
			pose.conformation().safely_append_polymer_residue_after_seqpos( *ala, pep_end, true );
			++pep_end;
		}
		for( Size i = 1; i <= n_prepend; ++i ){
			pose.conformation().safely_prepend_polymer_residue_before_seqpos( *ala, pep_begin, true );
			++pep_anchor;
			++pep_end;
		}

		for( Size i_omega = pep_begin; i_omega <= pep_end - 1; ++i_omega ){
			pose.set_omega( i_omega, 180.0 );
		}

		if( option[ pepspec::use_input_bb ] ){

			vector1< Size > ref_indices;
			for( Size ii = 1; ii <= n_ref_poses; ++ii ) ref_indices.push_back( ii );
			//from -n to +m
			int begin( pep_begin - pep_anchor );
			int end( pep_end - pep_anchor );
			numeric::random::random_permutation( ref_indices, RG ); //shuffle ref poses
			for( Size i_ref = 1; i_ref <= n_ref_poses; ++i_ref ){
				ref_index = ref_indices[ i_ref ];
				Pose ref_pose( ref_poses[ ref_index ] );
				Size ref_pep_anchor( ref_pep_anchors[ ref_index ] );
				Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pose.chain( ref_pep_anchor ) ) );
				Size ref_pep_end( ref_pose.conformation().chain_end( ref_pose.chain( ref_pep_anchor ) ) );
				int ref_begin = static_cast< int >( ref_pep_begin ) - static_cast< int >( ref_pep_anchor );
				int ref_end = static_cast< int >( ref_pep_end ) - static_cast< int >( ref_pep_anchor );

				//goto next if ref does not overlap
				if( begin < ref_begin || end > ref_end ){
					if( i_ref == n_ref_poses ) utility_exit_with_message( "N append/prepend > all ref poses!\n" );
					continue;
				}
				else{
					for( int ii = begin; ii <= end; ++ii ){
						//or use first avail then goto next res
						if( ii > begin ) pose.set_phi( pep_anchor + ii, ref_pose.phi( ref_pep_anchor + ii ) );
						pose.set_psi( pep_anchor + ii, ref_pose.psi( ref_pep_anchor + ii ) );
						if( ii < end ){
							pose.set_omega( pep_anchor + ii, ref_pose.omega( ref_pep_anchor + ii ) );
						}
					}
					break;
				}
			}
		}

		dump_pdb( pose, pdb_name );
		pdb_out_list_file << pdb_name << std::endl;

		/*
		   Real cst_distance_n( pose.residue( prot_anchor ).xyz("CA").distance( pose.residue( new_pep_anchor ).xyz("N") ) );
		   Real cst_distance_ca( pose.residue( prot_anchor ).xyz("CA").distance( pose.residue( new_pep_anchor ).xyz("CA") ) );
		   Real cst_distance_c( pose.residue( prot_anchor ).xyz("CA").distance( pose.residue( new_pep_anchor ).xyz("C") ) );
		   std::string cst_file_name_str( option[ out::file::o ]+"_"+string_of(peptide_loop)+".cst" );
		   char const *cst_file_name = cst_file_name_str.c_str();
		   std::fstream cst_file( cst_file_name, std::ios::out );
		   cst_file << "[ atompairs ]\n";
		   cst_file << "CA\t" + string_of( prot_anchor ) + "\tN\t" + string_of( new_pep_anchor ) + "\tHARMONIC\t" + string_of( cst_distance_n ) + "\t0.1\n";
		   cst_file << "CA\t" + string_of( prot_anchor ) + "\tCA\t" + string_of( new_pep_anchor ) + "\tHARMONIC\t" + string_of( cst_distance_ca ) + "\t0.1\n";
		   cst_file << "CA\t" + string_of( prot_anchor ) + "\tC\t" + string_of( new_pep_anchor ) + "\tHARMONIC\t" + string_of( cst_distance_c ) + "\t0.1\n";
		   cst_list_file << cst_file_name_str << std::endl;
		 */
	}


}


int main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init(argc, argv);

		run_pep_prep();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
