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

 #include <core/fragment/ConstantLengthFragSet.hh>
 #include <core/fragment/FragSet.hh>
 #include <core/fragment/Frame.hh>
 #include <core/fragment/picking_old/vall/util.hh>
 #include <core/fragment/picking_old/FragmentLibraryManager.hh>

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
 #include <protocols/simple_moves/FragmentMover.hh>

 #include <protocols/viewer/viewers.hh>

 #include <core/types.hh>

 #include <core/scoring/sasa.hh>

// #include <core/util/prof.hh> // profiling
// #include <core/util/CacheableData.hh> // profiling

// #include <core/sequence/SequenceMapping.hh>

 #include <core/chemical/AtomTypeSet.hh>

 #include <core/chemical/AA.hh>
 #include <core/conformation/Residue.hh>
 #include <core/conformation/ResidueMatcher.hh>
 #include <core/pack/rotamer_set/RotamerCouplings.hh>
 #include <core/chemical/ResidueTypeSet.hh>
 #include <core/chemical/ResidueTypeSelector.hh>
 #include <core/conformation/ResidueFactory.hh>
 #include <core/chemical/VariantType.hh>
 #include <core/chemical/util.hh>
 #include <core/chemical/ChemicalManager.hh>

 #include <core/scoring/rms_util.hh>
 #include <core/scoring/etable/Etable.hh>
 #include <core/scoring/ScoringManager.hh>
 #include <core/scoring/ScoreFunction.hh>
 #include <core/scoring/ScoreFunctionFactory.hh>
 #include <core/scoring/Ramachandran.hh>
 #include <core/scoring/hbonds/HBondSet.hh>
 #include <core/scoring/hbonds/hbonds.hh>
 #include <core/scoring/etable/count_pair/CountPairFactory.hh>
 #include <core/scoring/etable/count_pair/CountPairAll.hh>
 #include <core/scoring/etable/count_pair/CountPairFunction.hh>
 #include <core/scoring/etable/EtableEnergy.hh>
 #include <core/scoring/etable/Etable.hh>
 #include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>

 #include <core/pack/rotamer_trials.hh>
 #include <core/pack/pack_rotamers.hh>
 #include <core/pack/task/PackerTask.hh>
 #include <core/pack/task/TaskFactory.hh>
 #include <core/pack/task/operation/TaskOperations.hh>
 #include <protocols/toolbox/IGEdgeReweighters.hh>
 #include <core/pack/task/IGEdgeReweightContainer.hh>

 #include <core/kinematics/FoldTree.hh>
 #include <core/kinematics/MoveMap.hh>
 #include <core/kinematics/util.hh>

 #include <core/pose/Pose.hh>
 #include <core/pose/util.hh>
 #include <core/pose/PDBPoseMap.hh>
 #include <core/pose/PDBInfo.hh>
 #include <core/pose/metrics/CalculatorFactory.hh>

 #include <basic/options/util.hh>//option.hh>
// #include <basic/options/after_opts.hh>

 #include <basic/basic.hh>
 #include <basic/MetricValue.hh>

 #include <basic/database/open.hh>

#include <devel/init.hh>

 #include <utility/vector1.hh>
 #include <utility/file/file_sys_util.hh>

 #include <numeric/xyzVector.hh>
 #include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>

 #include <core/scoring/constraints/ConstraintSet.hh>
 #include <core/scoring/func/FlatHarmonicFunc.hh>
 #include <core/scoring/constraints/AtomPairConstraint.hh>
 #include <core/scoring/constraints/CoordinateConstraint.hh>
 #include <core/scoring/constraints/ConstraintIO.hh>

// // C++ headers
 #include <cstdlib>
 #include <fstream>
 #include <iostream>
 #include <string>

 #include <basic/Tracer.hh>

 #include <basic/options/keys/in.OptionKeys.gen.hh>
 #include <basic/options/keys/out.OptionKeys.gen.hh>
 #include <basic/options/keys/score.OptionKeys.gen.hh>
 #include <basic/options/keys/pepspec.OptionKeys.gen.hh>
 #include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <protocols/toolbox/pose_metric_calculators/DecomposeAndReweightEnergiesCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/ResidueDecompositionByChainCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>


using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
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

static thread_local basic::Tracer TR( "apps.pilot.chrisk/pepspec" );


Size prot_chain( 0 );
Size prot_begin( 0 );
Size prot_anchor( 0 );
Size prot_end( 0 );

Size pep_jump( 2 );
Size pep_chain( 0 );
Size pep_begin( 0 );
Size pep_anchor( 0 );
Size pep_end( 0 );

std::string input_seq;

//load homol coord cst data
//parse cst line as vector1[ 10 ] name1 res1 x y z x0 sd tol
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

class myMC{
public:
	myMC(
		Pose pose,
		Real score,
		Real temp_init
	);
	Pose last_pose;
	Real last_score;
	Real temp;
	bool accept_;
	bool accept();
	void roll(
		Pose & pose,
		Real score
	);
};

myMC::myMC(
	Pose pose,
	Real score,
	Real temp_init
)
{
	last_pose = pose;
	last_score = score;
	temp = temp_init;
	accept_ = true;
}

void
myMC::roll(
	Pose & pose,
	Real score
)
{
	accept_ = false;
	if( score <= last_score ) accept_ = true;
	else{
		Real dE = score - last_score;
		Real p = std::exp( -1 * dE / temp );
		if( numeric::random::rg().uniform() < p ) accept_ = true;
	}
	if( accept_ ){
		last_pose = pose;
		last_score = score;
	}
	else{
		pose = last_pose;
		score = last_score;
	}
}

bool
myMC::accept()
{
	return accept_;
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


/*
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

	if( option[ pepspec::homol_csts ].user() ) scorefxn->set_weight( coordinate_constraint, option[ OptionKeys::constraints::cst_weight ] );
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
//			if ( !options::option[ options::OptionKeys::out::file::output_virtual ]() &&
//				(int)atom.type() == (int)rsd.atom_type_set().n_atomtypes() ) continue;

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
	out << "ENDMDL\n\n\n";
	out << "###residue_total_energies###\n";
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Real residue_total_energy( pose.energies().residue_total_energies( i ).dot( scorefxn->weights() ) );
		out << string_of( i ) << " " << pose.residue( i ).name1() << " " << pose.energies().residue_total_energies( i ).weighted_string_of( scorefxn->weights() ) << "\t" << residue_total_energy << "\n";
	}
}
pose::Pose
my_boltzmann(
	pose::Pose pose,
	Real score
)
{

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

	TR << "sel flex_prot, ";
	for( Size i = 1; i <= pose.total_residue(); ++i ){
		if( is_nbr[ i ] ) TR << "resi " << string_of( i ) << " or ";
	}
	TR << std::endl;
}
*/

/// @details  This function will make a sequence mutation while trying to preserve the variants
void
make_sequence_change(
		Size const seqpos,
		AA const & new_aa,
		pose::Pose & pose
		)
{
	Size const which_his_variant = 1;
	conformation::Residue const & current_rsd( pose.residue( seqpos ) );
	if ( current_rsd.aa() == new_aa ) return; // already done

	ResidueTypeCOPs rsd_types
		( ResidueTypeSelector().set_aa( new_aa ).match_variants( current_rsd.type() ).select( current_rsd.residue_type_set() ) );

	Size rsd_types_index( 1 );
	std::string const errmsg
		( "make_sequence_change failed: new_aa= "+name_from_aa(new_aa)+" rsd_types.size()= "+string_of( rsd_types.size() ) );

	if ( new_aa == aa_his ) {
		if ( rsd_types.size() != 2 || which_his_variant > 2 ) utility_exit_with_message( errmsg );
		rsd_types_index = which_his_variant;
	} else if ( rsd_types.size() != 1 ) {
		utility_exit_with_message( errmsg );
	}

	conformation::ResidueOP new_rsd( ResidueFactory::create_residue( *(rsd_types[ rsd_types_index ] ),
				current_rsd, pose.conformation() ) );
	pose.replace_residue( seqpos, *new_rsd, false );
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
//				TR<< "fa_rep: " << string_of( clash ) << " at " << pose.residue( j ).name1() << string_of( j ) << "-" << pose.residue( seqpos ).name1() << string_of( seqpos ) << std::endl;
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
//				TR<< "fa_rep: " << string_of( clash ) << " at " << pose.residue( j ).name1() << string_of( j ) << "-" << pose.residue( seqpos ).name1() << string_of( seqpos ) << std::endl;
				is_clash = true;
				break;
			}
		}
		if( is_clash == true ) break;
	}
	return is_clash;
}

vector1< std::pair< Size, Size > >
get_clash_pairs(
		pose::Pose pose,
		vector1< bool > is_checked,
		scoring::ScoreFunctionOP const & scorefxn,
		Real const clash_threshold
	 )
{
	using namespace scoring;
	using namespace chemical;
	vector1< std::pair< Size, Size > > clash_pairs;

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
				clash_pairs.push_back( std::pair< Size, Size >( seqpos, j ) );
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
				clash_pairs.push_back( std::pair< Size, Size >( seqpos, j ) );
			}
		}
	}
	return clash_pairs;
}

std::string
pep_rmsd_analysis(
	pose::Pose pose
)
{

	pose::Pose ref_pose;
	std::string ref_name( option[ in::file::native ]() );
	import_pose::pose_from_pdb( ref_pose, ref_name );

	Size ref_pep_anchor_in( option[ pepspec::pep_anchor ] );
	if( option[ pepspec::native_pep_anchor ].user() ) ref_pep_anchor_in = option[ pepspec::native_pep_anchor ];
	std::string ref_pep_chain_in( option[ pepspec::pep_chain ] );
	if( option[ pepspec::native_pep_chain ].user() ) ref_pep_chain_in = option[ pepspec::native_pep_chain ];
	Size ref_pep_anchor( ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_in ) );
	Size ref_pep_chain( ref_pose.chain( ref_pep_anchor ) );
	Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pep_chain ) );
	Size ref_pep_end( ref_pose.conformation().chain_end( ref_pep_chain ) );
	Size this_pep_begin( pep_begin );
	Size this_pep_end( pep_end );

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

	Size pep_nterm( pep_anchor - this_pep_begin );
	Size ref_pep_nterm( ref_pep_anchor - ref_pep_begin );
	Size pep_cterm( this_pep_end - pep_anchor );
	Size ref_pep_cterm( ref_pep_end - ref_pep_anchor );

	Size nterm( pep_nterm );
	Size cterm( pep_cterm );
	if( pep_nterm < ref_pep_nterm ) ref_pep_begin += ref_pep_nterm - pep_nterm;
	else if( pep_nterm > ref_pep_nterm ){
		this_pep_begin += pep_nterm - ref_pep_nterm;
		nterm = ref_pep_nterm;
	}
	if( pep_cterm < ref_pep_cterm ) ref_pep_end -= ref_pep_cterm - pep_cterm;
	else if( pep_cterm > ref_pep_cterm ){
		this_pep_end -= pep_cterm - ref_pep_cterm;
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
		core::scoring::ScoreFunctionOP full_scorefxn(  get_score_function() );
		core::scoring::superimpose_pose( pose, ref_pose, atom_map );

	}
	Real total_rmsd( 0 );
	std::string rmsd_analysis( "" );
	for( int i_seq = 0; i_seq <= nterm + cterm; ++i_seq ){
		Size ref_pep_seqpos = ref_pep_begin + i_seq;
		Size pep_seqpos = this_pep_begin + i_seq;
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
	pose::Pose pose
)
{

	pose::Pose ref_pose;
	std::string ref_name( option[ in::file::native ]() );
	import_pose::pose_from_pdb( ref_pose, ref_name );

	Size ref_pep_anchor_in( option[ pepspec::pep_anchor ] );
	if( option[ pepspec::native_pep_anchor ].user() ) ref_pep_anchor_in = option[ pepspec::native_pep_anchor ];
	std::string ref_pep_chain_in( option[ pepspec::pep_chain ] );
	if( option[ pepspec::native_pep_chain ].user() ) ref_pep_chain_in = option[ pepspec::native_pep_chain ];
	Size ref_pep_anchor( ref_pose.pdb_info()->pdb2pose( ref_pep_chain_in[0], ref_pep_anchor_in ) );
	Size ref_pep_chain( ref_pose.chain( ref_pep_anchor ) );
	Size ref_pep_begin( ref_pose.conformation().chain_begin( ref_pep_chain ) );
	Size ref_pep_end( ref_pose.conformation().chain_end( ref_pep_chain ) );
	Size this_pep_begin( pep_begin );
	Size this_pep_end( pep_end );

	Size pep_nterm( pep_anchor - this_pep_begin );
	Size ref_pep_nterm( ref_pep_anchor - ref_pep_begin );
	Size pep_cterm( this_pep_end - pep_anchor );
	Size ref_pep_cterm( ref_pep_end - ref_pep_anchor );

	Size nterm( pep_nterm );
	Size cterm( pep_cterm );
	if( pep_nterm < ref_pep_nterm ) ref_pep_begin += ref_pep_nterm - pep_nterm;
	else if( pep_nterm > ref_pep_nterm ){
		this_pep_begin += pep_nterm - ref_pep_nterm;
		nterm = ref_pep_nterm;
	}
	if( pep_cterm < ref_pep_cterm ) ref_pep_end -= ref_pep_cterm - pep_cterm;
	else if( pep_cterm > ref_pep_cterm ){
		this_pep_end -= pep_cterm - ref_pep_cterm;
		cterm = ref_pep_cterm;
	}

	Real total_phi( 0 );
	Real total_psi( 0 );
	Real total_omega( 0 );
	std::string phipsi_analysis( "" );
	for( int i_seq = 0; i_seq <= nterm + cterm; ++i_seq ){
		Size ref_pep_seqpos = ref_pep_begin + i_seq;
		Size pep_seqpos = this_pep_begin + i_seq;
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
	pose::Pose & pose
)
{
	add_lower_terminus_type_to_pose_residue( pose, pep_begin );
	add_upper_terminus_type_to_pose_residue( pose, pep_end );
	pose.conformation().update_polymeric_connection( pep_begin );
	pose.conformation().update_polymeric_connection( pep_begin + 1 );
	pose.conformation().update_polymeric_connection( pep_end );
	pose.conformation().update_polymeric_connection( pep_end - 1 );

}

void add_pep_res(
	pose::Pose & pose,
	bool add_nterm,
	bool add_cterm
)
{

		ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
		ResidueOP vrt( ResidueFactory::create_residue( rsd_set.name_map( "GLY" ) ) );
//		ResidueOP vrt( ResidueFactory::create_residue( rsd_set.name_map( "VirtBB" ) ) );

		if( add_cterm ){
			pose.conformation().safely_append_polymer_residue_after_seqpos( *vrt, pep_end, true );
			pep_end = pep_end + 1;
			pose.set_omega( pep_end - 1, 180.0 );
			pose.conformation().update_polymeric_connection( pep_end );
			pose.conformation().update_polymeric_connection( pep_end - 1 );
			add_variant_type_to_pose_residue( pose, "VIRTUAL_BB", pep_end );
		}

		if( add_nterm ){
			pose.conformation().safely_prepend_polymer_residue_before_seqpos( *vrt, pep_begin, true );
			pep_end = pep_end + 1;
			pep_anchor = pep_anchor + 1;
			pose.set_omega( pep_begin, 180.0 );
			pose.conformation().update_polymeric_connection( pep_begin );
			pose.conformation().update_polymeric_connection( pep_begin + 1 );
			add_variant_type_to_pose_residue( pose, "VIRTUAL_BB", pep_begin );
		}

		//replace termini
		add_termini( pose );
}

void remove_pep_res(
	pose::Pose & pose,
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
		add_termini( pose );
}

void
initialize_peptide(
	pose::Pose & pose
)
{

	if( option[ pepspec::remove_input_bb ] ){
		//remove non-anchor positions
		while( pep_anchor != pep_begin ){
			pose.conformation().delete_residue_slow( pep_begin );
			pep_anchor = pep_anchor - 1;
			pep_end = pep_end - 1;
			pose.conformation().update_polymeric_connection( pep_begin );
			pose.conformation().update_polymeric_connection( pep_begin + 1 );
		}
		while( pep_anchor != pep_end ){
			pose.conformation().delete_residue_slow( pep_end );
			pep_end = pep_end - 1;
			pose.conformation().update_polymeric_connection( pep_end );
			pose.conformation().update_polymeric_connection( pep_end - 1 );
		}
	}

	Pose start_pose( pose );
	ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
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
	Residue pep_anchor_res = pose.residue( pep_anchor );

	//must maintain fold tree!
	pose.fold_tree( start_pose.fold_tree() );

	//if using input backbone and designing, mutate all res to cg_res_type
	std::string cg_res_type( option[ pepspec::cg_res_type ] );
	ResidueOP res( ResidueFactory::create_residue( rsd_set.name_map( cg_res_type ) ) );
	if( option[ pepspec::use_input_bb ] && !option[ pepspec::no_design ] && !option[ pepspec::input_seq ].user() ){
		for( Size mut_site = pep_begin; mut_site <= pep_end; mut_site++ ){
			if(mut_site==pep_anchor) continue;
			make_sequence_change( mut_site, chemical::AA( aa_from_oneletter_code( res->name1() ) ) , pose );
		}
	}

}

void
set_pep_csts(
	pose::Pose & pose
)
{
	using namespace core::scoring::constraints;
	for( Size i_cst = 1; i_cst <= pep_coord_csts.size(); ++i_cst ){
		pep_coord_cst pep_cst( pep_coord_csts[ i_cst ] );
		Size seqpos( pep_cst.pep_pos + pep_anchor );
		if( seqpos < pep_begin ) continue;
		if( seqpos > pep_end ) break;
		if( numeric::random::rg().uniform() > option[ pepspec::p_homol_csts ] ) continue;
		Vector pep_cst_vector;
		pep_cst_vector.x( pep_cst.x );
		pep_cst_vector.y( pep_cst.y );
		pep_cst_vector.z( pep_cst.z );
		ConstraintCOP this_cst( new CoordinateConstraint( AtomID( pose.residue( seqpos ).atom_index( pep_cst.atom_name ), seqpos ), AtomID( pose.residue( prot_anchor ).atom_index( "CA" ), prot_anchor ), pep_cst_vector, new FlatHarmonicFunc( pep_cst.x0, pep_cst.sd, pep_cst.tol ) ) );
//		ConstraintCOP cst( new CoordinateConstraint( AtomID( pose.residue( i ).atom_index( "CA" ), i ), AtomID( pose.residue( pep_anchor ).atom_index( "CA" ), pep_anchor ), pose.residue( i ).xyz( "CA" ), new FlatHarmonicFunc( 0.0, 0.1, 2.0 ) ) );
		pose.add_constraint( this_cst );
	}
}

/*
void
refine_fa_pep_bb(
	Pose & pose,
	vector1< bool > is_pep,
	scoring::ScoreFunctionOP cen_scorefxn
)
{

	kinematics::MoveMapOP mm ( new kinematics::MoveMap );
	mm->set_bb( is_pep );
	mm->set_chi( pep_anchor, true );

	rigid::RigidBodyPerturbMoverOP rb_mover = new rigid::RigidBodyPerturbMover( pep_jump, 0.1, 0.0 );
	rb_mover->apply( pose );

	protocols::simple_moves::SmallMoverOP rep_small_mover( new protocols::simple_moves::SmallMover( mm, 2.0, 10 ) );
	rep_small_mover->angle_max( 'H', 1.5 );
	rep_small_mover->angle_max( 'E', 1.5 );
	rep_small_mover->angle_max( 'L', 1.5 );

	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm, cen_scorefxn, "dfpmin", 0.001, true );

	RandomMoverOP rand_mover( new protocols::moves::RandomMover() );
	rand_mover->add_mover( rep_small_mover, 9 );
//	rand_mover->add_mover( rb_mover, 1 );
//	rand_mover->add_mover( min_mover, 1 );

	MonteCarloOP mc_rep ( new MonteCarlo( pose, *cen_scorefxn, 0.8 ) );
	TrialMoverOP rep_trial = new TrialMover( rand_mover, mc_rep );
	for( Size i = 1; i <= 1000; ++i ){
		rep_trial->apply( pose );
		mc_rep->boltzmann( pose );
	}
//	mc_rep->recover_low( pose );
}
*/

/// @brief helper code for fragments generation, copied from S.M.Lewis
core::fragment::FragSetCOP
make_frags(
	core::Size const start,
	core::Size const stop,
	std::string const & seq
)
{
	core::Size const frags_length( 1 ); //magic number: 3mer fragments!!
	core::fragment::FragSetOP fragset( new core::fragment::ConstantLengthFragSet( frags_length ) );
	std::string ss_string( frags_length, 'L' );
	core::fragment::FragDataOPs list;

	for( core::Size j = start; j <= stop - frags_length + 1; ++j ){
		std::string const seqsubstr( seq, j-1, frags_length ); //j-1 accounts for string [] from 0
		list =  core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( ss_string, seqsubstr, 200, false ); //magic number: 200 fragments per position (not duplicated - this will be like robetta server fragments)
		core::fragment::FrameOP frame;
		frame = new core::fragment::Frame( j );
		frame->add_fragment( list );
		fragset->add( frame );
	}
	return fragset;
}

//@brief helper code for fragments generation, copied from S.M.Lewis
core::fragment::FragSetCOP
make_1mer_frags(
	core::Size const seqpos_start,
	core::Size const seqpos_stop,
	std::string const & seq,
	Size const nfrags
)
{
	core::Size const frags_length( 1 );
	core::fragment::FragSetOP fragset( new core::fragment::ConstantLengthFragSet( frags_length ) );
	//80% L, 10% H, 10% S
	core::fragment::FragDataOPs list;

	for( core::Size j = seqpos_start; j <= seqpos_stop - frags_length + 1; ++j ){
		//get substr
		//get ss str
		std::string ss_string( "L" );
		if( option[ pepspec::ss_type ].user() ) ss_string = std::string( option[ pepspec::ss_type ] );
		else{
			if( numeric::random::rg().uniform() < 0.1 ) ss_string = "H";
			else if( numeric::random::rg().uniform() > 0.9 ) ss_string = "S";
		}
		//get fraglist from ss or ss and seq
		if( option[ pepspec::input_seq ].user() && seq[ j - 1 ] != 'X' ){
			std::string const seqsubstr( seq, j-1, frags_length ); //j-1 accounts for string [] from 0
			list =  core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( ss_string, seqsubstr, nfrags, true );
		}
		else list =  core::fragment::picking_old::vall::pick_fragments_by_ss( ss_string, nfrags, true );
		core::fragment::FrameOP frame;
		frame = new core::fragment::Frame( j );
		frame->add_fragment( list );
		fragset->add( frame );
	}
	return fragset;
}

void
gen_pep_bb_sequential(
	pose::Pose & pose,
	scoring::ScoreFunctionOP cen_scorefxn
)
{
	using namespace scoring;
	using namespace chemical;
	using namespace scoring::methods;
	using namespace scoring::etable;
	using namespace scoring::etable::count_pair;

	Size n_build_loop( option[ pepspec::n_build_loop ] );
	Real clash_cutoff( option[ pepspec::clash_cutoff ] );
	Size n_prepend( option[ pepspec::n_prepend ] );
	Size n_append( option[ pepspec::n_append ] );

	( *cen_scorefxn )( pose );
	Real min_score( pose.energies().total_energies()[ vdw ] );

	Pose replace_res_pose( pose );
	ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
	std::string cg_res_type( option[ pepspec::cg_res_type ] );
	ResidueOP res( ResidueFactory::create_residue( rsd_set.name_map( cg_res_type ) ) );
	Size break_loop( 0 );
	Size n_prepended( 0 );
	Size n_appended( 0 );
	while( ( n_prepended < n_prepend || n_appended < n_append ) ){
		Size this_prepend( n_prepend );
		Size this_append( n_append );
		if( option[ pepspec::gen_pep_bb_sequential ] || !option[ pepspec::homol_csts ].user() ){
			this_prepend = static_cast< int >( numeric::random::rg().uniform() * ( n_prepend - n_prepended + 1 ) );
			this_append = static_cast< int >( numeric::random::rg().uniform() * ( n_append - n_appended + 1 ) );
		}

		std::string this_input_seq( input_seq );
		for( Size ii = 1; ii <= this_prepend; ++ii ){
			pose.prepend_polymer_residue_before_seqpos( *res, pep_begin, true );
			pep_end = pep_end + 1;
			pep_anchor = pep_anchor + 1;
			pose.set_omega( pep_begin, 180.0 );
			pose.conformation().update_polymeric_connection( pep_begin );
			pose.conformation().update_polymeric_connection( pep_begin + 1 );
		}
		for( Size ii = 1; ii <= this_append; ++ii ){
			pose.append_polymer_residue_after_seqpos( *res, pep_end, true );
			pep_end = pep_end + 1;
			pose.set_omega( pep_end - 1, 180.0 );
			pose.conformation().update_polymeric_connection( pep_end );
			pose.conformation().update_polymeric_connection( pep_end - 1 );
		}
		//use input seq?
		this_input_seq.erase( n_prepended + this_prepend + 1, n_append - n_appended - this_append );
		this_input_seq.erase( 0, n_prepend - n_prepended - this_prepend );
		if( option[ pepspec::input_seq ].user() ){
			for( Size ii = pep_begin; ii <= pep_end; ++ii ){
				Size this_input_seq_pos( ii - pep_begin );
				char this_input_seq_aa( this_input_seq[ this_input_seq_pos ] );
				if( this_input_seq_aa != 'X' ){
					make_sequence_change( ii, chemical::AA( aa_from_oneletter_code( this_input_seq_aa ) ), pose );
				}
			}
		}

		pose.update_residue_neighbors();

		//turn on homol_csts
		if( option[ pepspec::homol_csts ].user() ){
			cen_scorefxn->set_weight( coordinate_constraint, option[ OptionKeys::constraints::cst_weight ] );
			set_pep_csts( pose );
		}

		//setup pep vectors
		Size nres_pep( pep_end - pep_begin + 1 );
		vector1< bool > is_pep( pose.total_residue(), false );
		for ( Size i = 1; i <= pose.total_residue(); ++i ) is_pep[ i ] = ( i >= pep_begin && i <= pep_end );
//		vector1< bool > is_insert( nres_pep, false );
		kinematics::MoveMapOP mm_frag ( new kinematics::MoveMap );
		for( Size ii = 1; ii <= this_prepend + 1; ++ii ){
			Size moveable_seqpos( pep_begin + ii - 1 );
			mm_frag->set_bb( moveable_seqpos, true );
//			is_insert[ ii ] = true;
		}
		for( Size ii = nres_pep; ii >= nres_pep - this_append; --ii ){
			Size moveable_seqpos( pep_begin + ii - 1 );
			mm_frag->set_bb( moveable_seqpos, true );
//			is_insert[ ii ] = true;
		}

		( *cen_scorefxn )( pose );
		Real best_score( pose.energies().total_energies()[ total_score ] );
		MonteCarloOP mc_frag ( new MonteCarlo( pose, *cen_scorefxn, 2.0 ) );

		//replace prot w/ A/G/P?
		if( option[ pepspec::no_cen_rottrials ] ){
			for( Size seqpos = prot_begin; seqpos <= prot_end; ++seqpos ){
				if( pose.residue( seqpos ).name1() != 'G' && pose.residue( seqpos ).name1() != 'P' ){
					make_sequence_change( seqpos, chemical::aa_ala, pose );
				}
			}
		}

		//phi/psi insertion loop
		if( !option[ pepspec::use_input_bb ] ){
			Size this_build_loop( static_cast< int >( n_build_loop * std::max( this_prepend, this_append ) / std::max( n_prepend, n_append ) ) );

			//get frags
			Size const nfrags( this_build_loop );
			core::fragment::FragSetCOP fragset( make_1mer_frags( pep_begin, pep_end, this_input_seq, nfrags ) );
			protocols::simple_moves::ClassicFragmentMover frag_mover( fragset, mm_frag );

			for ( Size build_loop_inner = 1; build_loop_inner <= this_build_loop; ++build_loop_inner ) {
/*
				//choose an insertion position = seqpos
				Size insert_peppos = 0;
				while( true ){
					insert_peppos = static_cast< int > ( nres_pep * numeric::random::rg().uniform() ) + 1;
					if( is_insert[ insert_peppos ] ) break;
				}
				Size insert_seqpos( pep_begin + insert_peppos - 1 );
*/
				frag_mover.apply( pose );
				//random perturb angles +/- 1 degree
//				pose.set_phi( insert_seqpos, pose.phi( insert_seqpos ) + 1 * ( 2 * numeric::random::rg().uniform() - 1 ) );
//				pose.set_psi( insert_seqpos, pose.psi( insert_seqpos ) + 1 * ( 2 * numeric::random::rg().uniform() - 1 ) );

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
		( *cen_scorefxn )( pose );

		bool this_has_clash( false );

		//check no_cen clash
			//break if has clash already
			if( has_clash( pose, is_pep, cen_scorefxn, clash_cutoff ) ){
				this_has_clash = true;
				//replace prot residues
				for( Size seqpos = prot_begin; seqpos <= prot_end; ++seqpos ){
					if( pose.residue( seqpos ).name1() != 'G' && pose.residue( seqpos ).name1() != 'P' ){
						pose.replace_residue( seqpos, replace_res_pose.residue( seqpos ), true );
					}
				}
			}
			//else replace rotamers and repack if necessary
			else if( option[ pepspec::no_cen_rottrials ] ){
				//ignore preexisting clashes
				vector1< bool > ignore_clash( pose.total_residue(), false );
				for( Size i = 1; i <= pose.total_residue() - 1; ++i ){
					vector1< bool > this_clash( pose.total_residue(), false );
					this_clash[ i ] = true;
					if( has_clash( pose, this_clash, cen_scorefxn, clash_cutoff ) ){
						ignore_clash[ i ] = true;
					}
				}
				//replace prot residues
				for( Size seqpos = prot_begin; seqpos <= prot_end; ++seqpos ){
					pose.replace_residue( seqpos, replace_res_pose.residue( seqpos ), true );
				}
				//find clashes and repack
				if( has_clash( pose, is_pep, cen_scorefxn, clash_cutoff ) ){
					vector1< bool > repack_this( pose.total_residue(), false );
					vector1< std::pair< Size, Size > > clash_pairs( get_clash_pairs( pose, is_pep, cen_scorefxn, clash_cutoff ) );
					//get all prot clash res and their nbrs
					for( Size i = 1; i <= clash_pairs.size(); ++i ){
						std::pair< Size, Size > clash_pair( clash_pairs[ i ] );
						Size this_clash( clash_pair.second );
						EnergyGraph const & energy_graph( pose.energies().energy_graph() );
						for ( graph::Graph::EdgeListConstIter
								ir  = energy_graph.get_node( this_clash )->const_edge_list_begin(),
								ire = energy_graph.get_node( this_clash )->const_edge_list_end();
								ir != ire; ++ir ){
							Size this_nbr( (*ir)->get_other_ind( this_clash ) );
							repack_this[ this_nbr ] = true;
						}
					}
					//and repack
					pack::task::TaskFactoryOP rp_task_factory( new pack::task::TaskFactory );
					pack::task::operation::RestrictResidueToRepackingOP restrict_to_repack_taskop( new pack::task::operation::RestrictResidueToRepacking() );
					pack::task::operation::PreventRepackingOP prevent_repack_taskop( new pack::task::operation::PreventRepacking() );
					for ( Size i=1; i<= pose.total_residue(); ++i ) {
						if ( repack_this[ i ] ) restrict_to_repack_taskop->include_residue( i );
						else prevent_repack_taskop->include_residue( i );
					}
					rp_task_factory->push_back( new pack::task::operation::IncludeCurrent() );
					rp_task_factory->push_back( restrict_to_repack_taskop );
					rp_task_factory->push_back( prevent_repack_taskop );
					protocols::simple_moves::RotamerTrialsMoverOP dz_rottrial ( new protocols::simple_moves::RotamerTrialsMover( cen_scorefxn, rp_task_factory ) );
					dz_rottrial->apply( pose );
					//check clashes one last time
					vector1< bool > check_clash( pose.total_residue(), false );
					for( Size i = 1; i <= pose.total_residue() - 1; ++i ) if( is_pep[ i ] || ( repack_this[ i ] && !ignore_clash[ i ] ) ) check_clash[ i ] = true;

					this_has_clash = has_clash( pose, check_clash, cen_scorefxn, clash_cutoff );
				}
			}
		//check cen clash
		else{
			Real clash_score( pose.energies().total_energies()[ vdw ] );
			this_has_clash = ( clash_score > min_score + 0.03 );
		}

		//if gen clash, step back a step
		if( break_loop < 5 && this_has_clash ){
			++break_loop;
			for( Size ii = 1; ii <= this_prepend; ++ii ){
				pose.conformation().delete_residue_slow( pep_begin );
				pep_anchor = pep_anchor - 1;
				pep_end = pep_end - 1;
				pose.conformation().update_polymeric_connection( pep_begin );
				pose.conformation().update_polymeric_connection( pep_begin + 1 );
			}
			for( Size ii = 1; ii <= this_append; ++ii ){
				pose.conformation().delete_residue_slow( pep_end );
				pep_end = pep_end - 1;
				pose.conformation().update_polymeric_connection( pep_end );
				pose.conformation().update_polymeric_connection( pep_end - 1 );
			}
		}
		else{
			n_prepended += this_prepend;
			n_appended += this_append;
		}
	}
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
	protocols::simple_moves::SmallMoverOP cg_small( new protocols::simple_moves::SmallMover( mm_move, 5.0, 1 ) );
	cg_small->angle_max( 'H', 10.0 );
	cg_small->angle_max( 'E', 10.0 );
	cg_small->angle_max( 'L', 10.0 );

	MonteCarloOP mc_cg ( new MonteCarlo( pose, *cen_scorefxn, 2.0  ) );
	Real best_score( pose.energies().total_energies().dot( cen_scorefxn->weights() ) );
	for( Size cgmove_loop = 1; cgmove_loop <= n_iter; cgmove_loop++ ){
		cg_small->apply( pose );
		if( mc_cg->boltzmann( pose ) ){
			Real test_score( pose.energies().total_energies().dot( cen_scorefxn->weights() ) );
			if( test_score <= best_score ){
				best_score = test_score;
				mc_cg->reset( pose );
			}
		}
	}
	mc_cg->recover_low( pose );
}

void
mutate_random_residue(
	Pose & pose,
	vector1< bool > is_mutable,
	ScoreFunctionOP soft_scorefxn,
	ScoreFunctionOP full_scorefxn
)
{
	//pick rand seqpos
	assert( is_mutable.size() == pose.total_residue() );
	bool mutate( false );
	Size seqpos( 0 );
	while( !mutate ){
		seqpos = static_cast< int >( numeric::random::rg().uniform() * is_mutable.size() + 1 );
		mutate = is_mutable[ seqpos ];
	}

	//pick rand aa and mutate
	make_sequence_change( seqpos, chemical::AA( static_cast< int > ( 20 * numeric::random::rg().uniform() + 1 ) ), pose );

	//rottrial seqpos only
	pack::task::TaskFactoryOP mut_task_factory( new pack::task::TaskFactory );
	{
		pack::task::operation::RestrictResidueToRepackingOP restrict_to_repack_taskop( new pack::task::operation::RestrictResidueToRepacking() );
		pack::task::operation::PreventRepackingOP prevent_repack_taskop( new pack::task::operation::PreventRepacking() );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( i == seqpos ) restrict_to_repack_taskop->include_residue( i );
			else prevent_repack_taskop->include_residue( i );
		}
		mut_task_factory->push_back( new pack::task::operation::InitializeFromCommandline() );
		mut_task_factory->push_back( restrict_to_repack_taskop );
		mut_task_factory->push_back( prevent_repack_taskop );
	}
	protocols::simple_moves::RotamerTrialsMoverOP mut_rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, mut_task_factory ) );
	mut_rottrial->apply( pose );

	//get seqpos nbrs from energy map
	vector1< bool > is_nbr( pose.total_residue(), false );
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( graph::Graph::EdgeListConstIter
					ir  = energy_graph.get_node( seqpos )->const_edge_list_begin(),
					ire = energy_graph.get_node( seqpos )->const_edge_list_end();
				ir != ire; ++ir ) {
		Size i_nbr( (*ir)->get_other_ind( seqpos ) );
		is_nbr[ i_nbr ] = true && ( i_nbr != pep_anchor );
	}

	//rottrial seqpos nbrs too
	pack::task::TaskFactoryOP rp_task_factory( new pack::task::TaskFactory );
	{
		pack::task::operation::RestrictResidueToRepackingOP restrict_to_repack_taskop( new pack::task::operation::RestrictResidueToRepacking() );
		pack::task::operation::PreventRepackingOP prevent_repack_taskop( new pack::task::operation::PreventRepacking() );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( ( i == seqpos || is_nbr[ i ] ) && i != pep_anchor ) restrict_to_repack_taskop->include_residue( i );
			else prevent_repack_taskop->include_residue( i );
		}
		rp_task_factory->push_back( new pack::task::operation::InitializeFromCommandline() );
		rp_task_factory->push_back( new pack::task::operation::IncludeCurrent() );
		rp_task_factory->push_back( restrict_to_repack_taskop );
		rp_task_factory->push_back( prevent_repack_taskop );
	}
	protocols::simple_moves::RotamerTrialsMoverOP rp_rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, rp_task_factory ) );
	rp_rottrial->apply( pose );

/*
	if( pose.residue( seqpos ).name3() != "PRO" && pose.residue( seqpos ).name3() != "GLY" && pose.residue( seqpos ).name3() != "ALA" ){
		MonteCarloOP mc( new MonteCarlo( pose, *soft_scorefxn, 1.0 ) );
		protocols::simple_moves::sidechain_moves::SidechainMoverOP sidechain_mover( new protocols::simple_moves::sidechain_moves::SidechainMover() );
		sidechain_mover->set_task_factory( rp_task_factory );
		sidechain_mover->set_prob_uniform( 0.05 );
		for( Size ii = 1; ii <= 20; ++ii ){
			sidechain_mover->apply( pose );
			mc->boltzmann( pose );
		}
		mc->recover_low( pose );
	}
*/

	kinematics::MoveMapOP mm_min ( new kinematics::MoveMap );
	mm_min->set_chi( seqpos );
	mm_min->set_chi( is_nbr );
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "dfpmin", 0.001, true );
	min_mover->apply( pose );
}

void
packmin_unbound_pep(
	pose::Pose & pose,
	scoring::ScoreFunctionOP full_scorefxn
)
{
	pose.update_residue_neighbors();
	( *full_scorefxn )( pose );
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().or_include_current( true );
	task->restrict_to_repacking();
	protocols::simple_moves::PackRotamersMoverOP pack( new protocols::simple_moves::PackRotamersMover( full_scorefxn, task, 1 ) );
	pack->apply( pose );
	kinematics::MoveMapOP mm ( new kinematics::MoveMap );
	mm->set_chi( true );
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm, full_scorefxn, "dfpmin", 0.001, true );
	min_mover->apply( pose );
}

Real
get_binding_score(
	Pose pose,
	Size pep_chain,
	ScoreFunctionOP full_scorefxn
)
{
	Pose pep_pose( pose.split_by_chain( pep_chain ) );
	Real total_score( pose.energies().total_energies().dot( full_scorefxn->weights() ) );
	packmin_unbound_pep( pep_pose, full_scorefxn );
	( *full_scorefxn )( pep_pose );
	Real pep_score = pep_pose.energies().total_energies().dot( full_scorefxn->weights() );
	Real bind_score = total_score - pep_score;
	return bind_score;
}


void
print_pep_analysis(
	std::string pdb_name,
	std::fstream & out_file,
	Pose pose,
	Real prot_score,
	ScoreFunctionOP full_scorefxn,
	bool dump_pdb
)
{
	//must remove constraints!
	( *full_scorefxn )( pose );
	std::string output_seq;
	for(Size i = pep_begin; i <= pep_end; i++){
		output_seq.append( 1, pose.residue( i ).name1() );
	}
	Real total_score( pose.energies().total_energies().dot( full_scorefxn->weights() ) );
	out_file << pdb_name + "\t" << output_seq << "\t" << pose.energies().total_energies().weighted_string_of( full_scorefxn->weights() );
	out_file<<"\ttotal_score:\t"<<total_score<<"\t";
	out_file<<"total-prot_score:\t"<<total_score - prot_score<<"\t";
	if( option[ pepspec::homol_csts ].user() ){
		full_scorefxn->set_weight( coordinate_constraint, option[ OptionKeys::constraints::cst_weight ] );
		( *full_scorefxn )( pose );
		out_file << "coordinate_constraint:\t" << pose.energies().total_energies()[ coordinate_constraint ] << "\t";
		out_file << "total_cst_score:\t" << pose.energies().total_energies().dot( full_scorefxn->weights() ) << "\t";
//		pose.constraint_set( 0 );
		full_scorefxn->set_weight( coordinate_constraint, 0 );
		( *full_scorefxn )( pose );
	}
	//score pep unbound state//
	if( option[ pepspec::binding_score ] ){
		Pose pep_pose( pose.split_by_chain( pep_chain ) );
		packmin_unbound_pep( pep_pose, full_scorefxn );
		( *full_scorefxn )( pep_pose );
		Real pep_score = pep_pose.energies().total_energies().dot( full_scorefxn->weights() );
		Real bind_score = total_score - pep_score;
		out_file<<"binding_score:\t"<<bind_score<<"\t";
		out_file<<"binding-prot_score:\t"<<bind_score - prot_score<<"\t";
//		TR << pdb_name + ".pep\t" << output_seq << "\t" << pep_pose.energies().total_energies().weighted_string_of( full_scorefxn->weights() ) + "\ttotal_score:\t" << pep_score << "\n";
	}
	if( option[ pepspec::upweight_interface ] ){
			basic::MetricValue< Real > weighted_interface_score;
			pose.metric("energy_decomposition", "weighted_total", weighted_interface_score);
			out_file << "reweighted_score:\t" + string_of( weighted_interface_score.value() ) + "\t";
	}
	basic::MetricValue< Real > interface_score;
	pose.metric( "interface_delta_energies", "weighted_total", interface_score );
	out_file << "interface_score:\t" + string_of( interface_score.value() ) + "\t";
	if( option[ pepspec::calc_sasa ] ){
			basic::MetricValue< Real > delta_sasa;
			pose.metric( "sasa_interface", "delta_sasa", delta_sasa );
			out_file << "interface_sasa:\t" + string_of( delta_sasa.value() ) + "\t";
	}
//	if( option[ pepspec::rmsd_analysis ] ) out_file << pep_rmsd_analysis( pose );
//	if( option[ pepspec::phipsi_analysis ] ) out_file << pep_phipsi_analysis( pose );
	out_file<<std::endl;
	if( dump_pdb ) pose.dump_scored_pdb( pdb_name, *full_scorefxn );

}

void
RunPepSpec()
{
	//load loop lengths//
	Size n_peptides( option[ pepspec::n_peptides ] );
	Size n_cgrelax_loop( option[ pepspec::n_cgrelax_loop ] );

	//data out
	std::string out_nametag( "data" );
	if( option[ out::file::o ].user() ) out_nametag = option[ out::file::o ];
	std::string pdb_dir( out_nametag + ".pdbs" );
	std::string pdb_name;
	utility::file::create_directory( pdb_dir );
	std::string out_file_name_str( out_nametag + ".spec" );
	char const *out_file_name = out_file_name_str.c_str();
	std::fstream out_file( out_file_name, std::ios::out );

	//load input structure pdbs/csts//
	pose::Pose pose;
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

	Size n_prepend( option[ pepspec::n_prepend ] );
	Size n_append( option[ pepspec::n_append ] );

	bool const save_low_pdbs( option[ pepspec::save_low_pdbs ] );
	bool const save_all_pdbs( option[ pepspec::save_all_pdbs ] );

	//load homol_csts?
	if( option[ pepspec::homol_csts ].user() ){
		//parse cst line as: atom_name pep_pos x y z x0 sd tol
		std::string cst_filename( option[ pepspec::homol_csts ] );
		std::ifstream cst_file( cst_filename.c_str() );
		std::string str_buffer;
		while( !getline( cst_file, str_buffer, '\t' ).eof() ) {
			pep_coord_cst pep_cst;

			pep_cst.atom_name = str_buffer;
			getline( cst_file, str_buffer, '\t' );
			{
				std::istringstream istr_buffer( str_buffer );
				istr_buffer >> pep_cst.pep_pos;
			}
			getline( cst_file, str_buffer, '\t' );
			{
				std::istringstream istr_buffer( str_buffer );
				istr_buffer >> pep_cst.x;
			}
			getline( cst_file, str_buffer, '\t' );
			{
				std::istringstream istr_buffer( str_buffer );
				istr_buffer >> pep_cst.y;
			}
			getline( cst_file, str_buffer, '\t' );
			{
				std::istringstream istr_buffer( str_buffer );
				istr_buffer >> pep_cst.z;
			}
			getline( cst_file, str_buffer, '\t' );
			{
				std::istringstream istr_buffer( str_buffer );
				istr_buffer >> pep_cst.x0;
			}
			getline( cst_file, str_buffer, '\t' );
			{
				std::istringstream istr_buffer( str_buffer );
				istr_buffer >> pep_cst.sd;
			}
			getline( cst_file, str_buffer, '\n' );
			{
				std::istringstream istr_buffer( str_buffer );
				istr_buffer >> pep_cst.tol;
			}
			pep_coord_csts.push_back( pep_cst );
		}
	}


	//define scoring functions//
	core::scoring::ScoreFunctionOP full_scorefxn( get_score_function() );
	core::scoring::ScoreFunctionOP soft_scorefxn( ScoreFunctionFactory::create_score_function( option[ pepspec::soft_wts ] ) );
	core::scoring::ScoreFunctionOP cen_scorefxn;
	if(  option[ pepspec::cen_wts ].user() ){
			cen_scorefxn =  ScoreFunctionFactory::create_score_function( option[ pepspec::cen_wts ] );
	}
	else{
			cen_scorefxn->set_weight( fa_rep, 0.5 );
	}
//	Real const max_interchain_wt( cen_scorefxn->get_weight( interchain_contact ) );

	//////////////fragment library/////////////////
/*
	VallLibrarian librarian;
	librarian.add_fragment_gen( new LengthGen( 1 ) );
	librarian.catalog( FragmentLibraryManager::get_instance()->get_Vall() );
*/


	// set up the pose_metric_calculators
	core::pose::metrics::CalculatorFactory & calculator_factory(core::pose::metrics::CalculatorFactory::Instance());
	protocols::toolbox::pose_metric_calculators::ResidueDecompositionByChainCalculatorOP res_decomp_calculator(new protocols::toolbox::pose_metric_calculators::ResidueDecompositionByChainCalculator());
	calculator_factory.register_calculator("residue_decomposition", res_decomp_calculator);
	protocols::toolbox::pose_metric_calculators::DecomposeAndReweightEnergiesCalculatorOP decomp_reweight_calculator( new protocols::toolbox::pose_metric_calculators::DecomposeAndReweightEnergiesCalculator( "residue_decomposition" ) );
	vector1< Real > reweight_vector( 1.0, 3 ); reweight_vector.push_back( 2.0 );
	decomp_reweight_calculator->master_weight_vector( reweight_vector );
	decomp_reweight_calculator->num_sets( ( Size ) 2 );
	calculator_factory.register_calculator("energy_decomposition", decomp_reweight_calculator);
	core::pose::metrics::PoseMetricCalculatorOP int_calculator( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( (Size)1, (Size)2 ) );
	calculator_factory.register_calculator( "interface", int_calculator );
	core::pose::metrics::PoseMetricCalculatorOP int_delta_energy_calculator( new core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator( "interface" ) );
	calculator_factory.register_calculator( "interface_delta_energies", int_delta_energy_calculator );
	core::pose::metrics::PoseMetricCalculatorOP int_sasa_calculator( new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator( (Size)1, (Size)2 ) );
	calculator_factory.register_calculator( "sasa_interface", int_sasa_calculator );

	int pose_index( 0 );
	if( option[ pepspec::run_sequential ] ) n_peptides = pdb_filenames.size();
	for( Size peptide_loop = 1; peptide_loop <= n_peptides; ++peptide_loop ){

		/*
		VallData vall( option[ pepspec::frag_file ] );
		TorsionFragmentLibrary lib;
		lib.resize( 20 );
		//lib indexed same as AA enum
		std::string ss_seq( "L" );
		if( option[ pepspec::ss_type ].user() ) ss_seq = option[ pepspec::ss_type ];
		for( Size i = 1; i <= 20; ++i ){
			char aa( oneletter_code_from_aa( chemical::AA( i ) ) );
			std::string frag_seq;
			frag_seq.push_back( aa );
			Real seq_wt( 0.1 );
			Real ss_wt( 0.1 );
			vall.get_frags( 1000, frag_seq, ss_seq, seq_wt, ss_wt, false, false, true, lib[ i ] );
		}
		*/

		//load random start pdb//
		if( option[ pepspec::run_sequential ] ) ++pose_index;
		else pose_index = static_cast< int >( numeric::random::rg().uniform() * pdb_filenames.size() + 1 );
		std::string pdb_filename( pdb_filenames[ pose_index ] );
		TR<<"Initializing "<< out_nametag + "_" + string_of( peptide_loop ) + " with " + pdb_filename << std::endl;
		import_pose::pose_from_pdb( pose, pdb_filename );

//		ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
		//convert user input to internal values//
		if( option[ pepspec::pep_anchor ].user() ){
			Size const pep_anchor_in( option[ pepspec::pep_anchor ] );
			std::string const pep_chain_in( option[ pepspec::pep_chain ] );
			pep_anchor = pose.pdb_info()->pdb2pose( pep_chain_in[ 0 ], pep_anchor_in );
			pep_chain = pose.chain( pep_anchor );
			pep_begin = pose.conformation().chain_begin( pep_chain );
			pep_end = pose.conformation().chain_end( pep_chain );
		}
		else{
			pep_chain = 2;
			pep_begin = pose.conformation().chain_begin( pep_chain );
			pep_end = pose.conformation().chain_end( pep_chain );
			pep_anchor = pep_begin + static_cast< int >( ( pose.conformation().chain_end( pep_chain ) -  pose.conformation().chain_begin( pep_chain ) ) * numeric::random::rg().uniform() );
		}

		if( pep_anchor == 0 ) utility_exit_with_message( "pep_chain / pep_anchor combo not found\n" );
		for( Size i = 1; i <= pose.conformation().num_chains(); ++i ){
			if( i == pep_chain ) continue;
			if( !( pose.residue( pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
			else{
				prot_chain = i;
				break;
			}
		}
		prot_begin = pose.conformation().chain_begin( prot_chain );
		prot_end =pose.conformation().chain_end( prot_chain );
		prot_anchor =prot_begin;

		Real cutoff( option[ pepspec::interface_cutoff ] );
		std::string output_seq;

		//gen fold tree//
		FoldTree f( pose.total_residue() );

		std::string pep_anchor_root( "CB" );
		if( pose.residue( pep_anchor ).name1() == 'G' ) pep_anchor_root = "CA";
		pep_jump = 2;
		if( prot_chain < pep_chain ){
			pep_jump = f.new_jump( prot_anchor, pep_anchor, pep_begin - 1 );
			f.set_jump_atoms( pep_jump, "CA", pep_anchor_root );
			if( pep_end != pose.total_residue() ) f.new_jump( prot_anchor, pep_end + 1, pep_end );
			if( prot_end + 1 != pep_begin ) f.new_jump( prot_anchor, prot_end + 1, prot_end );
		}
		else{
			pep_jump = f.new_jump( pep_anchor, prot_anchor, prot_begin - 1 );
			f.set_jump_atoms( pep_jump, pep_anchor_root, "CA" );
			if( prot_end != pose.total_residue() ) f.new_jump( prot_anchor, prot_end + 1, prot_end );
			if( pep_end + 1 != prot_begin ) f.new_jump( pep_anchor, pep_end + 1, pep_end );
		}
		f.reorder( prot_anchor );
		TR << f << "\n";

		//set foldtree in the pose//
		pose.fold_tree( f );

		pose::Pose start_pose( pose );
		Pose prot_pose = pose.split_by_chain( prot_chain );
		( *full_scorefxn )( prot_pose );
		packmin_unbound_pep( prot_pose, full_scorefxn );
		Real prot_score( prot_pose.energies().total_energies().dot( full_scorefxn->weights() ) );

		protocols::viewer::add_conformation_viewer( pose.conformation(), "pepspec_pose" );

		//init peptide
		initialize_peptide( pose );

		//load input seq?
		input_seq = std::string( n_prepend + n_append + 1, 'X' );
		if( option[ pepspec::input_seq ].user() ) input_seq = std::string( option[ pepspec::input_seq ] );

		bool add_nterm( true );
		if( pep_anchor == pep_begin ) add_nterm = false;
		bool add_cterm( true );
		if( pep_anchor == pep_end ) add_cterm = false;

		Residue pep_anchor_res = pose.residue( pep_anchor );

		//convert to CG residues//
		if( !option[ pepspec::no_cen ] ) core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
		( *cen_scorefxn )( pose );

		Real min_vdw( pose.energies().total_energies()[ vdw ] );

		gen_pep_bb_sequential( pose, cen_scorefxn );

		if( option[ pepspec::homol_csts ].user() ){
			cen_scorefxn->set_weight( coordinate_constraint, option[ OptionKeys::constraints::cst_weight ] );
			set_pep_csts( pose );
		}

		( *cen_scorefxn )( pose );
		if( peptide_loop > 1 && pose.energies().total_energies()[ vdw ] > min_vdw + 0.1 ){
			--peptide_loop;
			continue;
		}

		//define prot, pep, and anchor residues//
		vector1< bool > is_pep( pose.total_residue(), false ), is_mutable( pose.total_residue(), false ), is_prot( pose.total_residue(), false ), is_anchor( pose.total_residue(), false );
		vector1< Size > pep_res_vec;
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			is_pep[i] = ( i >= pep_begin && i <= pep_end );
			if( is_pep[ i ] ) pep_res_vec.push_back( i );
			is_prot[i] = ( i >= prot_begin && i <= prot_end );
			is_anchor[i] = ( ( i == pep_anchor || i == prot_anchor ) );
			is_mutable[i] = ( i >= pep_begin && i <= pep_end && !(i==pep_anchor) );
			if( option[ pepspec::input_seq ].user() && i >= pep_begin ) is_mutable[i] = ( input_seq[ i - pep_begin ] == 'X' && i != pep_anchor );
		}

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

			//small-small movemap//
			kinematics::MoveMapOP mm_move ( new kinematics::MoveMap );
			mm_move->set_bb( is_pep );

			//small CG backbone moves//
			Size n_cgperturb_iter( 100 );
			if( option[ pepspec::use_input_bb ] ) n_cgperturb_iter = option[ pepspec::n_build_loop ];
			if( cgrelax_loop > 1 || option[ pepspec::use_input_bb ] ){
				perturb_pep_bb( pose, mm_move, cen_scorefxn, n_cgperturb_iter );
			}

			//debug dump CG data and pdb//
			if( option[ pepspec::dump_cg_bb ] ){
				//remove buffer res
				std::string pdb_name( pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + ".pdb" );
				if( n_cgrelax_loop > 1 ) pdb_name = pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + "_" + string_of( cgrelax_loop ) + ".pdb";
				if( option[ pepspec::save_low_pdbs ] ) pose.dump_scored_pdb( pdb_name, *full_scorefxn );
				( *cen_scorefxn )( pose );
				out_file << pdb_name + "\t"<<pose.energies().total_energies().weighted_string_of( cen_scorefxn->weights() );
				out_file<<"\ttotal_score:\t"<<pose.energies().total_energies()[ total_score ]<<"\t";
//				if( option[ pepspec::rmsd_analysis ] ) out_file << pep_rmsd_analysis( pose );
//				if( option[ pepspec::phipsi_analysis ] ) out_file << pep_phipsi_analysis( pose );
				out_file<<std::endl;
				continue;
			}

			//switch back to fullatom
			if( option[ pepspec::add_buffer_res ] ){
				core::pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_BB", pep_begin );
				core::pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_BB", pep_end );
			}
			core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
			if( option[ pepspec::add_buffer_res ] ){
				core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_BB", pep_begin );
				core::pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_BB", pep_end );
			}

			//replace prot residues w/ original rotamers
			for(Size resnum = prot_begin; resnum <= prot_end; ++resnum){
				pose.replace_residue( resnum, start_pose.residue( resnum ), true );
			}
			pose.replace_residue( pep_anchor, pep_anchor_res, true );

			//replace termini
			add_termini( pose );
			//turn off csts//
			full_scorefxn->set_weight( coordinate_constraint, 0 );

			//randomize pep sequence//
			if( !( option[ pepspec::no_design ] || option[ pepspec::input_seq ].user() ) ){
				for( Size mut_site = pep_begin; mut_site <= pep_end; ++mut_site ){
					if( mut_site==pep_anchor ) continue;
					int resindex;
					resindex = static_cast< int > ( 20 * numeric::random::rg().uniform() + 1 );
					make_sequence_change( mut_site, chemical::AA(resindex), pose );
				}
			}


			//define neighbors
			vector1< bool > is_pep_nbr( pose.total_residue(), false );
			vector1< Size > nbr_res_vec;
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
								nbr_res_vec.push_back( i );
								break;
							}
						}
					}
				}
			}
			mm_min->set_chi( is_pep_nbr );

			//define movers//
			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( mm_min, full_scorefxn, "dfpmin", 0.001, true );

			MonteCarloOP mc_relax ( new MonteCarlo( pose, *full_scorefxn, 1.0 ) );

			Real bind_score( get_binding_score( pose, pep_chain, full_scorefxn ) );
			myMC mc_bind( pose, bind_score, 1.0 );

			//define design task and repack task
			pack::task::TaskFactoryOP dz_task_factory( new pack::task::TaskFactory );
			{
				pack::task::operation::RestrictResidueToRepackingOP restrict_to_repack_taskop( new pack::task::operation::RestrictResidueToRepacking() );
				pack::task::operation::PreventRepackingOP prevent_repack_taskop( new pack::task::operation::PreventRepacking() );
				for ( Size i=1; i<= pose.total_residue(); ++i ) {
					if ( is_pep[ i ] && i != pep_anchor ) {
						if( option[ pepspec::no_design ] || ( option[ pepspec::input_seq ].user() && input_seq[ i - pep_begin ] != 'X' ) ){
							restrict_to_repack_taskop->include_residue( i );
						}
					} else if ( is_pep_nbr[i] && is_prot[i] ) {
						restrict_to_repack_taskop->include_residue( i );
					} else {
						prevent_repack_taskop->include_residue( i );
					}
				}
				if( !option[ pepspec::diversify_pep_seqs ].user() ) dz_task_factory->push_back( new pack::task::operation::InitializeFromCommandline() );
				dz_task_factory->push_back( new pack::task::operation::IncludeCurrent() );
				dz_task_factory->push_back( restrict_to_repack_taskop );
				dz_task_factory->push_back( prevent_repack_taskop );
			}
			std::string wts_name( option[ score::weights ] ), soft_wts_name( option[ pepspec::soft_wts ] );
			//full_wts
			{
				pack::task::PackerTaskOP dz_task( dz_task_factory->create_task_and_apply_taskoperations( pose ));

				//upweight interface?
				if( option[ pepspec::upweight_interface ] ){
					Real upweight( 2.0 );
					core::pack::task::IGEdgeReweightContainerOP IGreweight = dz_task->set_IGEdgeReweights();
					core::pack::task::IGEdgeReweighterOP upweighter = new protocols::toolbox::ResidueGroupIGEdgeUpweighter( upweight, pep_res_vec, nbr_res_vec  );
					IGreweight->add_reweighter( upweighter );
				}

				protocols::simple_moves::PackRotamersMoverOP dz_pack( new protocols::simple_moves::PackRotamersMover( full_scorefxn, dz_task, 1 ) );
				protocols::simple_moves::RotamerTrialsMoverOP dz_rottrial ( new protocols::simple_moves::RotamerTrialsMover( full_scorefxn, dz_task_factory ) );
				SequenceMoverOP design_seq = new SequenceMover;
				design_seq->add_mover( dz_pack );
				design_seq->add_mover( dz_rottrial );
				design_seq->add_mover( min_mover );
				( *full_scorefxn )( pose );
				design_seq->apply( pose );
				if( option[ pepspec::binding_score ] ){
					bind_score = get_binding_score( pose, pep_chain, soft_scorefxn );
					mc_bind.roll( pose, bind_score );
				}
				else mc_relax->boltzmann( pose );
				pdb_name = pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + "_full";
				print_pep_analysis( pdb_name, out_file, pose, prot_score, full_scorefxn, save_all_pdbs );
			}
			//design again with soft_wts
			if( option[ pepspec::soft_wts ].user() ){
				pack::task::PackerTaskOP dz_task( dz_task_factory->create_task_and_apply_taskoperations( pose ));

				//upweight interface?
				if( option[ pepspec::upweight_interface ] ){
					Real upweight( 2.0 );
					core::pack::task::IGEdgeReweightContainerOP IGreweight = dz_task->set_IGEdgeReweights();
					core::pack::task::IGEdgeReweighterOP upweighter = new protocols::toolbox::ResidueGroupIGEdgeUpweighter( upweight, pep_res_vec, nbr_res_vec  );
					IGreweight->add_reweighter( upweighter );
				}

				protocols::simple_moves::PackRotamersMoverOP dz_pack( new protocols::simple_moves::PackRotamersMover( soft_scorefxn, dz_task, 1 ) );
				protocols::simple_moves::RotamerTrialsMoverOP dz_rottrial ( new protocols::simple_moves::RotamerTrialsMover( soft_scorefxn, dz_task_factory ) );
				SequenceMoverOP design_seq = new SequenceMover;
				design_seq->add_mover( dz_pack );
				design_seq->add_mover( dz_rottrial );
				design_seq->add_mover( min_mover );
				( *soft_scorefxn )( pose );
				design_seq->apply( pose );
				if( option[ pepspec::binding_score ] ){
					bind_score = get_binding_score( pose, pep_chain, soft_scorefxn );
					mc_bind.roll( pose, bind_score );
				}
				else mc_relax->boltzmann( pose );
				pdb_name = pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + "_soft";
				print_pep_analysis( pdb_name, out_file, pose, prot_score, soft_scorefxn, save_all_pdbs );
			}
			//now try MCM desin against binding score
			if( option[ pepspec::diversify_pep_seqs ] ){
				Size iter_per_res( option[ pepspec::diversify_lvl ] );
				for( Size ii = 1; ii <= iter_per_res * ( pep_end - pep_begin + 1 ); ++ii ){
					mutate_random_residue( pose, is_mutable, soft_scorefxn, full_scorefxn );
					if( option[ pepspec::binding_score ] ){
						bind_score = get_binding_score( pose, pep_chain, full_scorefxn );
						mc_bind.roll( pose, bind_score );
						if( mc_bind.accept() ){
							pdb_name = pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + "_mc" + string_of( ii );
							print_pep_analysis( pdb_name, out_file, pose, prot_score, full_scorefxn, save_all_pdbs );
						}
					}
					else if( mc_relax->boltzmann( pose ) ){
						pdb_name = pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + "_mc" + string_of( ii );
						print_pep_analysis( pdb_name, out_file, pose, prot_score, full_scorefxn, save_all_pdbs );
					}
				}
			}
			//recover best
			mc_relax->recover_low( pose );
			pdb_name = pdb_dir + "/" + out_nametag + "_" + string_of( peptide_loop ) + ".pdb";
			print_pep_analysis( pdb_name, out_file, pose, prot_score, full_scorefxn, save_low_pdbs );

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
		return -1;
	}

	return 0;
}


