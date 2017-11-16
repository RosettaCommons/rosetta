// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

//core library
#include <math.h>
#include <fstream>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.hh>

#include <core/chemical/util.hh>
#include <numeric/random/random.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/scoring/types.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/sasa.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/tools/make_map.hh>
#include <utility/string_util.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rtmin.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>

//protocols library (Movers)
#include <devel/init.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>
//#include <protocols/moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/scoring/dssp/Dssp.hh>

//calculator stuff
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/MetricValueGetter.hh>
#include <basic/MetricValue.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>


//utilities

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rot_anl.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

//local options
namespace basic { namespace options { namespace OptionKeys {
}}}//basic::options::OptionKeys

////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace basic;
using namespace id;
using namespace pose;
using namespace pack;
using namespace conformation;
using namespace chemical;
using namespace kinematics;
using namespace scoring;
using namespace options;
using namespace basic::options::OptionKeys;
using namespace optimization;
namespace OK = OptionKeys;
using utility::vector1;
using std::string;
using import_pose::pose_from_pdb;
// deprecated though
using namespace ObjexxFCL;
using basic::Warning;
using basic::Error;

static basic::Tracer TR( "chrisk" );


//local options
namespace chrisk
{
basic::options::StringOptionKey tag( "chrisk:tag" );
basic::options::StringOptionKey insert_pdb( "chrisk:insert_pdb" );
basic::options::RealOptionKey rt_tol( "chrisk:rt_tol" );
basic::options::RealOptionKey rmsd_tol( "chrisk:rmsd_tol" );
basic::options::IntegerOptionKey seqsep_tol( "chrisk:seqsep_tol" );
basic::options::IntegerOptionKey term_tol( "chrisk:term_tol" );
basic::options::IntegerOptionKey trim_tol( "chrisk:trim_tol" );
basic::options::IntegerOptionKey nterm_in( "chrisk:nterm_in" );
basic::options::IntegerOptionKey cterm_in( "chrisk:cterm_in" );
basic::options::BooleanOptionKey loops_only( "chrisk:loops_only" );
basic::options::RealOptionKey design_radius( "chrisk:design_radius" );
basic::options::RealOptionKey repack_radius( "chrisk:repack_radius" );
basic::options::IntegerOptionKey n_add_linker( "chrisk:n_add_linker" );
}

/*
// domain pose class is just a pose plus a vector of ResidueOP's
class DomainPose : public core::pose::Pose{
public:
DomainPose(
core::pose::Pose src_pose,
vector1< ResidueOP > src_map
) :
src_map_( src_map )
{}

ResidueOP src_res(
Size seqpos
){
return src_map_[ seqpos ];
}

vector1< ResidueOP > src_map()
{
return src_map_;
}

private:
vector1< ResidueOP > src_map_;

}
*/

RT
get_atomatom_rt(
	pose::Pose const & pose,
	Size seqpos1,
	string atom1,
	Size seqpos2,
	string atom2
){
	assert( seqpos1 < seqpos2 );
	Residue rsd1( pose.residue( seqpos1 ) );
	Residue rsd2( pose.residue( seqpos2 ) );
	//get atom ids
	core::id::AtomID id1( rsd1.atom_index( atom1 ), seqpos1 );
	core::id::AtomID id2( rsd2.atom_index( atom2 ), seqpos2 );
	//get stubs for Nterm and Cterm
	AtomTree const & at ( pose.conformation().atom_tree() );
	Stub stub1( at.atom( id1 ).get_stub() );
	Stub stub2( at.atom( id2 ).get_stub() );
	kinematics::RT rt( stub1, stub2 );
	return rt;
}

void
find_backbone_rts_in_pose(
	Pose const & pose,
	kinematics::RT const & rt,
	Real const rt_tol,
	Size const seqsep_tol,
	Size const term_tol,
	bool loops_only,
	vector1< std::pair< Size, Size > > & cutpairs )
{
	for ( Size i = term_tol;
			i <= pose.size() - seqsep_tol; ++i ) {
		if ( loops_only && pose.secstruct( i ) != 'L' ) continue;
		for ( Size j = i+1; j <= i + seqsep_tol; ++j ) {
			if ( j > pose.size() - term_tol ) break;
			if ( loops_only && pose.secstruct( j ) != 'L' ) continue;
			//is rt from res i to j within rt_tol of given rt?
			std::pair< Size, Size > cuts( i, j );
			kinematics::RT this_rt( get_atomatom_rt( pose, cuts.first, "CA", cuts.second, "CA" ) );
			if ( this_rt.distance_squared( rt ) > rt_tol ) continue;
			cutpairs.push_back( cuts );
			TR << "Cutpoint pair:\t" << i << "\t" << j << std::endl;
		}
	}
}

//align 2 noncontig stretches of residues on bb heavyatoms
//the two seqpos vectors are assumed to be aligned!
Real
calc_nonlocal_segment_bb_rmsd(
	Pose const & pose1,
	vector1< Size > const seqpos1s,
	Pose & pose2,
	vector1< Size > const seqpos2s
){
	assert( seqpos1s.size() == seqpos2s.size() );
	//create the atomid map
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose2, id::AtomID::BOGUS_ATOM_ID() ); // maps every atomid to bogus atom
	for ( Size iseq = 1; iseq <= seqpos1s.size(); ++iseq ) {
		Size seqpos1( seqpos1s[ iseq ] );
		Size seqpos2( seqpos2s[ iseq ] );
		//atom index 1-->4 is N,CA,C,O
		for ( Size iatm = 1; iatm <= 4; ++iatm ) {
			id::AtomID const id1( iatm, seqpos1 );
			id::AtomID const id2( iatm, seqpos2 );
			atom_map[ id2 ] = id1;
		}
	}
	return( superimpose_pose( pose2, pose1, atom_map ) );
}

//for a given set of residues in src_pose, find all residue sets that match in aln_pose
//generalized to use segments greater than length 1 in the future
void
find_nonlocal_segment_alignments_in_pose(
	Pose & pose,
	Pose const & src_pose,
	vector1< Size > const src_nonlocal_segment,
	Real const rmsd_tol,
	Size const seqsep_tol,
	Size const term_tol,
	bool loops_only,
	vector1< vector1< Size > > & aln_nonlocal_segments )
{
	for ( Size i = term_tol;
			i <= pose.size() - seqsep_tol; ++i ) {
		if ( loops_only && pose.secstruct( i ) != 'L' ) continue;
		for ( Size j = i + 1; j <= i + seqsep_tol; ++j ) {
			if ( j > pose.size() - term_tol ) break;
			if ( loops_only && pose.secstruct( j ) != 'L' ) continue;
			//is rmsd for these rsds within tol when aligned w/ rsds in src pair?
			vector1< Size > aln_nonlocal_segment;
			aln_nonlocal_segment.push_back( i );
			aln_nonlocal_segment.push_back( j );
			Real this_rmsd( calc_nonlocal_segment_bb_rmsd(
				src_pose, src_nonlocal_segment, pose, aln_nonlocal_segment ) );
			if ( this_rmsd > rmsd_tol ) continue;
			aln_nonlocal_segments.push_back( aln_nonlocal_segment );
			TR << "Cutpoint pair:\t" << aln_nonlocal_segment[ 1 ] << "\t"
				<< aln_nonlocal_segment[ 2 ] << "\trmsd:\t" << this_rmsd << std::endl;
		}
	}
}

//for a given cutpair in src_pose, find all cutpairs that match in aln_pose
//just an easy interface b/t std::pair and vector1
void
find_cutpair_alignments_in_pose(
	Pose & pose,
	Pose const & src_pose,
	std::pair< Size, Size > const src_cutpair,
	Real const rmsd_tol,
	Size const seqsep_tol,
	Size const term_tol,
	bool loops_only,
	vector1< std::pair< Size, Size > > & aln_cutpairs
){
	vector1< Size > src_nonlocal_segment;
	src_nonlocal_segment.push_back( src_cutpair.first );
	src_nonlocal_segment.push_back( src_cutpair.second );
	vector1< vector1< Size > > aln_nonlocal_segments;
	find_nonlocal_segment_alignments_in_pose( pose, src_pose, src_nonlocal_segment,
		rmsd_tol, seqsep_tol, term_tol, loops_only, aln_nonlocal_segments );
	for ( Size ipair = 1; ipair <= aln_nonlocal_segments.size(); ++ipair ) {
		//TODO: assuming is length 2 nonloc segment
		aln_cutpairs.push_back( std::pair< Size, Size >(
			aln_nonlocal_segments[ ipair ][ 1 ], aln_nonlocal_segments[ ipair ][ 2 ] ) );
	}
}

// setup secstruct via dssp (copied from rhiju)
void
setup_secstruct_dssp( pose::Pose & pose )
{
	scoring::dssp::Dssp dssp( pose );
	FArray1D_char dssp_secstruct( pose.size() );
	dssp.dssp_reduced( dssp_secstruct );
	for ( Size i = 1; i <= pose.size(); i++ ) {
		pose.set_secstruct(i,  dssp_secstruct(i) );
	}

}

void
copy_segment_bb(
	Pose & pose,
	Size length,
	Pose const & src_pose,
	Size pose_start,
	Size src_pose_start )
{
	for ( Size i = 1; i <= length; ++i ) {
		Size seqpos( pose_start + i - 1 );
		Size src_seqpos( src_pose_start + i - 1 );
		pose.set_phi( seqpos, src_pose.phi( src_seqpos ) );
		pose.set_psi( seqpos, src_pose.psi( src_seqpos ) );
		pose.set_omega( seqpos, src_pose.omega( src_seqpos ) );
	}
}

//return the jump that would connect 2 bonded polymer res
Jump
get_rsdrsd_jump(
	Pose pose,
	Size const seqpos1,
	Size const seqpos2
){
	assert( seqpos2 <= pose.size() );
	assert( seqpos1 <= seqpos2 );
	kinematics::FoldTree ft( pose.fold_tree() );
	ft.new_jump( seqpos1, seqpos2, seqpos2 - 1 );
	ft.set_jump_atoms( 1, "CA", "CA" );
	pose.fold_tree( ft );
	// TR << "JUMP\t" << pose.fold_tree().upstream_atom( 1 ) << "\tto\t" << pose.fold_tree().downstream_atom( 1 ) << std::endl;
	return pose.jump( 1 );
}

//remove terminal residues from pose
//nterm, cterm are not deleted
void
trim_pose_termini(
	Pose & pose,
	Size const nterm,
	Size const cterm
){
	assert( nterm <+ cterm );
	assert( cterm <= pose.size() );
	//delete Cterminus
	if ( cterm < pose.size() ) {
		Size const cterm_orig( pose.conformation().chain_end( pose.chain( cterm ) ) );
		for ( Size iseq = cterm_orig; iseq > cterm; --iseq ) {
			pose.delete_polymer_residue( iseq );
		}
	}
	//and Nterm
	if ( nterm > 1 ) {
		//need to root the foldtree at new nterm to delete upstream atoms
		kinematics::FoldTree ft_del( pose.size() );
		ft_del.clear();
		ft_del.add_edge( Edge( nterm, 1, Edge::PEPTIDE ) );
		ft_del.add_edge( Edge( nterm, pose.size(), Edge::PEPTIDE ) );
		ft_del.reorder( nterm );
		pose.fold_tree( ft_del );
		//delete n terminal
		Size const nterm_orig( pose.conformation().chain_begin( pose.chain( nterm ) ) );
		for ( Size iseq = nterm_orig; iseq <= nterm - nterm_orig; ++iseq ) {
			pose.delete_polymer_residue( pose.conformation().chain_begin( pose.chain( nterm ) ) );
		}
		//is this necessary??
		kinematics::FoldTree ft_reset( pose.size() );
		pose.fold_tree( ft_reset );
	}
}


//append pose to src_pose by jump from src_seqpos to ap_seqpos
void
append_pose_by_jump(
	Pose & src_pose,
	Size const src_seqpos,
	Pose const & ap_pose,
	Size const ap_seqpos
	//let's add functionality for jump atoms, new chain later
){
	//append first residue
	src_pose.append_residue_by_jump( ap_pose.residue( ap_seqpos ), src_seqpos, "CA", "CA" );
	Size new_jump_seqpos_downstr( src_pose.size() );
	//append Cterminal by bond
	for ( Size iseq = ap_seqpos + 1; iseq <= ap_pose.size(); ++iseq ) {
		src_pose.append_polymer_residue_after_seqpos(
			ap_pose.residue( iseq ), new_jump_seqpos_downstr + iseq - ap_seqpos - 1, false );
	}
	//prepend Nterminal by bond
	for ( Size iseq = ap_seqpos - 1; iseq >= 1; --iseq ) {
		src_pose.prepend_polymer_residue_before_seqpos(
			ap_pose.residue( iseq ), new_jump_seqpos_downstr - iseq - ap_seqpos + 1, false );
	}
}

//trim poses at cutpair, then append poses 2…N to pose 1 by jumps from pose1 cterm cutpoint
void
build_pose_by_jumps_from_domains(
	vector1< Pose > domain_poses,
	vector1< std::pair< Size, Size > > const & domain_cutpairs,
	Pose & construct_pose
){
	assert( domain_poses.size() == domain_cutpairs.size() );
	//start with first domain
	trim_pose_termini( domain_poses[ 1 ],
		domain_cutpairs[ 1 ].first, domain_cutpairs[ 1 ].second );
	construct_pose = domain_poses[ 1 ];
	Size domain_1_cterm( construct_pose.size() );
	//then append the rest
	for ( Size ipose = 2; ipose <= domain_poses.size(); ++ipose ) {
		trim_pose_termini( domain_poses[ ipose ],
			domain_cutpairs[ ipose ].first, domain_cutpairs[ ipose ].second );
		//TR << "appending pose " << ipose << " resi 1 to pose 1 resi " << domain_1_cterm << std::endl;
		append_pose_by_jump( construct_pose, domain_1_cterm, domain_poses[ ipose ], Size( 1 ) );
	}
}

//find heavyatom nbrs of all incl seqpositions
//assumes default false! only sets things to true
void
find_rsd_nbrs(
	Pose const & pose,
	vector1< bool > const & is_incl,
	Real const nbr_radius,
	vector1< bool > & is_nbr
){
	for ( Size seqpos1 = 1; seqpos1 <= pose.size(); ++seqpos1 ) {
		if ( !is_incl[ seqpos1 ] ) continue; //only look at true positions
		Residue const & rsd1( pose.residue( seqpos1 ) );
		for ( Size seqpos2 = 1; seqpos2 <= pose.size(); ++seqpos2 ) {
			if ( seqpos2 == seqpos1 ) { //is its own nbr
				is_nbr[ seqpos2 ] = true;
				continue;
			}
			Residue const & rsd2( pose.residue( seqpos2 ) );
			for ( Size iatom1 = 1; iatom1 <= rsd1.nheavyatoms(); ++iatom1 ) {
				if ( is_nbr[ seqpos2 ] ) break;
				for ( Size iatom2 = 1; iatom2 <= rsd2.nheavyatoms(); ++iatom2 ) {
					if ( rsd1.xyz( iatom1 ).distance( rsd2.xyz( iatom2 ) ) < nbr_radius ) {
						is_nbr[ seqpos2 ] = true;
						break;
					}
				}
			}
		}
	}
}

/*
given a two cutpairs, dump the insert pose (by trimming termini)
then print the blueprint
let the overlap belong to the insertee
if overlap sectructs are same, use that secstruct for linker, else use any
optional N extras linker residues on both ends
*/
void
dump_insert_pdb_and_remodel_blueprint(
	Pose const & cutpose,
	std::pair< Size, Size > const cut_cutpair,
	Pose inpose,
	std::pair< Size, Size > const in_cutpair,
	Real design_link_nbr_radius,
	Real repack_link_nbr_radius,
	Size n_add_linker
){
	Size n_remodel_linker( 2 ); //nres to allow remodeling incl overlap

	//does overlap secstruct match? check before modifying inpose
	bool is_matching_ss_link1(
		cutpose.secstruct( cut_cutpair.first ) == inpose.secstruct( in_cutpair.first ) );
	bool is_matching_ss_link2(
		cutpose.secstruct( cut_cutpair.second ) == inpose.secstruct( in_cutpair.second ) );

	//inpose cutpairs are now out of date!
	//make this basename.cut1_cut2.pdb
	std::string cut_pdbname( cutpose.pdb_info()->name() );
	Size cut_pdbnamestart( 0 );
	if ( cut_pdbname.find_last_of( "/" ) < ( cut_pdbname.size() - 1 ) ) cut_pdbnamestart = cut_pdbname.find_last_of( "/" ) + 1;
	cut_pdbname = cut_pdbname.substr( cut_pdbnamestart, cut_pdbname.size() - cut_pdbnamestart - 4 );

	//create a ScoreFunction from commandline options (default is score12)
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	//and add decompose bb_hbond energies option

	std::string in_name( "inst." + utility::to_string( in_cutpair.first ) +
		"-" + utility::to_string( in_cutpair.second ) );
	std::string cut_name( cut_pdbname + "." + utility::to_string( cut_cutpair.first ) +
		"-" + utility::to_string( cut_cutpair.second ) );
	std::string in_pdbname( in_name + ".pdb" );

	Size in_chain( inpose.chain( in_cutpair.first ) );

	//let the overlap belong to the insertee, since inserter's bb can't move in remodel
	//this now just trims to the end of chain instead of end of pose
	trim_pose_termini( inpose, in_cutpair.first + 1, in_cutpair.second - 1 );
	inpose.dump_pdb( in_pdbname );

	//need to go ahead and figure out secstruct, design requirements
	//define linkers in boolean vectors
	vector1< bool > is_link_cut( cutpose.size(), false );
	is_link_cut[ cut_cutpair.first ] = true;
	is_link_cut[ cut_cutpair.second ] = true;
	vector1< bool > is_link_in( inpose.size(), false );
	is_link_in[ 1 ] = true;
	is_link_in[ inpose.size() ] = true;
	//find nbrs of linkers
	vector1< bool > is_design_nbr_cut( cutpose.size(), false );
	vector1< bool > is_design_nbr_in( inpose.size(), false );
	vector1< bool > is_repack_nbr_cut( cutpose.size(), false );
	vector1< bool > is_repack_nbr_in( inpose.size(), false );
	find_rsd_nbrs( cutpose, is_link_cut, design_link_nbr_radius, is_design_nbr_cut );
	find_rsd_nbrs( cutpose, is_link_cut, repack_link_nbr_radius, is_repack_nbr_cut );
	find_rsd_nbrs( inpose, is_link_in, design_link_nbr_radius, is_design_nbr_in );
	find_rsd_nbrs( inpose, is_link_in, repack_link_nbr_radius, is_repack_nbr_in );

	//open the blueprint file for writing
	std::string tag( option[ chrisk::tag ] );
	std::string bpt_name( tag + "." + cut_name + "." + in_name + ".bpt" );
	char const *bpt_name_char = bpt_name.c_str();
	std::fstream bpt( bpt_name_char, std::ios::out );

	//for 1st half of cutpose
	for ( Size iseq = 1; iseq <= cut_cutpair.first; ++iseq ) {
		char const aa( cutpose.residue( iseq ).name1() );
		char ss( '.' ); //default for insertee
		if ( cut_cutpair.first - iseq < n_remodel_linker ) {
			if ( is_matching_ss_link1 ) ss = cutpose.secstruct( iseq );
			else ss = 'D';
		}
		//print required data
		bpt << iseq << " " << aa << " " << ss;
		//if is nbr, print design control
		//use remodel default for non-designed positions
		if ( is_design_nbr_cut[ iseq ] ) bpt << " " << "ALLAAxc";
		bpt << std::endl;
	}
	//for additional linker
	if ( !is_matching_ss_link1 ) {
		for ( Size iseq = 1; iseq <= n_add_linker; ++iseq ) {
			char const aa( 'x' );
			char ss( 'D' ); //default any ss for new linkers
			//if inpose, cutpose 1st linker match ss, then use that
			bpt << "0 " << aa << " " << ss << " ALLAAxc" << std::endl;
		}
	}
	//for insert
	for ( Size iseq = inpose.conformation().chain_begin( in_chain );
			iseq <= inpose.conformation().chain_end( in_chain ); ++iseq ) {
		char const aa( 'x' );
		char ss( 'I' ); //default for insertee
		//print required data
		bpt << "0 " << aa << " " << ss;
		//if is nbr, print design control
		//must set all insert's tasks manually
		if ( is_design_nbr_in[ iseq ] ) bpt << " ALLAAxc";
		else if ( is_repack_nbr_in[ iseq ] ) bpt << " NATAA";
		else bpt << " NATRO";
		bpt << std::endl;
	}
	//for additional linker
	if ( !is_matching_ss_link2 ) {
		for ( Size iseq = 1; iseq <= n_add_linker; ++iseq ) {
			char const aa( 'x' );
			char ss( 'D' ); //default any ss for new linkers
			//if inpose, cutpose 1st linker match ss, then use that
			if ( is_matching_ss_link2 ) ss = cutpose.secstruct( cut_cutpair.second );
			bpt << "0 " << aa << " " << ss << " ALLAAxc" << std::endl;
		}
	}
	//for 2nd half of cutpose
	for ( Size iseq = cut_cutpair.second; iseq <= cutpose.size(); ++iseq ) {
		char const aa( cutpose.residue( iseq ).name1() );
		char ss( '.' ); //default for insertee
		if ( iseq - cut_cutpair.second < n_remodel_linker ) {
			if ( is_matching_ss_link2 ) ss = cutpose.secstruct( iseq );
			else ss = 'D';
		}
		//print required data
		bpt << iseq << " " << aa << " " << ss;
		//if is nbr, print design control
		//use remodel default for non-designed positions
		if ( is_design_nbr_cut[ iseq ] ) bpt << " " << "ALLAAxc";
		bpt << std::endl;
	}
}


/*
// build a hybrid sequence, then just copy residues from the source poses
// the overlapping positions belong to the inserter, not the insertee
void
simple_domain_insertion(
Pose const & cutpose,
std::pair< Size, Size > const cutpose_cutpair,
Pose const & inpose,
std::pair< Size, Size > const inpose_cutpair,
Pose & newpose
){
std::string cutpose_seq( cutpose.sequence() );
std::string inpose_seq( inpose.sequence() );
//1st half of cutpose
std::string newpose_seq( cutpose_seq.substr( 1, cutpose_cutpair.first - 1 ) );
//inpose segment
newpose_seq += inpose_seq.substr( inpose_cutpair.first,
inpose_cutpair.second - inpose_cutpair.first + 1 );
//2nd half of cutpose
newpose_seq += cutpose_seq.substr( cutpose_cutpair.second + 1,
cutpose_seq.size() - cutpose_cutpair.second );
TR << "New pose sequence: " << newpose_seq << std::endl;

//now make the new pose
pose::make_pose_from_sequence( newpose, newpose_seq,
inpose.residue( 1 ).residue_type_set() );
//1st half of cutpose
copy_segment_bb( newpose, cutpose_cutpair.first,
cutpose, 1, 1 );
//inpose segment
copy_segment_bb( newpose, inpose_cutpair.second - inpose_cutpair.first + 1,
inpose, cutpose_cutpair.first, inpose_cutpair.first );
//2nd half of cutpose
copy_segment_bb( newpose, cutpose_seq.size() - cutpose_cutpair.second,
cutpose, cutpose_cutpair.first + inpose_cutpair.second - inpose_cutpair.first,
cutpose_cutpair.second + 1 );
}
*/

/*
multiplex these thigs:
-input proteins
-domain insertions

-circular permutants

we should just build sequences and copy backbone angles

what we need to do
there exist X cutpair for a given N/C termini (loops only?)
find cutpair - done
insert one domain into another w/ optional linkers
-build new seq
-copy correct phi/psi
there exist Y N/C termini
there exist Z circular permutants


*/

void
go(
	pose::Pose cutpose
){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// INPUT
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//tag is "." or "tag."
	// std::string tag( option[ chrisk::tag ] );
	// tag = tag + ".";
	// std::string const scname( tag + "pdb" );

	//init options
	Real const rmsd_tol( option[ chrisk::rmsd_tol ] ); //total rmsd deviation
	Real const rt_tol( option[ chrisk::rt_tol ] ); //total rt deviation
	Size const seqsep_tol( option[ chrisk::seqsep_tol ] ); //min dist b/t cutpair
	Size const term_tol( option[ chrisk::term_tol ] ); //termini buffer
	bool const loops_only( option[ chrisk::loops_only ] );
	Real repack_link_nbr_radius( option[ chrisk::repack_radius ] );
	Real design_link_nbr_radius( option[ chrisk::design_radius ] );
	Size n_add_linker( option[ chrisk::n_add_linker ] );

	pose::Pose inpose;
	std::string inpose_fname( option[ chrisk::insert_pdb ]() );
	core::import_pose::pose_from_file( inpose, inpose_fname , core::import_pose::PDB_file);

	//create a ScoreFunction from commandline options (default is score12)
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

	//silent-type output
	// io::silent::SilentFileData pose_silent_data;


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// /END INPUT
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup poses
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//setup secstruct
	setup_secstruct_dssp( cutpose );
	setup_secstruct_dssp( inpose );

	//just kick it to cen to speed things up for now
	// core::util::switch_to_residue_type_set( cutpose, core::chemical::CENTROID );
	// core::util::switch_to_residue_type_set( inpose, core::chemical::CENTROID );
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// /END setup poses
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// find possible insertion points
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// HEY! instead, let's scan thru all insertion point cutpoint *pairs*
	// (instead of just 1xN), then pick one pair of cutpair at random

	//define cutpair on inpose, where insertion has to satisfy
	Size nterm_in_start( 1 ); //nterm cutpoint for inserter
	Size nterm_in_end( Size( option[ chrisk::trim_tol ] ) + 1 ); //nterm cutpoint for inserter
	Size cterm_in_start( inpose.size() - Size( option[ chrisk::trim_tol ] ) ); //cterm cutpoint for inserter
	Size cterm_in_end( inpose.size() ); //cterm cutpoint for inserter
	//no looping if user-defined
	if ( option[ chrisk::nterm_in ].user() ) nterm_in_start = nterm_in_end = Size( option[ chrisk::nterm_in ] );
	if ( option[ chrisk::cterm_in ].user() ) cterm_in_start = cterm_in_end = Size( option[ chrisk::cterm_in ] );
	// assert( nterm_in < cterm_in );
	for ( Size nterm_in = nterm_in_start; nterm_in <= nterm_in_end; ++nterm_in ) {
		for ( Size cterm_in = cterm_in_start; cterm_in <= cterm_in_end; ++cterm_in ) {

			std::pair< Size, Size > inpose_cutpair( nterm_in, cterm_in );
			//find possible cutpair
			vector1< std::pair< Size, Size > > cutpose_cutpairs;
			//find matching cutpair based on RT?
			// kinematics::RT inpose_rt( get_atomatom_rt( inpose, inpose_cutpair.first, "CA", inpose_cutpair.second, "CA" ) );
			// find_backbone_rts_in_pose( cutpose, inpose_rt, rt_tol, seqsep_tol,
			//   term_tol, loops_only, cutpose_cutpairs );
			//or find matching cutpair based on bb rmsd?
			find_cutpair_alignments_in_pose( cutpose, inpose, inpose_cutpair, rmsd_tol, seqsep_tol,
				term_tol, loops_only, cutpose_cutpairs );

			//if( cutpose_cutpairs.size() < 1 ) utility_exit_with_message( "Failed to find insertion point within valid tolerance\n" );
			//this will loop through alternatives
			for ( Size icutpair = 1; icutpair <= cutpose_cutpairs.size(); ++icutpair ) {
				std::pair< Size, Size > cutpose_cutpair( cutpose_cutpairs[ icutpair ] );
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// /END find possible insertion points
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				id::AtomID_Map< id::AtomID > atm_map;
				for ( Size iatm = 1; iatm <= 4; ++iatm ) {
					atm_map[ id::AtomID( iatm, inpose_cutpair.first ) ] = id::AtomID( iatm, cutpose_cutpair.first );
					atm_map[ id::AtomID( iatm, inpose_cutpair.second ) ] = id::AtomID( iatm, cutpose_cutpair.second );
				}
				scoring::superimpose_pose( cutpose, inpose, atm_map );
				std::string cut_pdbname( cutpose.pdb_info()->name() );
				Size cut_pdbnamestart( 0 );
				if ( cut_pdbname.find_last_of( "/" ) < ( cut_pdbname.size() - 1 ) ) cut_pdbnamestart = cut_pdbname.find_last_of( "/" ) + 1;
				cut_pdbname = cut_pdbname.substr( cut_pdbnamestart, cut_pdbname.size() - cut_pdbnamestart - 4 );
				std::string cut_name( cut_pdbname + "." + utility::to_string( cutpose_cutpair.first ) +
					"-" + utility::to_string( cutpose_cutpair.second ) + ".pdb" );
				cutpose.dump_pdb( cut_name );

				dump_insert_pdb_and_remodel_blueprint( cutpose, cutpose_cutpair, inpose,
					inpose_cutpair, design_link_nbr_radius, repack_link_nbr_radius, n_add_linker );

				/*
				//put domains into a vector (in order)
				vector1< Pose > domain_poses;
				domain_poses.push_back( cutpose );
				domain_poses.push_back( inpose );
				domain_poses.push_back( cutpose );
				//and their cutpair (in same order!)
				vector1< std::pair< Size, Size > > domain_cutpair;
				domain_cutpair.push_back(
				std::pair< Size, Size >( Size( 1 ), cutpose_cutpair.first ) );
				domain_cutpair.push_back( inpose_cutpair );
				domain_cutpair.push_back(
				std::pair< Size, Size >( cutpose_cutpair.second, cutpose.size() ) );

				//assemble new construct by insertion
				//the new domain insertion construct
				Pose pose;
				build_pose_by_jumps_from_domains( domain_poses, domain_cutpair, pose );
				//jump 1 --> inserter, jump 2 --> other half of insertee
				pose.fold_tree().show( TR );
				//set inserter jump to be contiguous w/ insertee's backbone
				//can do this using alignment
				//remember, insertee is missing both its cutpair!
				//new jump goes from CA->CA?
				pose.set_jump( 1, get_rsdrsd_jump( cutpose, cutpose_cutpair.first - 1, cutpose_cutpair.first ) );
				pose.dump_pdb( "test_domain_assembly.pdb" );
				*/

			} //end cutpairs
		} //end cterm_in
	} //end nterm_in
} //end go

void
load_coords()
{
	using namespace core;

	import_pose::pose_stream::MetaPoseInputStream input =  import_pose::pose_stream::streams_from_cmd_line();
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	//iter over cutposes (insertees/hosts)
	while ( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		go( pose );
	}
}

void*
my_main( void*)
{
	load_coords();
	exit(0);

}

int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility;

	option.add( chrisk::tag, "nametag" ).def( "struct" );
	option.add( chrisk::insert_pdb, "pdb to insert into another" ).def( "struct" );
	option.add( chrisk::rt_tol, "RT dev tolerance" ).def( 1.0 );
	option.add( chrisk::rmsd_tol, "RT dev tolerance" ).def( 3.0 );
	option.add( chrisk::seqsep_tol, "max seq separation bt insert cutpair" ).def( 5 );
	option.add( chrisk::term_tol, "max seq separation bt termini and cutpair" ).def( 5 );
	option.add( chrisk::trim_tol, "max seq eliminated from inserter n,cterm" ).def( 3 );
	option.add( chrisk::nterm_in, "nterm_in residue for insertion" ).def( 1 );
	option.add( chrisk::cterm_in, "cterm_in residue for insertion" ).def( 1 );
	option.add( chrisk::loops_only, "insert into loops" ).def( true );
	option.add( chrisk::repack_radius, "radius around linker to repack" ).def( 7.0 );
	option.add( chrisk::design_radius, "radius around linker to design" ).def( 5.0 );
	option.add( chrisk::n_add_linker, "n additional linker res to add" ).def( 0 );

	devel::init(argc, argv);

	protocols::viewer::viewer_main( my_main );

}


