// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/util.cc
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Jeliazko Jeliazkov (jeliazkov@jhu.edu)

// Project Headers
#include <protocols/antibody/util.hh>

// Core Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/init_id_map.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <core/scoring/rms_util.tmpl.hh>


// Protocol Headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/docking/util.hh>
#include <protocols/interface/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Numeric Headers
#include <numeric/xyz.functions.hh>


// Basic Headers
#include <fstream>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_constants.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <cmath>

static basic::Tracer TR( "antibody.util" );



namespace protocols {
namespace antibody {

using namespace core;
using namespace protocols::antibody::clusters;
using utility::vector1;

utility::vector1<bool>
get_cdr_bool_from_tag(utility::tag::TagCOP tag, std::string const & name, bool include_cdr4 /* false */){
	core::Size total_cdrs;
	if ( include_cdr4 ) {
		total_cdrs = CDRNameEnum_proto_total;
	} else {
		total_cdrs = CDRNameEnum_total;
	}

	utility::vector1<bool> cdrs (total_cdrs, false);
	vector1<std::string> cdr_strings = utility::string_split_multi_delim(tag->getOption<std::string>(name), ":,'`~+*&|;. ");
	AntibodyEnumManager manager = AntibodyEnumManager();
	for ( core::Size i = 1; i <= cdr_strings.size(); ++i ) {
		CDRNameEnum cdr = manager.cdr_name_string_to_enum(cdr_strings[i]);
		cdrs[cdr] = true;
	}
	return cdrs;
}

void
attributes_for_get_cdr_bool_from_tag(utility::tag::AttributeList& attlist,
	std::string const& tagname, std::string const& Description) {
	using namespace utility::tag;

	std::string recognized_cdrs;
	AntibodyEnumManager manager = AntibodyEnumManager();
	for ( auto& cdr_def : manager.get_recognized_cdr_definitions() ) {
		recognized_cdrs += cdr_def + "|";
	}

	attlist + XMLSchemaAttribute(
		tagname, xs_string,
		Description + (!Description.empty() ? "\n" : "" ) +
		//"List of CDR regions (string) devided by one of the following characters: \":,'`~+*&|;. \""
		"List of CDR regions (string) devided by one of the following characters: \":,'`~+*|;. \""
		"Recognized CDRs: \"" + recognized_cdrs + "\"");
}

/// @brief Get a set of loops for a boolean vector of CDRNameEnums including any stem residues.
protocols::loops::LoopsOP
get_cdr_loops(
	AntibodyInfoCOP ab_info,
	core::pose::Pose const & pose,
	utility::vector1<bool> const & cdrs,
	core::Size stem_size /* 0 */ ) {

	debug_assert( cdrs.size() <= CDRNameEnum_proto_total );
	protocols::loops::LoopsOP cdr_loops( new protocols::loops::Loops() );
	for ( core::Size i = 1; i <= cdrs.size(); ++i ) {

		if ( ! cdrs[ i ] ) continue;
		auto cdr = static_cast<CDRNameEnum>( i );
		if ( ab_info->is_camelid() && ab_info->get_CDR_chain( cdr ) == "L" ) continue;

		cdr_loops->add_loop(ab_info->get_CDR_loop(cdr, pose, stem_size));
	}
	return cdr_loops;
}


core::pack::task::TaskFactoryOP setup_packer_task(pose::Pose & pose_in ) {
	using namespace pack::task;
	using namespace pack::task::operation;

	TR.Debug << "Setting Up Packer Task" << std::endl;

	core::pack::task::TaskFactoryOP tf( new TaskFactory );
	tf->clear();

	tf->push_back(utility::pointer::make_shared< OperateOnCertainResidues >(utility::pointer::make_shared< PreventRepackingRLT >(), utility::pointer::make_shared< ResidueLacksProperty >("PROTEIN") ));
	tf->push_back(utility::pointer::make_shared< InitializeFromCommandline >() );
	tf->push_back(utility::pointer::make_shared< IncludeCurrent >() );
	tf->push_back(utility::pointer::make_shared< RestrictToRepacking >() );
	tf->push_back(utility::pointer::make_shared< NoRepackDisulfides >() );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	pack::rotamer_set::UnboundRotamersOperationOP
		unboundrot( new pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot->initialize_from_command_line();

	operation::AppendRotamerSetOP unboundrot_operation( new operation::AppendRotamerSet( unboundrot ) );
	tf->push_back( unboundrot_operation );

	// adds scoring bonuses for the "unbound" rotamers, if any
	core::pack::dunbrack::load_unboundrot( pose_in ); ///FIXME: this is dangerous here..... JQX
	//TODO:
	//JQX: need to understand this pose_in!!!!

	TR.Debug << "Done: Setting Up Packer Task" << std::endl;

	return tf;

} // setup_packer_task

/*    void
dle_extreme_repack(
pose::Pose & pose_in,
int repack_cycles,
ObjexxFCL::FArray1D_bool & allow_repack,
bool rt_min,
bool rotamer_trials,
bool force_one_repack,
bool use_unbounds
)
{
using namespace pose;

// Exit if not fullatom
if( !pose_in.fullatom() ) {
std::cout << "Repack called in centroid mode" << std::endl;
std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
std::cout << "------NOT REPACKING-----------" << std::endl;
return;
}
// Saving parameters
bool initial_rot_trial = score_get_try_rotamers();
bool initial_min_rot = get_minimize_rot_flag();
score_set_minimize_rot( rt_min );
score_set_try_rotamers( rotamer_trials );
// initial allowed chi movement
FArray1D_bool old_chi_move( pose_in.size(), false );
for( int i = 1; i <= pose_in.size(); i++ ) {
// storing old
old_chi_move(i) = pose_in.get_allow_chi_move(i);
// setting new
pose_in.set_allow_chi_move( i, allow_repack(i) || old_chi_move(i) );
}
Score_weight_map weight_map( score12 );
Monte_carlo mc( pose_in, weight_map, 2.0 );
// repack idealized native
Pose start_native_pose;
start_native_pose = pose_in;
for(int i=1; i <= repack_cycles; i++) {
pose_in = start_native_pose;
if( use_unbounds )
dle_pack_with_unbound( pose_in, allow_repack, true  ); // include_current = true
else
pose_in.repack( allow_repack, true ); //include_current =  true
pose_in.score( weight_map );
score_set_minimize_rot( false );
score_set_try_rotamers( false );
if( force_one_repack && (i == 1) ) mc.reset( pose_in );
mc.boltzmann( pose_in );
score_set_minimize_rot( rt_min );
score_set_try_rotamers( rotamer_trials );
}
pose_in = mc.low_pose();
pose_in.score( weight_map );

// Restoring Globals
score_set_minimize_rot( initial_rot_trial );
score_set_try_rotamers( initial_min_rot );
pose_in.set_allow_chi_move( old_chi_move );

return;
}
*/


bool cutpoints_separation( core::pose::Pose & pose, AntibodyInfoOP & antibody_info ) {

	bool closed_cutpoints = true;

	for ( auto const & elem : *( antibody_info->get_AllCDRs_in_loopsop() ) ) {
		core::Size cutpoint   = elem.cut();
		//Real separation = 10.00; // an unlikely high number
		Real separation = cutpoint_separation( pose, cutpoint );

		if ( separation > 1.9 ) {
			closed_cutpoints = false;
			break;
		}
	}
	return( closed_cutpoints );
} // cutpoints_separation

Real cutpoint_separation(pose::Pose & pose_in, core::Size cutpoint) {

	core::Size const N ( 1 ); // N atom
	core::Size const C ( 3 ); // C atom

	// Coordinates of the C atom of cutpoint res and N atom of res cutpoint+1
	numeric::xyzVector_float peptide_C(pose_in.residue( cutpoint ).xyz( C )),
		peptide_N( pose_in.residue( cutpoint + 1 ).xyz( N ) );
	//   Real cutpoint_separation=distance(peptide_C, peptide_N);
	Real cutpoint_separation=peptide_C.distance(peptide_N);

	return( cutpoint_separation );
} // cutpoint_separation


Real global_loop_rmsd (const pose::Pose & pose_in, const pose::Pose & native_pose,loops::LoopsOP current_loop ) {
	if ( pose_in.size() != native_pose.size() ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "The pose sequence length does not match that of native_pose");
	}

	using namespace scoring;

	core::Size loop_start = (*current_loop)[1].start();
	core::Size loop_end = (*current_loop)[1].stop();

	using ObjexxFCL::FArray1D_bool;
	FArray1D_bool superpos_partner ( pose_in.size(), false );

	for ( core::Size i = loop_start; i <= loop_end; ++i ) superpos_partner(i) = true;

	using namespace core::scoring;
	Real rmsG = rmsd_no_super_subset( native_pose, pose_in, superpos_partner, is_protein_backbone_including_O );
	return ( rmsG );
}

void align_to_native( core::pose::Pose & pose,
	core::pose::Pose const & native_pose,
	AntibodyInfoOP const ab_info,
	AntibodyInfoOP const native_ab_info,
	std::string const & reqeust_chain) {


	std::string pose_seq        = pose.sequence();
	std::string native_pose_seq = native_pose.sequence();
	if ( pose_seq != native_pose_seq   ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, " the pose sequence does not match native_pose sequence ");
	}


	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, core::id::AtomID::BOGUS_ATOM_ID() );

	// loop over the L and H chains
	for ( core::Size i_chain=1; i_chain<=ab_info->get_AntibodyFrameworkInfo().size(); i_chain++ ) {
		vector1<FrameWork>        chain_frmwk =        ab_info->get_AntibodyFrameworkInfo()[i_chain];
		vector1<FrameWork> native_chain_frmwk = native_ab_info->get_AntibodyFrameworkInfo()[i_chain];

		if (  (chain_frmwk[1].chain_name == reqeust_chain ) || (reqeust_chain =="LH")  ) {
			// loop over the segments on the framework of one chain
			for ( core::Size j_seg=1; j_seg<=chain_frmwk.size(); j_seg++ ) { // for loop of the framework segments
				core::Size count=0;

				// loop over the residues on one segment on one framework of one chain
				for ( core::Size k_res=chain_frmwk[j_seg].start; k_res<= chain_frmwk[j_seg].stop; k_res++ ) {
					count++;
					core::Size res_counter = k_res;
					core::Size nat_counter = native_chain_frmwk[j_seg].start+count-1;
					//TR << "Matching Residue " << res_counter << " (" << pose.residue(res_counter).name() << ")" << " with  " << nat_counter << " (" << native_pose.residue(nat_counter).name() << ")" << std::endl;

					// loop over the backbone atoms including Oxygen
					for ( core::Size latm=1; latm <= 4; latm++ ) {
						core::id::AtomID const id1( latm, res_counter );
						core::id::AtomID const id2( latm, nat_counter );
						atom_map[ id1 ] = id2;
					}
				}

			}// end of the for loop

		}//end of the if statement
	}

	core::scoring::superimpose_pose( pose, native_pose, atom_map );

} // align_to_native()

vector1<bool>
select_epitope_residues(AntibodyInfoCOP ab_info, core::pose::Pose const & pose, core::Size const interface_distance) {
	vector1<bool> epitope(pose.size(), false);
	if ( ! ab_info->antigen_present() ) return epitope;

	std::string interface = ab_info->get_antibody_chain_string()+"_"+ab_info->get_antigen_chain_string();

	//I really should have just used a TF and operation.  Oh well.  Now we have an interface namespace...

	vector1<bool> interface_residues = protocols::interface::select_interface_residues(pose, interface, interface_distance);

	//Turn off L or H residues at the interface
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		char chain = core::pose::get_chain_from_chain_id(pose.chain(i), pose);
		if ( chain == 'L' || chain == 'H' ) {
			interface_residues[i] = false;
		}
	}
	return interface_residues;

}

void
check_fix_aho_cdr_numbering(AntibodyInfoCOP ab_info, CDRNameEnum cdr, core::pose::Pose & pose){
	if ( ab_info->get_current_AntibodyNumberingScheme() != "AHO_Scheme" ) {
		return;
	}

	core::Size max_length = ab_info->get_CDR_end_PDB_num(cdr) - ab_info->get_CDR_start_PDB_num(cdr);
	core::Size cdr_length = ab_info->get_CDR_length(cdr, pose);
	if ( cdr_length > max_length ) {
		TR << "adding insertion codes for long cdr loops." << std::endl;

		std::string const & alphabet( utility::ALPHANUMERICS );

		core::Size start_residues  = ceil(max_length/2.0);
		core::Size end_residues = floor(max_length/2.0);

		//TR <<cdr_length << " Ceil "<< start_residues << " Floor " << end_residues << std::endl;

		core::Size i_code_start = ab_info->get_CDR_start(cdr, pose) + start_residues +1;
		core::Size i_code_pdb = pose.pdb_info()->number(i_code_start -1);
		core::Size i_code_end = ab_info->get_CDR_end(cdr, pose) - end_residues;
		core::Size i = 0;

		for ( core::Size resnum = i_code_start; resnum <= i_code_end; ++resnum ) {
			if ( i > alphabet.size() ) {
				TR <<"Cannot have a CDR with insertion codes more than "<< alphabet.size() <<" skipping rest of fix.." << std::endl;
				return;
			}
			//TR << resnum <<" "<< alphabet[ i ] << std::endl;
			//TR << pose.pdb_info()->number(resnum) << " "<< i_code_pdb << std::endl;
			pose.pdb_info()->icode(resnum, alphabet[ i ]);
			pose.pdb_info()->number(resnum, i_code_pdb);
			i+=1;
		}

	} else {
		return;
	}

}

void
check_fix_aho_cdr_numbering(AntibodyInfoCOP ab_info, core::pose::Pose & pose){

	for ( auto const & cdr : ab_info->get_all_cdrs_present() ) {
		check_fix_aho_cdr_numbering(ab_info, cdr, pose);
	}

}

//JQX:
// Description (Jason contributed)
// assuming a loop 10,11,12,13,14,15, cut_point is 13/14
//                  or loop = Loop (10, 15, 13)
// 1. (Rosetta Default)
//           set_single_loop_fold_tree()
//           jump from 8->17 (-2 amd +2)
// 2. (Aroop's)
//            simple_one_loop_fold_tree()
//            jump from 9->16 (-1 and +1)
// 3. (Aroop's)
//            simple_fold_tree()
//            jump from 10->15, but you need manuly input
// 4. (Aroop's)
//            setup_simple_fold_tree()
//            jump from 10->15, but you need manuly input

void simple_one_loop_fold_tree(
	pose::Pose & pose_in,
	loops::Loop const & loop ) {
	using namespace kinematics;

	TR.Debug <<  "Setting up simple one loop fold tree" << std::endl;

	//setup fold tree for this loop
	FoldTree f;
	f.clear();
	core::Size nres = pose_in.size();
	core::Size jumppoint1 = loop.start() - 1;
	core::Size jumppoint2 = loop.stop() + 1;

	if ( jumppoint1 < 1 )   jumppoint1 = 1;
	if ( jumppoint2 > nres ) jumppoint2 = nres;

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, loop.cut(), Edge::PEPTIDE );
	f.add_edge( loop.cut() + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, 1 );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR.Debug <<  "Finished setting up simple one loop fold tree" << std::endl;

	return;
} // simple_one_loop_fold_tree

void simple_fold_tree(
	pose::Pose & pose_in,
	core::Size jumppoint1,
	core::Size cutpoint,
	core::Size jumppoint2 ) {
	using namespace kinematics;

	TR.Debug <<  "Setting up simple fold tree" << std::endl;

	//setup fold tree for this loop
	FoldTree f;
	f.clear();
	core::Size nres = pose_in.size();

	if ( jumppoint1 < 1 )   jumppoint1 = 1;
	if ( jumppoint2 > nres ) jumppoint2 = nres;

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
	f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, 1 );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR.Debug <<  "Finished setting up simple fold tree" << std::endl;

	return;
} // simple_fold_tree

std::string
setup_LH_A_foldtree(AntibodyInfoCOP ab_info, core::pose::Pose & pose){
	vector1< int > movable_jumps(1, 1);
	vector1< char > antigen_chains = ab_info->get_antigen_chains();
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	std::string dock_chains = "LH_"+antigen;
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	return dock_chains;

}

std::string
setup_A_LH_foldtree(AntibodyInfoCOP ab_info, core::pose::Pose & pose){
	vector1< int > movable_jumps(1, 1);
	vector1< char > antigen_chains = ab_info->get_antigen_chains();
	std::string antigen(antigen_chains.begin(), antigen_chains.end());
	std::string dock_chains = antigen+"_LH";
	protocols::docking::setup_foldtree(pose, dock_chains, movable_jumps);
	return dock_chains;
}


///////////////////////////////////////////////////////////////////////////
///
/// @brief tests if a loop has H3 like base charachteristics
///
/// @details Uses the Shirai rules to find out if the dihedral angle
///           formed by CA atoms of residues n-2,n-1,n and n+1 conform to a
///           kinked/extended structure in accordance with the sequence. If
///           there is a match, a true value is returned
///
/// @param[in] pose: full actual protein
///            loop_begin: seq numbered loop begin corresponding to pose
///            size: size of loop to compute loop_end
///
/// @global_read reads -command line flag -base stored in dle_ns
///              to determine to do the complete H3 filter check or just do
///              a prediction of the H3 base type based on the
///              aforementioned dihedral angle
///
/// @references Structural classification of CDR-H3 in antibodies
///             Hiroki Shirai, Akinori Kidera, Haruki Nakamura
///             FEBS Letters 399 (1996) 1-8
///
/// @author Aroop 02/04/2010
///
///////////////////////////////////////////////////////////////////////////

//TODO:
//JQX:
//work with Daisuke to put the L89 creteria into the code

bool CDR_H3_filter_legacy_code_with_old_rule(const pose::Pose & pose_in, loops::Loop & input_loop, bool is_camelid) {


	TR.Debug <<  "Checking Kink/Extended CDR H3 Base Angle" << std::endl;


	std::string const light_chain = "L";

	if ( is_camelid ) {
		return( true );
	}

	// Values read from plot in reference paper. Fig 1 on Page 3
	// Values adjusted to match data from antibody training set
	Real const kink_lower_bound = -10.00; // Shirai: 0
	Real const kink_upper_bound = 70.00; // Shirai: 70
	Real const extended_lower_bound = 125.00; // Shirai: ~180
	Real const extended_upper_bound = 185.00; // Shirai: ~180

	// Hydrogen Bond maximum value is 3.9 Angstroms - not used
	// Real const h_bond(3.9);
	// Salt Bridge maximum value is 2.0 Angstroms - not used
	// Real const s_bridge(4.0);

	// chop out the loop:
	//JQX: 2 residues before h3, one residue after h3. Matched Rosetta2!
	core::Size start(input_loop.start()-2);
	core::Size stop(input_loop.stop()+1);


	//bool is_kinked( false );  // unused ~Labonte
	//bool is_extended( false );  // unused ~Labonte
	bool is_H3( false );

	// extract 3 letter residue codes for the chopped loop
	std::vector <std::string> aa_name; // loop residue 3 letter codes
	//JQX: pay attention here!! It is vector, not vector1! too painful to compare to R2 code
	//     just make vector, so it can match R2 code easily
	for ( core::Size ii=start; ii<=stop; ii++ ) {
		aa_name.push_back(pose_in.residue(ii).name3() );
		//            TR<<pose_in.residue(ii).name1()<<std::endl;
	}

	core::Size const CA(2);   // CA atom position in full_coord array
	// base dihedral angle to determine kinked/extended conformation

	Real base_dihedral( numeric::dihedral_degrees(
		pose_in.residue( stop ).xyz( CA ),
		pose_in.residue( stop - 1).xyz( CA ),
		pose_in.residue( stop - 2).xyz( CA ),
		pose_in.residue( stop - 3).xyz( CA ) ) );

	//pose_in.dump_pdb("check_cter_dihedral.pdb");


	TR << "Base Dihedral: " << base_dihedral << std::endl;


	// JQX: the code in the below if statement was in Rosetta 2, but Aroop did not port it into R3
	//      Maybe he has already tested that, the performance was better if using extra
	//      sequence creteria. But for now, I still port the code in here. Maybe for some reason
	//      one still decides not to use the sequence rules.

	bool H3_base_only=false;

	if ( H3_base_only ) {
		std::string base;
		if ( (base_dihedral > kink_lower_bound) && (base_dihedral < kink_upper_bound) ) {
			base = "KINK";
			TR << "              " << base << std::endl;
		} else if ( (base_dihedral > extended_lower_bound) && (base_dihedral < extended_upper_bound) ) {
			base = "EXTENDED";
			TR << "              " << base << std::endl;
		} else {
			base = "NEUTRAL";
			TR << "              " << base << std::endl;
		}

		return( true );
	}


	// setting up pseudo-periodic range used in extended base computation
	if ( base_dihedral < kink_lower_bound ) {
		base_dihedral = base_dihedral + 360.00;
	}

	// Rule 1a for standard kink
	if ( (aa_name[aa_name.size()-3] != "ASP") && (aa_name[aa_name.size()-1] == "TRP") ) { //aa_name.size()-3 = n-1
		//aa_name.size()-1 = n+1
		if ( (base_dihedral > kink_lower_bound) && (base_dihedral < kink_upper_bound) ) {
			// std::cout << "KINK Found" << std::endl; // aroop_temp remove
			//is_kinked = true;  // unused ~Labonte
			is_H3 = true;
		}
	}

	// Rule 1b for standard extended form
	if (  ( aa_name[ aa_name.size() - 3 ] == "ASP" ) &&
			( ( aa_name[1] != "LYS" ) && ( aa_name[1] != "ARG" ) ) &&
			( is_H3 != true )    ) {   //aa_name[1] = 0 position

		if ( ( base_dihedral>extended_lower_bound) && (base_dihedral<extended_upper_bound) ) {
			// std::cout << "EXTENDED Found" << std::endl; // aroop_temp remove
			//is_extended = true;  // unused ~Labonte
			is_H3 = true;
		}

		if ( !is_H3 ) {
			// Rule 1b extension for special kinked form
			bool is_basic(false); // Special basic residue exception flag
			for ( core::Size ii = 2; ii <= core::Size(aa_name.size() - 5); ii++ ) {      //aa_name.size() - 5 = n-3
				if ( aa_name[ii] == "ARG" || aa_name[ii] == "LYS" ) {      //aa_name[2] =  0 position
					is_basic = true;
					break;
				}
			}

			if ( !is_basic ) {
				core::Size rosetta_number_of_L49 = pose_in.pdb_info()->pdb2pose(light_chain, 49 );
				std::string let3_code_L49 = pose_in.residue( rosetta_number_of_L49 ).name3();
				if ( let3_code_L49 == "ARG" || let3_code_L49 == "LYS" ) {
					is_basic = true;
				}
			}
			if ( is_basic && ( base_dihedral > kink_lower_bound ) &&
					( base_dihedral < kink_upper_bound ) ) {
				// aroop_temp remove
				// std::cout << "KINK (special 1b) Found" << std::endl;
				//is_kinked = true;  // unused ~Labonte
				is_H3 = true;
			}
		}
	}

	// Rule 1c for kinked form with salt bridge
	if ( ( aa_name[ aa_name.size() - 3 ] == "ASP") &&
			( ( aa_name[1] == "LYS") || ( aa_name[1] == "ARG" ) ) &&
			( (aa_name[0] != "LYS" ) && ( aa_name[0] != "ARG" ) ) &&
			( is_H3 != true) ) {
		if ( (base_dihedral > kink_lower_bound ) &&
				(base_dihedral < kink_upper_bound ) ) {
			// aroop_temp remove
			// std::cout << "KINK (w sb) Found" << std::endl;
			//is_kinked = true;  // unused ~Labonte
			is_H3 = true;
		}
		if ( !is_H3 ) {
			bool is_basic(false); // Special basic residue exception flag
			core::Size rosetta_number_of_L46 = pose_in.pdb_info()->pdb2pose( light_chain, 46 );
			std::string let3_code_L46 = pose_in.residue( rosetta_number_of_L46 ).name3();
			if ( let3_code_L46 == "ARG" || let3_code_L46 == "LYS" ) is_basic = true;
			if ( is_basic && (base_dihedral > extended_lower_bound ) &&
					( base_dihedral < extended_upper_bound ) ) {
				// aroop_temp remove
				// std::cout << "EXTENDED (special 1c) Found" << std::endl;
				//is_extended = true;  // unused ~Labonte
				is_H3 = true;
			}
		}
	}

	// Rule 1d for extened form with salt bridge
	if (  (  aa_name[ aa_name.size() - 3 ] == "ASP") &&
			( ( aa_name[1] == "LYS") || ( aa_name[1] == "ARG" ) ) &&
			( ( aa_name[0] == "LYS") || ( aa_name[0] == "ARG") ) &&
			( is_H3 != true )       ) {
		if ( ( base_dihedral > extended_lower_bound ) &&
				( base_dihedral < extended_upper_bound ) ) {
			// aroop_temp remove
			// std::cout << "EXTENDED (w sb) Found" << std::endl;
			//is_extended = true;  // unused ~Labonte
			is_H3 = true;
		}
	}

	TR.Debug <<  "Finished Checking Kink/Extended CDR H3 Base Angle: " << is_H3 << std::endl;

	return is_H3;
} // CDR_H3_filter

bool CDR_H3_cter_filter(const pose::Pose & pose_in, AntibodyInfoOP ab_info) {

	TR <<  "Checking Kink/Extended CDR H3 Base Angle" << std::endl;

	if ( ab_info->is_camelid() ) {
		return( true );
	}

	// Values read from plot in reference paper. Fig 1 on Page 3
	// Values adjusted to match data from antibody training set
	Real const kink_lower_bound = -10.00; // Shirai: 0
	Real const kink_upper_bound = 70.00; // Shirai: 70
	Real const extended_lower_bound = 125.00; // Shirai: ~180
	Real const extended_upper_bound = 185.00; // Shirai: ~180

	// Hydrogen Bond maximum value is 3.9 Angstroms - not used
	// Real const h_bond(3.9);
	// Salt Bridge maximum value is 2.0 Angstroms - not used
	// Real const s_bridge(4.0);

	// chop out the loop:
	//JQX: 2 residues before h3, one residue after h3. Matched Rosetta2!
	core::Size start(  ab_info->get_CDR_loop(h3, pose_in, Aroop).start()  -  2 );
	core::Size stop(  ab_info->get_CDR_loop(h3, pose_in, Aroop).stop()  +  1  );


	bool matched_kinked( false );
	bool matched_extended( false );


	// extract 3 letter residue codes for the chopped loop
	std::vector <std::string> aa_name; // loop residue 3 letter codes
	//JQX: pay attention here!! It is vector, not vector1! too painful to compare to R2 code
	//     just make vector, so it can match R2 code easily
	for ( core::Size ii=start; ii<=stop; ii++ ) {
		aa_name.push_back(pose_in.residue(ii).name3() );
		//            TR<<pose_in.residue(ii).name1()<<std::endl;
	}

	core::Size const CA(2);   // CA atom position in full_coord array
	// base dihedral angle to determine kinked/extended conformation

	Real base_dihedral( numeric::dihedral_degrees(
		pose_in.residue( stop ).xyz( CA ),
		pose_in.residue( stop - 1).xyz( CA ),
		pose_in.residue( stop - 2).xyz( CA ),
		pose_in.residue( stop - 3).xyz( CA ) ) );

	//pose_in.dump_pdb("check_cter_dihedral.pdb");

	TR << "Base Dihedral: " << base_dihedral << std::endl;

	// setting up pseudo-periodic range used in extended base computation
	if ( base_dihedral < kink_lower_bound ) {
		base_dihedral = base_dihedral + 360.00;
	}


	if ( (base_dihedral > kink_lower_bound ) && (base_dihedral < kink_upper_bound ) ) {
		if ( ab_info->get_H3_kink_type()==Kinked ) {
			matched_kinked = true;
		}
	}
	if ( (base_dihedral > extended_lower_bound ) && (base_dihedral < extended_upper_bound ) ) {
		if ( ab_info->get_H3_kink_type()==Extended ) {
			matched_extended = true;
		}
	}


	bool passed;
	if ( matched_kinked || matched_extended ) {
		passed=true;
	} else {
		passed=false;
	}

	TR <<  "Finished Checking Kink/Extended CDR H3 Base Angle: " << passed << std::endl;

	return passed;
}

bool
is_H3_rama_kinked(std::string const & rama){
	//TR <<"Rama: "<< rama << std::endl;
	std::string last_two = rama.substr(rama.length() - 2);

	//TR <<"Last two "<< last_two << std::endl;

	if ( last_two == "AB" || last_two == "DB" ) {
		//TR << "kinked. " << std::endl;
		return true;
	} else {
		return false;
	}
}

void
kink_constrain_antibody_H3( core::pose::Pose & pose, AntibodyInfoOP const antibody_info ) {
	// Get relevant residue numbers
	core::Size kink_begin = antibody_info->kink_begin( pose );
	kink_constrain_antibody_H3( pose, kink_begin );
}

void
kink_constrain_antibody_H3( core::pose::Pose & pose, core::Size kink_begin ) {

	using namespace core::scoring::constraints;

	TR << " Automatically setting kink constraint! " << std::endl;

	// Constraints operate on AtomIDs:
	core::Size CA( 2 ); // CA is atom 2; AtomIDs can only use numbers
	id::AtomID const atom_100x( CA, kink_begin );
	id::AtomID const atom_101( CA, kink_begin + 1 );
	id::AtomID const atom_102( CA, kink_begin + 2 );
	id::AtomID const atom_103( CA, kink_begin + 3 );

	// Yeah, yeah this is turrible. I'm trying to get stuff done over here, ok?
	// Eventually this will read from a database file
	// Generate functions
	Real alpha_x0 = 0.678; // radians
	Real alpha_sd = 0.41; // radians
	Real alpha_tol = 0.205; // radians
	scoring::func::FlatHarmonicFuncOP alpha_func( new scoring::func::FlatHarmonicFunc( alpha_x0, alpha_sd, alpha_tol ) );

	Real tau_x0 = 1.761; // radians
	Real tau_sd = 0.194; // radians
	Real tau_tol = 0.0972; // radians
	scoring::func::FlatHarmonicFuncOP tau_func( new scoring::func::FlatHarmonicFunc( tau_x0, tau_sd, tau_tol ) );

	// Instantiate constraints
	ConstraintOP tau_cst( new AngleConstraint( atom_100x, atom_101, atom_102, tau_func ) );
	ConstraintOP alpha_cst( new DihedralConstraint( atom_100x, atom_101, atom_102, atom_103, alpha_func ) );

	// Cache to pose
	pose.add_constraint( tau_cst );
	pose.add_constraint( alpha_cst );
}




void
qq_constrain_antibody( core::pose::Pose & pose, core::Size VH_qq_resi, core::Size VL_qq_resi ) {
	// borrowed from kink_constrain_antibody_H3...

	using namespace core::scoring::constraints;

	TR << " Automatically setting qq constraint. " << std::endl;

	// Constraints operate on AtomIDs:
	Size OE1( 8 ); // OE1 is atom 8; AtomIDs can only use numbers
	Size NE2( 9 ); // NE2 is atom 9; AtomIDs can only use numbers

	id::AtomID const VH_OE1( OE1, VH_qq_resi );
	id::AtomID const VL_OE1( OE1, VL_qq_resi );
	id::AtomID const VH_NE2( NE2, VH_qq_resi );
	id::AtomID const VL_NE2( NE2, VL_qq_resi );

	//    Separate parameters for two H-bonds
	//    Real qq_x0_1 = 2.924; // H bond distance (based on distribution in antibody database 01/09/2020)
	//    Real qq_sd_1 = 0.244; // std (based on distribution in antibody database 01/09/2020)
	//    Real qq_tol_1 = ( 0.5 * qq_sd_1 );
	//    Real qq_sd_1_adjust = ( qq_sd_1/ sqrt(2.0) );
	//    scoring::func::FlatHarmonicFuncOP qq_func1( new scoring::func::FlatHarmonicFunc( qq_x0_1, qq_sd_1_adjust, qq_tol_1) );
	//
	//    Real qq_x0_2 = 2.913; // H bond distance (based on distribution in antibody database 01/09/2020)
	//    Real qq_sd_2 = 0.238; // std (based on distribution in antibody database 01/09/2020)
	//    Real qq_tol_2 = ( 0.5 * qq_sd_2);
	//    Real qq_sd_2_adjust = ( qq_sd_2/ sqrt(2.0) );
	//    scoring::func::FlatHarmonicFuncOP qq_func2( new scoring::func::FlatHarmonicFunc( qq_x0_2, qq_sd_2_adjust, qq_tol_2) );
	//


	//  Joint parameters for two H-bonds
	Real qq_x0 = 2.91; // H bond distance (based on distribution in antibody database 01/09/2020)
	Real qq_sd = 0.23; // std (based on distribution in antibody database 01/09/2020)
	Real qq_tol = ( 0.5 * qq_sd );
	Real qq_sd_adjust = ( qq_sd/ sqrt(2.0) );
	scoring::func::FlatHarmonicFuncOP qq_func( new scoring::func::FlatHarmonicFunc( qq_x0, qq_sd_adjust, qq_tol) );

	// Instantiate constraints
	ConstraintOP qq_cst_1( new AtomPairConstraint( VH_OE1, VL_NE2, qq_func ) );
	ConstraintOP qq_cst_2( new AtomPairConstraint( VL_OE1, VH_NE2, qq_func ) );

	// Cache to pose
	pose.add_constraint( qq_cst_1 );
	pose.add_constraint( qq_cst_2 );
}





PDBLandmarkCOP
get_matching_landmark(
	AntibodyNumbering const & numbering,
	PDBLandmark const & query_landmark,
	AntibodyNumberingSchemeEnum const from_scheme,
	AntibodyNumberingSchemeEnum const to_scheme){

	utility::vector1<PDBLandmarkOP> current_landmarks = numbering.numbering_scheme_transform.at(from_scheme);
	utility::vector1<PDBLandmarkOP> new_landmarks = numbering.numbering_scheme_transform.at(to_scheme);

	for ( core::Size x = 1; x <= current_landmarks.size(); ++x ) {
		PDBLandmarkCOP landmark = current_landmarks[x];

		//TR<< "Attempt match:<"<< query_landmark.get_string()<<">"<< std::endl;

		if ( *landmark == query_landmark ) {

			PDBLandmarkOP new_landmark = new_landmarks[x];

			//TR << "Matched "<< "<" << query_landmark.get_string() << ">" << std::endl;
			//TR << "To " << "<" <<new_landmark->get_string() << ">" << std::endl;

			return new_landmark;

		} else {
			continue;
		}
	}

	PDBLandmarkOP empty_landmark( new PDBLandmark("X", 0, ' ') );
	return empty_landmark;

}

vector1< Real >
get_conserved_residue_list( char chain ) {

	vector1< core::Size > conserved_frh_residues{10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,36,37,38,39,40,41,42,43,44,45,46,47,48,49,66,69,70,71,72,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,103,104,105};

	vector1< core::Size > conserved_frl_residues{10,11,12,13,14,15,16,17,18,19,20,21,22,23,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,57,58,59,60,61,62,63,64,65,66,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,98,99,100};

	vector1< core::Size > conserved_residues;
	if ( chain == 'H' ) {
		conserved_residues = conserved_frh_residues;
	}

	if ( chain == 'L' ) {
		conserved_residues = conserved_frl_residues;
	}

	return conserved_residues;
}

} // namespace antibody
} // namespace protocols

