// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/domain_insertion/FusePosesNtoCMover.cc
/// @brief  cc file for FusePosesNtoCMover
/// @author Florian Richter, flosopher@gmail.com, february 2013


// Unit headers
#include <devel/domain_insertion/FusePosesNtoCMover.hh>
#include <devel/domain_insertion/FusePosesNtoCMoverCreator.hh>

// package headers

// Project headers
#include <basic/Tracer.hh>

#include <core/chemical/VariantType.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/sequence/ABEGOManager.hh>

#include <protocols/relax/FastRelax.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>


// C++ headers
#include <set>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace devel {
namespace domain_insertion {

static THREAD_LOCAL basic::Tracer tr( "devel.domain_insertion.FusePosesNtoCMover" );

FusePosesNtoCMover::FusePosesNtoCMover()
:
	fuse_pose_(/* NULL */), sfxn_(nullptr), sfxn_nocb_(nullptr),
	//mostly arbitrary numbers, but can be modified through RS tag
	superpose_window_(3), relax_window_(0), max_number_allowed_clashes_(10),
	hardsphere_clash_limit_(2.5), max_allowed_rms_(5.0),
	nterm_span_tag_(""),
	relax_mover_(/* NULL */), rs_specified_relax_mover_(false), debugmode_(false)
{
	chains_to_use_.push_back('A');
	chain_mappings_.clear();
	fuse_positions_.clear();
	fusion_regions_.clear();
	jumps_to_move_.clear();
}

FusePosesNtoCMover::FusePosesNtoCMover( FusePosesNtoCMover const & ) = default;


FusePosesNtoCMover::~FusePosesNtoCMover()= default;

protocols::moves::MoverOP
FusePosesNtoCMover::clone() const{
	return protocols::moves::MoverOP( new FusePosesNtoCMover( *this ) );
}


// XRW TEMP std::string
// XRW TEMP FusePosesNtoCMover::get_name() const
// XRW TEMP {
// XRW TEMP  return "FusePosesNtoCMover";
// XRW TEMP }


/// @details
/// strategy:
/// 1. loop over all possible combinations of candidate seqpos in the two poses
/// 2. for each pair
///    do the superposition, calculate rms and count clashes
/// 3. if both values are below the cutoffs, call the relax mover
/// 4. for every relaxed pose, report energy and rms over stitched region, output
/// relatively simple, aye?
void
FusePosesNtoCMover::apply( core::pose::Pose & pose )
{

	//0. some preparations
	this->setup_candidate_pdbpos( candidate_pdbpos_pose_, nterm_span_tag_ );
	this->setup_chain_mappings( pose );

	//1.
	// remember: fuse pose will be new c-terminus
	for ( Size i(1); i <= candidate_pdbpos_pose_.size(); ++i ) {

		for ( Size j(1); j <= candidate_pdbpos_fuse_pose_.size(); ++j ) {

			Size pose_pdbpos( candidate_pdbpos_pose_[i] ), fusepose_pdbpos( candidate_pdbpos_fuse_pose_[j] );

			std::string combo_identifier(utility::to_string(pose_pdbpos)+"_"+utility::to_string(fusepose_pdbpos) );

			core::id::AtomID_Map< core::id::AtomID > superpose_map( this->generate_superposition_map( pose, *fuse_pose_, pose_pdbpos, fusepose_pdbpos ) );

			//2.
			Real this_pair_rmsd( core::scoring::rms_at_all_corresponding_atoms( *fuse_pose_, pose, convert_AtomID_Map_to_std_map( superpose_map) ) );

			core::scoring::superimpose_pose( *fuse_pose_, pose, superpose_map );

			core::pose::PoseOP trunc_pose( this->truncate_pose_at_fusion_site( pose, pose_pdbpos, false ) );
			core::pose::PoseOP trunc_fusepose( this->truncate_pose_at_fusion_site( *fuse_pose_, fusepose_pdbpos, true ) );

			Size num_hardsphere_clashes( this->count_hardsphere_clashes( *trunc_pose, *trunc_fusepose ) );

			if ( debugmode_ ) {
				tr << "RMS for superposition on residues " << pose_pdbpos << " and " << fusepose_pdbpos << " is " << this_pair_rmsd << " and there are " << num_hardsphere_clashes << " hardsphere clashes." << std::endl;
				fuse_pose_->dump_pdb("fusepose_super_"+combo_identifier+"_"+utility::to_string( this_pair_rmsd)+"rms_"+utility::to_string(num_hardsphere_clashes)+"clash.pdb");

				trunc_pose->dump_pdb("truncated_pose.pdb");
				trunc_fusepose->dump_pdb("truncated_fusepose.pdb");
			}

			//if pose is too bad at this point, skip
			if ( (this_pair_rmsd > max_allowed_rms_ ) || ( num_hardsphere_clashes > max_number_allowed_clashes_ ) ) continue;

			//3.
			core::pose::PoseOP fused_pose_raw( new core::pose::Pose(*trunc_pose) );

			this->fuse_poses( *fused_pose_raw, *trunc_fusepose );
			core::pose::PoseOP fused_pose_rlx( new core::pose::Pose(*fused_pose_raw) );

			tr << "Scorefunction has chainbreak weight of " << sfxn_->get_weight( core::scoring::chainbreak ) << std::endl;
			(*sfxn_)( *fused_pose_raw );

			tr << "the raw fusion for " << combo_identifier << " has a chainbreak score of " << fused_pose_raw->energies().total_energies()[ core::scoring::chainbreak ] << std::endl;

			//now that we have the fused pose, we should relax it
			//this can be done through an rs specified mover, but if
			if ( rs_specified_relax_mover_ ) relax_mover_->apply( *fused_pose_rlx );

			//no mover was specified, so we'll custom make one
			else {
				//a special fold tree is probably necessary
				core::kinematics::FoldTree old_ft( fused_pose_rlx->fold_tree() );
				this->setup_relax_fold_tree( *fused_pose_rlx );
				relax_mover_ = generate_default_relax_mover( *fused_pose_rlx );
				relax_mover_->apply( *fused_pose_rlx );
				fused_pose_rlx->fold_tree( old_ft );

			}

			//4.
			this->analyze_fused_pose( *trunc_pose, *trunc_fusepose, *fused_pose_raw, *fused_pose_rlx, combo_identifier+"  "+utility::to_string( this_pair_rmsd) );
			//if( debugmode_ ){
			fused_pose_rlx->dump_pdb("fusedpose_rlx_super_"+utility::to_string(pose_pdbpos)+"_"+utility::to_string(fusepose_pdbpos)+".pdb");
			//}
		}// loop over candidate seqpos fuse pose

	} //loop over candidate seqpos pose

} //apply function

void
FusePosesNtoCMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{
	if ( tag->hasOption("fuse_pose") ) {
		fuse_pose_ = core::pose::PoseOP( new core::pose::Pose() );
		core::import_pose::pose_from_file( *fuse_pose_, tag->getOption<std::string>( "fuse_pose") , core::import_pose::PDB_file);
	} else utility_exit_with_message("FusePosesNtoCMover needs to be supplied with a fuse_pose.");

	if ( tag->hasOption("chains") ) {
		tr << "Chains specified in tag: ";
		chains_to_use_.clear();
		std::string chains( tag->getOption<std::string>("chains") );
		for ( std::string::const_iterator chain_it = chains.begin(); chain_it != chains.end(); ++chain_it ) {
			chains_to_use_.push_back( *chain_it );
			tr << *chain_it << " ";
		}
		tr << std::endl;
	}

	if ( tag->hasOption("nterm_span") ) {
		nterm_span_tag_ = tag->getOption<std::string>("nterm_span");
	} else utility_exit_with_message("FusePosesNtoCMover needs to be supplied with an nterm_span, i.e. the candidate residues for the N-terminal partner.");

	if ( tag->hasOption("cterm_span") ) {
		this->setup_candidate_pdbpos( candidate_pdbpos_fuse_pose_, tag->getOption<std::string>("cterm_span"));
	} else utility_exit_with_message("FusePosesNtoCMover needs to be supplied with an cterm_span, i.e. the candidate residues for the C-terminal partner.");

	if ( tag->hasOption("superpose_window") ) {
		superpose_window_ = tag->getOption<Size>( "superpose_window");
	}
	if ( tag->hasOption("relax_window") ) {
		relax_window_ = tag->getOption<Size>( "relax_window");
	}
	if ( tag->hasOption("max_number_allowed_clashes") ) {
		max_number_allowed_clashes_ = tag->getOption<Size>( "max_number_allowed_clashes");
	}
	if ( tag->hasOption("hardsphere_clash_limit") ) {
		hardsphere_clash_limit_ = tag->getOption<Size>( "hardsphere_clash_limit");
	}
	if ( tag->hasOption("max_allowed_rms") ) {
		max_allowed_rms_ = tag->getOption<Real>( "max_allowed_rms");
	}
	if ( tag->hasOption("debugmode") ) {
		debugmode_ = tag->getOption<bool>( "debugmode");
	}

	sfxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();

	{
		using namespace core::scoring;
		if ( sfxn_->has_zero_weight( chainbreak ) ) {
			tr << "Scorefunction specified has no chainbreak term turned on, setting the weight to 1.0";
			sfxn_->set_weight( chainbreak, 1.0 );
		}
		sfxn_nocb_ = sfxn_->clone();
		sfxn_nocb_->set_weight( chainbreak, 0.0);
	}
}

// note move this function to a separate mover
/// @details
/// the goal is to make a fold tree that insulates movement
/// to the fusion regions, without touching parts of the pose
/// that are not close to the fusion
/// this function only does something for dimeric poses
/// for monomeric poses nothind needs to be done,
/// and for trimeric fusion shit's too complicated
void
FusePosesNtoCMover::setup_relax_fold_tree(
	core::pose::Pose & pose
)
{
	if ( chain_mappings_.size() == 1 ) return;

	if ( chain_mappings_.size() > 2 ) utility_exit_with_message("Can't yet deal with fusions between poses that have more than two chains.");

	//ok, for two chains it's relatively easy
	//the chemical edge for the first chain is replaced by:
	//  a. a chemical edge till the beginning of fusion_regions (fusion_begin)
	//  b. a chemical edge from fusion_begin to the cut/fuse_position
	//  c. a jump from fusion_begin to cut+1, this is the only jump to move in default relax
	//  d. a chemical edge from cut+1 to the end of the cahin
	//
	//the chemical edge for the second chain is replaced by:
	//  e. a chemical edge till the cut
	//  f. a jump from the end of fusion_region1 chain to the end of fusion_region2 chain
	//  g. a reverse chemical edge from the end of the chain to cut+1

	core::kinematics::FoldTree const & old_fold_tree( pose.fold_tree() );
	core::kinematics::FoldTree new_ft;
	Size jump_count(old_fold_tree.num_jump() + 1);
	Size chain1_span_end(0);
	jumps_to_move_.clear();

	if ( debugmode_ ) {
		tr << "Incoming foldtree: " << std::endl << old_fold_tree << std::endl;
		tr.flush();
	}

	for ( auto const & e : old_fold_tree ) {

		//sanity check: we need to make sure that we don't have edges
		//originating or terminating in the fusion regions
		if ( ((e.start() >= (int) fusion_regions_[1].first) && (e.start() <= (int) fusion_regions_[1].second))
				|| ((e.stop() >= (int) fusion_regions_[1].first) && (e.stop() <= (int) fusion_regions_[1].second))
				|| ((e.start() >= (int) fusion_regions_[2].first) && (e.start() <= (int) fusion_regions_[2].second))
				|| ((e.stop() >= (int) fusion_regions_[2].first) && (e.stop() <= (int) fusion_regions_[2].second)) ) {
			utility_exit_with_message("Incipient fold tree has edges originating or terminating in the fusion regions, don't know what to do.");
		}
		//is this the chemical edge spanning the chain 1 fusion?
		if ( (e.start() < (int) fusion_regions_[1].first) && (e.stop() > (int) fusion_regions_[1].second) && !e.is_jump() ) {
			new_ft.add_edge( e.start(), fusion_regions_[1].first - 1, e.label() ); //a.
			new_ft.add_edge( fusion_regions_[1].first -1, fuse_positions_[1], e.label() ); //b.
			new_ft.add_edge( fusion_regions_[1].first-1, fuse_positions_[1] + 1, jump_count ); //c.
			jumps_to_move_.push_back( jump_count );
			++jump_count;
			new_ft.add_edge( fuse_positions_[1] + 1, e.stop(), e.label() ); //d.
			chain1_span_end = e.stop();
		} else if ( (e.start() < (int) fusion_regions_[2].first) && (e.stop() > (int) fusion_regions_[2].second) && !e.is_jump() ) {
			//or is this the chemical edge spanning the chain 2 fusion?
			new_ft.add_edge( e.start(), fuse_positions_[2], e.label() ); //e.
			if ( chain1_span_end ==0 ) utility_exit_with_message("Edge spanning a later chain was first detected in FoldTree.");
			new_ft.add_edge( chain1_span_end, e.stop(), jump_count ); //f.
			new_ft.add_edge( e.stop(), fuse_positions_[2] + 1, e.label() ); //g.
		} else {
			//otherwise we keep it
			new_ft.add_edge( e.start(), e.stop(), e.label() );
		}
	}

	tr << "The following fold tree was set up: " << std::endl << new_ft << std::endl;
	if ( !new_ft.check_fold_tree() ) utility_exit_with_message("Invalid fold tree after trying to set up for dimeric fusion relaxation");
	pose.fold_tree( new_ft );
}

//moving to sepaate mover end

protocols::moves::MoverOP
FusePosesNtoCMover::generate_default_relax_mover(
	core::pose::Pose const & //pose
) const
{
	using core::pack::task::operation::TaskOperationCOP;

	Size repeats( debugmode_ ? 1 : 5 );
	protocols::relax::FastRelaxOP frelax( new protocols::relax::FastRelax( sfxn_, repeats ) );

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_chi(true);
	mm->set_bb( false );
	mm->set_jump( false );
	for ( Size i = 1; i <= jumps_to_move_.size(); ++i ) mm->set_jump( jumps_to_move_[i], true );
	for ( Size i = 1; i <= fusion_regions_.size(); ++i ) {
		for ( Size j = fusion_regions_[i].first; j <= fusion_regions_[i].second; ++j ) mm->set_bb( j, true );
	}
	frelax->set_movemap( mm );

	core::pack::task::TaskFactoryOP taskf( new core::pack::task::TaskFactory() );
	taskf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline() ) );
	taskf->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking() ) );

	frelax->set_task_factory( taskf );

	return frelax;
}

void
FusePosesNtoCMover::setup_candidate_pdbpos(
	utility::vector1< Size > & pdbpos,
	std::string const tag_input
) const
{
	pdbpos.clear();

	utility::vector1< std::string > resvec( utility::string_split(tag_input, ',' ) );
	runtime_assert( resvec.size() < 3 );
	if ( resvec.size() == 1 ) {
		pdbpos.push_back(  core::Size( utility::string2int( resvec[1] ) ) );
		return;
	}
	Size lower_res( core::Size( utility::string2int( resvec[1] )));
	Size upper_res( core::Size( utility::string2int( resvec[2] )));
	runtime_assert( lower_res <= upper_res );

	for ( Size i = lower_res; i <= upper_res; ++i ) {
		pdbpos.push_back( i );
	}
}


core::id::AtomID_Map< core::id::AtomID >
FusePosesNtoCMover::generate_superposition_map(
	core::pose::Pose const & pose,
	core::pose::Pose const & fuse_pose,
	Size const pose_pdbpos,
	Size const fuse_pose_pdbpos
) const
{

	core::id::AtomID_Map< core::id::AtomID > atom_map( core::id::BOGUS_ATOM_ID );
	core::pose::initialize_atomid_map( atom_map, fuse_pose, core::id::BOGUS_ATOM_ID );

	//remember: pose is nterminal half, fuse_pose cterminal half
	for ( Size chain_it(1); chain_it <= chains_to_use_.size(); ++chain_it ) {

		char chain( chains_to_use_[ chain_it ] );
		Size pose_seqpos( pose.pdb_info()->pdb2pose( chain, pose_pdbpos ) );
		Size fuse_pose_seqpos( fuse_pose.pdb_info()->pdb2pose( chain, fuse_pose_pdbpos ) );

		//we have to figure out what the number of the chain is in the
		//internal numbering of the conformation objects, and what the
		//respective first and last residues are
		Size pose_chain_no( pose.chain( pose_seqpos) );
		Size fuse_pose_chain_no( fuse_pose.chain( fuse_pose_seqpos) );

		Size pose_chain_first_seqpos( pose.conformation().chain_begin( pose_chain_no ) );
		Size pose_chain_last_seqpos( pose.conformation().chain_end( pose_chain_no ) );

		Size fuse_pose_chain_first_seqpos( fuse_pose.conformation().chain_begin( fuse_pose_chain_no ) );
		Size fuse_pose_chain_last_seqpos( pose.conformation().chain_end( fuse_pose_chain_no ) );
		//chain numberings figured out

		//necessary because we're doing the fusion through superposition
		runtime_assert( (pose_seqpos != pose_chain_last_seqpos) || ( fuse_pose_seqpos != fuse_pose_chain_first_seqpos ) );

		//always superpose on the actual cut residues
		if ( fuse_pose_seqpos != fuse_pose_chain_first_seqpos ) {
			set_superpose_residues_in_atom_map( atom_map, fuse_pose, pose, fuse_pose_seqpos - 1, pose_seqpos, debugmode_ );
		}
		if ( pose_seqpos != pose_chain_last_seqpos ) {
			set_superpose_residues_in_atom_map( atom_map, fuse_pose, pose, fuse_pose_seqpos, pose_seqpos +1, debugmode_ );
		}

		//and if a window is requested, also add these residues to the
		//superposition
		for ( Size i=1; i <= superpose_window_; ++i ) {
			//the following two need to be integers because they can get negative..
			int pose_up_pos( pose_seqpos - i);
			int fuse_pose_up_pos( fuse_pose_seqpos - 1 - i);

			//but make sure we don't reach past the ends of the pose
			if ( (pose_up_pos >= (int) pose_chain_first_seqpos ) && (fuse_pose_up_pos >= (int) fuse_pose_chain_first_seqpos ) ) {
				//in the below call, the two ints get silenty type-converted to Size,
				//but this should be fine, since in the above if clause we make sure
				//that neither of them is negative
				set_superpose_residues_in_atom_map( atom_map, fuse_pose, pose, fuse_pose_up_pos, pose_up_pos, debugmode_ );
			}

			Size pose_down_pos(pose_seqpos + i + 1);
			Size fuse_pose_down_pos( fuse_pose_seqpos + i );
			if ( (pose_down_pos <= pose_chain_last_seqpos) && (fuse_pose_down_pos <= fuse_pose_chain_last_seqpos ) ) {
				set_superpose_residues_in_atom_map( atom_map, fuse_pose, pose, fuse_pose_down_pos, pose_down_pos, debugmode_ );
			}
		}// loop over window
	} //loop over chains
	return atom_map;
} //generate_superposition_map

FusePosesNtoCMover::Size
FusePosesNtoCMover::count_hardsphere_clashes(
	core::pose::Pose const &, //pose,
	core::pose::Pose const & //fuse_pose,
) const
{

	//stubbed out
	return 0; //nothing ever clashes ;)
}

/// @details
/// right now this assumes that the variant types have already been setup
/// stitches the poses together according to the previously defined chain
/// mappings. if fuse pose has additional chain (i.e. ligands ), these
/// will also be copied
void
FusePosesNtoCMover::fuse_poses(
	core::pose::Pose & pose,
	core::pose::Pose const & fuse_pose
)
{
	std::set< Size > fused_fusepose_chains;
	fuse_positions_.clear();
	fusion_regions_.clear();

	for ( Size i(1); i <= chain_mappings_.size(); ++i ) {

		Size fusepose_chain( chain_mappings_[i].first ), pose_chain( chain_mappings_[i].second );
		fused_fusepose_chains.insert( fusepose_chain);
		tr << "About to fuse n-terminal residue " << pose.residue_type( pose.conformation().chain_end( pose_chain ) ).name() << " to c-terminal residue " << fuse_pose.residue_type( fuse_pose.conformation().chain_begin( fusepose_chain ) ).name() << std::endl;

		//the following line should work because chain_mappings is sorted
		//according to the chain ordering of the pose
		fuse_positions_.push_back( pose.conformation().chain_end( pose_chain ) );

		for ( Size insres( fuse_pose.conformation().chain_begin( fusepose_chain ) ), chain_end(fuse_pose.conformation().chain_end( fusepose_chain ) ); insres <= chain_end; ++insres ) {

			pose.append_polymer_residue_after_seqpos( fuse_pose.residue( insres ), pose.conformation().chain_end( pose_chain ), false );
		}

	} //loop over chain mappings

	//now we need to loop over additional fuse pose chains and copy them over
	//random decision: append by jump to the end of the first fused chain
	for ( Size i(1); i <= fuse_pose.conformation().num_chains(); ++i ) {
		if ( fused_fusepose_chains.find( i ) != fused_fusepose_chains.end() ) continue;

		Size anchor_res = pose.conformation().chain_end( chain_mappings_[1].second );
		pose.append_residue_by_jump(
			fuse_pose.residue( fuse_pose.conformation().chain_begin( i ) ), anchor_res, "", "", true );

		for ( Size rescounter( fuse_pose.conformation().chain_begin( i ) + 1), chain_end( fuse_pose.conformation().chain_end( i ) ); rescounter <= chain_end; ++rescounter ) {
			//no idea if this shit is going to work
			pose.append_polymer_residue_after_seqpos( fuse_pose.residue( rescounter ), pose.size() /* hope this is correct */, false );
		}
	}
	if ( debugmode_ ) {
		tr << "The fused positions in pose numbering are: ";
		for ( Size i =1; i <= fuse_positions_.size(); ++i ) tr << fuse_positions_[i] << ", ";
		tr << std::endl;
	}

	//now we also set up the surrounding fusion regions, using dssp to do so
	core::scoring::dssp::Dssp ss_pose(pose);
	//ss_pose.insert_ss_into_pose( pose );
	//the 4 is very arbitrary, but since the fused pose is not put together yet
	//dssp is prolly not accurate at the actual fuseposition
	Size offset(4);

	for ( Size i = 1; i <= fuse_positions_.size(); ++i ) {
		char fuse_region_ss( ss_pose.get_dssp_secstruct( fuse_positions_[i] - offset ) );
		int nterm_end(0); Size cterm_end(0), counter(1);

		while ( (nterm_end == 0 ) || (cterm_end == 0 ) ) {
			if ( nterm_end == 0 ) {
				Size nterm_check_pos = fuse_positions_[i] - offset - counter;
				if ( nterm_check_pos < 1 ) nterm_end = 1;
				else if ( ss_pose.get_dssp_secstruct( nterm_check_pos ) != fuse_region_ss ) nterm_end = nterm_check_pos;
			}

			if ( cterm_end == 0 ) {
				Size cterm_check_pos = fuse_positions_[i] + offset + counter;
				if ( cterm_check_pos > pose.conformation().chain_end( fuse_positions_[i] ) ) cterm_end = pose.conformation().chain_end( fuse_positions_[i] );
				else if ( ss_pose.get_dssp_secstruct( cterm_check_pos ) != fuse_region_ss ) cterm_end = cterm_check_pos;
			}
			counter++;
		} //while nterm_end / cterm_end
		if ( relax_window_ != 0 ) {
			nterm_end - relax_window_ <= 0 ? nterm_end = 0 : nterm_end = nterm_end - relax_window_;
			cterm_end + relax_window_ > pose.conformation().chain_end( fuse_positions_[i] ) ? cterm_end = pose.conformation().chain_end( fuse_positions_[i] ) : cterm_end = cterm_end + relax_window_;
		}
		fusion_regions_.push_back( std::pair< Size, Size >( nterm_end, cterm_end ) );
		if ( debugmode_ ) {
			tr << "adding fusion region between " << nterm_end << " and " << cterm_end << std::endl;
			tr << "The residue types at the fusion point are " << pose.residue_type( fuse_positions_[i]).name() << " and " << pose.residue_type( fuse_positions_[i] + 1).name() << std::endl;
		}
	} //i loop over fusion positions

} //fuse_poses


/// @details
/// we'll also take care of inserting the proper chainbreak variants here
/// although maybe it's better to do this in the above fuse_poses function
core::pose::PoseOP
FusePosesNtoCMover::truncate_pose_at_fusion_site(
	core::pose::Pose const & pose,
	Size const fusion_pdbpos,
	bool nterm_truncation
) const
{

	core::pose::PoseOP to_return( new core::pose::Pose( pose ) );
	for ( Size chain_it(1); chain_it <= chains_to_use_.size(); ++chain_it ) {

		char chain( chains_to_use_[ chain_it ] );
		Size fusion_seqpos( to_return->pdb_info()->pdb2pose( chain, fusion_pdbpos ) );
		Size pose_chain_no( to_return->chain( fusion_seqpos) );

		Size last_truncate_pos( nterm_truncation ? to_return->conformation().chain_begin( pose_chain_no ) : to_return->conformation().chain_end( pose_chain_no ) );

		if ( fusion_seqpos == last_truncate_pos ) {

			if ( nterm_truncation ) {
				if ( to_return->residue_type( fusion_seqpos ).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) {
					core::conformation::remove_lower_terminus_type_from_conformation_residue( to_return->conformation(), fusion_seqpos );
				}
				add_variant_type_to_conformation_residue( to_return->conformation(), core::chemical::CUTPOINT_UPPER, fusion_seqpos );
			} else {
				if ( to_return->residue_type( fusion_seqpos ).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) ) {
					core::conformation::remove_upper_terminus_type_from_conformation_residue( to_return->conformation(), fusion_seqpos );
				}
				add_variant_type_to_conformation_residue( to_return->conformation(), core::chemical::CUTPOINT_LOWER, fusion_seqpos );
			}
			continue;
		}

		if ( nterm_truncation ) {
			to_return->conformation().delete_residue_range_slow( last_truncate_pos, fusion_seqpos - 1 );
			Size newfirstpos( to_return->conformation().chain_begin( pose_chain_no ) );

			if ( to_return->residue_type( newfirstpos ).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) {
				core::conformation::remove_lower_terminus_type_from_conformation_residue( to_return->conformation(), newfirstpos );
			}
			add_variant_type_to_conformation_residue( to_return->conformation(), core::chemical::CUTPOINT_UPPER, newfirstpos );
		} else {
			to_return->conformation().delete_residue_range_slow( fusion_seqpos + 1, last_truncate_pos);
			Size newlastpos( to_return->conformation().chain_end( pose_chain_no ) );

			if ( to_return->residue_type( newlastpos ).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) ) {
				core::conformation::remove_upper_terminus_type_from_conformation_residue( to_return->conformation(), newlastpos );
			}
			add_variant_type_to_conformation_residue( to_return->conformation(), core::chemical::CUTPOINT_LOWER, newlastpos );

		}

	} //loop over all chains
	return to_return;
}


/// @details
/// still unclear about how to do this best,
/// to start with, let's do the following
/// 1. energy diff between fusion_rlx and cterm+nterm
/// 2. rms between fusion_raw and fusion_rlx on the fusion regions
/// 3. ABEGO for the fusion regions
void
FusePosesNtoCMover::analyze_fused_pose(
	core::pose::Pose & nterm_half,
	core::pose::Pose & cterm_half,
	core::pose::Pose & fusion_raw,
	core::pose::Pose & fusion_rlx,
	std::string const identifier
) const
{

	std::string report_string("FUSEANALYSIS "+identifier+" ");

	//1.
	core::Real score_no_cb_fusion( (*sfxn_nocb_)( fusion_rlx ) ), score_no_cb_halves( (*sfxn_nocb_)(nterm_half) + (*sfxn_nocb_)(cterm_half) );
	core::Real scorediff( score_no_cb_fusion - score_no_cb_halves );
	core::Real cb_fusion( (*sfxn_)(fusion_rlx) - score_no_cb_fusion );

	//2.
	std::list< Size > residues; //rms code wants a list
	for ( Size i(1); i <= fusion_regions_.size(); ++i ) {
		for ( Size j( fusion_regions_[i].first); j <= fusion_regions_[i].second; ++j ) residues.push_back( j );
	}

	core::Real rms( core::scoring::CA_rmsd( fusion_raw, fusion_rlx, residues ) );

	//3.
	utility::vector1< std::string > complete_abegos( core::sequence::get_abego( fusion_rlx ) );
	core::scoring::dssp::Dssp ss_frlx(fusion_rlx);
	//ss_frlx.insert_ss_into_pose( fusion_rlx );
	utility::vector1< std::string > fusion_region_abegos;
	utility::vector1< std::string > fusion_region_dssp;
	for ( Size i(1); i <= fusion_regions_.size(); ++i ) {
		runtime_assert( fusion_regions_[i].first < fuse_positions_[i] );
		runtime_assert( fuse_positions_[i] < fusion_regions_[i].second );
		std::string this_reg_abego("");
		std::string this_reg_dssp("");
		for ( Size j( fusion_regions_[i].first); j <= fuse_positions_[i]; ++j ) {
			this_reg_abego += complete_abegos[j];
			this_reg_dssp += ss_frlx.get_dssp_secstruct( j );
		}
		this_reg_abego += "-"; //highlighting fusion
		this_reg_dssp += "-";
		for ( Size j( fuse_positions_[i] + 1); j <= fusion_regions_[i].second; ++j ) {
			this_reg_abego += complete_abegos[j];
			this_reg_dssp += ss_frlx.get_dssp_secstruct( j );
		}
		fusion_region_abegos.push_back( this_reg_abego );
		fusion_region_dssp.push_back( this_reg_dssp);
	}

	report_string += ( utility::to_string( scorediff ) + "   " + utility::to_string( cb_fusion ) + "   " + utility::to_string( rms ) + "   ");
	for ( Size i(1); i <= fusion_regions_.size(); ++i ) report_string += (fusion_region_abegos[i] + "  ");
	for ( Size i(1); i <= fusion_regions_.size(); ++i ) report_string += (fusion_region_dssp[i] + "  ");

	std::cout << report_string << std::endl;
}

void
FusePosesNtoCMover::setup_chain_mappings(
	core::pose::Pose const & pose )
{
	chain_mappings_.clear();
	for ( Size chain(1); chain <= chains_to_use_.size(); ++chain ) {

		for ( Size i = 1; i <= pose.conformation().num_chains() ; ++i ) {
			char pdb_chain_code( pose.pdb_info()->chain( pose.conformation().chain_begin( i ) ) );

			if ( pdb_chain_code == chains_to_use_[ chain ] ) {
				for ( Size j =1; j <= fuse_pose_->conformation().num_chains() ; ++j ) {
					char fusepdb_chain_code( fuse_pose_->pdb_info()->chain( fuse_pose_->conformation().chain_begin( j ) ) );
					if ( fusepdb_chain_code == chains_to_use_[chain] ) {
						chain_mappings_.push_back( std::pair<Size, Size>( j, i ) );
						break;
					}
				}//loop over fuse_pose chains
			}
			if ( chain_mappings_.size() == chain ) break; //in case this chain was found
		} //loop over pose chains
	}
	if ( chain_mappings_.size() != chains_to_use_.size() ) utility_exit_with_message("The desired chain numbers were not found in each of the two input poses.");

	//we sort the mappings in ascending order of pose chains
	std::sort( chain_mappings_.begin(), chain_mappings_.end(), custom_pair_compare );
}


/// @details
/// right now superimposes onto N, CA, C, but could
/// be changed
void
set_superpose_residues_in_atom_map(
	core::id::AtomID_Map< core::id::AtomID > & atom_map,
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size const seqpos1,
	core::Size const seqpos2,
	bool debugout
)
{
	if ( debugout ) tr << "superimposing pose1 residue " << seqpos1 << " and pose2 residue " << seqpos2 << std::endl;
	atom_map.set( core::id::AtomID( pose1.residue( seqpos1 ).atom_index("CA"), seqpos1 ), core::id::AtomID( pose2.residue( seqpos2 ).atom_index("CA"), seqpos2 ) );

	atom_map.set( core::id::AtomID( pose1.residue( seqpos1 ).atom_index("N"), seqpos1 ), core::id::AtomID( pose2.residue( seqpos2 ).atom_index("N"), seqpos2 ) );

	atom_map.set( core::id::AtomID( pose1.residue( seqpos1 ).atom_index("C"), seqpos1 ), core::id::AtomID( pose2.residue( seqpos2 ).atom_index("C"), seqpos2 ) );
}


// XRW TEMP std::string
// XRW TEMP FusePosesNtoCMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return FusePosesNtoCMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FusePosesNtoCMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FusePosesNtoCMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FusePosesNtoCMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "FusePosesNtoCMover";
// XRW TEMP }

std::string FusePosesNtoCMover::get_name() const {
	return mover_name();
}

std::string FusePosesNtoCMover::mover_name() {
	return "FusePosesNtoCMover";
}

void FusePosesNtoCMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction chain_string;
	chain_string.name( "chain_string" );
	chain_string.base_type( xs_string );
	chain_string.add_restriction( xsr_pattern, "[" + chr_chains_nonrepeated() + "]*" );
	xsd.add_top_level_element( chain_string );

	AttributeList attlist; // TO DO: add attributes to this list
	attlist
		+ XMLSchemaAttribute::required_attribute( "fuse_pose", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "chains", "chain_string", "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "nterm_span", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "cterm_span", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "superpose_window", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "relax_window", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "max_number_allowed_clashes", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "hardsphere_clash_limit", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "max_allowed_rms", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute( "debugmode", xsct_rosetta_bool, "XRW TO DO" );
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string FusePosesNtoCMoverCreator::keyname() const {
	return FusePosesNtoCMover::mover_name();
}

protocols::moves::MoverOP
FusePosesNtoCMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new FusePosesNtoCMover );
}

void FusePosesNtoCMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FusePosesNtoCMover::provide_xml_schema( xsd );
}



SetupCoiledCoilFoldTreeMover::SetupCoiledCoilFoldTreeMover()
: Mover(),
	chain2_cutpos_(0), add_chainbreak_variants_(false)
{}

SetupCoiledCoilFoldTreeMover::SetupCoiledCoilFoldTreeMover( SetupCoiledCoilFoldTreeMover const & ) = default;


SetupCoiledCoilFoldTreeMover::~SetupCoiledCoilFoldTreeMover() = default;

protocols::moves::MoverOP
SetupCoiledCoilFoldTreeMover::clone() const{
	return protocols::moves::MoverOP( new SetupCoiledCoilFoldTreeMover( *this ) );
}


// XRW TEMP std::string
// XRW TEMP SetupCoiledCoilFoldTreeMover::get_name() const
// XRW TEMP {
// XRW TEMP  return "SetupCoiledCoilFoldTreeMover";
// XRW TEMP }


/// @details
/// the chemical edge for the first chain is kept
///
///the chemical edge spanning the cutpos is replaced by
///  a. a chemical edge from the start of that edge till chain2_cutpos_
///  b. backwards chemical edge from the end of that edge till chain2_cutpos_ +1
///
///the jump from chain1 to chain2 is extended by
/// another jump from the end of chain 1 till the end of chain2


void
SetupCoiledCoilFoldTreeMover::apply( core::pose::Pose & pose )
{

	core::Size cut_seqpos( pose.pdb_info()->pdb2pose( coiled_coil_chains_[2], chain2_cutpos_ ) );

	//figure out what the pose numbers of the desired chains is
	utility::vector1< core::Size > pose_chains;
	for ( Size chain(1); chain <= coiled_coil_chains_.size(); ++chain ) {
		for ( Size i = 1; i <= pose.conformation().num_chains(); ++i ) {
			char pdb_chain_code( pose.pdb_info()->chain( pose.conformation().chain_begin( i ) ) );
			if ( pdb_chain_code == coiled_coil_chains_[ chain ] ) pose_chains.push_back( i );
		} //loop over pose chains
	}
	//figured out, now mod the fold tree

	core::kinematics::FoldTree const & old_fold_tree( pose.fold_tree() );
	core::kinematics::FoldTree new_ft;
	Size jump_count(old_fold_tree.num_jump() + 1);

	tr << "CoiledCoilFoldTreeMover Incoming foldtree: " << std::endl << old_fold_tree << std::endl << " and cut_seqpos is " << cut_seqpos << std::endl;
	tr.flush();
	bool relevant_jump_found(false);
	for ( auto const & e : old_fold_tree ) {

		std::cout << "looping edge with label " << e. label() << ", start " << e.start() << " and stop " << e.stop() << std::endl;

		if ( (e.label() == core::kinematics::Edge::PEPTIDE) && (e.start() < (int) cut_seqpos) && (e.stop() > (int) cut_seqpos) ) {
			std::cout << "MEEP" << std::endl;
			new_ft.add_edge( e.start(), cut_seqpos, e.label() );
			new_ft.add_edge( e.stop(), cut_seqpos + 1, e.label() );
		} else if ( (e.is_jump() ) && ( pose.chain( e.start() ) == (int) pose_chains[1]) && ( pose.chain( e.stop() ) == (int) pose_chains[2] ) ) {
			relevant_jump_found = true;
			new_ft.add_edge( pose.conformation().chain_begin( pose_chains[1] ), pose.conformation().chain_begin( pose_chains[2] ), e.label() );
			new_ft.add_edge( pose.conformation().chain_end( pose_chains[1] ), pose.conformation().chain_end( pose_chains[2] ), jump_count );
		} else {
			new_ft.add_edge( e.start(), e.stop(), e.label() );
		}
	} //loop over fold tree edges

	tr << "CoiledCoilFoldTreeMover new foldtree: " << std::endl << new_ft << std::endl;
	if ( !relevant_jump_found ) utility_exit_with_message("A jump connecting the two chains of the coiled coil was not found, dunno what to do.");
	if ( !new_ft.check_fold_tree() ) utility_exit_with_message("Invalid fold tree after trying to set up for dimeric coiled coil.");
	pose.fold_tree( new_ft );

	if ( add_chainbreak_variants_ ) {
		add_variant_type_to_conformation_residue( pose.conformation(), core::chemical::CUTPOINT_UPPER, cut_seqpos + 1 );
		add_variant_type_to_conformation_residue( pose.conformation(), core::chemical::CUTPOINT_LOWER, cut_seqpos );

	}
}

void
SetupCoiledCoilFoldTreeMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{

	if ( tag->hasOption("chains") ) {
		tr << "Chains specified in tag: ";
		coiled_coil_chains_.clear();
		std::string chains( tag->getOption<std::string>("chains") );
		for ( std::string::const_iterator chain_it = chains.begin(); chain_it != chains.end(); ++chain_it ) {
			coiled_coil_chains_.push_back( *chain_it );
			tr << *chain_it << " ";
		}
		tr << std::endl;
	}

	if ( tag->hasOption("chain2_cutpos") ) {
		chain2_cutpos_ = tag->getOption<core::Size>("chain2_cutpos", 0);
	}

	if ( tag->hasOption("add_chainbreak_variants") ) {
		add_chainbreak_variants_ = tag->getOption<bool>("add_chainbreak_variants", false);
	}
}


// XRW TEMP std::string
// XRW TEMP SetupCoiledCoilFoldTreeMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SetupCoiledCoilFoldTreeMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SetupCoiledCoilFoldTreeMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SetupCoiledCoilFoldTreeMover() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SetupCoiledCoilFoldTreeMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "SetupCoiledCoilFoldTreeMover";
// XRW TEMP }

std::string SetupCoiledCoilFoldTreeMover::get_name() const {
	return mover_name();
}

std::string SetupCoiledCoilFoldTreeMover::mover_name() {
	return "SetupCoiledCoilFoldTreeMover";
}

void SetupCoiledCoilFoldTreeMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;

	// AMW namespacing for safety because otherwise kill me
	XMLSchemaRestriction chain_string;
	chain_string.name( "SetupCoiledCoilFoldTreeMover_chain_string" );
	chain_string.base_type( xs_string );
	chain_string.add_restriction( xsr_pattern, "[" + chr_chains_nonrepeated() + "]*" );
	xsd.add_top_level_element( chain_string );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "chains", "SetupCoiledCoilFoldTreeMover_chain_string", "XRW TO DO" )
		+ XMLSchemaAttribute( "chain2_cutpos", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "add_chainbreak_variants", xsct_rosetta_bool, "XRW TO DO" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string SetupCoiledCoilFoldTreeMoverCreator::keyname() const {
	return SetupCoiledCoilFoldTreeMover::mover_name();
}

protocols::moves::MoverOP
SetupCoiledCoilFoldTreeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetupCoiledCoilFoldTreeMover );
}

void SetupCoiledCoilFoldTreeMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetupCoiledCoilFoldTreeMover::provide_xml_schema( xsd );
}



/// @details effin stooopid that this is necessary
std::map< core::id::AtomID, core::id::AtomID >
convert_AtomID_Map_to_std_map(
	core::id::AtomID_Map< core::id::AtomID > const & atom_map )
{

	std::map< core::id::AtomID, core::id::AtomID > to_return;

	for ( core::Size rescounter(1); rescounter <= atom_map.size(); ++rescounter ) {
		core::Size atoms_this_res( atom_map.n_atom( rescounter ) );
		for ( core::Size atcounter(1); atcounter <= atoms_this_res; ++atcounter ) {

			core::id::AtomID to_probe( atcounter, rescounter );
			core::id::AtomID mapped_id( atom_map[ to_probe ] );
			if ( mapped_id != core::id::BOGUS_ATOM_ID ) {
				to_return.insert( std::pair< core::id::AtomID, core::id::AtomID>( to_probe, mapped_id ) );
			}
		} //loop over atoms for res
	}//loop over residues
	return to_return;
}


}
}
