// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <protocols/motifs/MotifDnaPacker.hh>
#include <protocols/motifs/MotifDnaPackerCreator.hh>

#include <protocols/dna/typedefs.hh>
#include <protocols/dna/DnaInterfacePacker.hh>
#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/util.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/dna/WatsonCrickRotamerCouplings.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/motifs/motif_utils.hh>
#include <protocols/motifs/MotifSearch.hh>
#include <protocols/toolbox/rotamer_set_operations/SpecialRotamerRotSetOps.hh>
#include <protocols/motifs/Motif.hh> // REQUIRED FOR WINDOWS

//#include <devel/blab/opte/sidechain_relax.hh>
//#include <devel/blab/motif/MotifData.hh>
//#include <devel/blab/motif/util.hh>
//#include <devel/blab/motif/util_classes.hh>
//#include <devel/blab/motif/loop_rebuild.hh>
//#include <devel/blab/motif/loop_rebuild_movers.hh>
//#include <devel/blab/loops/util.hh>
//#include <devel/dna/util.hh>

#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
//#include <core/pack/task/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//#include <basic/tracer.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh> // lead_zero_string_of

// c++ headers
#include <vector> // for rot_to_pack
#include <iostream>
#include <sstream>


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace protocols {
namespace motifs {

using utility::vector1;
using utility::string_split;
using namespace core;
using namespace chemical;
using namespace conformation;
using namespace optimization;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace pose;
using namespace scoring;
using namespace protocols::dna;

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

using basic::t_warning;
using basic::t_info;
using basic::t_debug;
static THREAD_LOCAL basic::Tracer TR( "devel.motifs.MotifDnaPacker" );

std::string
MotifDnaPackerCreator::keyname() const
{
	return MotifDnaPackerCreator::mover_name();
}

protocols::moves::MoverOP
MotifDnaPackerCreator::create_mover() const
{
	return protocols::moves::MoverOP( new MotifDnaPacker );
}

std::string
MotifDnaPackerCreator::mover_name()
{
	return "MotifDnaPacker";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief lightweight default constructor
MotifDnaPacker::MotifDnaPacker()
: protocols::moves::Mover("MotifDnaPacker"),
	dna_packer_(/* 0 */),
	motif_search_(/* 0 */),
	scorefxn_(/* 0 */),
	pdboutput_(/* 0 */),
	targeted_dna_(/* 0 */),
	run_motifs_(false),
	expand_motifs_(false),
	aromatic_motifs_(false),
	minimize_dna_(false),
	special_rotweight_(-40.0),
	num_repacks_(5),
	flex_dna_sugar_(false),
	dna_design_(false)
{
	init_options();
}

/// @brief functional constructor
MotifDnaPacker::MotifDnaPacker(
	ScoreFunctionOP scorefxn_in,
	bool minimize,
	std::string filename_root
) : protocols::moves::Mover("MotifDnaPacker"),
	scorefxn_( scorefxn_in ),
	run_motifs_(false),
	expand_motifs_(false),
	aromatic_motifs_(false),
	minimize_dna_(false),
	special_rotweight_(-40.0),
	num_repacks_(5),
	flex_dna_sugar_(false),
	dna_design_(false)
{
	init_options();
	filename_root_ = filename_root;
	dna_packer_ = protocols::dna::DnaInterfacePackerOP( new DnaInterfacePacker( scorefxn_in, minimize, filename_root ) );
	protocols::motifs::MotifSearchOP motif_search( new protocols::motifs::MotifSearch );
	motif_search_ = motif_search;
	pdboutput_ = protocols::dna::PDBOutputOP( new PDBOutput );
}

/// @brief destructor
MotifDnaPacker::~MotifDnaPacker()= default;

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MotifDnaPacker::fresh_instance() const
{
	return protocols::moves::MoverOP( new MotifDnaPacker );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MotifDnaPacker::clone() const
{
	return protocols::moves::MoverOP( new MotifDnaPacker( *this ) );
}

std::string
MotifDnaPacker::get_name() const {
	return MotifDnaPackerCreator::mover_name();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief
/// @author sthyme
void
MotifDnaPacker::apply( Pose & pose )
{
	init_standard( pose );

	TaskFactoryOP taskfactory( new TaskFactory );
	taskfactory->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	if ( option[ OptionKeys::packing::resfile ].user() ) {
		taskfactory->push_back( TaskOperationCOP( new ReadResfile ) );
	}

	utility::vector1< core::Size > protein_positions, dna_positions, dna_design_positions, preventrepack;
	utility::vector1< utility::vector1< core::Size > > base_partners;
	core::Size nres( pose.size() );
	for ( Size p_index(1); p_index<=nres; ++p_index ) {
		if ( pose.residue_type( p_index ).is_DNA() ) {
			dna_positions.push_back( p_index );
		}
		if ( pose.residue_type( p_index ).is_protein() ) {
			protein_positions.push_back( p_index );
		} else {
			preventrepack.push_back( p_index );
		}
	}

	bool target_dna_empty( false );
	if ( ! targeted_dna_.empty() ) {
		dna_design_positions = protocols::motifs::defs2vector( pose, targeted_dna_ );
		core::scoring::dna::BasePartner const & partner( core::scoring::dna::retrieve_base_partner_from_pose( pose ) );
		utility::vector1< core::Size > complete_ddps;
		for ( Size i(1); i<=dna_design_positions.size(); ++i ) {
			complete_ddps.push_back(dna_design_positions[i]);
			complete_ddps.push_back(partner[dna_design_positions[i]]);
		}
		base_partners.push_back( complete_ddps );
	} else {
		target_dna_empty = true;
		core::scoring::dna::BasePartner const & partner( core::scoring::dna::retrieve_base_partner_from_pose( pose ) );
		utility::vector1< core::Size > complete_design_positions;
		for ( Size i(1); i<=( dna_positions.size() / 2 ); ++i ) {
			utility::vector1< core::Size > pair;
			if ( partner[dna_positions[i]] != 0 ) {
				pair.push_back( dna_positions[i] );
				pair.push_back( partner[dna_positions[i]] );
			}
			if ( ! pair.empty() ) {
				base_partners.push_back( pair );
			}
		}
		// NEED THE INPUTS HERE TO BE FROM THE OPTIONS? OR DEFAULT AND SPREAD TO ALL DNAIFs
		RestrictDesignToProteinDNAInterfaceOP rdtpdi( new RestrictDesignToProteinDNAInterface );
		protocols::dna::DnaInterfaceFinderOP complete_interface_( new protocols::dna::DnaInterfaceFinder( 10 * 10, 3.9 * 3.9, 6., true ) );
		complete_interface_->determine_protein_interface( pose, protein_positions, dna_positions );
		//// rdtpdi->set_reference_pose( starting_pose_ );
		rdtpdi->copy_interface( complete_interface_ );
		rdtpdi->set_forget_chains_and_interface( false );
		TaskFactoryOP taskfactory0( new TaskFactory( *taskfactory ) );
		taskfactory0->push_back( rdtpdi );
		PackerTaskOP ptask_ = taskfactory0->create_task_and_apply_taskoperations( pose );
		for ( core::Size j(1); j<=nres; ++j ) {
			if ( !pose.residue_type(j).is_protein() ) continue;
			ResidueLevelTask const & rtask( ptask_->residue_task(j) );
			if ( rtask.being_designed() ) {
				complete_design_positions.push_back(j);
			}
		}
		if ( run_motifs_ ) {
			motif_search_->run( pose, complete_design_positions );
		}
	}

	for ( Size bp(1); bp<=base_partners.size(); ++bp ) {
		pose = *starting_pose_;
		utility::vector1< core::Size > design_positions;
		std::map< core::Size, pack::rotamer_set::Rotamers > rotamer_map;
		std::map< core::Size, std::set< std::string > > types_map;
		std::set< core::Size > src_pos;
		std::list< std::string > info_lines;
		std::stringstream info_lines_str;
		// NEED THE INPUTS HERE TO BE FROM THE OPTIONS? OR DEFAULT AND SPREAD TO ALL DNAIFs
		protocols::dna::DnaInterfaceFinderOP interface_( new protocols::dna::DnaInterfaceFinder( 10 * 10, 3.9 * 3.9, 6., true ) );
		interface_->determine_protein_interface( pose, protein_positions, base_partners[bp] );
		RestrictDesignToProteinDNAInterfaceOP rdtpdi2( new RestrictDesignToProteinDNAInterface );
		// need or not??
		////// rdtpdi->set_reference_pose( starting_pose_ );
		rdtpdi2->copy_interface( interface_ );

		DnaChainsOP dna_chains( new DnaChains );
		find_basepairs( pose, *dna_chains );
		DnaDesignDefOPs targetdefs;
		for ( Size i(1); i<=(base_partners[bp].size()); ++i ) {
			if ( dna_chains->is_top(base_partners[bp][i]) ) {
				std::stringstream def;
				def << pose.pdb_info()->chain( base_partners[bp][i]) << "." << pose.pdb_info()->number(base_partners[bp][i]) << "." << pose.residue(base_partners[bp][i]).name3();
				DnaDesignDefOP target( new DnaDesignDef( def.str() ) );
				targetdefs.push_back( target );
			}
		}
		rdtpdi2->set_forget_chains_and_interface( false );
		TaskFactoryOP taskfactory1( new TaskFactory( *taskfactory ) );
		TaskFactoryOP taskfactory2( new TaskFactory( *taskfactory ) );

		targeted_dna_ = targetdefs;
		if ( ! dna_design_ ) {
			rdtpdi2->copy_targeted_dna( targetdefs );
		} else {
			utility::vector1< Size > restricted_DNA;
			for ( Size d(1) ; d <= dna_positions.size(); ++d ) {
				bool r_DNA( true );
				for ( Size i(1); i<=(base_partners[bp].size()); ++i ) {
					if ( dna_positions[d] == base_partners[bp][i] ) {
						r_DNA = false;
					}
				}
				if ( r_DNA ) {
					restricted_DNA.push_back( dna_positions[d] );
				}
			}
			OperateOnCertainResiduesOP fixedDNA( new OperateOnCertainResidues );
			fixedDNA->residue_indices( restricted_DNA );
			fixedDNA->op( ResLvlTaskOperationCOP( new PreventRepackingRLT ) );
			WatsonCrickRotamerCouplingsOP wrc( new WatsonCrickRotamerCouplings );
			taskfactory1->push_back( fixedDNA );
			taskfactory1->push_back( wrc );
			taskfactory2->push_back( fixedDNA );
			taskfactory2->push_back( wrc );
		}

		taskfactory1->push_back( rdtpdi2 );
		taskfactory2->push_back( rdtpdi2 );
		PackerTaskOP ptask_ = taskfactory2->create_task_and_apply_taskoperations( pose );
		for ( core::Size j(1); j<=nres; ++j ) {
			if ( !pose.residue_type(j).is_protein() ) continue;
			ResidueLevelTask const & rtask( ptask_->residue_task(j) );
			if ( rtask.being_designed() ) {
				design_positions.push_back(j);
				std::stringstream numberandchain;
				numberandchain << pose.pdb_info()->number( j ) << pose.pdb_info()->chain( j );
				info_lines_str << numberandchain.str() << ",";
			}
		}
		info_lines.push_back( info_lines_str.str() );

		if ( design_positions.empty() ) continue;

		std::string mot_name3( "NOMotifs" );
		std::stringstream filename;
		filename << filename_root_ << "_" << mot_name3;

		core::Real zero_special_rotweight( -0.0 );
		scorefxn_->set_weight( special_rot, ( zero_special_rotweight ) );

		dna_packer_->task_factory( taskfactory1 );
		dna_packer_->set_filename_root( filename.str() );
		dna_packer_->apply( pose );

		bool const overwrite_old_info(true);
		pdboutput_->add_info( "REMARK DESIGNED POSITIONS " + filename_root_ + ":", info_lines, !overwrite_old_info );
		pdboutput_->score_function( *scorefxn_ );
		(*pdboutput_)(pose, dna_packer_->pdbname() + ".pdb");

		dna_packer_->clear_initialization();

		if ( run_motifs_ && ( ! target_dna_empty ) ) {
			motif_search_->run( pose, design_positions );
			if ( option[ OptionKeys::motifs::quick_and_dirty ].user() ) break;
			run_motifs( pose, design_positions, src_pos, rotamer_map, types_map, info_lines, taskfactory1 );
		} else {
			if ( option[ OptionKeys::motifs::quick_and_dirty ].user() ) break;
			run_motifs( pose, design_positions, src_pos, rotamer_map, types_map, info_lines, taskfactory1 );
		}
		if ( run_motifs_ && expand_motifs_ ) {
			if ( option[ OptionKeys::motifs::quick_and_dirty ].user() ) break;
			expand_motifs( pose, design_positions, src_pos, rotamer_map, types_map, info_lines, taskfactory1 );
		}
		if ( run_motifs_ && aromatic_motifs_ ) {
			if ( option[ OptionKeys::motifs::quick_and_dirty ].user() ) break;
			aromatic_motifs( pose, design_positions, src_pos, rotamer_map, types_map, info_lines, taskfactory1 );
		}
	}
}

void
MotifDnaPacker::init_standard( Pose & pose )
{
	//protocols::motifs::make_dna_mutations( pose );
	starting_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new pose::Pose( pose ) ) );
	pdboutput_->reference_pose( *starting_pose_ );
	if ( option[ OptionKeys::dna::design::dna_defs ].user() ) {
		load_dna_design_defs_from_options( targeted_dna_ );
	} else if ( option[ OptionKeys::motifs::target_dna_defs ].user() ) {
		utility::vector1< std::string > str_def( option[ OptionKeys::motifs::target_dna_defs ]().vector() );
		load_dna_design_defs_from_strings( targeted_dna_, str_def );
	}
	dna_packer_->reference_pose( *starting_pose_ );
	// Don't need pdbs output without the dna minimization, want the pdboutput to all come from this code
	//dna_packer_->pdboutput( pdboutput_ );
}

/*void
MotifDnaPacker::minimize_dna( Pose & pose )
{

using namespace devel::blab::motif;

core::Size const nres( pose.size() );
MotifData & md( get_nonconst_motif_data( pose ) );
utility::vector1< core::Size > flexandcut = protocols::motifs::defs2vector( pose, targeted_dna_ );
sort( flexandcut.begin(), flexandcut.end() );
md.clear();
MotifData::Segment & flex_dna( md.segment("FLEX_DNA") );
MotifData::Segment & cutpoint_dna( md.segment("CUTPOINT_DNA") );
if ( flexandcut.empty() ) {
flex_dna = get_motif_data( pose ).segment( "FLEX_DNA" );
cutpoint_dna = get_motif_data( pose ).segment( "CUTPOINT_DNA" );
} else {
Size fac1( flexandcut[1] - 1 );
Size fac2( flexandcut[flexandcut.size()] + 1 );
if ( fac1 != 0 ) {
if ( (! pose.residue(fac1).is_terminus()) && pose.residue(fac1).is_DNA() ) {
flex_dna.push_back(fac1);
cutpoint_dna.push_back(fac1);
}
}
for ( Size i(1); i<=flexandcut.size(); ++i ) {
if ( ! pose.residue(flexandcut[i]).is_terminus() ) {
flex_dna.push_back(flexandcut[i]);
cutpoint_dna.push_back(flexandcut[i]);
}
}
if( ! (fac2>pose.size() ) ) {
if ( (! pose.residue( fac2 ).is_terminus()) && pose.residue( fac2 ).is_DNA() ) {
flex_dna.push_back( fac2 );
cutpoint_dna.push_back( fac2 );
}
}
}

// figure out interface positions
// new way to figure out interface in order to avoid having to include /apps/phil/interface.hh
utility::vector1< bool > is_interface_protein( nres, false );
core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory( *(dna_packer_->task_factory()) );
core::pack::task::PackerTaskOP ptask_ = tf->create_task_and_apply_taskoperations( pose );
for( core::Size j(1); j<=ptask_->total_residue(); ++j ) {
if ( !pose.residue_type(j).is_protein() ) continue;
ResidueLevelTask const & rtask( ptask_->residue_task(j) );
if( rtask.being_designed() || rtask.being_packed() ) {
is_interface_protein[j] = true;
}
}

utility::vector1< bool > is_flex_dna( nres, false );
for ( MotifData::Segment::const_iterator pos = flex_dna.begin(); pos != flex_dna.end(); ++pos ) {
is_flex_dna[ *pos ] = is_flex_dna[ scoring::dna::retrieve_base_partner_from_pose( pose )[ *pos ] ] = true;
}

/// choose a dna "loop", which just means choosing a random cutpoint in the dna, from the CUTPOINT_DNA segment
protocols::loops::Loop const dna_loop( devel::blab::motif::setup_dna_loop_with_random_cutpoint( pose ) );
protocols::loops::Loop const protein_loop;*/ // an empty loop
//setup_loop_rebuild_fold_tree( pose, protein_loop, dna_loop, false /*anchor_dna_jumps_in_backbone*/,
//               ( numeric::random::uniform() < 0.5 ) );
/* setup_loop_rebuild_cutpoint_variants( pose, protein_loop );
// which positions should be repacked/designed ?
utility::vector1< Size > const designable_positions( get_motif_data( pose ).segment("DESIGN_PROTEIN", true ) );
utility::vector1< bool > is_mutable( protocols::motifs::bools_from_sizes( nres, designable_positions ) ),
is_chimin_flexible( nres, false ), is_chimove_flexible( nres, false ), is_bb_flexible( nres, false ),
is_jump_flexible( nres, false );

for  ( core::Size i=1; i<= nres; ++i ) {
is_chimin_flexible[i] = is_flex_dna[i] || is_interface_protein[i];
is_chimove_flexible[i] = is_interface_protein[i] || is_mutable[i];
is_jump_flexible[i] = is_flex_dna[i];
is_bb_flexible[i] = is_flex_dna[i];
}

utility::vector1< core::Real > aa_bias;
core::Size outer_cycles( 10 );
bool const ramp_repulsive( false ), skip_mutation_moves_in_first_cycle( false ), dry_run( false );
devel::blab::opte::motifpacker_dna_minimize( *scorefxn_, aa_bias, is_bb_flexible, is_chimin_flexible,
is_chimove_flexible, is_mutable, is_jump_flexible,
outer_cycles, dry_run, ramp_repulsive,
skip_mutation_moves_in_first_cycle,*/
//pose, 0, false, flex_dna_sugar_/*update_pose_after_sequence_change_mover=0, vary_omega=false, flex_dna_sugar=false*/);
//}

void
MotifDnaPacker::init_options()
{
	if ( option[ OptionKeys::motifs::run_motifs ].user() ) run_motifs_ = true;
	if ( option[ OptionKeys::motifs::expand_motifs ].user() ) expand_motifs_ = true;
	if ( option[ OptionKeys::motifs::aromatic_motifs ].user() ) aromatic_motifs_ = true;
	if ( option[ OptionKeys::motifs::special_rotweight ].user() ) special_rotweight_ = option[ OptionKeys::motifs::special_rotweight ]();
	if ( option[ OptionKeys::motifs::num_repacks ].user() ) num_repacks_ = option[ OptionKeys::motifs::num_repacks ]();
	if ( option[ OptionKeys::motifs::minimize_dna ].user() ) minimize_dna_ = true;
	if ( option[ OptionKeys::motifs::flex_sugar ].user() ) flex_dna_sugar_ = true;
	//if( option[ OptionKeys::motifs::quick_and_dirty ].user() );
	if ( option[ OptionKeys::motifs::target_dna_defs ].user() ) dna_design_ = true;
}

void
MotifDnaPacker::run_motifs(
	Pose & pose,
	utility::vector1< core::Size > & design_positions,
	std::set< core::Size > & src_pos,
	std::map< core::Size, pack::rotamer_set::Rotamers > & rotamer_map,
	std::map< core::Size, std::set< std::string > > & types_map,
	std::list< std::string > & info_lines,
	pack::task::TaskFactoryOP taskfactory
)
{
	for ( core::Size i(1); i <= design_positions.size(); ++i ) {
		std::set< std::string > name3set;
		pack::rotamer_set::Rotamers motif_rotamers( motif_search_->bp_rotamers( design_positions[i] ) );
		src_pos.insert( design_positions[i] );
		pack::rotamer_set::Rotamers variant_rotamers;
		if ( ! motif_rotamers.empty() ) {
			for ( core::Size rot(1); rot <= motif_rotamers.size(); ++rot ) {
				conformation::ResidueOP variant_rot( pose::add_variant_type_to_residue( *(motif_rotamers[rot]), core::chemical::SPECIAL_ROT, pose ) );
				variant_rotamers.push_back( variant_rot );
				name3set.insert( (motif_rotamers[rot])->name3() );
			}
			rotamer_map[design_positions[i]] = variant_rotamers;
			if ( (! expand_motifs_) && (! aromatic_motifs_) ) {
				protocols::toolbox::rotamer_set_operations::SpecialRotamerRSOOP ms_rsoop( new protocols::toolbox::rotamer_set_operations::SpecialRotamerRSO( design_positions[i] ) );
				ms_rsoop->set_new_rots( variant_rotamers );
				taskfactory->push_back( TaskOperationCOP( new AppendRotamerSet( ms_rsoop ) ) );
			}
		}
		types_map[design_positions[i]] = name3set;
	}
	std::string tag("");
	if ( (! expand_motifs_) && (! aromatic_motifs_) ) {
		std::stringstream mot_name;
		bool special_rotweight_zero( true );
		if ( special_rotweight_zero ) {
			core::Real zero_special_rotweight( -0.0000000001 );
			scorefxn_->set_weight( special_rot, ( zero_special_rotweight ) );
			mot_name << "ZPW";
			std::string mot_name3( mot_name.str() );
			std::stringstream filename;
			filename << filename_root_ << "_" << mot_name3;
			dna_packer_->task_factory( taskfactory );
			dna_packer_->set_filename_root( filename.str() );
			dna_packer_->apply( pose );
			/*if ( minimize_dna_ ) {
			minimize_dna( pose );
			tag = "_min";
			}*/
			bool const overwrite_old_info(true);
			pdboutput_->add_info( "REMARK DESIGNED POSITIONS " + filename_root_ + ":", info_lines, !overwrite_old_info );
			pdboutput_->score_function( *scorefxn_ );
			(*pdboutput_)(pose, dna_packer_->pdbname() + tag + ".pdb");

			dna_packer_->clear_initialization();
		}
		pose = *starting_pose_;
		core::Real special_rotweight( special_rotweight_ );
		for ( core::Size trial(0); trial < num_repacks_; ++trial ) {
			core::Real special_rotweight2 = ( special_rotweight / 2 );
			scorefxn_->set_weight( special_rot, ( special_rotweight2 /*(trial+1)*/ ) );
			special_rotweight = special_rotweight2;
			std::stringstream filename;
			filename << filename_root_ << "_" << lead_zero_string_of( trial, 4 );
			dna_packer_->task_factory( taskfactory );
			dna_packer_->set_filename_root( filename.str() );
			dna_packer_->apply( pose );
			/*if ( minimize_dna_ ) {
			minimize_dna( pose );
			}*/
			bool const overwrite_old_info(true);
			pdboutput_->add_info( "REMARK DESIGNED POSITIONS " + filename_root_ + ":", info_lines, !overwrite_old_info );
			pdboutput_->score_function( *scorefxn_ );
			(*pdboutput_)(pose, dna_packer_->pdbname() + tag + ".pdb");

			dna_packer_->clear_initialization();
		}
	}
}

void
MotifDnaPacker::expand_motifs(
	Pose & pose,
	//mjo commenting out 'design_positions' because it is unused and causes a warning
	utility::vector1< core::Size > & /*design_positions*/,
	std::set< core::Size > & src_pos,
	std::map< core::Size, pack::rotamer_set::Rotamers > & rotamer_map,
	std::map< core::Size, std::set< std::string > > & types_map,
	std::list< std::string > & info_lines,
	pack::task::TaskFactoryOP taskfactory
)
{
	for ( std::map< core::Size, core::pack::rotamer_set::Rotamers >::const_iterator it( rotamer_map.begin() ),
			end( rotamer_map.end() ); it != end; ++it ) {
		std::set< std::string > name3s( types_map[it->first] );
		std::set< core::Size > current_pos( src_pos );
		current_pos.erase( it->first );
		for ( auto it2( name3s.begin() ),
				end2( name3s.end() ); it2 != end2; ++it2 ) {
			motif_expansion_inner_loop( pose, current_pos, it, it2, rotamer_map, info_lines, taskfactory );
		}
	}
}

void
MotifDnaPacker::motif_expansion_inner_loop(
	core::pose::Pose & pose,
	std::set< core::Size > current_pos,
	std::map< core::Size, core::pack::rotamer_set::Rotamers >::const_iterator it,
	std::set< std::string >::const_iterator it2,
	std::map< core::Size, core::pack::rotamer_set::Rotamers > & rotamer_map,
	std::list< std::string > & info_lines,
	pack::task::TaskFactoryOP taskfactory
) {
	std::string tag("");
	pose = *starting_pose_;
	TaskFactoryOP my_tf2;
	my_tf2 = TaskFactoryOP( new TaskFactory( *taskfactory ) );
	pack::rotamer_set::Rotamers restricted_rotamers;
	pack::rotamer_set::Rotamers src_rotamers( it->second );
	for ( core::Size rot2(1); rot2 <= src_rotamers.size(); ++rot2 ) {
		if ( src_rotamers[rot2]->name3() != *it2 ) continue;
		restricted_rotamers.push_back( src_rotamers[rot2] );
	}
	utility::vector1< bool > keep_aas( num_canonical_aas, false );
	keep_aas[ chemical::aa_from_name(*it2) ] = true;
	my_tf2->push_back( TaskOperationCOP( new RestrictAbsentCanonicalAAS( it->first, keep_aas ) ) );
	protocols::toolbox::rotamer_set_operations::SpecialRotamerRSOOP ms_rsoop( new protocols::toolbox::rotamer_set_operations::SpecialRotamerRSO( it->first ) );
	ms_rsoop->set_new_rots( restricted_rotamers );
	my_tf2->push_back( TaskOperationCOP( new AppendRotamerSet( ms_rsoop ) ) );

	for ( unsigned long current_po : current_pos ) {
		if ( rotamer_map[current_po].empty() ) continue;
		pack::rotamer_set::Rotamers other_rotamers;
		for ( core::Size it4(1); it4<=(rotamer_map[current_po]).size(); ++it4 ) {
			other_rotamers.push_back( rotamer_map[current_po][it4] );
		}
		protocols::toolbox::rotamer_set_operations::SpecialRotamerRSOOP ms_rsoop2( new protocols::toolbox::rotamer_set_operations::SpecialRotamerRSO( current_po ) );
		ms_rsoop2->set_new_rots( other_rotamers );
		my_tf2->push_back( TaskOperationCOP( new AppendRotamerSet( ms_rsoop2 ) ) );
	}

	std::stringstream mot_name;
	mot_name << *it2 << "_" << pose.pdb_info()->chain(it->first) << pose.pdb_info()->number(it->first);
	std::string mot_name2( mot_name.str() );

	bool special_rotweight_zero( true );
	if ( special_rotweight_zero ) {
		core::Real zero_special_rotweight( -0.0000000001 );
		scorefxn_->set_weight( special_rot, ( zero_special_rotweight ) );
		mot_name << "_ZPW";
		std::string mot_name3( mot_name.str() );
		std::stringstream filename;
		filename << filename_root_ << "_" << mot_name3;
		dna_packer_->task_factory( taskfactory );
		dna_packer_->set_filename_root( filename.str() );
		dna_packer_->apply( pose );
		/*if ( minimize_dna_ ) {
		minimize_dna( pose );
		tag = "_min";
		}*/
		bool const overwrite_old_info(true);
		pdboutput_->add_info( "REMARK DESIGNED POSITIONS " + filename_root_ + ":", info_lines, !overwrite_old_info );
		pdboutput_->score_function( *scorefxn_ );
		(*pdboutput_)(pose, dna_packer_->pdbname() + tag + ".pdb");

		dna_packer_->clear_initialization();
	}
	pose = *starting_pose_;
	core::Real special_rotweight( special_rotweight_ );
	for ( core::Size trial(0); trial < num_repacks_; ++trial ) {
		core::Real special_rotweight2 = ( special_rotweight / 2 );
		scorefxn_->set_weight( special_rot, ( special_rotweight2 /*(trial+1)*/ ) );
		special_rotweight = special_rotweight2;
		std::stringstream filename;
		filename << filename_root_ << "_" << mot_name2 << "_" << lead_zero_string_of( trial, 4 );
		dna_packer_->task_factory( my_tf2 );
		dna_packer_->set_filename_root( filename.str() );
		dna_packer_->apply( pose );
		/*if ( minimize_dna_ ) {
		minimize_dna( pose );
		}*/

		bool const overwrite_old_info(true);
		pdboutput_->add_info( "REMARK DESIGNED POSITIONS " + filename_root_ + ":", info_lines, !overwrite_old_info );
		pdboutput_->score_function( *scorefxn_ );
		(*pdboutput_)(pose, dna_packer_->pdbname() + tag + ".pdb");

		dna_packer_->clear_initialization();
	}
}

void
MotifDnaPacker::aromatic_motifs(
	Pose & pose,
	utility::vector1< core::Size > & design_positions,
	std::set< core::Size > & src_pos,
	std::map< core::Size, pack::rotamer_set::Rotamers > & rotamer_map,
	std::map< core::Size, std::set< std::string > > & types_map,
	std::list< std::string > & info_lines,
	pack::task::TaskFactoryOP taskfactory
)
{
	for ( std::map< core::Size, pack::rotamer_set::Rotamers >::const_iterator it( rotamer_map.begin() ),
			end( rotamer_map.end() ); it != end; ++it ) {
		std::set< std::string > name3s( types_map[it->first] );
		std::set< core::Size > current_pos( src_pos );
		current_pos.erase( it->first );
		bool aromatics( false );
		std::set< std::string > name3_arom;
		 for ( auto const & name3 : name3s ) {
			if ( name3 == "TYR" || name3 == "PHE" || name3 == "TRP" ) {
				aromatics = true;
				name3_arom.insert( name3 );
			}
		}
		if ( !aromatics ) {
			continue;
		}
		for ( auto it2( name3_arom.begin() ),
				end2( name3_arom.end() ); it2 != end2; ++it2 ) {
			motif_expansion_inner_loop( pose, current_pos, it, it2, rotamer_map, info_lines, taskfactory );
		}
	}
	aromatic_motifs_ = false;
	run_motifs( pose, design_positions, src_pos, rotamer_map, types_map, info_lines, taskfactory );
}

} // namespace motifs
} // namespace protocols
