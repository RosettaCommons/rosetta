// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/AntibodyDesignMover.cc
/// @brief Class that initially designs antibodies through grafting using an AntibodyDatabase + North_AHO numbering scheme
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/antibody/design/AntibodyDesignMover.hh>
#include <protocols/antibody/design/AntibodyDesignMoverCreator.hh>

#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/clusters/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/CDRDihedralConstraintMover.hh>

#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/GeneralAntibodyModeler.hh>
#include <protocols/antibody/design/AntibodySeqDesignTFCreator.hh>
#include <protocols/antibody/design/MutateFrameworkForCluster.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/design/AntibodyDesignEnumManager.hh>
#include <protocols/antibody/design/NativeAntibodySeq.hh>

// Core Includes
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.hh>

#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>

// Protocol Includes
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/util.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>

// Numeric Includes
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

//Utility
#include <math.h>
#include <basic/Tracer.hh>
#include <utility/py/PyAssert.hh>
#include <map>
#include <iterator>
#include <boost/algorithm/string.hpp>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.design.AntibodyDesignMover");

namespace protocols {
namespace antibody {
namespace design {
using namespace protocols::antibody;
using namespace protocols::grafting;
using namespace protocols::antibody::clusters;
using namespace protocols::antibody::constraints;
using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace core::scoring;
using basic::datacache::DataCache_CacheableData;
using core::Size;
using core::pose::Pose;
using core::pose::PoseOP;
using utility::vector1;

AntibodyDesignMover::AntibodyDesignMover():
	ab_info_(/* NULL */),
	graft_mover_(/* NULL */),
	scorefxn_(/* NULL */),
	scorefxn_min_(/* NULL */),
	paratope_epitope_cst_mover_(/* NULL */),
	paratope_cst_mover_(/* NULL */),
	cdr_dihedral_cst_mover_ (/* NULL */)

{
	design_enum_manager_ = AntibodyDesignEnumManagerOP( new AntibodyDesignEnumManager);
	read_command_line_options();
	set_defaults();

}

AntibodyDesignMover::AntibodyDesignMover( AntibodyInfoCOP ab_info ):
	graft_mover_(/* NULL */),
	scorefxn_(/* NULL */),
	scorefxn_min_(/* NULL */),
	paratope_epitope_cst_mover_(/* NULL */),
	paratope_cst_mover_(/* NULL */),
	cdr_dihedral_cst_mover_ (/* NULL */)

{
	ab_info_ = ab_info->clone();

	design_enum_manager_ = AntibodyDesignEnumManagerOP( new AntibodyDesignEnumManager);
	read_command_line_options();
	set_defaults();

}

AntibodyDesignMover::~AntibodyDesignMover(){}

void
AntibodyDesignMover::set_defaults(){
	///Conservative defaults.  Defaults here are also read and set from a database file.  To allow run-time manipulations and testing.
	overhang_ = 3;
	print_tracer_info_ = true;
	paratope_cdrs_.clear();
	paratope_cdrs_.resize(6, true);
	epitope_residues_.clear();
	cdr_set_options_.clear();
	cdr_graft_design_options_.clear();
	cdr_seq_design_options_.clear();
}

void
AntibodyDesignMover::read_command_line_options(){
	using namespace basic::options;
	AntibodyEnumManager manager = AntibodyEnumManager();

	outer_cycles_ = (option [OptionKeys::antibody::design::outer_cycle_rounds]());
	inner_cycles_ = option [OptionKeys::antibody::design::inner_cycle_rounds]();

	set_keep_top_designs(option [OptionKeys::antibody::design::top_designs]());
	set_dock_post_graft(option [OptionKeys::antibody::design::do_dock]());
	set_rb_min_post_graft(option [OptionKeys::antibody::design::do_rb_min]());
	set_dock_rounds(option [OptionKeys::antibody::design::dock_cycle_rounds]());
	//initial_perturb_ = basic::options::option [basic::options::OptionKeys::antibody::design::initial_perturb] ();
	benchmark_ = option [OptionKeys::antibody::design::random_start]();
	use_light_chain_type_ = option [OptionKeys::antibody::design::use_light_chain_type]();
	adapt_graft_ = option [OptionKeys::antibody::design::adapt_graft]();
	enable_adapt_graft_cartesian_ = option[ OptionKeys::antibody::design::enable_adapt_graft_cartesian]();

	idealize_graft_cdrs_ = option [OptionKeys::antibody::design::idealize_graft_cdrs]();
	add_log_to_pose_ = option [OptionKeys::antibody::design::add_graft_log_to_pdb]();

	outer_kt_ = option [OptionKeys::antibody::design::outer_kt]();
	inner_kt_ = option [OptionKeys::antibody::design::inner_kt]();

	//TR << "Design protocol: " << option[ OptionKeys::antibody::design::design_protocol]() << std::endl;
	design_protocol_ = design_enum_manager_->design_protocol_string_to_enum(option[ OptionKeys::antibody::design::design_protocol]());

	interface_dis_ = option [OptionKeys::antibody::design::interface_dis]();
	neighbor_dis_ = option [OptionKeys::antibody::design::neighbor_dis]();

	dock_min_dock_ = option[OptionKeys::antibody::design::dock_min_dock]();

	mutate_framework_for_cluster_ = option [OptionKeys::antibody::design::mutate_framework_for_cluster]();

	if ( basic::options::option [basic::options::OptionKeys::antibody::design::cdr_instructions].user() ) {
		instruction_file_ = basic::options::option [basic::options::OptionKeys::antibody::design::cdr_instructions]();
		TR << "Instructions file: " << instruction_file_ << std::endl;
	}

	//Epitope Constraint options
	use_epitope_constraints_ = option [OptionKeys::antibody::design::use_epitope_constraints]();
	if ( use_epitope_constraints_ && option [OptionKeys::antibody::design::epitope].user() ) {
		vector1<std::string> epitope_res_strings = option [ OptionKeys::antibody::design::epitope]();
		epitope_residues_ = protocols::antibody::design::get_pdb_numbering_from_strings(epitope_res_strings);
	}

	//Paratope Constraint options
	if ( option [OptionKeys::antibody::design::paratope].user() ) {
		paratope_cdrs_.clear();
		paratope_cdrs_.resize(6, false);
		vector1<std::string> cdrs = option [OptionKeys::antibody::design::paratope]();

		for ( core::Size i = 1; i <= cdrs.size(); ++i ) {
			CDRNameEnum cdr = manager.cdr_name_string_to_enum(cdrs[i]);
			paratope_cdrs_[cdr] = true;
		}
	} else {
		paratope_cdrs_.clear();
		paratope_cdrs_.resize(6, true);
	}

	//Design CDR Overrides
	if ( option [OptionKeys::antibody::design::design_cdrs].user() ) {
		utility::vector1<std::string> cdrs = option [OptionKeys::antibody::design::design_cdrs]();
		for ( core::Size i = 1; i <= cdrs.size(); ++i ) {
			CDRNameEnum cdr_enum = manager.cdr_name_string_to_enum(cdrs[i]);
			design_override_.push_back(cdr_enum);
		}
	}
}

protocols::moves::MoverOP
AntibodyDesignMover::clone() const {
	return protocols::moves::MoverOP( new AntibodyDesignMover(*this) );
}

protocols::moves::MoverOP
AntibodyDesignMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AntibodyDesignMover );
}

AntibodyDesignMover::AntibodyDesignMover( AntibodyDesignMover const & src ):
	protocols::moves::Mover( src ),
	cdr_set_options_( src.cdr_set_options_ ),
	cdr_graft_design_options_( src.cdr_graft_design_options_ ),
	cdr_seq_design_options_( src.cdr_seq_design_options_ ),
	cdr_set_( src.cdr_set_ ),
	cluster_based_CDRDBPose_indexes_( src.cluster_based_CDRDBPose_indexes_ ),
	length_based_CDRDBPose_indexes_( src.length_based_CDRDBPose_indexes_ ),
	overhang_(src.overhang_),
	outer_cycles_(src.outer_cycles_),
	inner_cycles_(src.inner_cycles_),
	dock_cycles_(src.dock_cycles_),
	num_top_designs_(src.num_top_designs_),
	interface_dis_(src.interface_dis_),
	neighbor_dis_(src.neighbor_dis_),
	top_scores_(src.top_scores_),
	total_permutations_(src.total_permutations_),
	dock_post_graft_(src.dock_post_graft_),
	rb_min_post_graft_(src.rb_min_post_graft_),
	paratope_cdrs_(src.paratope_cdrs_),
	epitope_residues_(src.epitope_residues_),
	adapt_graft_( src.adapt_graft_ ),
	enable_adapt_graft_cartesian_( src.enable_adapt_graft_cartesian_ ),
	benchmark_( src.benchmark_ ),
	use_light_chain_type_( src.use_light_chain_type_ ),
	use_epitope_constraints_( src.use_epitope_constraints_),
	print_tracer_info_( src.print_tracer_info_ ),
	idealize_graft_cdrs_( src.idealize_graft_cdrs_),
	add_log_to_pose_( src.add_log_to_pose_),
	graft_log_( src.graft_log_ ),
	accept_log_(src.accept_log_ ),
	outer_kt_(src.outer_kt_),
	inner_kt_(src.inner_kt_),
	enable_full_protocol_atom_pair_cst_(src.enable_full_protocol_atom_pair_cst_ ),
	design_protocol_( src.design_protocol_ ),
	design_override_(src.design_override_),
	cdrs_to_design_(src.cdrs_to_design_),
	instruction_file_(src.instruction_file_),
	dock_min_dock_(src.dock_min_dock_),
	stats_cutoff_( src.stats_cutoff_),
	mutate_framework_for_cluster_(src.mutate_framework_for_cluster_)
{
	using namespace protocols::grafting;
	using namespace protocols::simple_moves;
	
	if (src.scorefxn_) scorefxn_ = src.scorefxn_->clone();
	if (src.scorefxn_min_) scorefxn_min_ = src.scorefxn_min_->clone();
	if (src.scorefxn_cart_graft_) scorefxn_cart_graft_ = src.scorefxn_cart_graft_->clone();
	top_designs_.clear();
	for ( core::pose::PoseOP p : src.top_designs_){
		
		core::pose::PoseOP new_p = p->clone();
		top_designs_.push_back( new_p );
	}
	
	if ( src.ab_info_ ) ab_info_ = AntibodyInfoOP( new AntibodyInfo( *src.ab_info_));
	if ( src.design_enum_manager_ ) design_enum_manager_ = AntibodyDesignEnumManagerOP( new AntibodyDesignEnumManager( *src.design_enum_manager_ ));
	if ( src.seq_design_creator_ ) seq_design_creator_ = AntibodySeqDesignTFCreatorOP( new AntibodySeqDesignTFCreator( *src.seq_design_creator_));
	if ( src.graft_mover_ ) graft_mover_ = CCDEndsGraftMoverOP( new CCDEndsGraftMover( *src.graft_mover_));
	if ( src.anchored_graft_mover_) anchored_graft_mover_ = AnchoredGraftMoverOP( new AnchoredGraftMover( *src.anchored_graft_mover_));
	if ( src.framework_mutator_) MutateFrameworkForClusterOP( new MutateFrameworkForCluster( *src.framework_mutator_));
	if ( src.modeler_) modeler_ = GeneralAntibodyModelerOP( new GeneralAntibodyModeler( *src.modeler_ ));
	if ( src.cart_min_graft_) cart_min_graft_ = MinMoverOP( new MinMover( *src.cart_min_graft_ ));
	if ( src.mc_ ) mc_ = src.mc_->clone();
	if ( src.paratope_epitope_cst_mover_ ) paratope_epitope_cst_mover_ = ParatopeEpitopeSiteConstraintMoverOP( new constraints::ParatopeEpitopeSiteConstraintMover( *src.paratope_epitope_cst_mover_));
	if ( src.paratope_cst_mover_ ) paratope_cst_mover_ = ParatopeSiteConstraintMoverOP( new ParatopeSiteConstraintMover( *src.paratope_cst_mover_));
	if ( src.cdr_dihedral_cst_mover_ ) cdr_dihedral_cst_mover_ = CDRDihedralConstraintMoverOP( new CDRDihedralConstraintMover( *src.cdr_dihedral_cst_mover_));
	
}



void
AntibodyDesignMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const & ,
	moves::Movers_map const & ,
	Pose const &
){
	using namespace core::scoring;
	AntibodyEnumManager manager = AntibodyEnumManager();

	scorefxn_ = get_ab_design_global_scorefxn(tag, data);
	scorefxn_min_ = get_ab_design_min_scorefxn(tag, data);

	if ( tag->hasOption("design_cdrs") ) {
		utility::vector1<std::string> cdr_strings = utility::string_split_multi_delim(tag->getOption< std::string>("design_cdrs"), ":,'`~+*&|;.");
		for ( core::Size i = 1; i <= cdr_strings.size(); ++i ) {
			CDRNameEnum cdr_enum =manager.cdr_name_string_to_enum( cdr_strings[ i ] );
			design_override_.push_back(cdr_enum);
		}
	}

	//A little redundancy
	if ( tag->hasOption("instruction_file") ) {
		instruction_file_ = tag->getOption< std::string >("instruction_file");
	} else if ( tag->hasOption("instructions_file") ) {
		instruction_file_ = tag->getOption< std::string >("instructions_file");
	} else if ( tag->hasOption("cdr_instructions_file") ) {
		instruction_file_ = tag->getOption< std::string >("cdr_instructions_file");
	}


	design_protocol_ = design_enum_manager_->design_protocol_string_to_enum(tag->getOption< std::string >("design_protocol", "generalized_monte_carlo"));

	interface_dis_ = tag->getOption< core::Real >("interface_dis", interface_dis_);
	neighbor_dis_ = tag->getOption< core::Real >("neighbor_dis", neighbor_dis_);

	outer_cycles_ = tag->getOption< core::Size >("outer_cycles", outer_cycles_);
	inner_cycles_ = tag->getOption< core::Size >("inner_cycles", inner_cycles_);

	outer_kt_ = tag->getOption< core::Real >("outer_kt", outer_kt_);
	inner_kt_ = tag->getOption< core::Real >("inner_kt", inner_kt_);

	num_top_designs_ = tag->getOption< core::Size >("top_designs", num_top_designs_);

	dock_post_graft_ = tag->getOption< bool >("do_dock", dock_post_graft_);
	rb_min_post_graft_ = tag->getOption< bool >("do_rb_min", rb_min_post_graft_);
	dock_cycles_ = tag->getOption< core::Size >("dock_cycles", dock_cycles_);
	dock_min_dock_ = tag->getOption<bool>("dock_min_dock", dock_min_dock_);

	//initial_perturb_ = basic::options::option [basic::options::OptionKeys::antibody::design::initial_perturb] ();
	benchmark_ = tag->getOption< bool >("benchmark", benchmark_);

	use_light_chain_type_ = tag->getOption< bool >("use_light_chain_type", use_light_chain_type_);

	adapt_graft_ = tag->getOption< bool >("adapt_graft", adapt_graft_);
	enable_adapt_graft_cartesian_ = tag->getOption< bool >("enable_adapt_graft_cartesian", enable_adapt_graft_cartesian_);

	idealize_graft_cdrs_ = tag->getOption< bool >("idealize_graft_cdrs", idealize_graft_cdrs_);
	add_log_to_pose_ = tag->getOption< bool >("add_graft_log_to_pdb", add_log_to_pose_);

	mutate_framework_for_cluster_ = tag->getOption< bool >("mutate_framework_for_cluster", mutate_framework_for_cluster_);

	//Epitope Constraints
	tag->getOption<bool>("use_epitope_csts", use_epitope_constraints_);
	if ( use_epitope_constraints_ && tag->hasOption("epitope_residues") ) {
		vector1<std::string> epitope_res_strings = utility::string_split_multi_delim(tag->getOption<std::string>("epitope_residues"), ",+;");
		epitope_residues_ = protocols::antibody::design::get_pdb_numbering_from_strings(epitope_res_strings);
	}

	//Paratope Constraint options
	if ( tag->hasOption("paratope_cdrs") ) {
		paratope_cdrs_ = get_cdr_bool_from_tag(tag, "paratope_cdrs");
	} else {
		paratope_cdrs_.clear();
		paratope_cdrs_.resize(6, true);
	}
}


void
AntibodyDesignMover::set_scorefunction(ScoreFunctionOP scorefxn){
	scorefxn_=scorefxn->clone();
}

void
AntibodyDesignMover::set_scorefunction_min(core::scoring::ScoreFunctionOP min_scorefxn) {
	scorefxn_min_ = min_scorefxn->clone();
}

void
AntibodyDesignMover::set_keep_top_designs(core::Size top_designs){
	num_top_designs_ = top_designs;
}

void
AntibodyDesignMover::set_dock_post_graft(bool dock_post_graft){
	dock_post_graft_ = dock_post_graft;
}

void
AntibodyDesignMover::set_dock_rounds(core::Size dock_rounds){
	dock_cycles_ = dock_rounds;
}

void
AntibodyDesignMover::set_rb_min_post_graft(bool rb_min_post_graft){
	rb_min_post_graft_ = rb_min_post_graft;
}

void
AntibodyDesignMover::set_outer_cycles(core::Size graft_rounds){
	outer_cycles_=graft_rounds;
}

void
AntibodyDesignMover::set_inner_cycles(core::Size cycles) {
	inner_cycles_ = cycles;
}

void
AntibodyDesignMover::set_interface_dis(core::Real interface_dis){
	interface_dis_ = interface_dis;
}

void
AntibodyDesignMover::set_neighbor_dis(core::Real neighbor_dis){
	neighbor_dis_ = neighbor_dis;
}

void
AntibodyDesignMover::set_design_protocol(AntibodyDesignProtocolEnum design_protocol){
	design_protocol_ = design_protocol;
}

void
AntibodyDesignMover::set_use_epitope_constraints(bool use_epitope_csts){
	use_epitope_constraints_ = use_epitope_csts;
}

void
AntibodyDesignMover::set_paratope_cdrs(const vector1<bool>& cdrs) {
	if ( cdrs.size() != 6 ) {
		utility_exit_with_message("Cannot setup paratope cdrs - number passed does not equal total number of cdrs!");
	}
	paratope_cdrs_ = cdrs;

}

void
AntibodyDesignMover::set_epitope_residues(vector1<PDBNumbering> epitope_residues) {
	epitope_residues_ = epitope_residues;
}

void
AntibodyDesignMover::setup_native_clusters(core::pose::Pose & pose){
	ab_info_->setup_CDR_clusters(pose);
	ab_info_->get_CDR_cluster_set()->set_cacheable_cluster_data_to_pose(pose);
}

void
AntibodyDesignMover::setup_native_sequence(core::pose::Pose & pose){
	pose.data().set(core::pose::datacache::CacheableDataType::NATIVE_ANTIBODY_SEQ, DataCache_CacheableData::DataOP( new NativeAntibodySeq( pose, ab_info_) ));
}
//void
//AntibodyDesignMover::set_cdr_set(CDRDBPoseSet& cdr_set, core::Size overhang){
// cdr_set_ = cdr_set;
// overhang_ = overhang;
//}

void
AntibodyDesignMover::set_cdr_override(utility::vector1<CDRNameEnum> cdrs_to_design){
	design_override_ = cdrs_to_design;
}

void
AntibodyDesignMover::set_instruction_file(std::string instruction_file){
	instruction_file_ = instruction_file;
}

void
AntibodyDesignMover::initialize_cdr_set(core::pose::Pose const & pose){
	AntibodyDatabaseManager manager = AntibodyDatabaseManager(ab_info_);
	manager.use_light_chain_type(use_light_chain_type_);
	cdr_set_ = manager.load_cdr_poses(cdr_set_options_, pose);

}

void
AntibodyDesignMover::set_cdr_set_options(AntibodyCDRSetOptions cdr_set_options){
	cdr_set_options_ = cdr_set_options;
}

void
AntibodyDesignMover::set_graft_design_options(AntibodyCDRGraftDesignOptions graft_design_options){
	cdr_graft_design_options_ = graft_design_options;
}

void
AntibodyDesignMover::set_seq_design_options(AntibodyCDRSeqDesignOptions seq_design_options){
	cdr_seq_design_options_ = seq_design_options;
}

void
AntibodyDesignMover::setup_options_classes(){

	if ( cdr_set_options_.size() == 0 ) {
		cdr_set_options_  = protocols::antibody::design::get_cdr_set_options(instruction_file_);
	}
	if ( cdr_graft_design_options_.size() == 0 ) {
		cdr_graft_design_options_ = protocols::antibody::design::get_graft_design_options(instruction_file_);
	}
	if ( cdr_seq_design_options_.size() == 0 ) {
		cdr_seq_design_options_ = protocols::antibody::design::get_seq_design_options(instruction_file_);
	}

	if ( design_override_.size() > 0 ) {
		for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
			cdr_set_options_[i]->load(false);
			cdr_graft_design_options_[i]->design(false);
			cdr_seq_design_options_[i]->design(false);
		}
		for ( core::Size i = 1; i <= design_override_.size(); ++i ) {
			CDRNameEnum cdr_enum = design_override_[ i ];

			//Design override only works for l4 and h4 in regard to sequence design.
			if ( cdr_enum != l4 && cdr_enum != h4 ) {
				cdr_set_options_[cdr_enum]->load(true);
				cdr_graft_design_options_[cdr_enum]->design(true);
			}

			cdr_seq_design_options_[cdr_enum]->design(true);
		}
	}

	//Disable CDRs that are not present (aka camelid design)
	for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		if ( ! ab_info_->has_CDR( cdr ) ) {
			cdr_seq_design_options_[ cdr ]->design( false );
			cdr_graft_design_options_[ cdr ]->design( false );
		}
	}

	///Set it up so that if we design the cdr, we load CDRs from the database.  Only canonical for now.
	for ( core::Size i = 1; i <= 6; ++i ) {
		CDRGraftDesignOptionsOP options = cdr_graft_design_options_[i];
		cdr_set_options_[i]->load(options->design());
	}



}
void
AntibodyDesignMover::setup_cdr_pose_sampling_strategies() {
	typedef std::map< CDRNameEnum, utility::vector1< CDRDBPose > >::const_iterator it_type;

	for ( it_type it = cdr_set_.begin(); it != cdr_set_.end(); ++it ) {
		CDRNameEnum cdr = it->first;
		for ( core::Size index=1; index <= cdr_set_[ cdr ].size(); ++index ) {
			CDRClusterEnum cluster = cdr_set_[ cdr ][ index ].cluster;
			core::Size cluster_length = ab_info_->get_cluster_length(cluster);
			cluster_based_CDRDBPose_indexes_[ cdr ][ cluster ].push_back( index );
			length_based_CDRDBPose_indexes_[ cdr ][ cluster_length ][ cluster ].push_back( index );
		}
	}

}

void
AntibodyDesignMover::setup_default_graft_settings(){

	//CCDEndsGraftMover
	graft_mover_->set_cycles(50);
	graft_mover_->set_scaffold_flexibility(2, 2);
	graft_mover_->set_insert_flexibility(2, 2);
	graft_mover_->final_repack(true);
	graft_mover_->stop_at_closure(true);
	graft_mover_->neighbor_dis(neighbor_dis_);
	graft_mover_->set_fa_scorefunction(scorefxn_);
	graft_mover_->idealize_insert(idealize_graft_cdrs_);
	graft_mover_->copy_pdbinfo(true);
	//graft_mover_->set_skip_sampling(true);

	//AnchoredGraftMover
	anchored_graft_mover_->set_cycles(100);
	anchored_graft_mover_->set_scaffold_flexibility(2, 2);
	anchored_graft_mover_->set_insert_flexibility(2, 2);
	anchored_graft_mover_->final_repack(true);
	anchored_graft_mover_->stop_at_closure(true);
	//anchored_graft_mover_->conservative_smallmover(true);
	anchored_graft_mover_->neighbor_dis(neighbor_dis_);
	anchored_graft_mover_->set_fa_scorefunction(scorefxn_);
	anchored_graft_mover_->idealize_insert(idealize_graft_cdrs_);
	anchored_graft_mover_->copy_pdbinfo(true);

}

void
AntibodyDesignMover::setup_cart_minimizer(){
	cart_min_graft_ = simple_moves::MinMoverOP( new simple_moves::MinMover());
	cart_min_graft_->cartesian(true);
	cart_min_graft_->min_type("lbfgs_armijo_nonmonotone");
	cart_min_graft_->min_options()->max_iter(200);
	cart_min_graft_->tolerance(.01);
	cart_min_graft_->score_function(scorefxn_cart_graft_);
}

void
AntibodyDesignMover::setup_modeler(){

	modeler_ = GeneralAntibodyModelerOP( new GeneralAntibodyModeler(ab_info_) );
	modeler_->set_scorefunction(scorefxn_); //Note that modeler will setup its own docking scorefunctions with proper constraint weights.
	modeler_->set_scorefunction_min(scorefxn_min_);
	modeler_->interface_detection_dis(interface_dis_);
	modeler_->neighbor_detection_dis(neighbor_dis_);

	std::string ab_dock_chains = "A_" +ab_info_->get_antibody_chain_string();
	modeler_->ab_dock_chains(ab_dock_chains);
}

void
AntibodyDesignMover::setup_scorefxn() {
	using namespace basic::options;

	if ( !scorefxn_ ) {
		scorefxn_ = design::get_ab_design_global_scorefxn();
	}

	if ( !scorefxn_min_ ) {
		scorefxn_min_ = design::get_ab_design_min_scorefxn();
	}

	//Set atom-pair constraint to zero so that it does not affect cart
	scorefxn_cart_graft_ = scorefxn_->clone();
	scorefxn_cart_graft_->set_weight(atom_pair_constraint, 0);
	scorefxn_cart_graft_->set_weight_if_zero(cart_bonded, .5);
	scorefxn_cart_graft_->set_weight_if_zero(dihedral_constraint, option[ OptionKeys::antibody::design::dihedral_cst_weight]());
	scorefxn_cart_graft_->set_weight(pro_close, 0);
}

void
AntibodyDesignMover::setup_epitope_residues(const Pose & pose){

	//This is so we setup epitope at the very beginning of the protocol.
	if ( epitope_residues_.size() == 0 ) {
		vector1<bool> epitope = select_epitope_residues(ab_info_, pose, interface_dis_);
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( ! epitope[i] ) continue;

			PDBNumbering numbering;
			numbering.icode = pose.pdb_info()->icode(i);
			numbering.chain = pose.pdb_info()->chain(i);
			numbering.resnum = pose.pdb_info()->number(i);
			epitope_residues_.push_back(numbering);
		}
	}
}

void
AntibodyDesignMover::setup_paratope_epitope_constraints(Pose & pose){

	//Create classes if they don't yet exist
	if ( !paratope_epitope_cst_mover_ ) {
		paratope_epitope_cst_mover_ = constraints::ParatopeEpitopeSiteConstraintMoverOP( new ParatopeEpitopeSiteConstraintMover(ab_info_) );
	}
	if ( ! paratope_cst_mover_ ) {
		paratope_cst_mover_ = constraints::ParatopeSiteConstraintMoverOP( new ParatopeSiteConstraintMover(ab_info_) );
	}
	if ( ! cdr_dihedral_cst_mover_ ) {
		cdr_dihedral_cst_mover_ = constraints::CDRDihedralConstraintMoverOP( new CDRDihedralConstraintMover(ab_info_));
	}

	//The residues need to be regenerated if size is changing.
	// An observer could be used here, but it seems to be slower then re-initing everything.
	if ( use_epitope_constraints_ ) {
		paratope_epitope_cst_mover_->set_defaults(); //Clear everything
		paratope_epitope_cst_mover_->constrain_to_paratope_cdrs(paratope_cdrs_); //Regenerates paratope residues
		paratope_epitope_cst_mover_->constrain_to_epitope_residues(epitope_residues_, pose);
		paratope_epitope_cst_mover_->set_interface_distance(interface_dis_);
		paratope_epitope_cst_mover_->apply(pose);
	} else {
		paratope_cst_mover_->set_defaults(); //Clear everything
		paratope_cst_mover_->constrain_to_paratope_cdrs(paratope_cdrs_);
		paratope_cst_mover_->set_interface_distance(interface_dis_);
		paratope_cst_mover_->apply(pose);
	}
}

void
AntibodyDesignMover::finalize_setup(Pose & pose){

	seq_design_creator_ = AntibodySeqDesignTFCreatorOP( new AntibodySeqDesignTFCreator(ab_info_, cdr_seq_design_options_, false, 2 /* cdr_stem_size */) );
	framework_mutator_ = MutateFrameworkForClusterOP( new MutateFrameworkForCluster(ab_info_));

	framework_mutator_->set_scorefxn(scorefxn_);
	framework_mutator_->set_pack_shell(neighbor_dis_);
	framework_mutator_->keep_current(false);

	//Reinit CDRDBPoseSet.  No option to set your own CDRs at the moment.  Later maybe.
	initialize_cdr_set(pose);
	setup_cdr_pose_sampling_strategies();

	///////////Reinitialize.  Need to figure out how to do this properly via JD.
	top_designs_.clear();
	top_scores_.clear();
	graft_log_.clear();
	cdrs_to_design_.clear();

	////////// Setup Which CDRs we are grafting ////////////////////////////////

	for ( core::Size i=1; i<=CDRNameEnum_total; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if ( cdr_graft_design_options_[cdr]->design()  || cdr_seq_design_options_[ cdr ]->design() ) {
			if ( cdr_graft_design_options_[ cdr]->design() && cdr_set_[cdr].size() == 0 ) {
				//Exit so that we don't continue thinking we are graft designing some cdr when we really are not...
				utility_exit_with_message("CDR "  + ab_info_->get_CDR_name(cdr) + " set to graft design"
					" but no CDRs are in the set.  Please double check instructions settings.");
			}
			cdrs_to_design_.push_back(cdr);
		}
	}

	if ( cdrs_to_design_.size() == 0 ) {
		TR << "All CDRs fixed for low res graft designer or there are no CDRs in the set...." << std::endl;
		//Make sure our top designs has something in it for the antibody designer.
		check_for_top_designs(pose);
		return;
	}

	//////////// Setup GraftMovers and settings ////////////////////////////////
	graft_mover_ = protocols::grafting::CCDEndsGraftMoverOP( new CCDEndsGraftMover(
		ab_info_->get_CDR_start(cdrs_to_design_[1], pose)-1, ab_info_->get_CDR_end(cdrs_to_design_[1], pose)+1) );

	anchored_graft_mover_ = protocols::grafting::AnchoredGraftMoverOP( new AnchoredGraftMover(
		ab_info_->get_CDR_start(cdrs_to_design_[1], pose) - 1, ab_info_->get_CDR_end(cdrs_to_design_[1], pose)+1) );

	setup_default_graft_settings();

	total_permutations_ = 1;
	for ( core::Size i=1; i<=6; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if ( cdr_graft_design_options_[  cdr ]->design() ) {
			total_permutations_ *= cdr_set_[ cdr ].size();
		}
	}
	TR<< "///// Total CDRs in set /////"<<std::endl;
	for ( core::Size i = 1; i<=6; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		modeler_->cdr_overhang(cdr, 2); // 2 Residue overhang for CDR packing and CDR minimization to go with grafting.  North CDR definitions already go into the beta barrel by a few residues
		if ( cdr_graft_design_options_[ cdr ]->design() ) {
			TR << "/// "<<ab_info_->get_CDR_name(cdr)<<" "<<cdr_set_[cdr].size()<<std::endl;
		}
	}

	TR <<"Total possible CDR combinations: "<< total_permutations_ << std::endl;

}

bool
AntibodyDesignMover::apply_to_cdr(Pose & pose, CDRNameEnum cdr, core::Size index, bool min_post_graft /*true*/){
	using namespace protocols::docking;

	//Graft CDR into scaffold pose.  If unsuccessful, do not continue.
	// This will be refactored into a general graft mover soon.

	//Remove constraints to increase speed of insertion. Quicker to remove and re-add
	// Needs to be done every step due to MC - we need the constraints to match.  Presumably, the pose should have constraints
	// Just to make sure, we remove and re-add, which is pretty fast.
	pose.remove_constraints();
	core::pose::Pose original_pose = pose;
	std::string graft_id;

	try{

		if ( cdr_graft_design_options_[ cdr ]->design() ) {


			CDRDBPose & cdr_pose = cdr_set_[cdr][index];
			TR << "Grafting CDR from cluster " << ab_info_->get_cluster_name(cdr_pose.cluster) << " fragment "<< cdr_pose.pdb << std::endl;
			utility::to_string(index)+"_"+ab_info_->get_cluster_name(cdr_pose.cluster)+"_"+cdr_pose.pdb;
			core::Size start = ab_info_->get_CDR_start(cdr, pose)-1;
			core::Size end = ab_info_->get_CDR_end(cdr, pose)+1;
			graft_mover_->set_insert_region(start, end);
			anchored_graft_mover_->set_insert_region(start, end);
			Pose temp_pose = pose;
			std::pair<bool, core::Size> cb = run_graft(temp_pose, cdr, cdr_pose, graft_mover_);
			if ( cb.first && adapt_graft_ ) {
				TR << "Graft not closed. Adapting graft closure" << std::endl;
				Pose temp_pose = pose;
				cb = run_graft(temp_pose, cdr, cdr_pose, anchored_graft_mover_);
			}
			if ( ! cb.first ) {
				pose = temp_pose;
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Print data about graft closure to parse from output. Return if not closed.

			if ( print_tracer_info_ ) {
				CDRClusterEnum cluster = cdr_pose.cluster;
				std::string out = "DATA GRAFT_CLOSURE " + utility::to_string(!cb.first) + " " + pose.pdb_info()->name() + " "
					+ ab_info_->get_CDR_name(cdr) + " " + ab_info_->get_cluster_name(cluster) + " " + cdr_pose.pdb ;

				TR << out << std::endl;
				if ( add_log_to_pose_ ) {
					graft_log_.push_back(out);
				}
			}
			//Graft not successful.  Move on.
			if ( cb.first ) {
				return false;
			}

			////////// Update AntibodyInfo and DataCache ///////////////////////
			//
			//

			//Setup the clusters in AntibodyInfo before any minimization can presumably change the cluster.
			//This is primarily for SeqDesign during minimization.

			core::Size grafted_cdr_length = ab_info_->get_CDR_length(cdr, pose);

			CDRClusterOP cluster( new CDRCluster(
				pose,
				cdr,
				grafted_cdr_length,
				cdr_pose.cluster,
				ab_info_->get_CDR_start(cdr, pose),
				cdr_pose.distance) );


			// 1) Update AntibodyInfo with whatever info is in the datacache. This is important,
			//  as the MC object can revert to previous poses, making AntibodyInfo clusters info out of date.
			// 2) Update AbInfo with the new cluster.
			// 3) Update the datacache with this new cluster info.
			ab_info_->setup_CDR_clusters(pose, true /* setup data from any datacache */);
			ab_info_->set_CDR_cluster(cdr, cluster);
			ab_info_->get_CDR_cluster_set()->set_cacheable_cluster_data_to_pose(pose);
			set_native_cdr_sequence(ab_info_, cdr, pose);

			core::pose::add_comment(pose, "REMARK "+ab_info_->get_CDR_name( cdr )+"_origin", cdr_pose.pdb);

			////////// Constraints /////////////////////////////////////////////
			//
			//

			//Now that AntibodyInfo and the Cacheable data is updated.  Add the correct CDR constraints to all CDRs.
			for ( core::Size i = 1; i <= core::Size(ab_info_->get_total_num_CDRs()); ++i ) {
				CDRNameEnum cdr_in_pose = static_cast<CDRNameEnum>(i);
				cdr_dihedral_cst_mover_->set_cdr(cdr_in_pose);
				if ( cdr_in_pose == cdr ) {
					cdr_dihedral_cst_mover_->set_force_cluster(cdr_pose.cluster);
					cdr_dihedral_cst_mover_->apply(pose);
				} else {
					cdr_dihedral_cst_mover_->set_remove_any_set_forced_cluster();
					cdr_dihedral_cst_mover_->apply(pose);
				}
			}
		} else {
			// No Graft.  Make sure clusters are setup.
			graft_id = "no_graft";
			ab_info_->setup_CDR_clusters(pose, true /* setup data from any datacache */);
			ab_info_->get_CDR_cluster_set()->set_cacheable_cluster_data_to_pose(pose);


			////////// Constraints /////////////////////////////////////////////
			//
			//
			for ( core::Size i = 1; i <= core::Size(ab_info_->get_total_num_CDRs()); ++i ) {
				CDRNameEnum cdr_in_pose = static_cast<CDRNameEnum>(i);
				cdr_dihedral_cst_mover_->set_cdr(cdr_in_pose);
				cdr_dihedral_cst_mover_->set_remove_any_set_forced_cluster();
				cdr_dihedral_cst_mover_->apply(pose);
			}
		}

		///////// Framework Mutations for Cluster //////////////////////////////
		//
		// (Note: Too late to put in packer and extremely complicated to combine it)
		//
		if ( mutate_framework_for_cluster_ ) {
			framework_mutator_->set_cdr_only(cdr);
			framework_mutator_->apply(pose);
		}
		setup_paratope_epitope_constraints(pose);


		////////// Minimization ////////////////////////////////////////////////
		//
		//
		if ( min_post_graft ) {
			TR << "Beginning MinProtocol" << std::endl;
			protocols::moves::MonteCarlo mc = protocols::moves::MonteCarlo(pose, *scorefxn_, inner_kt_);
			for ( core::Size i = 1; i <= inner_cycles_; ++i ) {

				TR << "Inner round: "<< i << std::endl;
				run_optimization_cycle(pose, mc, cdr);
			}
			mc.recover_low(pose);
		}
		TR << "Cycle complete" << std::endl;
		return true;
	}
catch ( utility::excn::EXCN_Base& excn ) {
	std::cerr << "Exception : " << std::endl;
	excn.show( std::cerr );
	excn.show( TR );
	TR << "caught exception.  skipping graft. printing structures." << std::endl;

	//Attempt output of structures for debugging.
	using namespace protocols::jd2;

	TR << pose << std::endl;
	TR << pose.fold_tree() << std::endl;
	JobOP job = protocols::jd2::JobDistributor::get_instance()->current_job();

	//std::string graft_id = utility::to_string(index)+"_"+ab_info_->get_cluster_name(cdr_pose.cluster)+"_"+cdr_pose.pdb;
	//std::string cdr_prefix = "excn_cdr_out_"+graft_id;
	//std::string original_prefix = "excn_start_out_"+graft_id;
	std::string grafted_prefix = "excn_grafted_out_"+graft_id;
	TR << grafted_prefix << std::endl;

	//JobDistributor::get_instance()->job_outputter()->other_pose(job, *cdr_pose.pose, cdr_prefix);
	//JobDistributor::get_instance()->job_outputter()->other_pose(job, original_pose, original_prefix);
	JobDistributor::get_instance()->job_outputter()->other_pose(job, pose, grafted_prefix);

	pose = original_pose;
	return false;
}
}


std::pair<bool, core::Size>
AntibodyDesignMover::run_graft(core::pose::Pose& pose, CDRNameEnum const cdr, CDRDBPose & cdr_pose, AnchoredGraftMoverOP grafter){

	core::Size nter_flex = grafter->get_nterm_insert_flexibility();
	core::Size cter_flex = grafter->get_cterm_insert_flexibility();
	grafter->set_piece(*cdr_pose.pose, overhang_, overhang_);
	if ( cdr_pose.pose->size() - overhang_ - overhang_ <= 4 ) {
		grafter->set_insert_flexibility(1, 1);
	}
	grafter->apply(pose);
	grafter->set_insert_flexibility(nter_flex, cter_flex);

	//pose.pdb_info()->copy(*(cdr_pose.pose->pdb_info()), 1 + overhang_, cdr_pose.pose->size()-overhang_, grafter->start()+1);
	//pose.pdb_info()->obsolete(false);

	//Debugging cartesian+relax bug.  First, force removal of all cutpoint variants in the pose. Any cuts should be termini here.
	protocols::loops::remove_cutpoint_variants(pose, true);

	//Apply CDR constraints for only the grafted CDR
	cdr_dihedral_cst_mover_->set_cdr(cdr);
	cdr_dihedral_cst_mover_->set_force_cluster(cdr_pose.cluster);
	cdr_dihedral_cst_mover_->apply(pose);

	scorefxn_->score(pose); //Segfault prevention.

	modeler_->set_cdr_only(cdr, true);

	TR << "Checking Geometry" << std::endl;
	std::pair<bool, core::Size> cb = check_cb(pose, grafter->get_loops());

	if ( adapt_graft_ && cb.first && enable_adapt_graft_cartesian_ ) {
		TR << "optimizing graft with bb cartesian min" << std::endl;
		scorefxn_cart_graft_->score(pose);
		cart_min_graft_->set_movemap(modeler_->get_cdrs_movemap_with_overhang(pose, true /*min_bb*/, false /*min_sc*/, false /*min_neighbor_sc*/, false /*min_neighbor_bb*/));
		cart_min_graft_->apply(pose);
		//scorefxn_cart_graft_->show(TR, pose);

		protocols::loops::Loop loo = ab_info_->get_CDR_loop(cdr, pose, 4);
		cb = protocols::loops::has_severe_pep_bond_geom_issues(pose, loo, true, true, 1.5, 15, 15);
	}

	return cb;

}

void
AntibodyDesignMover::run_optimization_cycle(core::pose::Pose& pose, protocols::moves::MonteCarlo & mc, CDRNameEnum cdr) {


	//Setup Neighbor CDR minimization - instruction file/options controlled..
	// This should either be refactored into its own class, or the modeler should incorporate most of these.
	//
	scorefxn_->score( pose );
	CDRGraftDesignOptionsOP options = cdr_graft_design_options_[ cdr ];
	CDRSeqDesignOptionsOP seq_des_options = cdr_seq_design_options_[ cdr ];

	vector1<bool> cdrs_to_min( core::Size( CDRNameEnum_proto_total ), false );
	utility::vector1<CDRNameEnum> neighbor_min = options->neighbor_min();

	//Convert from vector of cdrs to vector1 bool.
	cdrs_to_min[ cdr ] = true;
	for ( core::Size i = 1; i <= neighbor_min.size(); ++i ) {
		TR <<"Add min neighbors : "<< ab_info_->get_CDR_name( neighbor_min[ i ] )<<std::endl;
		cdrs_to_min[ neighbor_min[ i ] ] = true;
	}

	modeler_->set_cdrs(cdrs_to_min);

	//Setup design task factory.
	TaskFactoryOP design_tf = seq_design_creator_->generate_tf_seq_design_graft_design( pose, cdr, cdrs_to_min );
	//TR <<"CDRDesignTF" << std::endl;
	//design_tf->create_task_and_apply_taskoperations(pose)->show(TR);

	//Setup dock design tf
	TaskFactoryOP dock_design_tf;

	if ( dock_post_graft_ ) {
		TR << "Setting up DockDesign TF" << std::endl;
		dock_design_tf = seq_design_creator_->generate_tf_seq_design( pose );
		seq_design_creator_->disable_design_for_non_designing_cdrs( dock_design_tf, pose );
		dock_design_tf->push_back(TaskOperationCOP( new protocols::toolbox::task_operations::RestrictToInterface( 1, interface_dis_) ));
	}
	//TR <<"DockDesign task" << std::endl;
	//dock_design_tf->create_task_and_apply_taskoperations(pose)->show(TR);


	//Dock - before min?
	if ( dock_post_graft_ ) {
		for ( core::Size i = 1; i<= dock_cycles_; i++ ) {

			TR << "Dock round "<< i <<std::endl;

			modeler_->set_task_factory(dock_design_tf);

			modeler_->dock_low_res(pose, true /*repack_interface*/ ); // This should change once the minimization does some neighbor detection.
			modeler_->minimize_interface(pose); //Seems to help
			modeler_->dock_high_res(pose, 3 /*first cycles*/, 10 /*second_cycles*/); //Normal DockMCM is 4/45. This should mainly just be quick to fix overlap from low_res dock.
			mc.boltzmann(pose); //Low res dock will definitely screw up the pose - so do MC after both low and high res
		}
	}

	if ( seq_des_options->design() ) {
		modeler_->set_task_factory(design_tf);
	}

	//Repack CDR?
	if ( options->mintype() == repack ) {
		modeler_->repack_cdrs(pose, options->min_neighbor_sc());
		mc.boltzmann(pose);
	} else if ( options->mintype() == relax ) {
		//Relax CDR?
		modeler_->relax_cdrs(pose, options->min_neighbor_sc(), false, options->min_rb());
		mc.boltzmann(pose);
	} else if ( options->mintype() == dualspace ) {
		modeler_->relax_cdrs(pose, options->min_neighbor_sc(), false, options->min_rb(), true /*dualspace*/);
		mc.boltzmann(pose);
	} else if ( options->mintype() == minimize ) {
		//Minimize CDR?

		modeler_->repack_cdrs(pose, options->min_neighbor_sc());
		mc.boltzmann(pose);
		modeler_->minimize_cdrs(pose, true, options->min_neighbor_sc(), options->min_rb());
		mc.boltzmann(pose);
	} else if ( options->mintype() == minimize_cartesian ) {

		modeler_->repack_cdrs(pose, options->min_neighbor_sc());
		mc.boltzmann(pose);
		modeler_->minimize_cdrs(pose, true, options->min_neighbor_sc(), options->min_rb(), true /*cartesian*/);
		mc.boltzmann(pose);
	} else if ( options->mintype() == backrub_protocol ) {
		modeler_->backrub_cdrs(pose, true, options->min_neighbor_sc());
		modeler_->repack_cdrs(pose, options->min_neighbor_sc());
		mc.boltzmann(pose);
	}
	//Interface Rigid Body minimize?
	if ( rb_min_post_graft_ ) {
		modeler_->minimize_interface(pose, true /* min_interface_sc */);
		mc.boltzmann(pose);
	}

	if ( dock_post_graft_ && dock_min_dock_ ) {
		for ( core::Size i = 1; i<= dock_cycles_; i++ ) {

			TR << "Dock round "<< i <<std::endl;

			modeler_->set_task_factory(dock_design_tf);

			//modeler_->dock_low_res(pose, true /*repack_interface*/ ); // This should change once the minimization does some neighbor detection.
			//modeler_->minimize_interface(pose); //Seems to help
			modeler_->dock_high_res(pose, 3 /*first cycles*/, 10 /*second_cycles*/); //Normal DockMCM is 4/45. This should mainly just be quick to fix overlap from low_res dock.
			mc.boltzmann(pose); //Low res dock will definitely screw up the pose - so do MC after both low and high res
		}
	}
	modeler_->reset_task_factory();
}

void
AntibodyDesignMover::check_for_top_designs(core::pose::Pose & pose){

	//Can be refactored to use utility::TopScoreSelector
	//From mc algorithm, you can have multiple poses that are equivalent...
	core::Real score = (*scorefxn_)(pose);
	vector1<core::Real>::iterator score_it = top_scores_.begin();
	vector1<core::pose::PoseOP>::iterator pose_it = top_designs_.begin();
	if ( top_scores_.size()==0 ) {
		top_scores_.push_back(score);
		top_designs_.push_back(core::pose::PoseOP( new Pose() ));
		*(top_designs_[top_designs_.size()]) = pose;
	} else {
		bool inserted = false;
		for ( core::Size i = 1; i<=top_scores_.size(); ++i ) {
			if ( score <= top_scores_[i] ) {
				top_scores_.insert(score_it+i-1, score);
				top_designs_.insert(pose_it+i-1, core::pose::PoseOP( new Pose() ));
				*(top_designs_[i]) = pose;
				inserted = true;
				break;
			}
		}
		if ( ! inserted && top_scores_.size() < num_top_designs_ ) {
			top_scores_.push_back(score);
			top_designs_.push_back(core::pose::PoseOP( new Pose() ));
			*(top_designs_[top_designs_.size()]) = pose;
		} else if ( inserted && top_scores_.size() > num_top_designs_ ) {
			top_scores_.pop_back();
			top_designs_.pop_back();
		}
	}
	//mc_->eval_lowest_score_pose(pose, false, true);
}

/// @brief Gets a list of vectors whose indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
vector1< vector1 < Size > >
AntibodyDesignMover::get_cdr_set_index_list(){
	vector1< vector1< Size > > index_list;
	vector1<core::Size> cdr_set_totals(6, 0);

	for ( core::SSize i=CDRNameEnum_start; i<=CDRNameEnum_total; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		cdr_set_totals[i] = cdr_set_[cdr].size();
	}

	vector1< vector1< Size > > empty_vector;
	get_all_graft_permutations(empty_vector, cdr_set_totals, 1);
	return index_list;

}

void
AntibodyDesignMover::setup_random_start_pose(core::pose::Pose& pose, vector1<CDRNameEnum>& cdrs_to_design){

	TR << "Randomizing starting CDRs" << std::endl;
	print_tracer_info_ = false;
	core::Size max_attempts = 1000;
	bool original_adaptation = adapt_graft_;
	for ( core::Size i = 1; i <= cdrs_to_design.size(); ++i ) {
		CDRNameEnum cdr = cdrs_to_design[i];

		if ( ! cdr_graft_design_options_[ cdr ]->design() ) { continue;}

		bool graft_closed = false;
		core::Size attempts = 0;

		core::pose::Pose temp_pose;
		while ( ! graft_closed ) {
			temp_pose = pose;
			core::Size cdr_index = numeric::random::rg().random_range(1, cdr_set_[cdr].size());
			graft_closed = this->apply_to_cdr(temp_pose, cdr, cdr_index, false); //No minimization
			++attempts;

			//Infinite loop prevention for stubborn CDRs
			if ( attempts == max_attempts ) {
				if ( ! adapt_graft_ ) {
					adapt_graft_ = true;
					attempts = 0;
				} else {
					utility_exit_with_message(
						"Could not setup randomized starting pose as graft for " + ab_info_->get_CDR_name(cdr)+" Max attempts reached!");
				}
			}
		}
		TR << "Attempts "<<attempts << " " << ab_info_->get_CDR_name(cdr) << std::endl;
		adapt_graft_ = original_adaptation;
		pose = temp_pose;
	}

	//JAB - removed 8/18/16.  I still have yet to use this pose for any benchmarking.
	/*
	protocols::jd2::JobOP job = protocols::jd2::JobDistributor::get_instance()->current_job();
	std::string prefix = "initial_benchmark_perturbation";
	protocols::jd2::JobDistributor::get_instance()->job_outputter()->other_pose(job, pose, prefix);
	*/


	print_tracer_info_ = true;
}

void
AntibodyDesignMover::run_basic_mc_algorithm(Pose & pose, vector1<CDRNameEnum>& cdrs_to_design, AntibodyDesignProtocolEnum mc_algorithm){

	using namespace utility;

	TR << "Running basic monte carlo algorithm " << std::endl;

	top_scores_.push_back((*scorefxn_)(pose));
	top_designs_.push_back(core::pose::PoseOP( new Pose() ));
	*(top_designs_[1]) = pose;
	mc_->set_last_accepted_pose(pose);

	//Setup Weights
	utility::vector1<core::Real> cdr_weights;
	numeric::random::WeightedSampler sampler;
	for ( core::Size i = 1; i <= cdrs_to_design.size(); ++i ) {
		core::Real weight = cdr_graft_design_options_[cdrs_to_design[ i ]]->weight();

		//CDRNameEnum cdr = static_cast<CDRNameEnum>(cdrs_to_design[ i ]);
		//TR << ab_info_->get_CDR_name(cdr) <<" Weight " << weight << std::endl;
		cdr_weights.push_back(weight);
	}
	sampler.weights(cdr_weights);

	//Convert to vector for index-based RNG.  Beautiful c++ code right here.  This would be so much better in c++ 11 or python
	//  Probably a better way, though using std::advance is still not ideal...
	std::map< CDRNameEnum, utility::vector1< CDRClusterEnum > > cluster_indexes;

	typedef std::map< CDRNameEnum,
		std::map< clusters::CDRClusterEnum,
		utility::vector1< core::Size > > >::iterator it_type;

	typedef std::map< clusters::CDRClusterEnum,
		utility::vector1< core::Size > >::iterator it_type_clusters;

	for ( it_type it = cluster_based_CDRDBPose_indexes_.begin(); it != cluster_based_CDRDBPose_indexes_.end(); ++it ) {

		//For Length indexing:
		for ( it_type_clusters it2 = cluster_based_CDRDBPose_indexes_[ it->first ].begin(); it2 != cluster_based_CDRDBPose_indexes_[ it->first ].end(); ++it2 ) {
			cluster_indexes[ it->first ].push_back( it2->first );
		}

		//For length and cluster indexing

	};

	//typedef std::map< clusters::CDRClusterEnum,
	//  utility::vector1< core::Size > > ::iterator cluster_based_CDRDBPose_indexes_it_type;

	//Choose random weighted CDR, graft in CDR, minimize
	for ( core::Size i = 1; i <= outer_cycles_; ++i ) {
		TR << "Outer round: " << i <<std::endl;


		CDRNameEnum cdr_type = cdrs_to_design[sampler.random_sample(numeric::random::rg())];
		core::Size cdr_index = 0;

		/// Only get a cdr_index if we are graft designing.
		if ( cdr_graft_design_options_[ cdr_type ]->design() ) {
			if ( mc_algorithm == generalized_monte_carlo ) {
				cdr_index = numeric::random::rg().random_range(1, cdr_set_[cdr_type].size());
			} else if ( mc_algorithm == even_cluster_monte_carlo ) {


				core::Size cluster_index = numeric::random::rg().random_range(1, cluster_indexes[ cdr_type ].size());
				CDRClusterEnum cluster = cluster_indexes[ cdr_type ][cluster_index];
				core::Size pre_cdr_index = numeric::random::rg().random_range(1, cluster_based_CDRDBPose_indexes_[ cdr_type ][ cluster ].size());
				cdr_index = cluster_based_CDRDBPose_indexes_[ cdr_type ][ cluster ][ pre_cdr_index];

				//cluster_based_CDRDBPose_indexes_it_type clus_it = cluster_based_CDRDBPose_indexes_.begin();
				//std::advance(clus_it, numeric::random::rg().random_range(0, cluster_based_CDRDBPose_indexes_.size()));
				//CDRClusterEnum cluster = clus_it->first;

			} else if ( mc_algorithm == even_length_cluster_monte_carlo ) {

				//Get a random length
				//Index needs start from 0 for std::advance
				std::map< core::Size, std::map< CDRClusterEnum, utility::vector1< core::Size > > >::iterator len_it = length_based_CDRDBPose_indexes_[ cdr_type ].begin();
				std::advance(len_it, numeric::random::rg().random_range(0, length_based_CDRDBPose_indexes_[ cdr_type ].size() - 1));
				core::Size length = len_it->first;

				//Get a random cluster within that length
				//Index needs start from 0 for std::advance
				std::map<CDRClusterEnum, utility::vector1< core::Size > >::iterator clus_it = length_based_CDRDBPose_indexes_[ cdr_type ][ length ].begin();
				std::advance(clus_it, numeric::random::rg().random_range(0, length_based_CDRDBPose_indexes_[ cdr_type ][ length ].size() - 1));
				CDRClusterEnum cluster = clus_it->first;

				//Get a random CDR within that cluster
				core::Size pre_cdr_index = numeric::random::rg().random_range(1, length_based_CDRDBPose_indexes_[ cdr_type ][ length ][ cluster ].size());

				cdr_index = length_based_CDRDBPose_indexes_[ cdr_type ][ length ][ cluster ][ pre_cdr_index ];

			}
		}

		bool successful = apply_to_cdr(pose, cdr_type, cdr_index);
		if ( successful ) {
			check_for_top_designs(pose);

			core::Real energy = scorefxn_->score(pose);
			bool accepted = mc_->boltzmann(pose);
			std::string out = to_string(i)+" "+to_string(energy)+" "+to_string(accepted);
			accept_log_.push_back(out);

			out = "FINAL "+to_string(i)+" "+to_string(scorefxn_->score(pose));
			accept_log_.push_back(out);
		} else {
			pose = mc_->last_accepted_pose();
		}
	}
	mc_->recover_low(pose);
	std::string out = "FINAL END "+to_string(scorefxn_->score(pose));
	accept_log_.push_back(out);
	mc_->show_counters();
}


void
AntibodyDesignMover::run_deterministic_graft_algorithm(core::pose::Pose & pose,vector1<CDRNameEnum>& cdrs_to_design){
	//Note:  Not feasible with >= 4 CDRs to try at each position. Should not be used with docking on.
	//Temporary fix to deterministically graft one CDR.  Needs to be fixed correctly soon.
	CDRNameEnum cdr = cdrs_to_design[1];
	for ( core::Size i = 1; i <= cdr_set_[cdr].size(); ++i ) {
		TR << "Graft round: " << i << std::endl;
		core::pose::Pose trial_pose = pose;
		bool graft_successful = apply_to_cdr(trial_pose, cdr, i);
		if ( ! graft_successful ) {
			TR << "Graft unsuccessful.  Skipping combination..." << std::endl;
			continue;
		}
		//Only check for top designs after all combinations have been grafted.
		check_for_top_designs(trial_pose);
		mc_->eval_lowest_score_pose(trial_pose, false, true);
	}
	mc_->recover_low(pose);
	mc_->show_counters();
}

void
AntibodyDesignMover::apply(core::pose::Pose & pose){

	//if (! protocols::antibody::clusters::check_if_pose_renumbered_for_clusters(pose)){
	// utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.");
	//}

	if ( !ab_info_ ) {
		ab_info_ = AntibodyInfoOP(new AntibodyInfo(pose, AHO_Scheme, North));
	}

	ab_info_->show(std::cout);


	//////////////////Create Instances.  Make sure everything is ready to begin.
	setup_scorefxn();
	setup_native_clusters(pose);
	setup_native_sequence(pose);
	setup_options_classes();
	setup_cart_minimizer();
	setup_epitope_residues(pose);
	setup_paratope_epitope_constraints(pose);
	setup_modeler();

	finalize_setup(pose);
	show(std::cout);

	if ( cdrs_to_design_.size() == 0 ) {
		return;
	}

	scorefxn_->show(pose);
	core::Real native_score = (*scorefxn_)(pose);

	///Energy Log:
	std::string out = "-1 "+utility::to_string(native_score)+" NA";
	accept_log_.push_back(out);
	if ( benchmark_ ) {
		this->setup_random_start_pose(pose, cdrs_to_design_);
		out = "0 "+utility::to_string(scorefxn_->score(pose))+" NA";
		accept_log_.push_back(out);
	}

	//////////// Setup Monte carlo ////////////////////////////////////////////
	mc_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo(pose, *scorefxn_, outer_kt_) );
	core::Real benchmark_start_score = scorefxn_->score(pose);

	//////////// Run Main Algorithms ///////////////////////////////////////////
	if ( design_protocol_ == deterministic_graft && cdrs_to_design_.size() == 1 ) {
		run_deterministic_graft_algorithm(pose, cdrs_to_design_);
	} else if ( design_protocol_ == deterministic_graft && cdrs_to_design_.size() > 1 ) {
		utility_exit_with_message("Cannot currently run the deterministic protocol for more than 1 "
			"grafting CDR.  Set CDRs using the instructions file or the option -design_cdrs");
	} else {
		run_basic_mc_algorithm(pose, cdrs_to_design_, design_protocol_);
	}

	//} else {
	// utility_exit_with_message("Design Protocol not understood: "+design_protocol_to_string( design_protocol_ )+" See -design_protocol option for available protocols");
	//}


	//////// Print Score Information ///////////////////////////////////////////
	TR << "Native Pose: " << native_score << std::endl;
	if ( benchmark_ ) {
		TR << "Benchmarked Pose: "<< benchmark_start_score << std::endl;
	}
	TR << "Final Pose: " << (*scorefxn_)(pose) << std::endl;

	for ( core::Size i = 2; i<= top_scores_.size(); ++i ) {
		TR << "Top Ensemble " << i << " : " << (*scorefxn_)(*top_designs_[i]) << std::endl;
	}

	//////// Add Logs for Benchmarking /////////////////////////////////////////
	if ( add_log_to_pose_ ) {
		for ( core::Size i = 1; i <= graft_log_.size(); ++i ) {
			core::pose::add_comment(pose, "GRAFT LOG "+utility::to_string(i), graft_log_[i]);

			for ( core::Size x = 1; x <= top_designs_.size(); ++x ) {
				core::pose::add_comment(*top_designs_[x], "GRAFT LOG "+utility::to_string(i), graft_log_[i]);
			}
		}
		for ( core::Size i = 1; i <= accept_log_.size(); ++i ) {
			core::pose::add_comment(pose, "ACCEPT LOG "+utility::to_string(i), accept_log_[i]);
			for ( core::Size x = 1; x <= top_designs_.size(); ++x ) {
				core::pose::add_comment(*top_designs_[x], "ACCEPT LOG "+utility::to_string(i), accept_log_[i]);
			}
		}
	}

	//Needs to be called for RosettaScripts integration - should be called for each 'get_other_pose' or whatever it is.
	add_cluster_comments_to_pose( pose, ab_info_ );
	check_fix_aho_cdr_numbering(ab_info_, pose);
}

////////////////////////////////////////////// Boiler Plate ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// XRW TEMP std::string
// XRW TEMP AntibodyDesignMover::get_name() const{
// XRW TEMP  return "AntibodyDesignMover";
// XRW TEMP }

void
AntibodyDesignMover::show(std::ostream & output) const{

	utility::vector1<std::string> strategies(SeqDesignStrategyEnum_total, "");
	strategies[ seq_design_profiles ] = "Profile-Based";
	strategies[ seq_design_conservative ] = "Conservative";
	strategies[ seq_design_basic ] = "Basic";
	strategies[ seq_design_none ] = "No Design";

	output << "////////////////////////////////////////////////////////////////////////////////" <<std::endl;
	output << "///                   Rosetta Antibody Design Settings                       ///" << std::endl;
	output << "///                                                                          ///" << std::endl;
	output << "// Design Protocol: " << design_enum_manager_->design_protocol_enum_to_string( design_protocol_ ) << std::endl;
	output << "//  Outer Cycles: "<< outer_cycles_ << std::endl;
	output << "//  Inner Cycles: "<< inner_cycles_ << std::endl;
	output << "//  Dock? " << std::boolalpha << dock_post_graft_ << std::endl;
	output << std::endl;
	for ( core::Size i=1; i<=core::Size( ab_info_->get_total_num_CDRs() ); ++i ) {

		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);


		CDRSetOptionsCOP set_options = cdr_set_options_[cdr];
		CDRGraftDesignOptionsCOP graft_options = cdr_graft_design_options_[cdr];
		CDRSeqDesignOptionsCOP seq_options = cdr_seq_design_options_[cdr];

		if ( (! graft_options->design() ) && (! seq_options->design() ) ) {
			continue;
		}

		output << "///////////////////////////////////// " << ab_info_->get_CDR_name(cdr) << " ///////////////////////////////////////"<< std::endl;
		//output << "///                                                                          ///" << std::endl;
		output << "///  Weight: "<< graft_options->weight() << std::endl;
		output << "/// "<< std::endl;
		output << "///  Graft? " << std::boolalpha << graft_options->design() << std::endl;
		//output << std::endl;
		output << "///  SeqDesign? " << std::boolalpha << seq_options->design() << std::endl;
		if ( seq_options->design() ) {
			output<< "///   Design Strategy: " << strategies[ seq_options->design_strategy() ] << std::endl;
			output<< "///   Fallback Strategy: " << strategies[ seq_options->fallback_strategy() ] << std::endl;
			//output<< "///   Fallback? " << std::boolalpha << seq_options->fallback() << std::endl;

		}

		output << "/// "<< std::endl;
		if ( graft_options->neighbor_min().size() >=1 ) {
			utility::vector1<CDRNameEnum> neighbor_cdrs = graft_options->neighbor_min();
			utility::vector1<std::string> neighbor_str;
			for ( core::Size i = 1; i <= neighbor_cdrs.size(); ++i ) {
				neighbor_str.push_back(ab_info_->get_CDR_name(neighbor_cdrs[ i ]));
			}
			output<< "///  Min Neighbors:  " << utility::to_string(neighbor_str) << std::endl;
		}

		if ( ! graft_options->design() ) continue;
		output<< "///  Current Clusters only? " << std::boolalpha << set_options->include_only_current_cluster() << std::endl;
		output<< "///  Center Clusters only? " << std::boolalpha << set_options->include_only_center_clusters() << std::endl;
		if ( cdr!=h3 ) {
			output << "////" << std::endl;
			output << "////// Length Types: " << std::endl;
			output << "////" << std::endl;
			output << "///  1 " << std::boolalpha << set_options->length_type()[1] << std::endl;
			output << "///  2 " << std::boolalpha << set_options->length_type()[2] << std::endl;
			output << "///  3 " << std::boolalpha << set_options->length_type()[3] << std::endl;
		}
		output << "////" << std::endl;
		output << "////// Lengths: " << std::endl;
		output << "////" << std::endl;
		output << "/// Min " << set_options->min_length() << std::endl;
		output << "/// Max " << set_options->max_length() << std::endl;
		output << "///" << std::endl;

		//Include/Leave out options:

		vector1<std::string> exclude_clusters;
		vector1<std::string> include_clusters;
		vector1<CDRClusterEnum> exclude_clusters_enum = set_options->exclude_clusters();
		vector1<CDRClusterEnum> include_clusters_enum = set_options->include_only_clusters();

		for ( core::Size i = 1; i <= exclude_clusters_enum.size(); ++i ) {
			exclude_clusters.push_back(ab_info_->get_cluster_name(exclude_clusters_enum[i]));
		}
		for ( core::Size i = 1; i <= include_clusters_enum.size(); ++i ) {
			include_clusters.push_back(ab_info_->get_cluster_name(include_clusters_enum[i]));
		}

		output << "////// Include Only: " << std::endl;
		print_str_vec("PDBs", set_options->include_only_pdbs(), output);
		print_str_vec("Clusters", include_clusters, output);
		print_str_vec("Species", set_options->include_only_species(), output);
		print_str_vec("Germlines", set_options->include_only_germlines(), output);

		output<<"////// Exclude: " << std::endl;
		print_str_vec("PDBs", set_options->exclude_pdbs(), output);
		print_str_vec("Clusters", exclude_clusters, output);
		print_str_vec("Species", set_options->exclude_species(), output);
		print_str_vec("Germlines", set_options->exclude_germlines(), output);
		output<< std::endl;
	}
	output << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////" <<std::endl;
}

void
AntibodyDesignMover::print_str_vec(std::string const name, utility::vector1<std::string> const & vec, std::ostream & output) const {
	if ( vec.size()>=1 ) {
		output<< "/// "<< name <<" "<<utility::to_string(vec)<<std::endl;
	}
}

std::string AntibodyDesignMover::get_name() const {
	return mover_name();
}

std::string AntibodyDesignMover::mover_name() {
	return "AntibodyDesignMover";
}

void AntibodyDesignMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace basic::options;

	AttributeList attlist;

	attributes_for_get_ab_design_min_scorefxn(attlist);
	attributes_for_get_ab_design_global_scorefxn(attlist);

	XMLSchemaRestriction ABcdr_enum;
	ABcdr_enum.name("ABcdr_definitions");
	ABcdr_enum.base_type(xs_string);
	AntibodyEnumManager ABenum_manager;
	for ( auto& cdr_def : ABenum_manager.get_recognized_cdr_definitions() ) {
		ABcdr_enum.add_restriction(xsr_enumeration, cdr_def);
	}
	xsd.add_top_level_element(ABcdr_enum);

	attlist + XMLSchemaAttribute(
		"design_cdrs", "ABcdr_definitions",
		"CDR regions to be designed");

	attlist + XMLSchemaAttribute(
		"instruction_file", xs_string,
		"Path to the CDR instruction file (see application documentation for format)");
	attlist + XMLSchemaAttribute(
		"instructions_file", xs_string,
		"used if instruction_file attribute is not specified");
	attlist + XMLSchemaAttribute(
		"cdr_instructions_file", xs_string,
		"used if instructions_file attribute is not specified");

	XMLSchemaRestriction ABdesign_enum;
	ABdesign_enum.name("ABdesign_protocols");
	ABdesign_enum.base_type(xs_string);
	AntibodyDesignEnumManager ABdesign_enum_manager;
	for ( auto& protocol : ABdesign_enum_manager.get_recognized_design_protocols() ) {
		ABdesign_enum.add_restriction(xsr_enumeration, protocol);
	}
	xsd.add_top_level_element(ABdesign_enum);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"design_protocol", "ABdesign_protocols",
		"Sets the design protocol (see app documentation for more information)",
		"GENERALIZED_MONTE_CARLO");

	attlist + XMLSchemaAttribute(
		"interface_dis", xsct_real,
		"Set the interface detection distance");

	attlist + XMLSchemaAttribute(
		"neighbor_dis", xsct_real,
		"Set the neighbor detection distance");

	attlist + XMLSchemaAttribute(
		"outer_cycles", xsct_non_negative_integer,
		"Set the number of outer cycles");

	attlist + XMLSchemaAttribute(
		"inner_cycles", xsct_non_negative_integer,
		"Set the number of inner (minimization) cycles");

	attlist + XMLSchemaAttribute(
		"outer_kt", xsct_real,
		"Temperature to use for outer cycle");

	attlist + XMLSchemaAttribute(
		"inner_kt", xsct_real,
		"Temperature to use for inner cycle");

	attlist + XMLSchemaAttribute(
		"top_designs", xsct_non_negative_integer,
		"Number of top designs to keep");

	attlist + XMLSchemaAttribute(
		"do_dock", xsct_rosetta_bool,
		"Run RosettaDock during the inner cycles? Significantly increases run time");

	attlist + XMLSchemaAttribute(
		"do_rb_min", xsct_rosetta_bool,
		"Perform rigid body minimization during the inner cycles?");

	attlist + XMLSchemaAttribute(
		"dock_cycles", xsct_non_negative_integer,
		"Change the number of time the dock protocol is run");

	attlist + XMLSchemaAttribute(
		"dock_min_dock", xsct_rosetta_bool,
		"Perform a quick high resolution dock after minimization?");

	attlist + XMLSchemaAttribute(
		"benchmark", xsct_rosetta_bool,
		"Start with random CDRs from the antibody design database for any undergoing GraftDesign");

	attlist + XMLSchemaAttribute(
		"use_light_chain_type", xsct_rosetta_bool,
		"Type of light chain to use when selecting CDR set");

	attlist + XMLSchemaAttribute(
		"adapt_graft", xsct_rosetta_bool,
		"Adapt graft closure?");

	attlist + XMLSchemaAttribute(
		"enable_adapt_graft_cartesian", xsct_rosetta_bool,
		"Use cartesian minimization in adapt graft");

	attlist + XMLSchemaAttribute(
		"idealize_graft_cdrs", xsct_rosetta_bool,
		"Idealize graft before inserting?");

	attlist + XMLSchemaAttribute(
		"add_graft_log_to_pdb", xsct_rosetta_bool,
		"Add a full grafting log to the pose including origin PDBs, clusters, etc.");

	attlist + XMLSchemaAttribute(
		"mutate_framework_for_cluster", xsct_rosetta_bool,
		"Should we add framework mutations for the specified cdr?");

	attlist + XMLSchemaAttribute(
		"use_epitope_csts", xsct_rosetta_bool,
		"Use the ParatopeEpitopeSiteConstraintMover during design instead of just the ParatopeSiteConstraintMover");

	attlist + XMLSchemaAttribute(
		"epitope_residues", xs_string,
		"Use these residues as the epitope residues (auto-detects by default )");

	attributes_for_get_cdr_bool_from_tag(attlist, "paratope_cdrs", "Specifically set the paratope as these cdrs.");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Part of the RosettaAntibody and RosettaAntibodyDesign (RAbD) Framework. Runs the full Rosetta Antibody Design protocol.",
		attlist );
}

std::string AntibodyDesignMoverCreator::keyname() const {
	return AntibodyDesignMover::mover_name();
}

protocols::moves::MoverOP
AntibodyDesignMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AntibodyDesignMover );
}

void AntibodyDesignMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AntibodyDesignMover::provide_xml_schema( xsd );
}


} //design
} //antibody
} //protocols
