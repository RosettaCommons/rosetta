// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#include <protocols/comparative_modeling/hybridize/HybridizeProtocol.hh>
#include <protocols/comparative_modeling/hybridize/FoldTreeHybridize.hh>
#include <protocols/comparative_modeling/hybridize/CartesianHybridize.hh>
#include <protocols/comparative_modeling/hybridize/TemplateHistory.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>

#include <core/sequence/util.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <numeric/random/WeightedSampler.hh>
#include <ObjexxFCL/format.hh>

// evaluation
#include <core/scoring/rms_util.hh>
//#include <protocols/comparative_modeling/Align_RmsdEvaluator.hh>
#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/loops/util.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <string>

static basic::Tracer TR( "protocols.comparative_modeling.hybridize.HybridizeProtocol" );
static numeric::random::RandomGenerator RG(541938);

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace kinematics;
using namespace sequence;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;

HybridizeProtocol::HybridizeProtocol() :
	template_weights_sum_(0)
{
	check_options();
	
	//read templates
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	read_template_structures( option[cm::hybridize::template_list]() );
	
	//read fragments
	if ( option[ in::file::frag9 ].user() ) {
		using namespace core::fragment;
		fragments9_ = new ConstantLengthFragSet( 9 );
		fragments9_ = FragmentIO().read_data( option[ in::file::frag9 ]() );
		for (core::fragment::FrameIterator i = fragments9_->begin(); i != fragments9_->end(); ++i) {
			core::Size position = (*i)->start();
			//library_[position] = **i;
		}
	}
	if ( option[ in::file::frag3 ].user() ) {
		using namespace core::fragment;
		fragments3_ = new ConstantLengthFragSet( 3 );
		fragments3_ = FragmentIO().read_data( option[ in::file::frag3 ]() );
		for (core::fragment::FrameIterator i = fragments9_->begin(); i != fragments9_->end(); ++i) {
			core::Size position = (*i)->start();
			//library_[position] = **i;
		}
	}
	
	// native
	if ( option[ in::file::native ].user() ) {
		native_ = new core::pose::Pose;
		core::import_pose::pose_from_pdb( *native_, option[ in::file::native ]() );
	}

}
	
core::Real HybridizeProtocol::get_gdtmm( core::pose::Pose &pose ) {
	core::Real gdtmm = 0;
	if ( native_ && native_->total_residue() > 0) {
		if ( !aln_ ) {
			core::sequence::SequenceOP model_seq ( new core::sequence::Sequence( pose.sequence(),  "model",  1 ) );
			core::sequence::SequenceOP native_seq( new core::sequence::Sequence( native_->sequence(), "native", 1 ) );
			aln_ = new core::sequence::SequenceAlignment;
			*aln_ = align_naive(model_seq,native_seq);
		}
		
		int n_atoms;
		ObjexxFCL::FArray2D< core::Real > p1a, p2a;
		protocols::comparative_modeling::gather_coords( pose, *native_, *aln_, n_atoms, p1a, p2a );
		
		core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
		gdtmm = core::scoring::xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
	}
	return gdtmm;
}
	
	
HybridizeProtocol::~HybridizeProtocol(){}

void HybridizeProtocol::add_template(std::string template_fn, std::string cst_fn, core::Real weight, core::Size cluster_id)
{
	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	core::pose::PoseOP template_pose = new core::pose::Pose();
	core::import_pose::pose_from_pdb( *template_pose, *residue_set, template_fn );

	// add secondary structure information to the template pose
	core::scoring::dssp::Dssp dssp_obj( *template_pose );
	dssp_obj.insert_ss_into_pose( *template_pose );
	
	// find ss chunks in template
	protocols::loops::Loops chunks = protocols::loops::extract_secondary_structure_chunks(*template_pose); 
	
	//break templates into contigs
	protocols::loops::Loops contigs = protocols::loops::extract_continuous_chunks(*template_pose); 
	
	template_fn_.push_back(template_fn);
	templates_.push_back(template_pose);
	template_cst_fn_.push_back(cst_fn);
	template_weights_.push_back(weight);
	template_clusterID_.push_back(cluster_id);
	template_chunks_.push_back(chunks);
	template_contigs_.push_back(contigs);
}

void HybridizeProtocol::read_template_structures(utility::file::FileName template_list)
{
	utility::io::izstream f_stream( template_list );
	std::string line;
	while (!f_stream.eof()) {
		getline(f_stream, line);
        
		if (line.size() == 0) continue;
		
		std::istringstream str_stream(line);
		std::string template_fn;
		std::string cst_fn;
		core::Size cluster_id;
		core::Real weight;
		str_stream >> template_fn >> cst_fn >> cluster_id >> weight;
		TR << template_fn << " " << cst_fn << " " << cluster_id << " " << weight << std::endl;
		add_template(template_fn, cst_fn, weight, cluster_id);
	}
	f_stream.close();
}

void HybridizeProtocol::read_template_structures(utility::vector1 < utility::file::FileName > const & template_filenames)
{
	templates_.clear();
	templates_.resize(template_filenames.size());

	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	for (core::Size i_ref=1; i_ref<= template_filenames.size(); ++i_ref) {
		templates_[i_ref] = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *(templates_[i_ref]), *residue_set, template_filenames[i_ref].name() );
		
		core::scoring::dssp::Dssp dssp_obj( *templates_[i_ref] );
		dssp_obj.insert_ss_into_pose( *templates_[i_ref] );
	}
}

void HybridizeProtocol::check_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ cm::hybridize::template_list ].user() ) {
		utility_exit_with_message("Error! Need the -cm::hybridize::template_list for input template structures.");
	}
}
	
void HybridizeProtocol::pick_starting_template(core::Size & initial_template_index,
							core::Size & initial_template_index_icluster,
							utility::vector1 < core::Size > template_index_icluster,
							utility::vector1 < core::pose::PoseOP > & templates_icluster,
							utility::vector1 < core::Real > & weights_icluster,
							utility::vector1 < protocols::loops::Loops > & template_chunks_icluster,
							utility::vector1 < protocols::loops::Loops > & template_contigs_icluster )
{
	template_index_icluster.clear();
	templates_icluster.clear();
	weights_icluster.clear();
	template_chunks_icluster.clear();
	template_contigs_icluster.clear();
	
	numeric::random::WeightedSampler weighted_sampler;
	weighted_sampler.weights(template_weights_);

	initial_template_index = weighted_sampler.random_sample(RG);
	Size cluster_id = template_clusterID_[initial_template_index];
	
	for (Size i_template = 1; i_template <= template_clusterID_.size(); ++i_template) {
		if ( cluster_id == template_clusterID_[i_template] ) {
			template_index_icluster.push_back(i_template);
			templates_icluster.push_back(templates_[i_template]);
			weights_icluster.push_back(template_weights_[i_template]);
			template_chunks_icluster.push_back(template_chunks_[i_template]);
			template_contigs_icluster.push_back(template_contigs_[i_template]);
			
			if (i_template == initial_template_index) {
				initial_template_index_icluster = template_index_icluster.size();
			}
		}
	}
}
	
void HybridizeProtocol::apply( core::pose::Pose & pose )
{
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;
	
	// pick starting template
	core::Size initial_template_index;
	core::Size initial_template_index_icluster; // index in the cluster
	utility::vector1 < core::Size > template_index_icluster; // index look back
	utility::vector1 < core::pose::PoseOP > templates_icluster;
	utility::vector1 < core::Real > weights_icluster;
	utility::vector1 < protocols::loops::Loops > template_chunks_icluster;
	utility::vector1 < protocols::loops::Loops > template_contigs_icluster;
	pick_starting_template(initial_template_index, initial_template_index_icluster, template_index_icluster, templates_icluster, weights_icluster, template_chunks_icluster, template_contigs_icluster);

    using namespace ObjexxFCL::fmt;
    TR << "Using initial template: " << I(4,initial_template_index) << " " << template_fn_[initial_template_index] << std::endl;

	// initialize template history
	TemplateHistoryOP history = new TemplateHistory(pose);
	history->setall( initial_template_index );
	pose.data().set( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY, history );
	
	// apply constraints
	core::scoring::constraints::ConstraintSetOP constraint_set;
	std::string cst_fn = template_cst_fn_[initial_template_index];
	if (!cst_fn.empty()) {
		constraint_set = ConstraintIO::get_instance()->read_constraints_new( cst_fn, new ConstraintSet, pose );
	}
	pose.constraint_set( constraint_set );
	
	// fold tree hybridize
	utility::vector1 < core::fragment::FragSetOP > frag_libs;
	frag_libs.push_back(fragments3_);
	frag_libs.push_back(fragments9_);
	core::scoring::ScoreFunctionOP scorefxn_stage1 = 
		core::scoring::ScoreFunctionFactory::create_score_function(option[cm::hybridize::stage1_weights](), option[cm::hybridize::stage1_patch]());
	FoldTreeHybridizeOP ft_hybridize(
		new FoldTreeHybridize(initial_template_index_icluster, templates_icluster, template_chunks_icluster, template_contigs_icluster,  frag_libs) ) ;
	ft_hybridize->set_scorefunction(scorefxn_stage1);
	ft_hybridize->apply(pose);
	
	//write gdtmm to output
	using namespace ObjexxFCL::fmt;
	core::Real gdtmm = get_gdtmm(pose);
	core::pose::setPoseExtraScores( pose, "GDTMM_after_stage1", gdtmm);
	TR << "GDTMM_after_stage1" << F(8,3,gdtmm) << std::endl;

	// cartesian fragment hybridize
	pose.constraint_set( constraint_set );  //reset constraints
	CartesianHybridizeOP cart_hybridize (
		new CartesianHybridize( templates_icluster, weights_icluster,template_chunks_icluster,template_contigs_icluster, fragments9_, fragments3_ ) );
	core::scoring::ScoreFunctionOP scorefxn_stage2 =
	core::scoring::ScoreFunctionFactory::create_score_function(option[cm::hybridize::stage2_weights](), option[cm::hybridize::stage2_patch]());
	cart_hybridize->set_scorefunction(scorefxn_stage2);
	cart_hybridize->apply(pose);
    //fragment_history = cart_hybridize->get_fragment_history(); // fragment indexed within the cluster

	// get fragment history
	runtime_assert( pose.data().has( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
	history = *( static_cast< TemplateHistory* >( pose.data().get_ptr( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY )() ));

	//write gdtmm to output
	gdtmm = get_gdtmm(pose);
	core::pose::setPoseExtraScores( pose, "GDTMM_after_stage2", gdtmm);
	TR << "GDTMM_after_stage2" << F(8,3,gdtmm) << std::endl;
	
	//look back and apply constraints
	TR << "History :";
    for (int i=1; i<= history->size(); ++i ) {
		TR << I(4,i);
	}
	TR << std::endl;
	TR << "History :";
	for (int i=1; i<= history->size(); ++i ) {
		TR << I(4, history->get(i));
	}
	TR << std::endl;

	(*scorefxn_stage2)(pose);

	// optional relax
}


protocols::moves::MoverOP HybridizeProtocol::clone() const { return new HybridizeProtocol( *this ); }
protocols::moves::MoverOP HybridizeProtocol::fresh_instance() const { return new HybridizeProtocol; }

std::string
HybridizeProtocol::get_name() const {
	return "HybridizeProtocol";
}

} // hybridize 
} // comparative_modeling 
} // protocols

