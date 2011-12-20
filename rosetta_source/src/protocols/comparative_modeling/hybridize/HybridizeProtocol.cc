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

#include <protocols/nonlocal/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

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
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>

// utility
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// evaluation
#include <core/scoring/rms_util.hh>
#include <protocols/evaluation/Align_RmsdEvaluator.hh>
#include <protocols/comparative_modeling/coord_util.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.comparative_modeling.hybridize.HybridizeProtocol" );

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace sequence;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;

HybridizeProtocol::HybridizeProtocol()
{
	//read templates
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	read_template_structures( option[cm::hybridize::templates]() );

	//break templates into chunks
	template_chunks_.clear();
	template_chunks_.resize(templates_.size());
	for (core::Size i_template=1; i_template<=templates_.size(); ++i_template) {
		// find ss chunks in template
		template_chunks_[i_template] = protocols::nonlocal::extract_secondary_structure_chunks(*templates_[i_template]); // add split by residue numbering in this function -ys
	}
	
	//break templates into contigs
	template_contigs_.clear();
	template_contigs_.resize(templates_.size());
	for (core::Size i_template=1; i_template<=templates_.size(); ++i_template) {
		// find ss chunks in template
		// template_contigs_[i_template] = extract_continuous_chunks(*templates_[i_template]); 
	}
	
	//read fragments
	if ( option[ OptionKeys::in::file::frag9 ].user() ) {
		using namespace core::fragment;
		fragments9_ = new ConstantLengthFragSet( 9 );
		fragments9_ = FragmentIO().read_data( option[ OptionKeys::in::file::frag9 ]() );
		for (core::fragment::FrameIterator i = fragments9_->begin(); i != fragments9_->end(); ++i) {
			core::Size position = (*i)->start();
			//library_[position] = **i;
		}
	}
	if ( option[ OptionKeys::in::file::frag3 ].user() ) {
		using namespace core::fragment;
		fragments3_ = new ConstantLengthFragSet( 3 );
		fragments3_ = FragmentIO().read_data( option[ OptionKeys::in::file::frag3 ]() );
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
	
void HybridizeProtocol::apply( core::pose::Pose & pose )
{
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

