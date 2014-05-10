// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentLoopInserter.cc
///
/// @brief
/// @author Tim Jacobs

//Unit
#include <devel/loop_creation/FragmentLoopInserter.hh>
#include <devel/loop_creation/FragmentLoopInserterCreator.hh>

//Numeric
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/BBTorsionAndAnglesSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

//protocols
#include <protocols/loops/Loop.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

//utility
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>


namespace devel {
namespace loop_creation {

static basic::Tracer TR( "protocols.loophash.FragmentLoopInserter" );
static numeric::random::RandomGenerator RG(7235646);

//****CREATOR METHODS****//
std::string
FragmentLoopInserterCreator::keyname() const
{
	return FragmentLoopInserterCreator::mover_name();
}

protocols::moves::MoverOP
FragmentLoopInserterCreator::create_mover() const {
	return new FragmentLoopInserter;
}

std::string
FragmentLoopInserterCreator::mover_name()
{
	return "FragmentLoopInserter";
}
//****END CREATOR METHODS****//

FragmentLoopInserter::FragmentLoopInserter():
num_flanking_residues_to_match_(3),
modify_flanking_regions_(false)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace utility::file;
	using namespace basic::options;
	using namespace core::fragment;

//	utility::vector1<core::Size> frag_sizes( option[ OptionKeys::loops::frag_sizes ].value() );
	FileVectorOption frag_files( option[ OptionKeys::loops::frag_files ] );

//	if( frag_sizes.size() != frag_files.size() ){
//		utility_exit_with_message( "You must specify as many fragment sizes as fragment file names " );
//	}

	for ( Size i = 1; i <= frag_files.size(); ++i ) {
//		Size const frag_size = Size(frag_sizes[i]);

//		FragSetOP frag_lib_op ( new ConstantLengthFragSet( frag_size ) );
//		frag_lib_op = (FragmentIO().read_data( frag_files[i] ));
		frag_sets_.push_back(FragmentIO().read_data( frag_files[i] ));
	}
}

protocols::moves::MoverOP
FragmentLoopInserter::clone() const {
	return( protocols::moves::MoverOP( new FragmentLoopInserter( *this ) ) );
}
protocols::moves::MoverOP
FragmentLoopInserter::fresh_instance() const {
	return protocols::moves::MoverOP( new FragmentLoopInserter );
}

std::string
FragmentLoopInserter::get_name() const {
	return "FragmentLoopInserter";
}

void
FragmentLoopInserter::apply(
	core::pose::Pose & pose
){
	using namespace core;

	//Sanity checks
	if(loop_anchor()<=0 || loop_anchor()>pose.total_residue()-1)
	{
		std::stringstream err;
		err << "Loop anchor " << loop_anchor() << " is invalid" << std::endl;
		utility_exit_with_message(err.str());
	}

	//Remove variants from anchor positions
	pose::remove_upper_terminus_type_from_pose_residue(pose, loop_anchor());
	pose::remove_lower_terminus_type_from_pose_residue(pose, loop_anchor()+1);

	//Check to see if we already have fragments for this anchor, if we do, build a random one,
	//if we don't, try to find some.
	utility::vector1<fragment::FragDataCOP> & low_rms_frags=anchor_frags_[loop_anchor()];

	if(low_rms_frags.size()==0){
		find_loop_fragments(pose);
	}

	build_fragment_loop(pose, low_rms_frags[low_rms_frags.size()]);
	low_rms_frags.pop_back();
}

void
FragmentLoopInserter::find_loop_fragments(
	core::pose::Pose & pose
){
	using namespace core;

	utility::vector1<fragment::FragDataCOP> & low_rms_frags=anchor_frags_[loop_anchor()];

	utility::vector1< numeric::xyzVector<Real> > pose_coords_to_match =
		get_pose_coords_to_match(pose);

	utility::vector1< numeric::xyzVector<Real> > frag_coords_to_match(num_flanking_residues_to_match_*2);
	for(Size frag_sets_index=1; frag_sets_index<=frag_sets_.size(); ++frag_sets_index)
	{
		fragment::FragSetOP cur_frag_set=frag_sets_[frag_sets_index];
		for(fragment::ConstFrameIterator it = cur_frag_set->begin(); it!=cur_frag_set->end(); ++it)
		{
			for(Size i = 1; i<=it->nr_frags(); ++i)
			{
				fragment::AnnotatedFragData const & cur_frag =
					dynamic_cast< const fragment::AnnotatedFragData & > (it->fragment(i));

				//Some fragments randomly don't have coords, should probably use torsion rmsds
				fragment::BBTorsionSRFDCOP test_residue =
					dynamic_cast< const fragment::BBTorsionSRFD* > (cur_frag.get_residue(1)());
				if(test_residue->x()==0 && test_residue->y()==0 && test_residue->z()==0)
				{
					continue;
				}


				//get coordinates for the helical residues on the beginning of the bridge fragment
				for(Size j=1; j<=num_flanking_residues_to_match_; j++)
				{
					Size cur_resnum = j;
					fragment::BBTorsionSRFDCOP fragment_residue =
						dynamic_cast< const fragment::BBTorsionSRFD* > (cur_frag.get_residue(cur_resnum)());

					frag_coords_to_match[j] = numeric::xyzVector<Real>(fragment_residue->x(),
						fragment_residue->y(), fragment_residue->z());
				}

				//get coordinates for the helical residues on the end of the bridge fragment
				for(Size j=1; j<=num_flanking_residues_to_match_; j++)
				{
					Size cur_resnum = cur_frag.size()-(num_flanking_residues_to_match_-j);
					fragment::BBTorsionSRFDCOP fragment_residue =
						dynamic_cast< const fragment::BBTorsionSRFD* > (cur_frag.get_residue(cur_resnum)());

					frag_coords_to_match[num_flanking_residues_to_match_+j] =
						numeric::xyzVector<Real>(fragment_residue->x(), fragment_residue->y(), fragment_residue->z());
				}


				Real rms = numeric::model_quality::calc_rms(
					frag_coords_to_match, pose_coords_to_match );

				if(rms<max_rms_){
					low_rms_frags.push_back(it->fragment_ptr(i));

//					pose::Pose temp_pose;
//					it->fragment_as_pose(i,temp_pose,chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ));

//					core::conformation::Residue new_res = temp_pose.residue(4);
//					core::conformation::Residue const & anchor_res = pose.residue(14);
//					core::conformation::orient_residue_for_ideal_bond(
//						new_res, new_res.lower_connect(), anchor_res, anchor_res.upper_connect(), pose.conformation());
//					pose.conformation().append_residue_by_jump(new_res, 14);

//					std::stringstream filename;
//					filename << "fragment_pdb_" << low_rms_frags.size() << ".pdb";
					//temp_pose.dump_pdb(filename.str());

//					TR << "RMSD" << rms << std::endl;

//					TR << "Fragment coords: " << std::endl;
//					for(Size k=1; k<=frag_coords_to_match.size(); ++k)
//						TR << frag_coords_to_match[k] << std::endl;
				}
			}
		}
	}
	TR << "Found " << low_rms_frags.size() << " fragments for anchor " << loop_anchor()
		<< " meeting RMS threshold of: " << max_rms_ << std::endl;
	if(low_rms_frags.size()==0){
		utility_exit_with_message("No low-rmsd loop fragments found. Try increasing max rmsd");
	}
	numeric::random::random_permutation(low_rms_frags,RG);
}

void
FragmentLoopInserter::build_fragment_loop(
	core::pose::Pose & pose,
	core::fragment::FragDataCOP fragment
){
	using namespace core;

	Size loop_size = fragment->size() - (num_flanking_residues_to_match_*2);

	Size fragment_begin = loop_anchor()-num_flanking_residues_to_match_+1;
	Size fragment_end = loop_anchor()+num_flanking_residues_to_match_;

	Size n_term_append_size=numeric::random::random_range(1, loop_size-1);
	Size c_term_append_size=loop_size-n_term_append_size;
	TR.Debug << "Attaching " << n_term_append_size << " residues to N-termini and "
		<< c_term_append_size << " residues to C-termini" << std::endl;

	//append idealized residues to the pose
	//TODO: use the rotamer id from the backbone db extra data
	chemical::ResidueTypeSetCAP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	for(Size i=0; i<n_term_append_size; ++i){

		Size append_seqpos = loop_anchor()+i;

		char aa_char = fragment->sequence(i+num_flanking_residues_to_match_+1);
		chemical::AA aa = chemical::aa_from_oneletter_code(aa_char);
		chemical::ResidueTypeCOPs const & rsd_type_list( restype_set->aa_map( aa ) );
		if(rsd_type_list.empty()){
			utility_exit_with_message("Could not find residue type for AA: " + std::string(1, aa_char));
		}
		chemical::ResidueType const & rsd_type( *(rsd_type_list[ 1 ]) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );

		TR.Debug << "Attaching residue " << aa_char << " after seqpos " << append_seqpos << std::endl;
		pose.append_polymer_residue_after_seqpos(*new_rsd, append_seqpos, true);
		++fragment_end;

//		std::stringstream filename;
//		filename << "testing_symm_" << i << ".pdb";
//		pose.dump_pdb(filename.str());
	}

	Size prepend_seqpos = loop_anchor()+loop_size-c_term_append_size+1;
	for(Size i=loop_size-1; i>=loop_size-c_term_append_size; --i){

		char aa_char = fragment->sequence(i+num_flanking_residues_to_match_+1);
		chemical::AA aa = chemical::aa_from_oneletter_code(aa_char);
		chemical::ResidueTypeCOPs const & rsd_type_list( restype_set->aa_map( aa ) );
		if(rsd_type_list.empty()){
			utility_exit_with_message("Could not find residue type for AA: " + std::string(1, aa_char));
		}
		chemical::ResidueType const & rsd_type( *(rsd_type_list[ 1 ]) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );

		TR.Debug << "Attaching residue " << aa_char << " before seqpos " << prepend_seqpos << std::endl;
		pose.prepend_polymer_residue_before_seqpos(*new_rsd, prepend_seqpos, true);
		++fragment_end;
	}
//	std::exit(1);
//	pose.dump_pdb("before_torsion_apply.pdb");

	for(Size i=1; i<=fragment->size(); ++i)
	{
		TR << "Fragment residue " << i << ":" << std::endl;
		fragment->get_residue(i)->show(TR);
		TR << std::endl;
	}

	//TEMP
	kinematics::FoldTree ft;

	ft.add_edge(1, fragment_begin-1, kinematics::Edge::PEPTIDE);

	ft.add_edge(1, fragment_begin, 1);
	ft.add_edge(fragment_begin, loop_anchor()+n_term_append_size, kinematics::Edge::PEPTIDE);

	ft.add_edge(1, fragment_end, 2);
	ft.add_edge(fragment_end, loop_anchor()+n_term_append_size+1, kinematics::Edge::PEPTIDE);

	ft.add_edge(1, fragment_end+1, 3);
	ft.add_edge(fragment_end+1, pose.total_residue(), kinematics::Edge::PEPTIDE);

	if(!ft.check_fold_tree())
	{
		utility_exit_with_message("LoophashLoopInserter made a bad fold tree. File a bug!");
	}
	pose.fold_tree(ft);

	fragment->apply(pose, fragment_begin, fragment_end);

//	pack_=true;
//	design_=true;
//	if(pack_ || design_)
//	{
//		protocols::simple_moves::PackRotamersMoverOP repack_mover =
//			new protocols::simple_moves::PackRotamersMover;
//		core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
//
//		core::pack::task::TaskFactoryOP task_factory = new pack::task::TaskFactory;
//
//		core::pack::task::operation::RestrictResidueToRepackingOP repack_res_task =
//			new core::pack::task::operation::RestrictResidueToRepacking();
//
//		core::pack::task::operation::PreventRepackingOP no_repack_res_task =
//			new core::pack::task::operation::PreventRepacking();
//
//		for(core::Size i=1; i<=pose.total_residue(); ++i)
//		{
//			if(i >= fragment_begin && i<= fragment_end)
//			{
//				repack_res_task->include_residue(i);
//			}
//			else
//			{
//				no_repack_res_task->include_residue(i);
//			}
//		}
//		task_factory->push_back(repack_res_task);
//		task_factory->push_back(no_repack_res_task);
//		repack_mover->task_factory(task_factory);
//		repack_mover->apply(pose);
//	}

	protocols::loops::Loop loop(
		loop_anchor()+1,
		loop_anchor()+loop_size,
		loop_anchor()+n_term_append_size
	);

	created_loop_=loop;
	modified_range(loop.start(), loop.stop());
}

utility::vector1< numeric::xyzVector<core::Real> >
FragmentLoopInserter::get_pose_coords_to_match(
	core::pose::Pose const & pose
){
	core::Size fragment_begin = loop_anchor()-num_flanking_residues_to_match_+1;
	core::Size fragment_end = loop_anchor()+num_flanking_residues_to_match_;

	utility::vector1< numeric::xyzVector<core::Real> > pose_coords_to_match(num_flanking_residues_to_match_*2);

	//get coordinates for the helical residues before the jump
	for(core::Size j=0; j<num_flanking_residues_to_match_; ++j)
	{
		pose_coords_to_match[j+1]=
			pose.residue(fragment_begin+j).atom("CA").xyz();
	}

	//get coordinates for the helical residues after the jump
	for(int j=num_flanking_residues_to_match_-1; j>=0; j--)
	{
		numeric::xyzVector<core::Real> res_xyz(pose.residue(fragment_end-j).atom("CA").xyz());

		pose_coords_to_match[num_flanking_residues_to_match_*2-j]=
			pose.residue(fragment_end-j).atom("CA").xyz();
	}

	TR << "Pose coords: " << std::endl;
	for(core::Size k=1; k<=pose_coords_to_match.size(); ++k)
		TR << "Pose coords: " << pose_coords_to_match[k] << std::endl;

	return pose_coords_to_match;
}

void
FragmentLoopInserter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const &
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	parse_loop_anchor(tag);

	//Maximum RMSD of torsion angles to flanking residues
	if(tag->hasOption("max_rms")){
		max_rms_ =
			tag->getOption<core::Real>("max_rms", 2.0);
	}
	else{
		utility_exit_with_message("You must specify the maximum rmsd of flanking regions using max_rms");
	}
}

} //loop creation
} //devel
