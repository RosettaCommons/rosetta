// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/SeqprofConsensusOperation.cc
/// @brief set designable residues to those observed in sequence profile
/// @author Florian Richter, floric@u.washington.edu, april 2011


// Unit Headers
#include <protocols/toolbox/task_operations/SeqprofConsensusOperation.hh>
#include <protocols/toolbox/task_operations/SeqprofConsensusOperationCreator.hh>
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToAlignedSegments.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/sequence/SequenceProfile.hh>

#include <core/io/ddg/PositionDdGInfo.hh>

// Utility Headers
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/vector0.hh>
#include <boost/foreach.hpp>


static THREAD_LOCAL basic::Tracer tr( "protocols.toolbox.task_operations.SeqprofConsensusOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

core::pack::task::operation::TaskOperationOP
SeqprofConsensusOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new SeqprofConsensusOperation );
}


/// @brief default constructor
SeqprofConsensusOperation::SeqprofConsensusOperation():
	TaskOperation(),
	seqprof_(/* NULL */),
	min_aa_probability_(0.0),
	prob_larger_current_(true),
	ignore_pose_profile_length_mismatch_(false),
	convert_scores_to_probabilities_( true ),
	restrict_to_aligned_segments_( /* NULL */ ),
	conservation_cutoff_aligned_segments_( -100000 ),
	protein_interface_design_( /* NULL */ ),
	conservation_cutoff_protein_interface_design_( -100000 ),
	debug_( false ),
	keep_native_(false),
	chain_num_(1),
	restrict_to_repacking_(true)

{
	if ( basic::options::option[ basic::options::OptionKeys::in::file::pssm ].user() ) {
		seqprof_filename_ = basic::options::option[ basic::options::OptionKeys::in::file::pssm ][1];
	} else {
		seqprof_filename_ = "";
	}
	if ( utility::file::file_exists( seqprof_filename_ ) ) {
		core::sequence::SequenceProfileOP seqprof( new core::sequence::SequenceProfile( seqprof_filename_ ) );
		seqprof->convert_profile_to_probs(); // was previously implicit in from-filename constructor
		seqprof_ = seqprof;
	}
}


/// @brief destructor
SeqprofConsensusOperation::~SeqprofConsensusOperation() {}

/// @brief clone
core::pack::task::operation::TaskOperationOP
SeqprofConsensusOperation::clone() const {
	return core::pack::task::operation::TaskOperationOP( new SeqprofConsensusOperation( *this ) );
}

/// @brief all AA that have a higher probability in the seqprofile
/// than the native residue are allowed. probability also
/// needs to be higher than min_aa_probability_
/// @details NOTE ON SYMMETRIC POSE BEHAVIOR:
/// pssm files are usually for one chain only, therefore
/// this task operation will only set the residue behavior for
/// the first chain/asymetric unit.
/// it could be possible to handle the symmetry setup here, i.e.
/// set up the residue level task for every symmetric copy, but
/// it's prolly better to let the symmetry machinery deal with that
/// mode of packer task symmetrization should be intersection
void
SeqprofConsensusOperation::apply( Pose const & pose, PackerTask & task ) const
{
	using namespace core::scoring::constraints;
	using namespace core::sequence;

	SequenceProfileOP seqprof = seqprof_;
	if ( !seqprof ) {
		seqprof = SequenceProfileOP( new core::sequence::SequenceProfile );
		tr<<"Sequence profile was not set until now. Attempting to read sequence profile from the pose's sequenceprofile constraints..."<<std::endl;

		core::pose::PoseOP chain;
		ConstraintCOPs constraints;
		if ( chain_num_ ) {
			chain = pose.split_by_chain( chain_num_ );
			constraints = chain->constraint_set()->get_all_constraints();
			tr<<"total number of residues in chain:"<<chain->total_residue()<<std::endl;
			tr<<"Total number of constraints in pose: "<<constraints.size()<<std::endl;
		} else {
			constraints = pose.constraint_set()->get_all_constraints();
			tr<<"total number of residues in pose:"<<pose.total_residue()<<std::endl;
			tr<<"Total number of constraints in pose: "<<constraints.size()<<std::endl;
		}

		core::Size cst_num( 0 );
		BOOST_FOREACH ( ConstraintCOP const c, constraints ) {
			if ( c->type() == "SequenceProfile" ) {
				SequenceProfileConstraintCOP seqprof_cst( utility::pointer::dynamic_pointer_cast< core::scoring::constraints::SequenceProfileConstraint const > ( c ) );
				runtime_assert( seqprof_cst != 0 );
				core::Size const seqpos( seqprof_cst->seqpos() );
				SequenceProfileCOP seqprof_pos( seqprof_cst->sequence_profile() );
				seqprof->prof_row( seqprof_pos->profile()[ seqpos ], seqpos );
				if ( basic::options::option[ basic::options::OptionKeys::out::file::use_occurrence_data ].value() ) {

					seqprof->probabilty_row( seqprof_pos->occurrence_data()[ seqpos ], seqpos );
					tr<<"Testing occurence_data for pos "<<seqpos<<std::endl;
					tr<<seqprof_pos->occurrence_data()[ seqpos ]<<std::endl;

				}
				cst_num++;
			}
		}
		tr<<"Added "<<cst_num<<" sequence profile positions to seqprof, taken from the pose's SequenceProfile constraints"<<std::endl;
	}
	if ( !seqprof ) {
		utility_exit_with_message("No sequence profile set. option -in:file:pssm not specified? no filename in tag specified? Sequence profile constraints not added to pose by other movers/filters?");
	}

	core::Size asymmetric_unit_res( pose.total_residue() );
	core::Size last_res( 0 );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		last_res = asymmetric_unit_res <= seqprof->profile().size() ? pose.total_residue() : seqprof->profile().size();
		tr<<" Pose is SYMMETRIC!!!"<<std::endl;
		core::conformation::symmetry::SymmetricConformation const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		asymmetric_unit_res = SymmConf.Symmetry_Info()->num_independent_residues();
		//  task.request_symmetrize_by_intersection();
	} else {
		tr<<"Pose is NOT--SYMMETRIC!!!"<<std::endl;
		last_res = asymmetric_unit_res <= seqprof->profile().size() ? pose.total_residue() : seqprof->profile().size() - 1 /*seqprof has size n+1 compared to its real contents; heaven knows why...*/;
	}
	tr<< "the size of sequence profile is: "<<last_res<<std::endl;
	/// following paragraph determines where PIDO and RestrictToAlignedInterface are defined.
	/// These are used in the following loop to restrict conservation profiles differently
	utility::vector1< core::Size > designable_interface, designable_aligned_segments;
	designable_interface.clear(); designable_aligned_segments.clear();
	if ( protein_interface_design() != NULL ) {
		core::pack::task::TaskFactoryOP temp_tf( new core::pack::task::TaskFactory );
		temp_tf->push_back( protein_interface_design_ );
		designable_interface = protocols::rosetta_scripts::residue_packer_states( pose, temp_tf, true/*designable*/, false/*packable*/ );
	}
	if ( restrict_to_aligned_segments() != NULL ) {
		core::pack::task::TaskFactoryOP temp_tf( new core::pack::task::TaskFactory );
		temp_tf->push_back( restrict_to_aligned_segments_ );
		designable_aligned_segments = protocols::rosetta_scripts::residue_packer_states( pose, temp_tf, true/*designable*/, false/*packable*/ );
	}
	if ( keep_native_ ) {
		tr<<"Adding native identities to allowed identities"<<std::endl;
	}
	tr<<"Allowing the following identities:"<<std::endl;

	core::Size const resi_begin = ( chain_num_ == 0 ? 1 : pose.conformation().chain_begin( chain_num_ ) );
	core::Size const resi_end   = ( chain_num_ == 0 ? pose.total_residue() : pose.conformation().chain_end( chain_num_ ) );

	runtime_assert( (seqprof->profile()).size()>=resi_end - resi_begin );

	for ( core::Size i = resi_begin; i <= resi_end; ++i ) {

		core::Real position_min_prob = min_aa_probability_;
		if ( protein_interface_design() != NULL && std::find( designable_interface.begin(), designable_interface.end(), i ) != designable_interface.end() ) {
			position_min_prob = conservation_cutoff_protein_interface_design();
		} else if ( restrict_to_aligned_segments() != NULL && std::find( designable_aligned_segments.begin(), designable_aligned_segments.end(), i ) != designable_aligned_segments.end() ) {
			position_min_prob = conservation_cutoff_aligned_segments();
		}

		if ( debug() ) {
			tr<<"At position "<<i<<", min_probability is: "<<position_min_prob<<std::endl;
		}

		if ( !pose.residue_type( i ).is_protein() ) continue;
		//std::cout << "SCO at pos " << i << " allows the following residues: ";
		utility::vector1< Real > const & pos_profile( (seqprof->profile())[ i - resi_begin + 1 ] );
		//actual aa proballities and not pssm scores

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		runtime_assert( pose.residue_type( i ).aa() <= (int) pos_profile.size() );
		core::Real current_prob( pos_profile[ pose.residue_type( i ).aa() ] );
		core::Real max_prob( -1000.0 ); // at this position, what is the maximal probability for a residue? these identities should always be allowed in design
		for ( core::Size aa = core::chemical::aa_ala; aa <= core::chemical::num_canonical_aas; ++aa ) {
			if ( max_prob <= pos_profile[ aa ] ) {
				max_prob = pos_profile[ aa ];
			}
		}
		for ( core::Size aa = core::chemical::aa_ala; aa <= core::chemical::num_canonical_aas; ++aa ) {
			core::Real prob( pos_profile[ aa ] );
			if ( prob >= position_min_prob || prob >= max_prob ) {
				if ( prob_larger_current_ ) {
					if ( prob >= current_prob ) keep_aas[ aa ] = true;
				} else keep_aas[ aa ] = true;
				//std::cout << " " << static_cast<core::chemical::AA>(aa) << " prob=" << prob << ", ";
			}
			if ( keep_native_ ) {
				keep_aas[ pose.residue_type(i).aa() ] = true;
			}
			if ( keep_aas[ aa ] && debug() ) {
				tr<<core::chemical::oneletter_code_from_aa( static_cast< core::chemical::AA >( aa ) );
			}
		}
		if ( debug() ) {
			tr<<std::endl;
		}
		//std::cout << " native " << pose.residue_type(i).aa() << " prob=" << native_prob << "." << std::endl;

		task.nonconst_residue_task(i).restrict_absent_canonical_aas( keep_aas );

	} //loop over all residues for which profile information exists

	if ( basic::options::option[ basic::options::OptionKeys::out::file::use_occurrence_data ].value() ) {
		for ( core::Size i = resi_begin; i <= resi_end; ++i ) {
			//for all non interface reisdues we allow only reidues that acctualy appear in native proteins.
			if ( protein_interface_design() != NULL && std::find( designable_interface.begin(), designable_interface.end(), i ) != designable_interface.end() ) {
				tr <<"Residue "<<i<<" is in the interface, not applying occurrence data"<<std::endl;
				continue;
			}
			utility::vector1< Real > const & pos_probability( (seqprof->occurrence_data())[ i - resi_begin + 1 ] );
			utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );

			// tr<<"The size for the probabilty vector is:"<<pos_probability.size()<<std::endl;

			if ( !pose.residue_type( i ).is_protein() ) continue;
			runtime_assert( pose.residue_type( i ).aa() <= (int) pos_probability.size() );
			bool has_probabilty=false;
			for ( core::Size aa = core::chemical::aa_ala; aa <= core::chemical::num_canonical_aas; ++aa ) {
				if ( pos_probability[ aa ] >0 ) {
					has_probabilty=true;
				}
			}
			if ( has_probabilty ) {
				for ( core::Size aa = core::chemical::aa_ala; aa <= core::chemical::num_canonical_aas; ++aa ) {
					core::Real prob( pos_probability[ aa ] );
					if ( prob >0 ) {
						keep_aas[ aa ] = true;
					}
					if ( keep_aas[ aa ] && debug() ) {
						tr<<core::chemical::oneletter_code_from_aa( static_cast< core::chemical::AA >( aa ) );
					}
				}
			} else { //if profile doesn't have probability data we revert to using pssm data
				core::Real position_min_prob = min_aa_probability_;
				if ( restrict_to_aligned_segments() != NULL && std::find( designable_aligned_segments.begin(), designable_aligned_segments.end(), i ) != designable_aligned_segments.end() ) {
					position_min_prob = conservation_cutoff_aligned_segments();
				}
				utility::vector1< Real > const & pos_profile( (seqprof->profile())[ i - resi_begin + 1 ] );

				for ( core::Size aa = core::chemical::aa_ala; aa <= core::chemical::num_canonical_aas; ++aa ) {
					core::Real prob( pos_profile[ aa ] );
					if ( prob >= position_min_prob  ) {
						keep_aas[ aa ] = true;
					}
					if ( keep_aas[ aa ] && debug() ) {
						tr<<core::chemical::oneletter_code_from_aa( static_cast< core::chemical::AA >( aa ) );
					}
				}
			}
			if ( debug() ) {
				tr<<std::endl;
			}
			//std::cout << " native " << pose.residue_type(i).aa() << " prob=" << native_prob << "." << std::endl;

			task.nonconst_residue_task(i).restrict_absent_canonical_aas( keep_aas );

		} //loop over all residues for which profile information exists
	}
	bool prot_res_without_profile_information_exist(false);
	for ( core::Size i = last_res + 1; i <= asymmetric_unit_res; ++i ) {
		if ( restrict_to_repacking_ ) { // Gideon Oct14, For cases where we want residues that don't have sequence constraints to design
			task.nonconst_residue_task(i).restrict_to_repacking();
		}
		if ( pose.residue_type( i ).is_protein() ) prot_res_without_profile_information_exist = true;
	}

	if ( prot_res_without_profile_information_exist ) {
		if ( ignore_pose_profile_length_mismatch_ ) tr<< tr.bgRed << "WARNING WARNING: the passed in pose is longer than the sequence profile specified. Double check whether the used sequence profile is correct. Setting every excess pose residue to repacking."<<tr.Reset<<std::endl;

		else utility_exit_with_message("The passed in pose is longer than the sequence profile specified. Double check whether the used sequence profile is correct.");
	}
} // apply

void
SeqprofConsensusOperation::parse_tag( TagCOP tag , DataMap & datamap )
{
	restrict_to_repacking_=tag->getOption< bool >("restrict_to_repacking", 1  );
	convert_scores_to_probabilities( tag->getOption< bool >("convert_scores_to_probabilities", 1  ) );
	if ( tag->hasOption("filename") ) {
		seqprof_filename_ = tag->getOption< String >( "filename" );
		tr<<"Loading seqprof from a file named: "<<seqprof_filename_<<std::endl;
		core::sequence::SequenceProfileOP seqprof( new core::sequence::SequenceProfile( seqprof_filename_ ) );
		if ( convert_scores_to_probabilities() ) {
			seqprof->convert_profile_to_probs(); // was previously implicit in from-filename constructor
		}
		seqprof_ = seqprof;
	} else {
		tr<<"Seqprof not loaded. Expecting another mover/filter to provide a sequence profile..."<<std::endl;
	}
	if ( tag->hasOption("min_aa_probability") ) min_aa_probability_ = tag->getOption< Real >("min_aa_probability" );
	if ( tag->hasOption("keep_native") ) keep_native_ = tag->getOption< Real >("keep_native",false );//Include native aa identity in allowed identities even if below min_aa_prob.
	chain_num_=tag->getOption< Size >("chain_num",1 );
	if ( tag->hasOption("probability_larger_than_current") ) prob_larger_current_ = tag->getOption< bool >("probability_larger_than_current");

	if ( tag->hasOption("ignore_pose_profile_length_mismatch") ) ignore_pose_profile_length_mismatch_ = tag->getOption< bool >("ignore_pose_profile_length_mismatch");

	utility::vector1< TagCOP > const sub_tags( tag->getTags() );
	BOOST_FOREACH ( TagCOP const sub_tag, sub_tags ) {
		if ( sub_tag->getName() == "RestrictToAlignedSegments" ) {
			restrict_to_aligned_segments_ = RestrictToAlignedSegmentsOperationOP( new RestrictToAlignedSegmentsOperation );
			tr<<"Within SeqprofConsensus I'm now reading a RestrictToAlignedSegments operation..."<<std::endl;
			restrict_to_aligned_segments_->parse_tag( sub_tag, datamap );
			conservation_cutoff_aligned_segments( tag->getOption< core::Real >( "conservation_cutoff_aligned_segments" ) );
			tr<<"conservation cutoff for aligned segments: "<<conservation_cutoff_aligned_segments()<<std::endl;
		} else if ( sub_tag->getName() == "ProteinInterfaceDesign" ) {
			protein_interface_design_ = ProteinInterfaceDesignOperationOP( new ProteinInterfaceDesignOperation );
			tr<<"Within SeqprofConsensus I'm now reading a ProteinInterfaceDesign operation..."<<std::endl;
			protein_interface_design_->parse_tag( sub_tag, datamap );
			conservation_cutoff_protein_interface_design( tag->getOption< core::Real >( "conservation_cutoff_protein_interface_design" ) );
			tr<<"conservation cutoff for protein interface design: "<<conservation_cutoff_protein_interface_design()<<std::endl;
		} else {
			utility_exit_with_message( "SeqprofConsensus subtag not recognized: " + sub_tag->getName() );
		}
	}//foreach
	if ( !(conservation_cutoff_protein_interface_design() <= conservation_cutoff_aligned_segments() ) ) {
		tr<<tr.bgRed<<"WARNING! conservation_cutoff_protein_interface_design() > conservation_cutoff_aligned_segments() "<<tr.Reset<<std::endl;
	}

	if ( !( conservation_cutoff_aligned_segments() <= min_aa_probability_ ) ) {
		tr<<tr.bgRed<<"WARNING!conservation_cutoff_aligned_segments()>min_aa_probability_ )"<<tr.Reset<<std::endl;
	}
	debug( tag->getOption< bool >( "debug", false ) );
}

core::sequence::SequenceProfileCOP
SeqprofConsensusOperation::seqprof() const
{
	return seqprof_;
}

void
SeqprofConsensusOperation::set_seqprof( core::sequence::SequenceProfileOP seqprof, bool reweight )
{
	if ( reweight ) {
		core::sequence::SequenceProfileOP reweightedprof( new core::sequence::SequenceProfile( *seqprof ) );
		reweightedprof->convert_profile_to_probs(); // was previously implicit in from-filename constructor
		seqprof_ = reweightedprof;
	} else {
		seqprof_ = seqprof;
	}
}

core::pack::task::operation::TaskOperationOP
RestrictConservedLowDdgOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictConservedLowDdgOperation );
}

void
SeqprofConsensusOperation::restrict_to_aligned_segments( RestrictToAlignedSegmentsOperationOP rtas ){ restrict_to_aligned_segments_ = rtas; }

RestrictToAlignedSegmentsOperationOP
SeqprofConsensusOperation::restrict_to_aligned_segments() const{ return restrict_to_aligned_segments_; }

ProteinInterfaceDesignOperationOP
SeqprofConsensusOperation::protein_interface_design() const{ return protein_interface_design_; }

void
SeqprofConsensusOperation::protein_interface_design( ProteinInterfaceDesignOperationOP pido ){
	protein_interface_design_ = pido;
}

RestrictConservedLowDdgOperation::RestrictConservedLowDdgOperation()
: Parent(),
	ddG_predictions_filename_(basic::options::option[ basic::options::OptionKeys::in::file::ddg_predictions_file ].value()),
	conservation_cutoff_(0.6),
	ddG_cutoff_(1.5),
	verbose_(false)
{
	position_ddGs_.clear();
	if ( utility::file::file_exists( ddG_predictions_filename_ ) ) {
		position_ddGs_ = core::io::PositionDdGInfo::read_ddg_predictions_file( ddG_predictions_filename_ );
	}
}

RestrictConservedLowDdgOperation::~RestrictConservedLowDdgOperation()
{}

core::pack::task::operation::TaskOperationOP
RestrictConservedLowDdgOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictConservedLowDdgOperation( *this ) );
}

void
RestrictConservedLowDdgOperation::parse_tag( TagCOP tag , DataMap & datamap )
{
	Parent::parse_tag( tag, datamap );
	if ( tag->hasOption("ddG_filename") ) {
		ddG_predictions_filename_ = tag->getOption< std::string >("ddG_filename" );
		position_ddGs_ = core::io::PositionDdGInfo::read_ddg_predictions_file( ddG_predictions_filename_ );
	}

	if ( tag->hasOption("conservation_cutoff") ) conservation_cutoff_ = tag->getOption< Real >("conservation_cutoff" );
	if ( tag->hasOption("ddG_cutoff") ) ddG_cutoff_ = tag->getOption< Real >("ddG_cutoff" );
	if ( tag->hasOption("verbose") ) verbose_ = tag->getOption< bool >("verbose" );
}

void
RestrictConservedLowDdgOperation::apply(
	Pose const & pose,
	PackerTask & task
) const
{
	if ( !this->seqprof() ) utility_exit_with_message("No sequence profile set. option -in:file:pssm not specified? no filename in tag specified?");

	if ( position_ddGs_.size() == 0 ) utility_exit_with_message("No ddG infos were read in. option -in:file:ddg_predictions_file not specified? no filename in tag specified?");

	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {

		if ( !pose.residue_type( i ).is_protein() ) continue;
		core::chemical::AA seqprof_wt_aa( this->seqprof_wt_aa( i ) );

		if ( position_untouchable( i, seqprof_wt_aa ) ) {
			if ( seqprof_wt_aa == pose.residue_type(i).aa() ) task.nonconst_residue_task( i ).restrict_to_repacking();
			else {
				utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
				keep_aas[ seqprof_wt_aa ] = true;
				keep_aas[ pose.residue_type(i).aa() ] = true;
				task.nonconst_residue_task(i).restrict_absent_canonical_aas( keep_aas );
			}
		} // if untouchable
	} //loop over all residues
}

bool
RestrictConservedLowDdgOperation::position_untouchable(
	core::Size seqpos,
	core::chemical::AA seqprof_wt
) const
{

	//note: first we deal with the alanine special case
	//obviousluy there is no ddG associated with mutating Ala to Ala,
	//so for alanine residues we return true through the conservation
	//criterion alone, the rationale being that conserved alanines
	//are probably important structurally
	if ( seqprof_wt == core::chemical::aa_ala ) {
		if ( (seqprof()->profile())[ seqpos ][ seqprof_wt ] > conservation_cutoff_ ) return true;
		else return false;
	}

	std::map< core::Size, core::io::PositionDdGInfo::PositionDdGInfoOP >::const_iterator posddg_it( position_ddGs_.find( seqpos ) );

	//note: if the position wasn't found, this could mean that the predictions file
	//was incomplete or that the ddG protocol couldn't calculate a proper ddG,
	//which is the case for disulfides and some modified residue types such as
	//phospho-ser etc. so let's spit out a warning and just apply the conservation_cutoff_
	if ( posddg_it == position_ddGs_.end() ) {

		tr << "Warning: no ddG information read for sequence position " << seqpos << ". This could either mean that the ddG predictions input file is incomplete or that the original PDB had a disulfide cys or other modified residue at this position. Decision whether residue is untouchable will be made based on sequence conservation alone." << std::endl;
		if ( (seqprof()->profile())[ seqpos ][ seqprof_wt ] > conservation_cutoff_ ) return true;
		else return false;
	}

	core::io::PositionDdGInfo::PositionDdGInfo const & pos_ddg( *(posddg_it->second) );
	if ( seqprof_wt != pos_ddg.wt_aa() ) utility_exit_with_message ("The wildtype aa for position "+utility::to_string( seqpos ) + " is different in the ddG file and the pssm file. Something's unclean somewhere." );
	std::map< core::chemical::AA, core::Real >::const_iterator ala_it( pos_ddg.mutation_ddGs().find( core::chemical::aa_ala ) );
	if ( ala_it ==  pos_ddg.mutation_ddGs().end() ) utility_exit_with_message("The ddG of mutating to Ala was not found for position "+utility::to_string( seqpos )+" in file "+ddG_predictions_filename_ + ".");

	if ( (ala_it->second > ddG_cutoff_) && ( (seqprof()->profile())[ seqpos ][ seqprof_wt ] > conservation_cutoff_) ) {
		if ( verbose_ ) tr << "Pos " << seqprof_wt << seqpos << " has ddG_cutoff of " << ala_it->second << " and profile frequency of " << (seqprof()->profile())[ seqpos ][ seqprof_wt ] << ", considered untouchable." << std::endl;
		return true;
	}
	return false;
}

core::chemical::AA
RestrictConservedLowDdgOperation::seqprof_wt_aa( core::Size seqpos ) const
{
	return core::chemical::aa_from_oneletter_code( (*(this->seqprof()))[seqpos] );
}

core::Real
RestrictConservedLowDdgOperation::position_ala_ddG( core::Size seqpos ) const
{
	if ( seqprof_wt_aa( seqpos ) == core::chemical::aa_ala ) return 0.0;

	std::map< core::Size, core::io::PositionDdGInfo::PositionDdGInfoOP >::const_iterator posddg_it( position_ddGs_.find( seqpos ) );
	if ( posddg_it == position_ddGs_.end() ) utility_exit_with_message("no ddg information read for sequence position "+ utility::to_string( seqpos ) );
	core::io::PositionDdGInfo::PositionDdGInfo const & pos_ddg( *(posddg_it->second) );
	std::map< core::chemical::AA, core::Real >::const_iterator ala_it( pos_ddg.mutation_ddGs().find( core::chemical::aa_ala ) );
	if ( ala_it ==  pos_ddg.mutation_ddGs().end() ) utility_exit_with_message("The ddG of mutating to Ala was not found for position "+utility::to_string( seqpos )+" in file "+ddG_predictions_filename_ + ".");
	return ala_it->second;
}

} // TaskOperations
} // toolbox
} // protocols

