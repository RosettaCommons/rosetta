// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/JointSequenceOperation.cc
/// @brief set designable residues to those observed in a set of structures
/// @author Rocco Moretti, rmoretti@u.washington.edu


// Unit Headers
#include <protocols/toolbox/task_operations/JointSequenceOperation.hh>
#include <protocols/toolbox/task_operations/JointSequenceOperationCreator.hh>

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <protocols/jd2/util.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <utility/vector1.hh>
#include <string>

#include <utility/vector0.hh>


static basic::Tracer TR( "protocols.toolbox.task_operations.JointSequenceOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

TaskOperationOP
JointSequenceOperationCreator::create_task_operation() const
{
	return TaskOperationOP( new JointSequenceOperation );
}

void JointSequenceOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	JointSequenceOperation::provide_xml_schema( xsd );
}

std::string JointSequenceOperationCreator::keyname() const
{
	return JointSequenceOperation::keyname();
}

/// @brief default constructor
JointSequenceOperation::JointSequenceOperation():
	TaskOperation(),
	use_current_pose_(true),
	use_natro_(false),
	use_fasta_(false),
	ubr_(/* 0 */),
	chain_(0)
{
}

/// @brief destructor
JointSequenceOperation::~JointSequenceOperation() {}

/// @brief clone
TaskOperationOP
JointSequenceOperation::clone() const {
	return TaskOperationOP( new JointSequenceOperation( *this ) );
}

/// @brief all AA that have a higher probability in the seqprofile
/// than the native residue are allowed. probability also
/// needs to be higher than min_aa_probability_
void
JointSequenceOperation::apply( Pose const & pose, PackerTask & task ) const
{

	core::conformation::symmetry::SymmetryInfoCOP syminfo = NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) syminfo = core::pose::symmetry::symmetry_info(pose);

	core::Size start;
	core::Size end;

	if ( chain_ > 0 ) {
		//core::pose::Pose new_pose;
		//new_pose = pose.split_by_chain( chain_ );
		//pose = new_pose;
		start = pose.conformation().chain_begin( chain_ );
		end =  pose.conformation().chain_end( chain_ );
	} else {
		start = 1;
		end = pose.size();
	}

	//TR << "CHAIN: " << chain_ << std::endl;
	//TR << "START: " << start << std::endl;
	//TR << "END: " << end << std::endl;

	core::Size seq_length = end - start + 1;


	// Iter through native sequence
	for ( std::vector<core::sequence::SequenceOP>::const_iterator iter(sequences_.begin()), end = sequences_.end(); iter != end; ++iter ) {
		//TR << "it " << **iter << " " << (*iter)->length() << " vs " << seq_length << std::endl;
		if ( (*iter)->length() != seq_length ) {
			std::string name("current pdb");
			if ( pose.pdb_info() ) {
				name = pose.pdb_info()->name();
			}
			TR.Warning << "considered sequence " << (*iter)->id() << " contains a different number of residues than " << name << std::endl;
		}
	}

	// Inter through current pose sequence
	for ( core::Size ii = start; ii <= end; ++ii ) {
		if ( core::pose::symmetry::is_symmetric(pose) && !syminfo->chi_is_independent(ii) ) continue;
		if ( !pose.residue_type( ii ).is_protein() ) continue;

		int na_ii = ii - start + 1;
		//NK120531 -- These lines commented out bc they are useless and crash for symmetric systems.
		//NK120531    std::vector<core::sequence::SequenceOP>::const_iterator iter(sequences_.begin());
		//NK120531    char aa( (*(*iter))[ na_ii ] ); // aa char of native pose

		utility::vector1< bool > allowed(core::chemical::num_canonical_aas, false);

		if ( use_current_pose_ ) {
			if ( pose.aa(ii) <= static_cast<int>(allowed.size()) ) allowed[ pose.aa(ii) ] = true;
		}
		for ( std::vector<core::sequence::SequenceOP>::const_iterator iter(sequences_.begin()), end = sequences_.end(); iter != end; ++iter ) {
			//if ( ii > (*iter)->length() ) continue; // ignore short references
			char aa( (*(*iter))[ na_ii ] );
			if ( core::chemical::oneletter_code_specifies_aa(aa) ) {
				if ( core::chemical::aa_from_oneletter_code(aa)<=static_cast<int>(allowed.size()) ) {
					//TR << "JointSeq: " << pose.aa(ii) << ii << " " << core::chemical::aa_from_oneletter_code(aa) << na_ii << std::endl;
					allowed[ core::chemical::aa_from_oneletter_code(aa)  ] = true;
				}
			}
		}
		task.nonconst_residue_task(ii).restrict_absent_canonical_aas( allowed );
	} //loop over all residues

	if ( use_natro_ && ubr_ ) {
		task.append_rotamerset_operation( ubr_ );
	}
} // apply

void
JointSequenceOperation::parse_tag( TagCOP tag , DataMap & )
{
	use_current_pose( tag->getOption< bool >( "use_current", true ) );
	use_natro( tag->getOption< bool >( "use_natro", false ) );

	// specifiy a chain, if 0 use all
	chain_ = tag->getOption< core::Size >("chain", 0 );

	// if a pdb file is supplied under native, which chain should be used to extract the sequence
	use_chain_ = tag->getOption < core::Size > ("use_chain" , 0 );

	if ( tag->getOption< bool >( "use_native", false ) ) {
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			add_native_pdb( basic::options::option[ basic::options::OptionKeys::in::file::native ] );
		} else {
			utility_exit_with_message("Native PDB not specified on command line.");
		}
	}

	if ( tag->getOption< bool >( "use_starting_as_native", false ) ) {
		core::pose::PoseCOP pose (protocols::jd2::get_current_jobs_starting_pose());
		if ( pose ) {
			add_native_pose( pose );
		} else {
			TR.Warning << "Unable to obtain a starting pose for use_starting_as_native" << std::endl;
		}
	}

	if ( tag->hasOption("filename") ) {
		add_pdb( tag->getOption< String >( "filename" ) );
	}
	if ( tag->hasOption("native") ) {
		add_native_pdb( tag->getOption< String >( "native" ) );
	}

	if ( tag->getOption< bool >( "use_fasta", false ) ) {
		// Use native FASTA sequence
		if ( basic::options::option[ basic::options::OptionKeys::in::file::fasta ].user() ) {
			std::string fasta_file( core::sequence::read_fasta_file_str( basic::options::option[ basic::options::OptionKeys::in::file::fasta ]()[1] )[1] );
			add_native_fasta( fasta_file );
		} else {
			utility_exit_with_message("Native FASTA file not specified on command line.");
		}
	}
}

void JointSequenceOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::attribute_w_default(
		"use_current", xsct_rosetta_bool,
		"Use residue identities from the current structure (input pose to apply() of the taskoperation)",
		"true"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"use_natro", xsct_rosetta_bool,
		"the task operation also adds the rotamers from the native structures (use_native/native) in the rotamer library.",
		"false"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"chain", xsct_non_negative_integer,
		"to which chain to apply, 0 is all chains",
		"0"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"use_chain", xsct_non_negative_integer,
		"given an additional input pdb, such as through in:file:native, which chain should the sequence be derived from. 0 is all chains.",
		"0"  )
		// AMW: This relies on the presence of a command-line option.
		+ XMLSchemaAttribute::attribute_w_default(
		"use_native", xsct_rosetta_bool,
		"Use residue identities from the structure listed with -in:file:native",
		"false"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"use_starting_as_native", xsct_rosetta_bool,
		"use starting pose as native pose",
		"false"  )
		+ XMLSchemaAttribute(
		"filename", xs_string ,
		"Use residue identities from the listed file" )
		+ XMLSchemaAttribute( "native", xs_string ,
		"use this file path as the native pdb" )
		// AMW: This relies on the presence of a command-line option.
		+ XMLSchemaAttribute::attribute_w_default(
		"use_fasta", xsct_rosetta_bool,
		"Use residue identities from a native sequence given by a FASTA file (specify the path to the FASTA file with the -in:file:fasta flag at the command line)",
		"false"  );

	task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"Prohibit designing to residue identities that aren't found at that position in any of the listed structures:");
}

/// @brief Add the sequence from the given filename to the set of allowed aas.
void
JointSequenceOperation::add_pdb( std::string filename )
{

	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
	add_pose( pose );
}

/// @brief Add the sequence from the given pose to the set of allowed aas.
void
JointSequenceOperation::add_pose( Pose const & pose )
{
	std::string name("unknown");
	core::pose::Pose new_pose;

	if ( use_chain_ > 0 ) {
		TR << "taking only chain " << use_chain_ << std::endl;
		new_pose = *pose.split_by_chain( use_chain_ ) ;
	} else {
		new_pose = pose;
	}

	if ( new_pose.pdb_info() ) {
		name = new_pose.pdb_info()->name();
	}
	sequences_.push_back( core::sequence::SequenceOP( new core::sequence::Sequence(new_pose.sequence(), name) ) );

}

/// @brief Add the sequence from the given filename to the set of allowed aas
/// and add the rotamers to the set of possible rotamers (if use_natro_ is set)
void
JointSequenceOperation::add_native_pdb( std::string filename ) {
	core::pose::PoseOP poseop( new core::pose::Pose );
	core::import_pose::pose_from_file( *poseop, filename , core::import_pose::PDB_file);
	add_native_pose( poseop );
}

/// @brief Add the sequence from the given pose to the set of allowed aas
/// and add the rotamers to the set of possible rotamers
void
JointSequenceOperation::add_native_pose( core::pose::PoseCOP posecop ){ // PoseCOP needed by UnboundRot, unfortunately
	if ( use_natro_ ) { // Deliberate check now to avoid keeping native poses around if we're never going to need them.
		ubr_->add_pose(posecop);
	}
	add_pose(*posecop);
}

/// @brief Should the current pose (pose supplied to apply) be used in addition to the other ones?
void
JointSequenceOperation::use_current_pose( bool ucp )
{
	use_current_pose_ = ucp;
}

/// @brief Should the rotamers for the native poses be used?
void
JointSequenceOperation::use_natro( bool unr ) {
	if ( !use_natro_ && unr ) {
		ubr_ = core::pack::rotamer_set::UnboundRotamersOperationOP( new core::pack::rotamer_set::UnboundRotamersOperation );
	}
	if ( use_natro_ && !unr ) {
		ubr_.reset(); // Allow owning pointer to garbage collect as necessary.
	}
	use_natro_ = unr;
}

/// @brief which chain should be considered
void
JointSequenceOperation::set_chain( core::Size chain){
	chain_ = chain;
}

/// @brief Read native as a FASTA file
void
JointSequenceOperation::add_native_fasta( std::string fasta_file ) {

	//core::sequence::SequenceOPSequenceOP native_sequence = new Sequence(fasta_file,"fasta",1);
	std::string name("unknown");
	sequences_.push_back( core::sequence::SequenceOP( new core::sequence::Sequence(fasta_file, name) ) );

	//string query_sequence ( read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence() );
	//core::pose::PoseOP poseop(new core::pose::Pose);
	//core::import_pose::pose_from_file( *poseop, filename , core::import_pose::PDB_file);
	//add_native_pose( poseop );

	//using core::sequence::read_fasta_file;
	//string fasta_sequence ( read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence() );
}

} // TaskOperations
} // toolbox
} // protocols
