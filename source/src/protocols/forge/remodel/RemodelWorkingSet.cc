// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelMover.cc
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

// Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

// project headers
#include <core/chemical/ResidueType.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/OrderedFragSet.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
// for switch typeset

// for resfile command map
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>


#include <protocols/forge/remodel/RemodelWorkingSet.hh>

// for yab managers
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/methods/util.hh>

// C++ Headers


using namespace core;
using namespace basic::options;

namespace protocols {
namespace forge {
namespace remodel {

static THREAD_LOCAL basic::Tracer TR( "protocols.forge.remodel.RemodelWorkingSet" );

RemodelWorkingSet::RemodelWorkingSet()
{
	hasInsertion = false;
	buildDisulfide = false;
}

///
/// @brief
/// copy constructor.
///
RemodelWorkingSet::RemodelWorkingSet(RemodelWorkingSet const & rval):
	loops( rval.loops ),
	translate_index( rval.translate_index ),
	begin( rval.begin ),
	end (rval.end ),
	copy_begin( rval.copy_begin ),
	copy_end (rval.copy_end ),
	src_begin (rval.src_begin),
	src_end (rval.src_end),
	manager( rval.manager ),
	task( rval.task ),
	rvjump_pose( rval.rvjump_pose)
{
	sequence = rval.sequence;
	ss = rval.ss;
	abego = rval.abego;
	hasInsertion = rval.hasInsertion;
}

RemodelWorkingSet & RemodelWorkingSet::operator= ( RemodelWorkingSet const & rval ){

	if ( this != & rval ) {
		loops =  rval.loops ;
		sequence = rval.sequence;
		ss = rval.ss;
		abego = rval.abego;
		translate_index = rval.translate_index;
		begin = rval.begin ;
		end = rval.end;
		copy_begin = rval.copy_begin ;
		copy_end = rval.copy_end ;
		src_begin = rval.src_begin;
		src_end = rval.src_end;
		hasInsertion = rval.hasInsertion;
		manager = rval.manager;
		task = rval.task ;
		rvjump_pose = rval.rvjump_pose;
	}
	return *this;
}

///
/// @brief
/// checks value of option -remodel::generic_aa
///
///
void RemodelWorkingSet::workingSetGen( pose::Pose const & input_pose, protocols::forge::remodel::RemodelData const & data ) {

	using namespace basic::options;
	using namespace protocols;

	//core::util::switch_to_residue_type_set( input_pose, core::chemical::CENTROID );

	// find rebuilding segments
	Size model_length = data.sequence.size();
	//bool NtermExt = false;  // unused ~Labonte
	bool CtermExt = false;
	//bool length_changed = false;

	// copy ss/seq from RemodelData so it can be passed elsewhere later
	ss = data.dssp_updated_ss;
	this->abego = data.abego;

	// this is purely experimental for matching fragment set
	if ( option[OptionKeys::remodel::repeat_structure].user() ) {
		ss.append( ss );
		this->abego.append(this->abego);
	} else {
		//ss.append("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD");
	}

	sequence = data.sequence;

	// find N term extension, if any
	std::string Xs = "xX";
	int first_ext;
	first_ext = (int)data.sequence.find_first_of(Xs);
	if ( first_ext == 0 ) {
		//NtermExt = true;  // unused ~Labonte
		TR << "workingSetGen(): N-terminal is extended" << std::endl;
	}
	// find C term extension, if any
	int last_ext;
	last_ext = (int)data.sequence.find_last_of(Xs);
	if ( last_ext == static_cast<int>(model_length -1) ) {
		TR << "workingSetGen(): C-terminal is extended" << std::endl;
		// if C-term extension, add extra degenerate ss type to get fragments past last res.
		CtermExt = true;
	}

	if ( data.blueprint[ model_length - 1 ].index != data.blueprint[ model_length - 1 ].original_index ) {
		//length_changed = true;  set but never used ~Labonte
		TR << "workingSetGen(): length change found. last blueprint line index: " << data.blueprint[ model_length - 1 ].index
			<< ", last blueprint line original_index: " << data.blueprint[ model_length - 1 ].original_index << std::endl;
	}

	//this is needed for manipulating denovo cases, as they are coded as Cterm
	//extensions.  Affects SegmentRebuld selections.
	if ( option[OptionKeys::remodel::repeat_structure].user() && CtermExt ) {
		model_length = model_length * 2;
	}

	// find all the indices.
	// identify truncation
	// for all the positions that are not extensions "x" or "X", put them in the
	// temp_for_truncation vector -- this corresponds to all the regions need to
	// be copied from the original pdb with deletions included.  fragment_pdb
	// requires a boolean vector that correspond to the length of the original
	// pdb, so initializes a "keep" vector of that size and set values to false.
	// As we iterate over the template_for_truncation, all the positions seen by
	// this iterative step are kept.  Because the truncated_pose is renumbered, we
	// also initializes a translate_index map to link the truncated index to the
	// original index.

	TR << "Adding lines to temp_for_truncation vector." << std::endl;

	std::vector< forge::remodel::LineObject > temp_for_truncation; // collection of positions to copy from original pdb
	for ( int i = 0, ie=(int)data.blueprint.size(); i < ie ; i++ ) {  // loop to extract positions to keep
		if ( data.blueprint[i].resname != "x" && data.blueprint[i].resname != "X" ) {
			temp_for_truncation.push_back( data.blueprint[i] );
		}
	}

	TR << "temp_for_truncation.size():" << temp_for_truncation.size() << std::endl;
	TR << "Setting up translate_index." << std::endl;

	// std::map<int,int> translate_index is now a class member variable
	utility::vector1< bool > keep( input_pose.total_residue(), false ); // by default, don't keep anything? this vector isn't used anyway so whatever.
	for ( int ii = 0, ie = (int)temp_for_truncation.size(); ii < ie; ii++ ) { // loop to update keep vector according to what's found

		keep[ temp_for_truncation[ ii ].original_index ] = true;
		if ( temp_for_truncation[ ii ].original_index != 0 ) { // correct for the use of "0" in marking extensions, original index should start from 1
			//TR << "temp_for_truncation[ ii ].original_index: " << temp_for_truncation[ ii ].original_index
			// << ", translate_index[ " << temp_for_truncation[ ii ].original_index << " ] set to " << ii + 1 << std::endl;
			translate_index[ temp_for_truncation[ ii ].original_index ] = ii + 1; // i is zero based, we want one based translation
		}
	}

	// soooo, temp_for_truncation, keep, and translate_index got created and then are not used anywhere else. wtf?
	// translate_index may get used somewhere since it's a class member variable, but the other two?

	std::vector< forge::remodel::LineObject > temp_for_copy;  // blueprint lines that have ss specified, what a horrible name. referring to this as 'lines_residues_to_remodel' in comments.
	std::vector< forge::remodel::LineObject > temp;       // blueprint lines that need to be copied b/c they have '.'???
	std::vector< forge::remodel::Segment > segmentStorageVector;
	//std::vector< forge::remodel::Segment > segment_to_copyVector;   // not used
	//std::vector< forge::remodel::Segment > segment_to_copyNewIndex; // not used

	for ( int i = 0, ie = (int)data.blueprint.size(); i < ie; i++ ) {
		if ( data.blueprint[i].sstype != "." ) { // first find the segments to be remodelled
			temp.push_back( data.blueprint[i] );
		} else if ( data.blueprint[i].sstype == "." ) { // parts to be copied
			temp_for_copy.push_back(data.blueprint[i]); // lines_residues_to_remodel
		} else {
			TR << "workingSetGen(): assignment error" << std::endl;
		}
	}

	// repeat structure loop over a second time; merge sections and update index
	if ( option[OptionKeys::remodel::repeat_structure].user() ) {
		//need to know the original index of the last element, for building
		//extensions or deletions across jxn points
		LineObject lastLO = data.blueprint.back();

		bool denovo = true;
		//have to loop to identify denovo case
		for ( int i = 0, ie = (int)data.blueprint.size(); i < ie; i++ ) {
			if ( data.blueprint[i].sstype == "." ) { //if anywhere hits this assignment, not de novo
				denovo = false;
			}
		}

		for ( int i = 0, ie = (int)data.blueprint.size(); i < ie; i++ ) {
			if ( data.blueprint[i].sstype != "." ) { // first find the segments to be remodeled
				LineObject LO = data.blueprint[i];
				//update indeces
				LO.index = LO.index + (int)data.blueprint.size();
				if ( LO.original_index != 0 ) { //in de novo case, the extension uses 0, don't increment.
					LO.original_index = LO.original_index + (int)lastLO.original_index;
				} else if ( LO.original_index == 0 && !denovo ) {
					LO.original_index = LO.original_index + (int)data.blueprint.size();
				} else {
					// de novo case don't increment
				}

				//TR << "LO object second time " << LO.index << " " << LO.original_index << std::endl;
				temp.push_back(LO);
			} else if ( data.blueprint[i].sstype == "." ) { // parts to be copied
				//not needed here
			} else {
				std::cout << "assignment error" << std::endl;
			}
		}
	}

	//save the first non-rebuilt position for potentially rooting a tree.
	if ( !temp_for_copy.empty() ) { // lines_residues_to_remodel
		safe_root_ = temp_for_copy.front().index; // lines_residues_to_remodel
	} else {
		safe_root_ = 1;
	}

	TR << "temp_for_copy (lines_residues_to_remodel): [ ";
	for ( Size ii=0; ii < temp_for_copy.size(); ++ii ) {
		TR << temp_for_copy[ ii ].original_index << "-" << temp_for_copy[ ii ].index << ", ";
	}
	TR << "]" << std::endl;

	// break up 'temp' into small segments
	protocols::forge::remodel::Segment segment;
	for ( int ii = 0, ie = (int)temp.size() - 1;  ii < ie ; ii++ ) { // compare the (i)-th and (i+1)-th element to find contiguous segments
		forge::remodel::LineObject first = temp[ ii ];
		forge::remodel::LineObject next  = temp[ ii+1 ];

		// compare the (i)-th and (i+1)-th element to find contiguous segments
		if ( next.index == ( first.index + 1 ) ) {
			segment.residues.push_back( first.index );

			// if reaching the end of the last segment
			if ( ii + 1 == (int)temp.size() - 1 ) {
				//TR << "next:" << next.index << std::endl;
				segment.residues.push_back( next.index );
				segmentStorageVector.push_back( segment );
				segment.residues.clear();
			}

		} else if ( next.index != first.index + 1 && (ii+1) == (int)temp.size() - 1 ) { // if there's a loner in the end by itself
			segment.residues.push_back(first.index);
			segmentStorageVector.push_back(segment);
			segment.residues.clear();

			segment.residues.push_back(next.index);
			segmentStorageVector.push_back(segment);
			segment.residues.clear();

		} else {
			segment.residues.push_back(first.index);
			// if reaching the end of the last segment
			if ( ii + 1 == (int)temp.size() - 1 ) {
				segment.residues.push_back(next.index);
				segmentStorageVector.push_back(segment);
				segment.residues.clear();
			}
			// TR << first.index << std::endl;
			segmentStorageVector.push_back( segment );
			segment.residues.clear();
		}
	}

	TR << "segmentStorageVector: [ ";
	for ( Size ii=0; ii < segmentStorageVector.size(); ++ii ) {
		TR << "segment " << ii << ": [ ";
		for ( Size jj=0; jj < segmentStorageVector[ ii ].residues.size(); ++jj ) {
			TR << segmentStorageVector[ ii ].residues[ jj ] << ", ";
		}
		TR << "], ";
	}
	TR << "]" << std::endl;

	using protocols::forge::build::BuildManager;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::build::SegmentInsert;
	using protocols::forge::build::SegmentInsertOP;
	using protocols::forge::build::BuildInstructionOP;
	typedef std::string String;
	typedef core::Size Size;

	using core::fragment::ConstantLengthFragSet;
	using core::fragment::ConstantLengthFragSetOP;
	using core::fragment::OrderedFragSetOP;
	using core::fragment::OrderedFragSet;
	using core::fragment::Frame;
	using core::fragment::FrameOP;
	using build::BuildInstructionOP;

	SegmentInsertOP segIns;

	// first get insertion index

	Size insertStartIndex = 0;
	Size insertEndIndex = 0;
	TR << "data.dssp_updated_ss: " << data.dssp_updated_ss << std::endl;
	if ( option[OptionKeys::remodel::domainFusion::insert_segment_from_pdb].user() ) {
		TR << "Processing insertion SS info..." << std::endl;
		insertStartIndex = data.dssp_updated_ss.find_first_of("I");
		insertEndIndex = data.dssp_updated_ss.find_last_of("I");
		TR << "Found insertion with insertStartIndex: " << insertStartIndex << " and insertEndIndex: " << insertEndIndex << std::endl;
	}

	// set generic aa type before assigning manager tasks.

	// for now only allow one letter code
	String build_aa_type =option[OptionKeys::remodel::generic_aa]; //defaults to VAL

	runtime_assert (build_aa_type.size() == 1);

	if ( option[OptionKeys::remodel::use_blueprint_sequence].user() ) {
		for ( int i = 0; i < (int)data.blueprint.size(); i++ ) {
			if ( data.blueprint[i].resname.compare("x") == 0  || data.blueprint[i].resname.compare("X") == 0 ) {
				aa.append(build_aa_type);
			} else {
				aa.append( data.blueprint[i].resname );
			}
		}
		//  runtime_assert( aa.size() == data.dssp_updated_ss.size());

	} else {
		if ( build_aa_type.compare("A") != 0 ) {
			//build the aa string to be the same length as dssp updated ss
			//use that length because repeat structures are bigger than blueprint
			for ( int i = 1; i<= (int)data.dssp_updated_ss.size(); i++ ) {
				aa.append(build_aa_type);
			}
		}
	}


	//debug
	std::cout << "AA for build: " << aa << std::endl;

	//BuildManager manager;

	// find the begin and end index
	for ( int i = 0, ie = (int)segmentStorageVector.size(); i < ie ; i++ ) {
		begin.push_back( segmentStorageVector[i].residues.front() );
		end.push_back( segmentStorageVector[i].residues.back() );

		core::Size idFront = segmentStorageVector[i].residues.front();
		core::Size idBack = segmentStorageVector[i].residues.back();
		core::Size seg_size = (int)data.blueprint.size();
		//core::Size rep_number =option[OptionKeys::remodel::repeat_structure];
		std::string DSSP = data.dssp_updated_ss;

		// Sachko: Changed Size (unsigned) to int (signed) to make this logic work properly.
		// Labonte (earlier comment): I don't know what happens when one assigns -1 to a Size, but this needs to be fixed.
		int head = -1, tail = -1, headNew = -1, tailNew = -1; //safety, init to negative values

		//use temp_For_copy to identify if it's de novo build; not empty means it's a loop case.
		if ( option[OptionKeys::remodel::repeat_structure].user() && !temp_for_copy.empty() ) {  // lines_residues_to_remodel
			//duplicate length of dssp and aastring
			DSSP += DSSP;
			aa += aa;

			if ( idFront > seg_size && idBack > seg_size ) { //ignore this type of assigment in repeat structures
				continue;
				/*
				if (idBack == 2*seg_size){ // can't hit the final residue in repeat unit pose
				idBack = idBack-1;
				}
				head = data.blueprint[ idFront-seg_size-1 ].original_index + seg_size;
				tail = data.blueprint[ idBack-seg_size-1 ].original_index + seg_size;
				headNew = data.blueprint[ idFront-seg_size-1 ].index + seg_size;
				tailNew = data.blueprint[ idBack-seg_size-1 ].index + seg_size;
				*/
			} else if ( idFront <= seg_size && idBack > seg_size ) { //spanning across repeat, adjust the tail
				if ( idBack == 2*seg_size ) { // can't hit the final residue in repeat unit pose
					idBack = idBack-1;
				}
				head = data.blueprint[ idFront-1 ].original_index;
				tail = data.blueprint[ idBack-seg_size-1 ].original_index + data.blueprint.back().original_index;
				if ( data.blueprint.back().original_index == 0 ) { //an extension
					tail = data.blueprint[ idBack-seg_size-1 ].original_index + seg_size;
				}
				headNew = data.blueprint[ idFront-1 ].index;
				tailNew = data.blueprint[ idBack-seg_size-1 ].index + seg_size;
			} else { //normal build in the first segment
				head = data.blueprint[ idFront-1 ].original_index;
				tail = data.blueprint[ idBack-1 ].original_index;
				headNew = data.blueprint[ idFront-1 ].index;
				tailNew = data.blueprint[ idBack-1 ].index;
			}
		} else if ( option[OptionKeys::remodel::repeat_structure].user() && temp_for_copy.empty() ) { //de novo case
			//std::cout << "idFront " << idFront << " idBack " << idBack << std::endl;
			Size range_limit = data.blueprint.size();
			if ( idBack >= range_limit ) { // can't assign beyond blueprint definition, as indices are missing as extensions
				idFront = 1;
				idBack = range_limit;
				head = data.blueprint[ idFront-1 ].original_index;
				tail = data.blueprint[ idBack-1 ].original_index;
				headNew = data.blueprint[ idFront-1 ].index;
				tailNew = data.blueprint[ idBack-1 ].index;
			}
		} else {
			head = data.blueprint[ idFront-1 ].original_index;
			tail = data.blueprint[ idBack-1 ].original_index;
			headNew = data.blueprint[ idFront-1 ].index;
			tailNew = data.blueprint[ idBack-1 ].index;
		}
		// Sachko on 02/14/2013
		// "tail" is set to 0 if that position did not exist in the first place, that is to insert a new atom.
		// So, here is a quick & dirty, hack, for now.
		//tail = tail<=head? data.blueprint[ idBack ].original_index : tail;
		//assert(tail>head);

		int gap = idBack - idFront +1;
		//int gap = segmentStorageVector[i].residues.back() - segmentStorageVector[i].residues.front() + 1;

		//debug
		//TR << "dssp size: " << data.dssp_updated_ss.size() << std::endl;
		TR << "head " << head << ":" << headNew << ", tail " << tail << ":" << tailNew << ", gap: " << gap
			<<  ", dssp_updated_ss.size(): " << data.dssp_updated_ss.size() << ", insert_ss: " << data.dssp_updated_ss.substr( headNew-1, gap ) << std::endl; // head-1 because dssp_updated_ss is 0 based std::string

		loops.add_loop( segmentStorageVector[i].residues.front(), segmentStorageVector[i].residues.back(), segmentStorageVector[i].residues.front()+1, 0, 0 );

		// process regions containing insertion
		if ( headNew <= static_cast<int>(insertStartIndex) && tailNew >= static_cast<int>(insertEndIndex) &&
				( (insertEndIndex - insertStartIndex) != 0 ) ) {
			TR << "segment contains insertion, skip normal SegmentRebuild instructions, use SegmentInsert instructions instead" << std::endl;

			String beforeInsert = data.dssp_updated_ss.substr( headNew-1, insertStartIndex-head + 1 ); // if we subtract 'head' and the first position was an extension, you get a really long string here
			//std::string beforeInsert = data.dssp_updated_ss.substr( headNew-1, insertStartIndex-headNew + 1 );
			String afterInsert = data.dssp_updated_ss.substr( insertEndIndex+1, tailNew-insertEndIndex - 1 );
			TR << "beforeInsert: " << beforeInsert << std::endl;
			TR << "afterInsert: " << afterInsert << std::endl;
			String blank;
			/*
			for (Size i=1; i<= data.insertionSS.size(); i++){
			blank.append("^");
			} // can't append extra ^ characters because insert_SS_string gets stored in SegmentInsert and used by the manager
			*/
			blank.append("^");
			std::string insert_SS_string = beforeInsert + blank + afterInsert;
			TR << "insert_SS_string: " << insert_SS_string << std::endl;

			using protocols::forge::build::SegmentInsertConnectionScheme::N; // default N2C insertion

			protocols::forge::build::SegmentInsertConnectionScheme::Enum connection_scheme = N; // default N2C insertion
			segIns = SegmentInsertOP( new SegmentInsert( Interval(head,tail), insert_SS_string , data.insertPose, false /*omega at junction*/, connection_scheme ) );
			manager.add( segIns );
			continue;
		}

		if ( head == 0 && segmentStorageVector[i].residues.front() == 1 ) { // N-term extension
			TR << "N-terminal extension found" << std::endl;
			manager.add( BuildInstructionOP( new SegmentRebuild( Interval(1,tail),  data.dssp_updated_ss.substr( headNew-1, gap ), aa.substr( headNew-1,gap )) ) );
		} else if ( tail == 0 && segmentStorageVector[i].residues.back() == static_cast<int>(model_length) ) {
			TR << "C-terminal extension found" << std::endl;
			gap = (int)data.blueprint.size()-segmentStorageVector[i].residues.front()+1;
			manager.add( BuildInstructionOP( new SegmentRebuild( Interval(head,input_pose.total_residue()), DSSP.substr( segmentStorageVector[i].residues.front()-1, gap ), aa.substr( segmentStorageVector[i].residues.front()-1, gap )) ) );
		} else if ( head != 0 && headNew == 1 && segmentStorageVector[i].residues.front() == 1 ) { // N-term deletion
			TR << "debug: N-term deletion" << std::endl;
			this->manager.add( BuildInstructionOP( new SegmentRebuild( Interval(1,tail),  DSSP.substr( headNew-1, gap ), aa.substr( headNew-1,gap )) ) );
		} else if ( tail != static_cast<int>(input_pose.total_residue()) && tailNew == static_cast<int>(model_length) && headNew == 1 &&
				segmentStorageVector[i].residues.back() == static_cast<int>(model_length) ) { // C-term deletion
			gap = (int)data.blueprint.size()-segmentStorageVector[i].residues.front()+1;
			TR << "debug: C-term deletion" << std::endl;
			this->manager.add( BuildInstructionOP( new SegmentRebuild( Interval(head,input_pose.total_residue()), DSSP.substr( segmentStorageVector[i].residues.front()-1, gap ), aa.substr( segmentStorageVector[i].residues.front()-1, gap )) ) );
		} else {
			TR << "normal rebuild" << std::endl;
			manager.add( BuildInstructionOP( new SegmentRebuild( Interval(head, tail), DSSP.substr( headNew-1, gap ), aa.substr( headNew-1, gap )) ) );
		}
	}

	return;
}

///
/// @brief
/// Takes in a pose and remodel data objects and constructs a packer task for that pose. Uses the resfile string
/// that was created while reading in the blueprint file.
/// What happens if the resfile string is "" ?
///
void RemodelWorkingSet::manualPackerTaskGen( pose::Pose const & built_pose, protocols::forge::remodel::RemodelData const & data ) {

	//task = pack::task::TaskFactory::create_packer_task( built_pose );
	core::pack::task::TaskFactoryOP tf = protocols::forge::methods::remodel_generic_taskfactory();

	//if need more operations added, put them here.

	// create the real task
	task = tf->create_task_and_apply_taskoperations( built_pose );

	//pose::PDBPoseMap map(built_pose.pdb_info()->pdb2pose());
	//TR << map.find(' ',1,' ') << "PDBPosemap" << std::endl;
	core::pack::task::parse_resfile_string( built_pose, *task, data.parsed_string_for_resfile );

}

} //namespace remodel
} //namespace forge
} //namespace protocols
