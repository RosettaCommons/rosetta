// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 13011 $
//  $Date: 2007-02-21 17:17:13 -0800 (Wed, 21 Feb 2007) $
//  $Author: possu $

// Rosetta Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
//test
// AUTO-REMOVED #include <core/pose/PDBPoseMap.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <protocols/forge/remodel/RemodelWorkingSet.hh>
//#include <devel/remodel/helpMenu.hh>

//for DSSP
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>

//fragset
//#include <core/fragment/OrderedFragSet.hh>

#include <core/chemical/ResidueType.hh>
#include <core/kinematics/FoldTree.hh>
 // for switch typeset

// for yab managers
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/methods/util.hh>
//#include <devel/fragment/picking/vall/util.hh> // pick_fragment_by_ss
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/build/SegmentInsert.hh>

// for resfile command map
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/TaskFactory.hh>

/*
//yab headers
#include "AtomPoint.hh"
#include "BoundingBox.hh"
#include "epigraft_functions.hh"
#include "Octree.hh"
#include "rootstock_types.hh"
#include "ccd_functions.hh"
*/
// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
//#include <ObjexxFCL/FArray3D.hh>
//#include <ObjexxFCL/FArray4D.hh>
//#include <ObjexxFCL/FArray5D.hh>
//#include <ObjexxFCL/format.hh>
//#include <ObjexxFCL/string.functions.hh>


// C++ Headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
// AUTO-REMOVED #include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>


// Utility Headers
// AUTO-REMOVED #include <utility/basic_sys_util.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/ocstream.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>
#include <utility/vector1.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/OrderedFragSet.fwd.hh>
#include <core/kinematics/Jump.hh>



//////////////// REMODEL

namespace protocols{
namespace forge{
namespace remodel{

static basic::Tracer TR("REMODELw");

/// @brief copy constructor
WorkingRemodelSet::WorkingRemodelSet(WorkingRemodelSet const & rval):
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
	hasInsertion= hasInsertion;
}

WorkingRemodelSet & WorkingRemodelSet::operator = ( WorkingRemodelSet const & rval ){
	if (this != & rval) {
		loops =  rval.loops ;
		sequence = rval.sequence;
		ss = rval.ss;
		translate_index = rval.translate_index;
		begin = rval.begin ;
		end =rval.end;
		copy_begin= rval.copy_begin ;
		copy_end =rval.copy_end ;
		src_begin =rval.src_begin;
		src_end =rval.src_end;
		hasInsertion= hasInsertion;
		manager = rval.manager;
		task =  rval.task ;
		rvjump_pose =  rval.rvjump_pose;
	}
	return *this;
}





void
protocols::forge::remodel::WorkingRemodelSet::workingSetGen(
	core::pose::Pose const & input_pose,
	protocols::forge::remodel::RemodelData const & data
)      {
	using namespace basic::options;
	//core::util::switch_to_residue_type_set( input_pose, core::chemical::CENTROID );

	//find rebuilding segments
	int model_length = (int)data.sequence.size();
	bool NtermExt = false;
	bool CtermExt = false;
	bool length_changed = false;

	// copy ss/seq from RemodelData so it can be passed elsewhere later
	this->ss = data.dssp_updated_ss;

	//this is purely experimental for matching fragment set
	if (option[ OptionKeys::remodel::repeat_structure].user()){
		this->ss.append(this->ss);
	} else {
		this->ss.append("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD");
	}

	this->sequence = data.sequence;

	//find N term extension, if any
	std::string Xs = "xX";
	int first_ext;
	first_ext = (int)data.sequence.find_first_of(Xs);
	if (first_ext == 0){
		NtermExt = true;
		std::cout << "N-terminal is extended" << std::endl;
	}
	//find C term extension, if any
	int last_ext;
	last_ext = (int)data.sequence.find_last_of(Xs);
	if (last_ext == model_length -1) {
		std::cout << "C-terminal is extended" << std::endl;
		// if C-term extension, add extra degenerate ss type to get fragments past last res.

		CtermExt = true;
	}

	if (data.blueprint[model_length-1].index != data.blueprint[model_length-1].original_index) {
		length_changed = true;
		std::cout << "length change found" << std::endl;
	}

	//this is needed for manipulating denovo cases, as they are coded as Cterm
	//extensions.  Affects SegmentRebuld selections.
	if (option[ OptionKeys::remodel::repeat_structure].user() && CtermExt){
		model_length = model_length * 2;
	}


// find all the indices.
	//identify truncation
	// for all the positions that are not extensions "x" or "X", put them in the
	// temp_for_truncation vector -- this corresponds to all the regions need to
	// be copied from the original pdb with deletions included.  fragment_pdb
	// requires a boolean vector that correspond to the length of the original
	// pdb, so initializes a "keep" vector of that size and set values to false.
	// As we iterate over the template_for_truncation, all the positions seen by
	// this iterative step are kept.  Because the truncated_pose is renumbered, we
	// also initializes a translate_index map to link the truncated index to the
	// original index.

	std::vector<protocols::forge::remodel::LineObject> temp_for_truncation; // collection of positions to copy from  original pdb
	//std::cout << input_pose.total_residue() << std::endl;
	utility::vector1<bool> keep(input_pose.total_residue(),false);
	for (int i = 0, ie=(int)data.blueprint.size(); i < ie ; i++){  // loop to extract positions to keep
		if (data.blueprint[i].resname != "x" && data.blueprint[i].resname != "X"){
			temp_for_truncation.push_back(data.blueprint[i]);
	//	std::cout << data.blueprint[i].resname << data.blueprint[i].original_index << std::endl;
		}
	}
	//std::map<int,int> translate_index; // moved to object data
	for (int i = 0, ie = (int)temp_for_truncation.size(); i < ie; i++){ // loop to update keep vector according to what's found
		keep[temp_for_truncation[i].original_index] = true;
	//	std::cout << temp_for_truncation[i].original_index << std::endl;
		if (temp_for_truncation[i].original_index != 0){ // correct for the use of "0" in marking extensions, original index should start from 1
			translate_index[temp_for_truncation[i].original_index] = i+1; // one based translation
		}
	}

	std::vector<protocols::forge::remodel::LineObject> temp_for_copy;
	std::vector<protocols::forge::remodel::LineObject> temp;
	std::vector<protocols::forge::remodel::Segment> segmentStorageVector;
	std::vector<protocols::forge::remodel::Segment> segment_to_copyVector;
	std::vector<protocols::forge::remodel::Segment> segment_to_copyNewIndex;
	for (int i = 0, ie = (int)data.blueprint.size(); i < ie; i++){
		if (data.blueprint[i].sstype != ".") { // first find the segments to be remodeled
			temp.push_back(data.blueprint[i]);
		}
		else if (data.blueprint[i].sstype == "."){ // parts to be copied
			temp_for_copy.push_back(data.blueprint[i]);
		}
		else {
			std::cout << "assignment error" << std::endl;
		}
	}

	if (option[ OptionKeys::remodel::repeat_structure].user()){ // repeat structure loop over a second time; merge sections and update index
		for (int i = 0, ie = (int)data.blueprint.size(); i < ie; i++){
			if (data.blueprint[i].sstype != ".") { // first find the segments to be remodeled
				LineObject LO = data.blueprint[i];
				//update indeces
				LO.index = LO.index + (int)data.blueprint.size();
				if (LO.original_index != 0){ //in de novo case, the extension uses 0, don't increment.
					LO.original_index = LO.original_index + (int)data.blueprint.size();
				}
				//TR << "LO object second time " << LO.index << " " << LO.original_index << std::endl;
				temp.push_back(LO);
			}
			else if (data.blueprint[i].sstype == "."){ // parts to be copied
				//not needed here
			}
			else {
				std::cout << "assignment error" << std::endl;
			}
		}
	}

	//save the first non-rebuilt position for potentially rooting a tree.
	if ( !temp_for_copy.empty() ){
		safe_root_ = temp_for_copy.front().index;
	}else {
		safe_root_ = 1;
		}


	//break up temp into small segments
	protocols::forge::remodel::Segment segment;
	for (int i = 0, ie = (int)temp.size()-1;  i < ie ; i++) { // compare the (i)-th and (i+1)-th element to find contiguous segments
		protocols::forge::remodel::LineObject first = temp[i];
		protocols::forge::remodel::LineObject next  = temp[i+1];
		if (next.index == (first.index+1)) {
			segment.residues.push_back(first.index);
	//		std::cout << first.index ;
			if (i+1 == (int)temp.size()-1){ // if reaching the end of the last segment
	//			std::cout << "next:" << next.index << std::endl;
				segment.residues.push_back(next.index);
				segmentStorageVector.push_back(segment);
				segment.residues.clear();
			}
		}
		else if (next.index != first.index+1 && (i+1) == (int)temp.size()-1 ){ //if there's a loner in the end by itself
        segment.residues.push_back(first.index);
        segmentStorageVector.push_back(segment);
        segment.residues.clear();
        segment.residues.push_back(next.index);
        segmentStorageVector.push_back(segment);
        segment.residues.clear();
    }
		else {
			segment.residues.push_back(first.index);
			if (i+1 == (int)temp.size()-1){ // if reaching the end of the last segment
				segment.residues.push_back(next.index);
				segmentStorageVector.push_back(segment);
				segment.residues.clear();
			}
	//		std::cout << first.index << std::endl;
			segmentStorageVector.push_back(segment);
			segment.residues.clear();
		}
	}

	//test only

	using protocols::forge::build::BuildManager;
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::build::SegmentInsert;
	using protocols::forge::build::SegmentInsertOP;
	typedef std::string String;
	typedef core::Size Size;

	using core::fragment::ConstantLengthFragSet;
	using core::fragment::ConstantLengthFragSetOP;
	using core::fragment::OrderedFragSetOP;
	using core::fragment::OrderedFragSet;
	using core::fragment::Frame;
	using core::fragment::FrameOP;
	using namespace basic::options;

	SegmentInsertOP segIns;

	//first get insertion index

	Size insertStartIndex =0;
	Size insertEndIndex =0;
	if (option[OptionKeys::remodel::domainFusion::insert_segment_from_pdb].user()){
		TR << "Processing insertion SS info..." << std::endl;
		insertStartIndex = data.dssp_updated_ss.find_first_of("I");
		insertEndIndex = data.dssp_updated_ss.find_last_of("I");
		TR << "debug: insertStartIndex: " << insertStartIndex << " insertEndIndex: " << insertEndIndex << std::endl;
	//process ss_string
	//	String beforeInsert = data.dssp_updated_ss.substr(0, insertStartIndex);
	//	String afterInsert = data.dssp_updated_ss.substr(insertEndIndex+1);
	//	data.ss_string = beforeInsert + data.insertionSS + afterInsert;
	//	TR << ss_string << std::endl;

	}

	//set generic aa type before assigning manager tasks.

	//for now only allow one letter code
	String build_aa_type = option[OptionKeys::remodel::generic_aa]; //devaults to VAL

	runtime_assert (build_aa_type.size() == 1);

	if (option[OptionKeys::remodel::use_blueprint_sequence].user()){
		for (int i = 0; i < (int)data.blueprint.size(); i++){
			if ( data.blueprint[i].resname.compare("x") == 0  || data.blueprint[i].resname.compare("X") == 0 ){
				aa.append(build_aa_type);
			}
			else {
				aa.append( data.blueprint[i].resname );
			}
		}
	//		runtime_assert( aa.size() == data.dssp_updated_ss.size());

		if (option[OptionKeys::remodel::repeat_structure].user()){
			String monomer_seq = aa;
			Size copies = option[OptionKeys::remodel::repeat_structure];
			while (copies > 1){ //first copy already made
				aa.append( monomer_seq );
				copies--;
			}
			//runtime_assert( aa.size() == data.dssp_updated_ss.size());
		}
	}	else {
		if ( build_aa_type.compare("A") != 0){
						//build the aa string to be the same length as dssp updated ss
						//use that length because repeat structures are bigger than blueprint
						for (int i = 1; i<= (int)data.dssp_updated_ss.size(); i++){
							aa.append(build_aa_type);
						}
		}
	}


	//debug
	std::cout << "AA for build: " << aa << std::endl;

	//BuildManager manager;

	// find the begin and end index
	for (int i = 0, ie = (int)segmentStorageVector.size(); i < ie ; i++){
		this->begin.push_back(segmentStorageVector[i].residues.front());
		this->end.push_back(segmentStorageVector[i].residues.back());

		core::Size idFront = segmentStorageVector[i].residues.front();
		core::Size idBack = segmentStorageVector[i].residues.back();
		core::Size seg_size = (int)data.blueprint.size();
		core::Size rep_number = option[ OptionKeys::remodel::repeat_structure];
		std::string DSSP = data.dssp_updated_ss;

		core::Size head = -1, tail = -1, headNew = -1, tailNew = -1; //safety, init to negative values

		if (option[ OptionKeys::remodel::repeat_structure].user() && input_pose.total_residue() == seg_size * rep_number){ //repeat and the blueprint do not match input pdb
			//duplicate length of dssp and aastring
			DSSP += DSSP;
			aa += aa;

			if (idFront > seg_size && idBack > seg_size){ //ignore this type of assigment in repeat structures
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
			}
			else if (idFront <= seg_size && idBack > seg_size){ //spanning across repeat, adjust the tail
				if (idBack == 2*seg_size){ // can't hit the final residue in repeat unit pose
					idBack = idBack-1;
				}
				head = data.blueprint[ idFront-1 ].original_index;
				tail = data.blueprint[ idBack-seg_size-1 ].original_index + seg_size;
				headNew = data.blueprint[ idFront-1 ].index;
				tailNew = data.blueprint[ idBack-seg_size-1 ].index + seg_size;
			} else { //normal build in the first segment
				head = data.blueprint[ idFront-1 ].original_index;
				tail = data.blueprint[ idBack-1 ].original_index;
				headNew = data.blueprint[ idFront-1 ].index;
				tailNew = data.blueprint[ idBack-1 ].index;
			}
		} else if ( option[ OptionKeys::remodel::repeat_structure].user() && temp_for_copy.empty()) { //de novo case
			//std::cout << "idFront " << idFront << " idBack " << idBack << std::endl;
			Size range_limit = data.blueprint.size();
			if (idBack >= range_limit){ // can't assign beyond blueprint definition, as indices are missing as extensions
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

		int gap = idBack - idFront +1;

		//debug
		//TR << "dssp size: " << data.dssp_updated_ss.size() << std::endl;
		TR << "head " << head << ":" << headNew << " tail " << tail << ":" << tailNew << " gap " << gap <<  " ss " << DSSP.size() << " " << DSSP.substr( headNew-1, gap ) << std::endl; // head -1 because dssp_updated_ss is 0 based std::string

		this->loops.add_loop(segmentStorageVector[i].residues.front(), segmentStorageVector[i].residues.back(), segmentStorageVector[i].residues.front()+1, 0, 0);

	  // process regions containing insertion
		if ( headNew <= insertStartIndex && tailNew >= insertEndIndex && ((insertEndIndex-insertStartIndex) != 0)){
			TR << "segment contain insertion, skip normal SegmentRebuild instructions, use SegmentInsert instructions instead" << std::endl;
		String beforeInsert = DSSP.substr(headNew-1, insertStartIndex-head+1);
		String afterInsert = DSSP.substr(insertEndIndex+1, tailNew-insertEndIndex-1);
		TR << "DEBUG beforeInsert: " << beforeInsert << std::endl;
		TR << "DEBUG afterInsert: " << afterInsert << std::endl;
		String blank;
		/*
		for (Size i=1; i<= data.insertionSS.size(); i++){
			blank.append("^");
		}*/

			blank.append("^");
		std::string insert_SS_string = beforeInsert + blank + afterInsert;
		TR << "DEBUG insert_SS_string: " << insert_SS_string << std::endl;

			using protocols::forge::build::SegmentInsertConnectionScheme::N;//default N2C insertion

			protocols::forge::build::SegmentInsertConnectionScheme::Enum connection_scheme = N;//default N2C insertion
			segIns = new SegmentInsert( Interval(head,tail), insert_SS_string , data.insertPose, false /*omega at junction*/, connection_scheme);
			this->manager.add(segIns);
			continue;
		}

		if (head == 0 && segmentStorageVector[i].residues.front() == 1 ){ // N-term extension
			TR << "debug: N-term ext" << std::endl;
			this->manager.add( new SegmentRebuild( Interval(1,tail),  DSSP.substr( headNew-1, gap ), aa.substr( headNew-1,gap )) );
		}
		else if (tail ==0 && segmentStorageVector[i].residues.back() == model_length){
		  TR << "debug: C-term ext" << std::endl;
			gap = (int)data.blueprint.size()-segmentStorageVector[i].residues.front()+1;
			this->manager.add( new SegmentRebuild( Interval(head,input_pose.total_residue()), DSSP.substr( segmentStorageVector[i].residues.front()-1, gap ), aa.substr( segmentStorageVector[i].residues.front()-1, gap )) );
		}
		else {
			TR << "debug: normal rebuild" << std::endl;
			this->manager.add( new SegmentRebuild( Interval(head, tail),  DSSP.substr( headNew-1, gap ), aa.substr( headNew-1, gap )) );
		}
	}

/*
//	ConstantLengthFragSetOP frag9( new ConstantLengthFragSet( 9 ) );
//	ConstantLengthFragSetOP frag3( new ConstantLengthFragSet( 3 ) );
	OrderedFragSetOP frag1;
	OrderedFragSetOP fragSet ( new OrderedFragSet );

	for ( protocols::loops::Loops::iterator itr = this->loops.v_begin(); itr != this->loops.v_end(); itr++){

		// setup regions
		// Pick fragments.  For now just use the 9-mer, 3-mer
		// breakdown to get things working.  This will be changed
		// to full-mer/variable length very soon.
		core::Size length = (*itr).size();
		//std::cout << "length: " << length << std::endl;
				for ( core::Size j = 0, je = length; j <= je; ++j ) {
				TR << "picking " << 200 << " 9-mers for position " << ( (*itr).start() + j ) << std::endl;
				String ss_sub = ss.substr( (*itr).start() + j - 1, 9 );
				FrameOP frame = new Frame( (*itr).start() + j, 9 );
				frame->add_fragment( core::fragment::picking::vall::pick_fragments_by_ss( ss_sub, 200 ) );
				fragSet->add( frame );
			}


		//pick the matching length fragment


		for ( core::Size j = 0, je = length; j <= je; ++j ) {
			TR << "picking " << 200 << " matching-mers for position " << ( (*itr).start() + j ) << std::endl;
			String ss_sub = ss.substr( (*itr).start() + j - 1, length );
			FrameOP frame = new Frame( (*itr).start() + j, length );
			frame->add_fragment( core::fragment::picking::vall::pick_fragments_by_ss( ss_sub, 200 ) );
			fragSet->add( frame );
		}



			for ( core::Size j = 0, je = length; j <= je; ++j ) {
				TR << "picking " << 200 << " 3-mers for position " << ( (*itr).start() + j ) << std::endl;
				String ss_sub = ss.substr( (*itr).start() + j - 1, 3 );
				FrameOP frame = new Frame( (*itr).start() + j, 3 );
				frame->add_fragment( core::fragment::picking::vall::pick_fragments_by_ss( ss_sub, 200 ) );
				fragSet->add( frame );
			}


		// make 1-mers from 3-mers
		frag1 = protocols::forge::components::smallmer_from_largemer( fragSet->begin(), fragSet->end(), 1 );
	}


	// Init VLB.  Be aware this is a bootstrap implementation, even
	// remotely sane results are not guaranteed. To get things pinned
	// down with the proper implementation and benchmarked is going to
	// take some time.
#if defined GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( model_pose.conformation(), "Remodel Test" );
#endif
	manager.modify(model_pose);
	//TR<< model_pose.psi(242)<< std::endl;
	//vlb.apply( model_pose );

	// setup loop building protocol
	protocols::loops::LoopMover_Perturb_QuickCCD_Moves loop_mover( this->loops, false );
	loop_mover.set_strict_loops(true);
	//loop_mover.add_fragments( frag9 );
	loop_mover.add_fragments( fragSet );
	//loop_mover.add_fragments( frag3 );
	loop_mover.add_fragments( frag1 );

	// Run loop modeling.  The loop movers return the original fold tree
	// after apply().  Do we want that...?  There's also no good way to
	// check that the loop mover actually finished with a closed loop,
	// which is troubling.  For now we work around this by post-evaluating
	// the chainbreak at the original cutpoint.
	loop_mover.apply( model_pose );

	model_pose.dump_pdb("test_vlb.pdb");
*/
	return;
}

void
protocols::forge::remodel::WorkingRemodelSet::manualPackerTaskGen(core::pose::Pose const & built_pose, protocols::forge::remodel::RemodelData const & data)
{
//	this->task = core::pack::task::TaskFactory::create_packer_task( built_pose );

  core::pack::task::TaskFactoryOP TF = protocols::forge::methods::remodel_generic_taskfactory();

  //if need more operations added, put them here.
  //

  //create the real task
  this->task = TF->create_task_and_apply_taskoperations( built_pose );

	//core::pose::PDBPoseMap map(built_pose.pdb_info()->pdb2pose());
	//TR << map.find(' ',1,' ') << "PDBPosemap" << std::endl;
	core::pack::task::parse_resfile_string( built_pose, *this->task, data.parsed_string_for_resfile );
}

} //namespace remodel
} //namespace forge
} //namespace protocols
