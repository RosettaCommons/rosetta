// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/matdes/StoreQuasiSymmetricTaskMover.cc
/// @brief The StoreQuasiSymmetricTaskMover accepts a set of TaskOperations and ensures that
/// equivalent positions that are designable in one quasisymmetric subunit are designable in all.
/// In addition, this is where the RotamerLinks that ensure that these positions maintain the
/// same AA identity are applied.
/// @author Neil King (neilking@uw.edu)
/// @edited Yang Hsia (yhsia@uw.edu)

// Unit Headers
#include <devel/matdes/StoreQuasiSymmetricTaskMover.hh>
#include <devel/matdes/StoreQuasiSymmetricTaskMoverCreator.hh>

//project headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/make_symmetric_task.hh>
#include <protocols/toolbox/task_operations/STMStoredTask.hh>
#include <basic/Tracer.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <ObjexxFCL/format.hh>

#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer TR( "devel.matdes.StoreQuasiSymmetricTaskMover" );

namespace devel {
namespace matdes {

// @brief default constructor
StoreQuasiSymmetricTaskMover::StoreQuasiSymmetricTaskMover() {}

// @brief destructor
StoreQuasiSymmetricTaskMover::~StoreQuasiSymmetricTaskMover() = default;

// @brief getters
char StoreQuasiSymmetricTaskMover::quasi_symm_comp() const { return quasi_symm_comp_[0]; }
core::Size StoreQuasiSymmetricTaskMover::num_quasi_repeats() const { return num_quasi_repeats_; }
core::Size StoreQuasiSymmetricTaskMover::offset_resis() const { return offset_resis_; }

// @brief setters
void StoreQuasiSymmetricTaskMover::quasi_symm_comp( std::string const quasi_symm_comp ) { quasi_symm_comp_ = quasi_symm_comp; }
void StoreQuasiSymmetricTaskMover::num_quasi_repeats( core::Size const num_quasi_repeats ) { num_quasi_repeats_ = num_quasi_repeats; }
void StoreQuasiSymmetricTaskMover::offset_resis( core::Size const offset_resis ) { offset_resis_ = offset_resis; }

void
StoreQuasiSymmetricTaskMover::apply( core::pose::Pose & pose )
{

	using namespace core;
	using namespace basic;
	using namespace pose;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace scoring;
	using namespace utility;

	// Create the raw PackerTask using the TaskFactory created in parse_my_tag
	runtime_assert( task_factory_ != nullptr );
	core::pack::task::PackerTaskOP raw_task = task_factory_->create_task_and_apply_taskoperations( pose );
	//core::pack::task::PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations( pose );   //raw_task vs task???

	// Get the designable positions from the raw task.
	// Also get the number of residues in the quasisymmetric subunits and in the
	// asymmetric unit as a whole.
	std::set< core::Size > design_pos;
	core::Size num_resis_quasi_subunits = 0;
	core::Size num_indy_resis = offset_resis();
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	vector1<bool> indy_resis = sym_info->independent_residues();
	core::Size ir= 1;

	//offset residues that will NOT be included in quasi-component (needs to be first in the chain)
	if ( int(offset_resis()) > 0 ) {
		ir= ( int(offset_resis()) + 1 );
		TR << "Offsetting quasi-equivalent residues by: " << offset_resis() << std::endl;
		TR << "Starting quasi-equivalent residue (ir)= " << ir <<std::endl;
	}

	//for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
	while ( ir<=sym_info->num_total_residues_without_pseudo() ) {                //loop through all residues (without_pseudo = without pseudo-atoms; all residues in pose)
		if ( indy_resis[ir] ) {                                  //if residue is a independent residue
			num_indy_resis++;                                    //+1 total independent residue count (indy resi = residue in a master subunit)
			if ( ( pose.residue(ir).is_protein() )
					&& ( get_component_of_residue(pose,ir) == quasi_symm_comp() ) ) {         //if residue is a protein AND if component is a designated quasi comp
				num_resis_quasi_subunits++;                               //+1 to TOTAL resi count of quasi subunits
				TR.Debug << "Residue " << ir << " is on a quasi component, total is now: " << num_resis_quasi_subunits << std::endl;
			}
			if ( raw_task->being_packed( ir ) ) {                          //check of the current residue is a designable residue
				design_pos.insert( ir );                               //if so, then add to design_pos list
			}
		}
		ir++;
	}
	TR << "num_indy_resis = " << num_indy_resis << std::endl;                  //number of residues that are in the master subunits (A and B combined)
	TR << "num_resis_quasi_subunits = " << num_resis_quasi_subunits << std::endl;        //number of residues that are in the quasi-equivalent subunit (A or B)

	// Iterate through the designable positions and make sure that if a position is designable in one
	// quasisymmetric subunit, it is also designable in the other
	std::set< core::Size >::iterator pos;
	for ( pos=design_pos.begin(); pos!=design_pos.end(); ++pos ) {               //ITERATOR, it goes through the design_pos vector
		if ( indy_resis[*pos] == 0 ) continue;                          //Make sure we're only looking in the asymmetric unit
		if ( get_component_of_residue(pose,*pos) != quasi_symm_comp() ) continue;         //Only look in the quasisymmetric subunits

		core::Size pos_minus_Xsubunit= 0;
		core::Size pos_plus_Xsubunit= 0;

		for ( Size repeat_mult=1; repeat_mult < num_quasi_repeats(); repeat_mult++ ) {
			pos_minus_Xsubunit = (int)*pos - ( repeat_mult * ( num_resis_quasi_subunits/num_quasi_repeats() ) );
			pos_plus_Xsubunit = (int)*pos + ( repeat_mult * ( num_resis_quasi_subunits/num_quasi_repeats() ) );

			if ( int(pos_minus_Xsubunit) > int(offset_resis())                                     //pos_minus is not a offset residue (or less than position 0)
					&& indy_resis[ pos_minus_Xsubunit ]                        //pos_minus is an independent residue
					&& get_component_of_residue(pose,pos_minus_Xsubunit) == quasi_symm_comp() ) {   //pos_minus is in a quasi component
				if ( design_pos.find( pos_minus_Xsubunit ) == design_pos.end() ) {          //pos_minus is NOT already in the design_pos list
					design_pos.insert( pos_minus_Xsubunit );
					TR.Debug << "Position (-" << repeat_mult << ") " << pos_minus_Xsubunit << " inserted on behalf of position " << *pos
						<< ": sele chain " << get_component_of_residue(pose,*pos) << " and resi " << pos_minus_Xsubunit << "+" << *pos
						<< std::endl;
				} else {
					TR.Debug << "Position (-" << repeat_mult << ") " << pos_minus_Xsubunit << " is already a design_pos, checked behalf of position " << *pos << std::endl;
				}
			} else if ( pos_plus_Xsubunit <= num_indy_resis                      //pos_plus is an independent residue
					&& indy_resis[ pos_plus_Xsubunit ]                      //pos_plus is an independent residue
					&& get_component_of_residue(pose,pos_plus_Xsubunit) == quasi_symm_comp() ) {  //pos_plus is in a quasi-component
				if ( design_pos.find( pos_plus_Xsubunit ) == design_pos.end() ) {           //pos_plus is NOT already in the design_pos list
					design_pos.insert( pos_plus_Xsubunit );
					TR.Debug << "Position (+" << repeat_mult << ") " << pos_plus_Xsubunit << " inserted on behalf of position " << *pos
						<< ": sele chain " << get_component_of_residue(pose,*pos) << " and resi " << pos_plus_Xsubunit << "+" << *pos
						<< std::endl;
				} else {
					TR.Debug << "Position (+" << repeat_mult << ") " << pos_plus_Xsubunit << " is already a design_pos, checked behalf of position " << *pos << std::endl;
				}
			}
		}
		//int pos_minus_1subunit = (int)*pos-((int)num_resis_quasi_subunits/int(num_quasi_repeats()));
		//TR << "pos_minus_1subunit: " << pos_minus_1subunit << " = " << *pos << " - ( " << num_resis_quasi_subunits << " / " << num_quasi_repeats() << " )" << std::endl;
		//int pos_plus_1subunit = (int)*pos+((int)num_resis_quasi_subunits/int(num_quasi_repeats()));
		//TR << "pos_plus_1subunit: " << pos_plus_1subunit << " = " << *pos << " + ( " << num_resis_quasi_subunits << " / " << num_quasi_repeats() << " )" << std::endl;
		//
		//if ( pos_minus_1subunit > int(offset_resis()) && indy_resis[ pos_minus_1subunit ] && get_component_of_residue(pose,pos_minus_1subunit) == quasi_symm_comp() ) {  // pos_minus_1sub is not a offset residue AND pos_minus_1sub is a indy resi AND component is a quasi component
		//  if ( design_pos.find( pos_minus_1subunit ) == design_pos.end() ) {
		//    design_pos.insert( pos_minus_1subunit );
		//    TR << "Position (-) " << pos_minus_1subunit << " inserted on behalf of position " << *pos << ": sele chain " << get_component_of_residue(pose,*pos) << " and resi " << pos_minus_1subunit << "+" << *pos << std::endl;
		//  }
		//} else if ( pos_plus_1subunit <= num_indy_resis && indy_resis[ pos_plus_1subunit ] && get_component_of_residue(pose,pos_plus_1subunit) == quasi_symm_comp() ) {
		//  if ( design_pos.find( pos_plus_1subunit ) == design_pos.end() ) {
		//    design_pos.insert( pos_plus_1subunit );
		//    TR << "Position (+) " << pos_plus_1subunit << " inserted on behalf of position " << *pos  << ": sele chain " << get_component_of_residue(pose,*pos) << " and resi " << pos_plus_1subunit << "+" << *pos << std::endl;
		//  }
	}

	// Create a new task; this is the one we will store.
	// Prevent_repacking at all positions that are not design positions
	// and apply rotamer links to all designable positions.
	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
	core::pack::rotamer_set::RotamerLinksOP links( new core::pack::rotamer_set::RotamerLinks );
	links->resize( num_indy_resis );
	std::string output = "select quasi_design_pos, resi ";
	core::Size ir_plus_Xsubunit = 0;

	//figure out what the residue number of the last residue in the first quasi repeat is
	core::Size first_quasi_repeat_last_resi = 0;
	if ( quasi_symm_comp() ==  'A' ) {
		first_quasi_repeat_last_resi = core::Size( offset_resis() + ( num_resis_quasi_subunits/num_quasi_repeats() ) );
	} else if ( quasi_symm_comp() == 'B' ) {
		first_quasi_repeat_last_resi = core::Size( ( num_resis_quasi_subunits/num_quasi_repeats() ) + ( num_indy_resis - num_resis_quasi_subunits ) );
	} else {
		TR << "WARNING!! " << quasi_symm_comp() << " is in an incompatible component! (not A or B)." << std::endl;
	}

	TR << "num_total_residues_without_pseudo: " << sym_info->num_total_residues_without_pseudo() << std::endl;

	for ( Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++ ) {         //loop through all residues
		if ( design_pos.find(ir) != design_pos.end() ) {                      //residue is found in design_pos (note: only master subunit resi)
			output += ObjexxFCL::string_of(ir)+"+";  //appends design positions to output string
			if ( pose.residue(ir).is_protein()                           //residue is a protein
					&& get_component_of_residue(pose,ir) != quasi_symm_comp() ) {           //residue is NOT a quasi-residue
				links->set_equiv( ir, ir );                              //For non-quasisymmetrical subunits, set up a dummy RotamerLink to this resi (itself)
				TR.Debug << "Residue: " << ir << " set with dummy RotamerLink. (non-quasi)" << std::endl;
			} else if ( pose.residue(ir).is_protein()                         //residue is a protein
					&& get_component_of_residue(pose,ir) == quasi_symm_comp()           //residue is a quasi-residue
					&& int(ir) <= int(offset_resis()) ) {                     //residue is a offset residue
				links->set_equiv( ir, ir );                              //For offset residues, set up a dummy RotamerLink to this resi
				TR.Debug << "Residue: " << ir << " set with dummy RotamerLink. (offset residue)" << std::endl;
			} else if ( pose.residue(ir).is_protein()                  //residue is a protein
					//generate a vector(list) of equivalent residues for quasi, then apply RotamerLinks
					//      else if ( pose.residue(ir).is_protein()                                //residue is a protein
					//           && ( get_component_of_residue(pose,ir) == quasi_symm_comp() )                //residue is in a quasi comp
					//           && ( int(ir) <= int( offset_resis() + ( num_resis_quasi_subunits/num_quasi_repeats() ) ) ) //residue is in the first quasi repeat
					//           && ( int(ir) > int(offset_resis()) ) ) {                          //residue is not a offset residue
					&& ( get_component_of_residue(pose,ir) == quasi_symm_comp() )  //residue is in a quasi comp
					&& ( int(ir) <= int(first_quasi_repeat_last_resi) )         //residue is in the first quasi repeat
					&& ( int(ir) > int(offset_resis()) ) ) {            //residue is not a offset residue

				//TR.Debug << "Residue: " << ir << " IS GETTING INTO THIS LOOP." << std::endl;

				utility::vector1< core::Size > list; list.push_back( ir );                    //generate list for linking
				for ( Size repeat_mult=1; repeat_mult < num_quasi_repeats(); repeat_mult++ ) {
					ir_plus_Xsubunit = ir + ( repeat_mult * ( num_resis_quasi_subunits/num_quasi_repeats() ) );
					if ( ( indy_resis[ir_plus_Xsubunit] )                             //ir_plus is a independent residue
							&& ( get_component_of_residue(pose,ir_plus_Xsubunit) == quasi_symm_comp() ) ) {     //ir_plus is a quasi-residue
						list.push_back( ir_plus_Xsubunit );
					}
				}

				links->set_equiv( ir, list );                                   //link residues based on list
				for ( Size repeat_mult=1; repeat_mult < num_quasi_repeats(); repeat_mult++ ) {
					ir_plus_Xsubunit = ir + ( repeat_mult * ( num_resis_quasi_subunits/num_quasi_repeats() ) );
					if ( ( indy_resis[ir_plus_Xsubunit] )                             //ir_plus is a independent residue
							&& ( get_component_of_residue(pose,ir_plus_Xsubunit) == quasi_symm_comp() ) ) {     //ir_plus is a quasi-residue
						links->set_equiv( ir_plus_Xsubunit, list );                        //note: this loop links each residue AND it's respective +X resis to the same list.
					}
				}
				TR << "RotamerLinks set for residue " << ir << ": sele chain " << get_component_of_residue(pose,ir) << " and resi ";
				for ( core::Size const & i : list ) {
					TR << i << "+";
				} TR << std::endl;
			} else if ( pose.residue(ir).is_protein()                 //residue is a protein
					//residues in linked but is NOT the first quasi subunit
					&& ( get_component_of_residue(pose,ir) == quasi_symm_comp() ) //residue is in a quasi comp
					&& ( int(ir) > int(first_quasi_repeat_last_resi) )      //residue is NOT in the first quasi repeat
					&& ( int(ir) > int(offset_resis()) ) ) {           //residue is not a offset residue

				TR.Debug << "RotamerLinks ALREADY set for residue " << ir << " from a previous quasi residue." << std::endl;
			} else {
				//this should catch design residues that for some reason are not being linked to something
				TR.Warning << "This residue: " << ir << " has NOT been assigned to anything. Something is WRONG." << std::endl;
			}

		} else {
			//all non-design residues are prevented from repacking
			task->nonconst_residue_task(ir).prevent_repacking();
			if ( indy_resis[ir] ) {
				TR.Debug << "Residue: " << ir << " is not a design_pos. It will not be designed." << std::endl;
			}
		}
	}
	task->rotamer_links( links ); //send to rotamer_links the final residue link list
	TR << output << std::endl;  //prints final list of design position string

	// Store the task
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // Does this need to be fixed or omitted?
	}

	// prints current pose residue information
	task->show( TR.Debug );
	TR.Debug.flush();

	// If the pose doesn't have STM_STORED_TASK data, put a blank STMStoredTask in there.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) {
		//ORIGINAL// protocols::toolbox::task_operations::STMStoredTaskOP blank_tasks = new protocols::toolbox::task_operations::STMStoredTask();
		protocols::toolbox::task_operations::STMStoredTaskOP blank_tasks( new protocols::toolbox::task_operations::STMStoredTask() );
		pose.data().set( core::pose::datacache::CacheableDataType::STM_STORED_TASKS, blank_tasks );
	}
	// Grab a reference to the data
	//ORIGINAL//   protocols::toolbox::task_operations::STMStoredTask & stored_tasks = *(                           static_cast< protocols::toolbox::task_operations::STMStoredTask* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS )() ) );
	protocols::toolbox::task_operations::STMStoredTask & stored_tasks = *( utility::pointer::static_pointer_cast< protocols::toolbox::task_operations::STMStoredTask > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) );
	//devel::matdes::STMStoredTask & stored_tasks = *( static_cast< devel::matdes::STMStoredTask* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS )() ) );

	// If you haven't set overwrite to true and your task name already exists, fail. Otherwise, put the task you've made into the data cache.
	if ( overwrite_ || !stored_tasks.has_task(task_name_) ) {
		stored_tasks.set_task( task, task_name_ );
	} else {
		utility_exit_with_message("A stored task with the name " + task_name_ + " already exists; you must set overwrite flag to true to overwrite." );
	}
}

void
StoreQuasiSymmetricTaskMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data_map );
	task_name_ = tag->getOption< std::string >( "task_name" );
	overwrite_ = tag->getOption< bool >( "overwrite", false );
	quasi_symm_comp( tag->getOption< std::string >( "quasi_symm_comp", "B" ) );
	num_quasi_repeats( tag->getOption< core::Size >( "num_quasi_repeats", 2 ) );
	offset_resis( tag->getOption< core::Size >( "offset_resis", 0 ) );
}

// @brief Identification
// XRW TEMP std::string StoreQuasiSymmetricTaskMoverCreator::keyname() const { return StoreQuasiSymmetricTaskMover::mover_name(); }
// XRW TEMP std::string StoreQuasiSymmetricTaskMover::mover_name() { return "StoreQuasiSymmetricTaskMover"; }
// XRW TEMP std::string StoreQuasiSymmetricTaskMover::get_name() const { return "StoreQuasiSymmetricTaskMover"; }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP StoreQuasiSymmetricTaskMoverCreator::create_mover() const {
// XRW TEMP  //return new StoreQuasiSymmetricTaskMover;
// XRW TEMP  return protocols::moves::MoverOP( new StoreQuasiSymmetricTaskMover );
// XRW TEMP }

protocols::moves::MoverOP
StoreQuasiSymmetricTaskMover::clone() const {
	//return new StoreQuasiSymmetricTaskMover( *this );
	return protocols::moves::MoverOP( new StoreQuasiSymmetricTaskMover( *this ) );
}

protocols::moves::MoverOP
StoreQuasiSymmetricTaskMover::fresh_instance() const {
	//return new StoreQuasiSymmetricTaskMover;
	return protocols::moves::MoverOP( new StoreQuasiSymmetricTaskMover );
}

std::string StoreQuasiSymmetricTaskMover::get_name() const {
	return mover_name();
}

std::string StoreQuasiSymmetricTaskMover::mover_name() {
	return "StoreQuasiSymmetricTaskMover";
}

void StoreQuasiSymmetricTaskMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "task_name", xsct_pose_cached_task_operation, "The name of the stored task to be cached in the Pose object's datacache.  The object can later be retrieved using this name." )
		+ XMLSchemaAttribute::attribute_w_default( "overwrite", xsct_rosetta_bool, "Will overwrite old stored tasks with the same name?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "quasi_symm_comp", xs_string, "Which component (if multi-component, A or B), is going to be quasi-equivalent." , "B" )
		+ XMLSchemaAttribute::attribute_w_default( "num_quasi_repeats", xsct_non_negative_integer, "How many subunits your quasi-equivalent building block consists of. For example, a trimer would be 3.", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "offset_resis", xsct_non_negative_integer, "If your building block is non-quasi-equivalent domains, this number denotes the number of residues to skip. Currently needs to be resis 1-x, the skipped portion must be on the N-terminus.", "0" ) ;

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This mover creates a stored task that links selected residues with RotamerLinks. The residues will remain identical in identity when the packer is called, but their rotamers are free to be packed differently. This mover was designed to take a shot at the quasi-equivalent design problem, where identical residues need to satisfy multiple interfaces at the same time. It is essentially a multi-state design problem.", attlist );
}

std::string StoreQuasiSymmetricTaskMoverCreator::keyname() const {
	return StoreQuasiSymmetricTaskMover::mover_name();
}

protocols::moves::MoverOP
StoreQuasiSymmetricTaskMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new StoreQuasiSymmetricTaskMover );
}

void StoreQuasiSymmetricTaskMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StoreQuasiSymmetricTaskMover::provide_xml_schema( xsd );
}


} // matdes
} // devel

