// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/matdes/StoreQuasiSymmetricTaskMover.cc
/// @brief The StoreQuasiSymmetricTaskMover accepts a set of TaskOperations and ensures that
/// equivalent positions that are designable in one quasisymmetric subunit are designable in all.
/// In addition, this is where the RotamerLinks that ensure that these positions maintain the
/// same AA identity are applied.
/// @author Neil King (neilking@uw.edu)

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
#include <devel/matdes/STMStoredTask.hh>
#include <basic/Tracer.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <ObjexxFCL/format.hh>

#include <utility/tag/Tag.hh>

// C++ Headers

// ObjexxFCL Headers

static thread_local basic::Tracer TR( "devel.matdes.StoreQuasiSymmetricTaskMover" );

namespace devel {
namespace matdes {

// @brief default constructor
StoreQuasiSymmetricTaskMover::StoreQuasiSymmetricTaskMover() {}

// @brief destructor
StoreQuasiSymmetricTaskMover::~StoreQuasiSymmetricTaskMover() {}

// @brief getters
char StoreQuasiSymmetricTaskMover::quasi_symm_comp() const { return quasi_symm_comp_[0]; }
core::Size StoreQuasiSymmetricTaskMover::num_quasi_repeats() const { return num_quasi_repeats_; }

// @brief setters
void StoreQuasiSymmetricTaskMover::quasi_symm_comp( std::string const quasi_symm_comp ) { quasi_symm_comp_ = quasi_symm_comp; }
void StoreQuasiSymmetricTaskMover::num_quasi_repeats( core::Size const num_quasi_repeats ) { num_quasi_repeats_ = num_quasi_repeats; }

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
	runtime_assert( task_factory_ );
	core::pack::task::PackerTaskOP raw_task = task_factory_->create_task_and_apply_taskoperations( pose );

	// Get the designable positions from the raw task.
	// Also get the number of residues in the quasisymmetric subunits and in the
  // asymmetric unit as a whole.
	std::set< core::Size > design_pos;
  core::Size num_resis_quasi_subunits = 0;
  core::Size num_indy_resis = 0;
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	vector1<bool> indy_resis = sym_info->independent_residues();
  for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
    if ( indy_resis[ir] ) {
      num_indy_resis++;
      if ( ( pose.residue(ir).is_protein() ) && ( get_component_of_residue(pose,ir) == quasi_symm_comp() ) ) {
      num_resis_quasi_subunits++;
      }
			if ( raw_task->being_packed( ir ) ) {
				design_pos.insert( ir );
			} 
    }
  }
	//std::cout << "num_indy_resis = " << num_indy_resis << std::endl;
  //std::cout<< "num_resis_quasi_subunits = " << num_resis_quasi_subunits << std::endl;

	// Iterate through the designable positions and make sure that if a position is designable in one
	// quasisymmetric subunit, it is also designable in the other
  std::set< core::Size >::iterator pos;
  for ( pos=design_pos.begin(); pos!=design_pos.end(); ++pos ) {
    if ( indy_resis[*pos] == 0 ) continue; // Make sure we're only looking in the asymmetric unit
    if ( get_component_of_residue(pose,*pos) != quasi_symm_comp() ) continue; // Only look in the quasisymmetric subunits
    int pos_minus_1subunit = (int)*pos-((int)num_resis_quasi_subunits/int(num_quasi_repeats()));
    int pos_plus_1subunit = (int)*pos+((int)num_resis_quasi_subunits/int(num_quasi_repeats()));
    if ( pos_minus_1subunit > 0 && indy_resis[ pos_minus_1subunit ] && get_component_of_residue(pose,pos_minus_1subunit) == quasi_symm_comp() ) {
      if ( design_pos.find( pos_minus_1subunit ) == design_pos.end() ) {
        design_pos.insert( pos_minus_1subunit );
        TR << "Position " << pos_minus_1subunit << " inserted on behalf of position " << *pos << std::endl;
      }
    } else if ( pos_plus_1subunit <= num_indy_resis && indy_resis[ pos_plus_1subunit ] && get_component_of_residue(pose,pos_plus_1subunit) == quasi_symm_comp() ) {
      if ( design_pos.find( pos_plus_1subunit ) == design_pos.end() ) {
        design_pos.insert( pos_plus_1subunit );
        TR << "Position " << pos_plus_1subunit << " inserted on behalf of position " << *pos << std::endl;
      }
		}
	}

	// Create a new task; this is the one we will store.
  // Prevent_repacking at all positions that are not design positions
  // and apply rotamer links to all designable positions.
	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
  core::pack::rotamer_set::RotamerLinksOP links( new core::pack::rotamer_set::RotamerLinks );
  links->resize( num_indy_resis );
  std::string output = "Final design_pos ";
  for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
    if (design_pos.find(ir) != design_pos.end()) {
      output += ObjexxFCL::string_of(ir)+"+";
      core::Size ir_plus_1subunit = ir+(num_resis_quasi_subunits/num_quasi_repeats());
      if ( pose.residue(ir).is_protein() && get_component_of_residue(pose,ir) != quasi_symm_comp() ) {
        links->set_equiv( ir, ir ); // For non-quasisymmetrical subunits, set up a dummy RotamerLink to this resi
      } else if ( pose.residue(ir).is_protein()
        && ( get_component_of_residue(pose,ir) == quasi_symm_comp() )
        && ( indy_resis[ir_plus_1subunit] )
        && ( get_component_of_residue(pose,ir_plus_1subunit) == quasi_symm_comp() ) )
      {
        utility::vector1< core::Size > list; list.push_back( ir ); list.push_back( ir_plus_1subunit );
        links->set_equiv( ir, list ); links->set_equiv( ir_plus_1subunit, list );
        TR.Debug << "Link set between residues " << ir << " and " << ir_plus_1subunit << std::endl;
      }
    } else {
      TR.Debug << "resi " << ir << " will not be designed" << std::endl;
      task->nonconst_residue_task(ir).prevent_repacking();
    }
  }
  task->rotamer_links( links );
  TR << output << std::endl;

	// Store the task
	if (core::pose::symmetry::is_symmetric(pose))
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // Does this need to be fixed or omitted?
	// If the pose doesn't have STM_STORED_TASK data, put a blank STMStoredTask in there.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) {
		devel::matdes::STMStoredTaskOP blank_tasks = new devel::matdes::STMStoredTask();
		pose.data().set( core::pose::datacache::CacheableDataType::STM_STORED_TASKS, blank_tasks );
	}
	// Grab a reference to the data
	devel::matdes::STMStoredTask & stored_tasks = *( static_cast< devel::matdes::STMStoredTask* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS )() ) );
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

}

// @brief Identification
std::string StoreQuasiSymmetricTaskMoverCreator::keyname() const { return StoreQuasiSymmetricTaskMoverCreator::mover_name(); }
std::string StoreQuasiSymmetricTaskMoverCreator::mover_name() { return "StoreQuasiSymmetricTaskMover"; }
std::string StoreQuasiSymmetricTaskMover::get_name() const { return "StoreQuasiSymmetricTaskMover"; }

protocols::moves::MoverOP
StoreQuasiSymmetricTaskMoverCreator::create_mover() const {
	return new StoreQuasiSymmetricTaskMover;
}

protocols::moves::MoverOP
StoreQuasiSymmetricTaskMover::clone() const {
	return new StoreQuasiSymmetricTaskMover( *this );
}

protocols::moves::MoverOP
StoreQuasiSymmetricTaskMover::fresh_instance() const {
	return new StoreQuasiSymmetricTaskMover;
}

} // matdes
} // devel

