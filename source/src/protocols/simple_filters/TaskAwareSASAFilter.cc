// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/TaskAwareSASAFilter.cc
/// @brief  Calculates SASA for a set of residues defined by TaskOperations
/// @author Neil King (neilking@u.washington.edu)

// Unit Headers
#include <protocols/simple_filters/TaskAwareSASAFilter.hh>
#include <protocols/simple_filters/TaskAwareSASAFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/sasa.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.TaskAwareSASAFilter" );

namespace protocols {
namespace simple_filters {

// @brief default constructor
TaskAwareSASAFilter::TaskAwareSASAFilter():
	task_factory_( /* NULL */ ),
	threshold_( 0 ),
	designable_only_( false ),
	sc_only_( false ),
	probe_radius_( 2.2 ),
	jump_id_( 0 )
{}

// @brief constructor with arguments
TaskAwareSASAFilter::TaskAwareSASAFilter( core::pack::task::TaskFactoryOP task_factory, core::Real const t, bool const d, bool const s, core::Real const r, core::Size const j ):
	task_factory_( task_factory ),
	threshold_( t ),
	designable_only_( d ),
	sc_only_( s ),
	probe_radius_( r ),
	jump_id_( j )
{}

// @brief copy constructor
TaskAwareSASAFilter::TaskAwareSASAFilter( TaskAwareSASAFilter const & rval ):
	Super( rval ),
	task_factory_( rval.task_factory_ ),
	threshold_( rval.threshold_ ),
	designable_only_( rval.designable_only_),
	sc_only_( rval.sc_only_ ),
	probe_radius_( rval.probe_radius_ ),
	jump_id_( rval.jump_id_ )
{}

// @brief destructor
TaskAwareSASAFilter::~TaskAwareSASAFilter() {}

protocols::filters::FilterOP
TaskAwareSASAFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new TaskAwareSASAFilter() );
}

protocols::filters::FilterOP
TaskAwareSASAFilter::clone() const{
	return protocols::filters::FilterOP( new TaskAwareSASAFilter( *this ) );
}

// @brief getters
core::pack::task::TaskFactoryOP TaskAwareSASAFilter::task_factory() const { return task_factory_; }
core::Real TaskAwareSASAFilter::threshold() const { return threshold_; }
bool TaskAwareSASAFilter::designable_only() const { return designable_only_; }
bool TaskAwareSASAFilter::sc_only() const { return sc_only_; }
core::Real TaskAwareSASAFilter::probe_radius() const { return probe_radius_; }
core::Size TaskAwareSASAFilter::jump_id() const { return jump_id_; }

// @brief setters
void TaskAwareSASAFilter::task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
void TaskAwareSASAFilter::threshold( core::Real const t ) { threshold_ = t; }
void TaskAwareSASAFilter::designable_only( bool const d ) { designable_only_ = d; }
void TaskAwareSASAFilter::sc_only( bool const s ) { sc_only_ = s; }
void TaskAwareSASAFilter::probe_radius( core::Real const r ) { probe_radius_ = r; }
void TaskAwareSASAFilter::jump_id( core::Size const j ) { jump_id_ = j; }

/// @brief
core::Real TaskAwareSASAFilter::compute( Pose const & p, bool const verbose ) const
{

	core::pose::Pose pose = p;

	runtime_assert( task_factory() != 0 );
	core::pack::task::PackerTaskCOP packer_task( task_factory()->create_task_and_apply_taskoperations( pose ) );

	// If a jump has been provided by the user, separate the pose by that jump.
	if ( jump_id() != 0 ) {
		int sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num(pose, jump_id() );
		protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( pose, sym_aware_jump_id ) );
		translate->step_size( 1000.0 );
		translate->apply( pose );
	}

	pose.dump_pdb("test.pdb");

	// Calculate SASA for each of the selected residues and spit that out to the log file.
	// Also add it to the total value.
	core::Real combined_sasa = 0;
	utility::vector1<Real> rsd_sasa(pose.n_residue(),0.0);
	core::id::AtomID_Map<Real> atom_sasa;
	core::id::AtomID_Map<bool> atom_mask;
	core::pose::initialize_atomid_map(atom_sasa,pose,0.0);
	core::pose::initialize_atomid_map(atom_mask,pose,false);
	for ( core::Size resi=1; resi<=pose.total_residue(); ++resi ) {
		core::conformation::Residue const rsd( pose.residue( resi ) );
		if ( rsd.type().name() == "VRT" ) { continue; }
		if ( packer_task->being_designed( resi ) || ( !designable_only() && packer_task->being_packed( resi ) ) ) {
			for ( core::Size atomi = 1; atomi <= pose.residue(resi).nheavyatoms(); atomi++ ) {
				atom_mask[core::id::AtomID(atomi,resi)] = true;
			}
		}
	}
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius(), false);//, atom_mask );
	utility::vector1<Real> resi_sasa(pose.n_residue(),0.0);

	for ( core::Size resi=1; resi<=pose.total_residue(); ++resi ) {
		core::conformation::Residue const rsd( pose.residue( resi ) );
		if ( rsd.type().name() == "VRT" ) { continue; }
		if ( packer_task->being_designed( resi ) || ( !designable_only() && packer_task->being_packed( resi ) ) ) {
			if ( sc_only() ) {
				// Use CA as the side chain for Glys
				if ( pose.residue(resi).name3()=="GLY" ) resi_sasa[resi] += atom_sasa[core::id::AtomID(2,resi)];
				for ( Size j=5; j <= pose.residue(resi).nheavyatoms(); j++ ) {
					resi_sasa[resi] += atom_sasa[core::id::AtomID(j,resi)];
				}
			} else {
				resi_sasa[resi] += rsd_sasa[resi];
			}
			if ( verbose ) { TR << "SASA of " << pose.residue(resi).name3() << resi << " is " << resi_sasa[resi] << std::endl; }
			combined_sasa += resi_sasa[resi];
		}
	}
	return( combined_sasa );

} // compute

// @brief returns true if the set of residues defined by the TaskOperations have a combined
// SASA higher than the threshold, false otherwise.
// Default threshold is 0 so that all structures will pass if a threshold is not explicitly
// specified.
bool TaskAwareSASAFilter::apply( Pose const & pose ) const
{
	// Get the combined sasa from the compute function and filter
	core::Real const combined_sasa( compute( pose, true ) );
	return( combined_sasa >= threshold() );
}

/// @brief parse xml
void
TaskAwareSASAFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	threshold( tag->getOption< core::Real >( "threshold", 0 ) );
	designable_only( tag->getOption< bool >( "designable_only", false ) );
	sc_only( tag->getOption< bool >( "sc_only", false ) );
	probe_radius( tag->getOption< core::Real >( "probe_radius", 2.2 ) );
	jump_id( tag->getOption< core::Size >( "jump", 0 ) );
}

core::Real
TaskAwareSASAFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose, false ) );
}

void
TaskAwareSASAFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out << "TaskAwareSASAFilter returns " << compute( pose, false ) << std::endl;
}

protocols::filters::FilterOP
TaskAwareSASAFilterCreator::create_filter() const { return protocols::filters::FilterOP( new TaskAwareSASAFilter ); }

std::string
TaskAwareSASAFilterCreator::keyname() const { return "TaskAwareSASA"; }


} // simple_filters
} // protocols
