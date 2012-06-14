// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/SaveResfileToDiskFilter.cc
/// @brief  Calculates AverageDegree within the context of an unbound oligomer
/// @author Neil King (neilking@u.washington.edu)

// Unit Headers
#include <devel/matdes/SaveResfileToDiskFilter.hh>
#include <devel/matdes/SaveResfileToDiskFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/conformation/Conformation.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <fstream>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static basic::Tracer TR("devel.matdes.SaveResfileToDiskFilter");

namespace devel {
namespace matdes {

// @brief default constructor
SaveResfileToDiskFilter::SaveResfileToDiskFilter() {}

// @brief constructor with arguments
SaveResfileToDiskFilter::SaveResfileToDiskFilter( core::pack::task::TaskFactoryOP task_factory, utility::vector1<core::Size> r, bool d, std::string n, std::string g ):
	task_factory_( task_factory ),
	selected_resis_( r ),
	designable_only_( d ),
	resfile_name_( n ),
	resfile_general_property_( g )
{}

// @brief copy constructor
SaveResfileToDiskFilter::SaveResfileToDiskFilter( SaveResfileToDiskFilter const & rval ):
	Super( rval ),
	task_factory_( rval.task_factory_ ),
	selected_resis_( rval.selected_resis_ ),
	designable_only_( rval.designable_only_ ),
	resfile_name_( rval.resfile_name_ ),
	resfile_general_property_( rval.resfile_general_property_ )
{}

protocols::filters::FilterOP
SaveResfileToDiskFilter::fresh_instance() const{
  return new SaveResfileToDiskFilter();
}

protocols::filters::FilterOP
SaveResfileToDiskFilter::clone() const{
  return new SaveResfileToDiskFilter( *this );
}

// @brief getters
core::pack::task::TaskFactoryOP SaveResfileToDiskFilter::task_factory() const { return task_factory_; }
utility::vector1<core::Size> SaveResfileToDiskFilter::selected_resis() const { return selected_resis_; }
bool SaveResfileToDiskFilter::designable_only() const { return designable_only_; }
std::string SaveResfileToDiskFilter::resfile_name() const { return resfile_name_; }
std::string SaveResfileToDiskFilter::resfile_general_property() const { return resfile_general_property_; }

// @brief setters
void SaveResfileToDiskFilter::task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
void SaveResfileToDiskFilter::selected_resis( utility::vector1<core::Size> const r ) { selected_resis_ = r; }
void SaveResfileToDiskFilter::designable_only( bool const d ) { designable_only_ = d; }
void SaveResfileToDiskFilter::resfile_name( std::string const n ) { resfile_name_ = n; }
void SaveResfileToDiskFilter::resfile_general_property( std::string const g ) { resfile_general_property_ = g; }

/// @brief Applies the TaskOperations specified in the xml, and then selects either
// the repackable or designable residues depending on what the user specifies
utility::vector1< core::Size > SaveResfileToDiskFilter::select_residues( Pose const & pose ) const
{

	utility::vector1< core::Size > selected_residues;

	// Prepare the PackerTask
	runtime_assert( task_factory() );
  core::pack::task::PackerTaskCOP task( task_factory()->create_task_and_apply_taskoperations( pose ) );

	// Find out which residues are packable or designable
	for( core::Size resi = 1; resi <= pose.total_residue(); ++resi ) {
		if ( designable_only_ ) {
    	if( task->residue_task( resi ).being_designed() && pose.residue(resi).is_protein() )
      	selected_residues.push_back( resi );
		} else {
      if( task->residue_task( resi ).being_packed() && pose.residue(resi).is_protein() )
        selected_residues.push_back( resi );
		}
  }
  if( selected_residues.empty() )
    TR.Warning << "WARNING: No residues were selected by your TaskOperations." << std::endl;

  return selected_residues;

} // selectResidues

/// @brief Write the resfile to disk
void
SaveResfileToDiskFilter::write_resfile( Pose const & pose, utility::vector1< core::Size > const & selected_residues ) const
{
	std::string resfile_to_write = resfile_name();

	if (resfile_to_write == "") {
		resfile_to_write = protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".resfile";
	}

  runtime_assert( resfile_to_write != "" );
  runtime_assert( !selected_residues.empty() );
	std::ofstream resfile;
  resfile.open( resfile_to_write.c_str(), std::ios::out );
  resfile << resfile_general_property() << "\nstart\n";
	for ( core::Size i=1; i<=selected_residues.size(); i++ ) {
		resfile << selected_residues[i] << '\t' << pose.pdb_info()->chain(selected_residues[i]) << " PIKAA " << pose.residue(selected_residues[i]).name1() << '\n';
  }
  resfile.close();

}

// @ brief
bool SaveResfileToDiskFilter::apply( Pose const & pose ) const
{

  utility::vector1< core::Size > selected_residues = select_residues( pose );
	if ( !selected_residues.empty() )
		write_resfile( pose, selected_residues );

	return 1;

} // apply

/// @brief parse xml
void
SaveResfileToDiskFilter::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
  task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	designable_only( tag->getOption< bool >( "designable_only", false ) );
	resfile_name( tag->getOption< std::string >( "resfile_name", "" ) );
	resfile_general_property( tag->getOption< std::string >( "resfile_general_property", "NATAA" ) );
}
/*
core::Real
SaveResfileToDiskFilter::report_sm( core::pose::Pose const & pose ) const
{
  return( compute( pose ) );
} 

void
SaveResfileToDiskFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
  out << "SaveResfileToDiskFilter returns " << compute( pose ) << std::endl;
}
*/
protocols::filters::FilterOP
SaveResfileToDiskFilterCreator::create_filter() const { return new SaveResfileToDiskFilter; }

std::string
SaveResfileToDiskFilterCreator::keyname() const { return "SaveResfileToDisk"; }


} // matdes
} // devel
