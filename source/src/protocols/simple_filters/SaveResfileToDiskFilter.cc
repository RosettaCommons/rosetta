// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/SaveResfileToDiskFilter.cc
/// @brief  Outputs a resfile based on a given set of taskoperations
/// @author Neil King (neilking@u.washington.edu)

// Unit Headers
#include <protocols/simple_filters/SaveResfileToDiskFilter.hh>
#include <protocols/simple_filters/SaveResfileToDiskFilterCreator.hh>

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
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <fstream>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static THREAD_LOCAL basic::Tracer TR( "protocols.matdes.SaveResfileToDiskFilter" );

namespace protocols {
namespace simple_filters {

// @brief default constructor
SaveResfileToDiskFilter::SaveResfileToDiskFilter():
	designable_only_( false ),
	renumber_pdb_( false )
{}

// @brief constructor with arguments
SaveResfileToDiskFilter::SaveResfileToDiskFilter( core::pack::task::TaskFactoryOP task_factory, utility::vector1<core::Size> const & r, bool const d, std::string const & n, std::string const & s, std::string const & p, std::string const & g, std::string const &  srp ):
	task_factory_( task_factory ),
	selected_resis_( r ),
	designable_only_( d ),
	renumber_pdb_( false ),
	resfile_name_( n ),
	resfile_suffix_( s ),
	resfile_prefix_( p ),
	resfile_general_property_( g ),
	selected_resis_property_( srp )
{}

// @brief copy constructor
SaveResfileToDiskFilter::SaveResfileToDiskFilter( SaveResfileToDiskFilter const & rval ):
	Super( rval ),
	task_factory_( rval.task_factory_ ),
	selected_resis_( rval.selected_resis_ ),
	designable_only_( rval.designable_only_ ),
	renumber_pdb_( rval.renumber_pdb_ ),
	resfile_name_( rval.resfile_name_ ),
	resfile_suffix_( rval.resfile_suffix_ ),
	resfile_prefix_( rval.resfile_prefix_ ),
	resfile_general_property_( rval.resfile_general_property_ ),
	selected_resis_property_( rval.selected_resis_property_ )
{}

// @brief destructor
SaveResfileToDiskFilter::~SaveResfileToDiskFilter() {}

protocols::filters::FilterOP
SaveResfileToDiskFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new SaveResfileToDiskFilter() );
}

protocols::filters::FilterOP
SaveResfileToDiskFilter::clone() const{
	return protocols::filters::FilterOP( new SaveResfileToDiskFilter( *this ) );
}

// @brief getters
core::pack::task::TaskFactoryOP SaveResfileToDiskFilter::task_factory() const { return task_factory_; }
utility::vector1<core::Size> SaveResfileToDiskFilter::selected_resis() const { return selected_resis_; }
bool SaveResfileToDiskFilter::designable_only() const { return designable_only_; }
bool SaveResfileToDiskFilter::renumber_pdb() const { return renumber_pdb_; }
std::string SaveResfileToDiskFilter::resfile_name() const { return resfile_name_; }
std::string SaveResfileToDiskFilter::resfile_suffix() const { return resfile_suffix_; }
std::string SaveResfileToDiskFilter::resfile_prefix() const { return resfile_prefix_; }
std::string SaveResfileToDiskFilter::resfile_general_property() const { return resfile_general_property_; }
std::string SaveResfileToDiskFilter::selected_resis_property() const { return selected_resis_property_; }

// @brief setters
void SaveResfileToDiskFilter::task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
void SaveResfileToDiskFilter::selected_resis( utility::vector1<core::Size> const r ) { selected_resis_ = r; }
void SaveResfileToDiskFilter::designable_only( bool const d ) { designable_only_ = d; }
void SaveResfileToDiskFilter::renumber_pdb( bool const p ) { renumber_pdb_ = p; }
void SaveResfileToDiskFilter::resfile_name( std::string const & n ) { resfile_name_ = n; }
void SaveResfileToDiskFilter::resfile_suffix( std::string const & s ) { resfile_suffix_ = s; }
void SaveResfileToDiskFilter::resfile_prefix( std::string const & p ) { resfile_prefix_ = p; }
void SaveResfileToDiskFilter::resfile_general_property( std::string const & g ) { resfile_general_property_ = g; }
void SaveResfileToDiskFilter::selected_resis_property( std::string const & srp ) { selected_resis_property_ = srp; }

/// @brief Applies the TaskOperations specified in the xml, and then selects either
// the repackable or designable residues depending on what the user specifies
utility::vector1< core::Size > SaveResfileToDiskFilter::select_residues( Pose const & pose ) const
{

	utility::vector1< core::Size > selected_residues;

	// Prepare the PackerTask
	runtime_assert( task_factory() != 0 );
	core::pack::task::PackerTaskCOP task( task_factory()->create_task_and_apply_taskoperations( pose ) );

	// Find out which residues are packable or designable
	for ( core::Size resi = 1; resi <= pose.total_residue(); ++resi ) {
		if ( designable_only_ ) {
			if ( task->residue_task( resi ).being_designed() && pose.residue(resi).is_protein() ) {
				selected_residues.push_back( resi );
			}
		} else {
			if ( task->residue_task( resi ).being_packed() && pose.residue(resi).is_protein() ) {
				selected_residues.push_back( resi );
			}
		}
	}
	if ( selected_residues.empty() ) {
		TR.Warning << "WARNING: No residues were selected by your TaskOperations." << std::endl;
	}

	return selected_residues;

} // selectResidues

/// @brief Write the resfile to disk
void
SaveResfileToDiskFilter::write_resfile( Pose const & pose, utility::vector1< core::Size > const & selected_residues ) const
{
	std::string resfile_to_write = resfile_name();

	if ( resfile_suffix() == "" && resfile_prefix() == "" && resfile_to_write == "" ) {
		resfile_to_write = protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".resfile";
	}

	if ( (resfile_suffix() != "" || resfile_prefix() != "") && resfile_to_write == "" ) {
		resfile_to_write = resfile_prefix() + protocols::jd2::JobDistributor::get_instance()->current_output_name() + resfile_suffix() + ".resfile";
	}

	runtime_assert( resfile_to_write != "" );
	runtime_assert( !selected_residues.empty() );
	std::ofstream resfile;
	resfile.open( resfile_to_write.c_str(), std::ios::out );
	resfile << resfile_general_property() << "\nstart\n";
	for ( core::Size i=1; i<=selected_residues.size(); i++ ) {
		if ( renumber_pdb() ) {
			resfile << selected_residues[i] << '\t' << pose.pdb_info()->chain(selected_residues[i]);
		} else {
			resfile << pose.pdb_info()->number( selected_residues[i] ) << pose.pdb_info()->icode( selected_residues[i] ) << '\t' << pose.pdb_info()->chain(selected_residues[i]);
		}
		if ( selected_resis_property() != "" ) {
			resfile << " " << selected_resis_property() << '\n';
		} else {
			resfile << " PIKAA " << pose.residue(selected_residues[i]).name1() << '\n';
		}
	}
	resfile.close();

}

// @ brief
bool SaveResfileToDiskFilter::apply( Pose const & pose ) const
{

	utility::vector1< core::Size > selected_residues = select_residues( pose );
	if ( !selected_residues.empty() ) {
		write_resfile( pose, selected_residues );
	}

	return 1;

} // apply

/// @brief parse xml
void
SaveResfileToDiskFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	designable_only( tag->getOption< bool >( "designable_only", false ) );
	renumber_pdb( tag->getOption< bool >( "renumber_pdb", false ) );
	resfile_suffix( tag->getOption< std::string >( "resfile_suffix", "" ) );
	resfile_prefix( tag->getOption< std::string >( "resfile_prefix", "" ) );
	resfile_name( tag->getOption< std::string >( "resfile_name", "" ) );
	resfile_general_property( tag->getOption< std::string >( "resfile_general_property", "NATAA" ) );
	selected_resis_property( tag->getOption< std::string >( "selected_resis_property", "" ) );
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
SaveResfileToDiskFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SaveResfileToDiskFilter ); }

std::string
SaveResfileToDiskFilterCreator::keyname() const { return "SaveResfileToDisk"; }


} // simple_filters
} // protocols
