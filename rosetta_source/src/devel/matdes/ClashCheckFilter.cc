// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/ClashCheckFilter.cc
/// @brief  Checks to see if any residues defined by the TaskOperations have clashing CA/CBs between interfaces (building block aware).
/// @author Jacob Bale (balej@u.washington.edu)

// Unit Headers
#include <devel/matdes/ClashCheckFilter.hh>
#include <devel/matdes/ClashCheckFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/sasa.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static basic::Tracer TR("devel.matdes.ClashCheckFilter");

namespace devel {
namespace matdes {

// @brief default constructor
ClashCheckFilter::ClashCheckFilter():
  task_factory_( NULL ),
  clash_dist_( 3.5 ),
	sym_dof_names_( "" ),
	nsub_bblock_( 1 )
{}

// @brief constructor with arguments
ClashCheckFilter::ClashCheckFilter( core::pack::task::TaskFactoryOP task_factory, core::Real const c, std::string const s, core::Size const n ):
	task_factory_( task_factory ),
	clash_dist_( c ),
	sym_dof_names_( s ),
	nsub_bblock_( n )
{}

// @brief copy constructor
ClashCheckFilter::ClashCheckFilter( ClashCheckFilter const & rval ):
	Super( rval ),
	task_factory_( rval.task_factory_ ),
	clash_dist_( rval.clash_dist_ ),
	sym_dof_names_( rval.sym_dof_names_ ),
	nsub_bblock_( rval.nsub_bblock_ )
{}

// @brief destructor
ClashCheckFilter::~ClashCheckFilter() {}

protocols::filters::FilterOP
ClashCheckFilter::fresh_instance() const{
  return new ClashCheckFilter();
}

protocols::filters::FilterOP
ClashCheckFilter::clone() const{
  return new ClashCheckFilter( *this );
}

// @brief getters
core::pack::task::TaskFactoryOP ClashCheckFilter::task_factory() const { return task_factory_; }
core::Real ClashCheckFilter::clash_dist() const { return clash_dist_; }
std::string ClashCheckFilter::sym_dof_names() const { return sym_dof_names_; }
core::Size ClashCheckFilter::nsub_bblock() const { return nsub_bblock_; }

// @brief setters
void ClashCheckFilter::task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
void ClashCheckFilter::clash_dist( core::Real const c ) { clash_dist_ = c; }
void ClashCheckFilter::sym_dof_names( std::string const s ) { sym_dof_names_ = s; }
void ClashCheckFilter::nsub_bblock( core::Size const n ) { nsub_bblock_ = n; }

/// @brief
core::Real ClashCheckFilter::compute( Pose const & pose ) const
{
	using namespace core;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace utility;
	typedef vector1<Size> Sizes;

	SymmetryInfoCOP sym_info = symmetry_info(pose);
	vector1<bool> indy_resis = sym_info->independent_residues();
	Real const clash_dist_sq = clash_dist_ * clash_dist_;
	vector1<std::string> sym_dof_name_list = string_split( sym_dof_names_ , ',' );

	Sizes intra_subs1, intra_subs2;

	if( sym_dof_name_list.size() == 2) {
	intra_subs1 = get_jump_name_to_subunits(pose,sym_dof_name_list[1]);
	intra_subs2 = get_jump_name_to_subunits(pose,sym_dof_name_list[2]);
	}

	runtime_assert( task_factory() );
  core::pack::task::PackerTaskCOP packer_task( task_factory()->create_task_and_apply_taskoperations( pose ) );

	utility::vector1<Real> clashing_pos;
	std::string select_clashing_pos("select clashing_pos, resi ");

	for(Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		if(sym_info->subunit_index(ir) != 1) continue;
		std::string atom_i = (pose.residue(ir).name3() == "GLY") ? "CA" : "CB";
		for(Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
			std::string atom_j = (pose.residue(jr).name3() == "GLY") ? "CA" : "CB";
			//If one component, then check for clashes between all residues in primary subunit and subunits with indices > nsub_bb
			if( sym_dof_names_ == "" ) {
      	if ( sym_info->subunit_index(jr) <= nsub_bblock_ ) continue;
			}
			//If two component, then check for clashes between all residues in primary subunitA and other building blocks, and all resis in primary subB and other building blocks. 
			else if( sym_dof_name_list.size() == 2 ) {
				Sizes const & isubs( get_component_of_residue(pose,ir)=='A'?intra_subs1:intra_subs2);
				if(get_component_of_residue(pose,ir)==get_component_of_residue(pose,jr)&&find(isubs.begin(),isubs.end(),sym_info->subunit_index(jr))!=isubs.end()) continue;
			} else {
				utility_exit_with_message("Clash check filter currently only works for 1 or 2 component symmetries");
			}	
			if(pose.residue(ir).xyz(atom_i).distance_squared(pose.residue(jr).xyz(atom_j)) <= clash_dist_sq) {
				clashing_pos.push_back(ir);
				select_clashing_pos.append(ObjexxFCL::string_of(ir) + "+");   
				break;
			}
		}
	}
	TR << select_clashing_pos << std::endl;
  return( clashing_pos.size() );

} // compute

// @brief returns true if the set of residues defined by the TaskOperations have a no clashes. False otherwise.
bool ClashCheckFilter::apply( Pose const & pose ) const
{
	// Get the number of clashes from the compute function and filter
  core::Real const clashes( compute( pose ) );

	TR<<"# clashing residues: "<< clashes <<". ";
	if( clashes > 0 ){
		TR<<"failing."<<std::endl;
		return false;
	}
	else {
		TR<<"passing."<<std::endl;
		return true;
	}

}

/// @brief parse xml
void
ClashCheckFilter::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
  task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
  clash_dist( tag->getOption< core::Real >( "clash_dist", 3.5 ) );
	sym_dof_names_ = tag->getOption< std::string >( "sym_dof_names", "" );
  nsub_bblock_ = tag->getOption<core::Size>("nsub_bblock", 1);
}

core::Real
ClashCheckFilter::report_sm( core::pose::Pose const & pose ) const
{
  return( compute( pose ) );
} 

void
ClashCheckFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
  out << "ClashCheckFilter returns " << compute( pose ) << std::endl;
}

protocols::filters::FilterOP
ClashCheckFilterCreator::create_filter() const { return new ClashCheckFilter; }

std::string
ClashCheckFilterCreator::keyname() const { return "ClashCheck"; }


} // matdes
} // devel
