// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/matdes/ClashCheckFilter.cc
/// @brief  Checks to see if any residues defined by the TaskOperations have clashing CA/CBs between interfaces (building block aware).
/// @author Jacob Bale (balej@u.washington.edu)

// Unit Headers
#include <protocols/matdes/ClashCheckFilter.hh>
#include <protocols/matdes/ClashCheckFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/sasa.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static basic::Tracer TR("protocols.matdes.ClashCheckFilter");

namespace protocols {
namespace matdes {

// @brief default constructor
ClashCheckFilter::ClashCheckFilter():
  task_factory_( NULL ),
  clash_dist_( 3.5 ),
	sym_dof_names_( "" ),
	nsub_bblock_( 1 ),
	threshold_( 0 ),
	verbose_( 0 ),
	write_( 0 )
{}

// @brief constructor with arguments
ClashCheckFilter::ClashCheckFilter( core::pack::task::TaskFactoryOP task_factory, core::Real const c, std::string const s, core::Size const n, core::Size const t, bool const v, bool const w ):
	task_factory_( task_factory ),
	clash_dist_( c ),
	sym_dof_names_( s ),
	nsub_bblock_( n ),
	threshold_( t ),
	verbose_( v ),
	write_( w )
{}

// @brief copy constructor
ClashCheckFilter::ClashCheckFilter( ClashCheckFilter const & rval ):
	Super( rval ),
	task_factory_( rval.task_factory_ ),
	clash_dist_( rval.clash_dist_ ),
	sym_dof_names_( rval.sym_dof_names_ ),
	nsub_bblock_( rval.nsub_bblock_ ),
	threshold_( rval.threshold_ ),
	verbose_( rval.verbose_ ),
	write_( rval.write_ )
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
core::Size ClashCheckFilter::threshold() const { return threshold_; }
bool ClashCheckFilter::verbose() const { return verbose_; }
bool ClashCheckFilter::write() const { return write_; }

// @brief setters
void ClashCheckFilter::task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }
void ClashCheckFilter::clash_dist( core::Real const c ) { clash_dist_ = c; }
void ClashCheckFilter::sym_dof_names( std::string const s ) { sym_dof_names_ = s; }
void ClashCheckFilter::nsub_bblock( core::Size const n ) { nsub_bblock_ = n; }
void ClashCheckFilter::threshold( core::Size const t ) { threshold_ = t; }
void ClashCheckFilter::verbose( bool const v ) { verbose_ = v; }
void ClashCheckFilter::write( bool const w ) { write_ = w; }

/// @brief
core::Size ClashCheckFilter::compute( Pose const & pose, bool const & v, bool const & w ) const
{
	using namespace core;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using core::id::AtomID;
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
	std::string select_clashing_pos("select clashing_pos, (");
  Size itype = 5;
  Size jtype = 5;
	bool clash = false;

	for(Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		if(sym_info->subunit_index(ir) != 1) continue;
		clash = false;
		for(Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
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
			if (chemical::name_from_aa(pose.aa(ir)) == "GLY") { // Use CAs instead of CBs for GLYs
				itype = 4;
			} else {
				itype = 5;
			}
      if (chemical::name_from_aa(pose.aa(jr)) == "GLY") {
        jtype = 4;
      } else {
        jtype = 5;
      }
  //    if (pose.xyz(AtomID(itype,ir)).distance_squared(pose.xyz(AtomID(jtype,jr))) <= contact_dist_sq) { // Check if the CBs/CAs are within contact dist
      for (Size ia=1; ia<=itype; ia++) {
        for (Size ja=1; ja<=jtype; ja++) {
          if (pose.xyz(AtomID(ia,ir)).distance_squared(pose.xyz(AtomID(ja,jr))) <= clash_dist_sq) { // Test for clashes
            if ( (((ia == 1) && (ja == 4)) || ((ia == 4) && (ja == 1))) && (pose.xyz(AtomID(ia,ir)).distance_squared(pose.xyz(AtomID(ja,jr))) >= 6.76) ) { // But don't count bb-bb h-bonds as clashes
              continue;
            } else {
							clash = true;
							clashing_pos.push_back(ir);
							if ( w ) { 
							write_to_pdb( pose, pose.residue(ir).name3(), ir, pose.residue(ir).atom_name(ia) ); 
							}
							core::Size output_resi = ir;
							if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
								output_resi = pose.pdb_info()->number( ir );
							}
							select_clashing_pos.append( " resi " + ObjexxFCL::string_of(output_resi) + " and name " + pose.residue(ir).atom_name(ia) + "+" );   
          	}
        	}		
					if ( clash ) break;
				}
				if ( clash ) break;
      }
//			}
			if (clash) break;
		}
	}
	if ( v ) {
		select_clashing_pos.erase(select_clashing_pos.end()-2,select_clashing_pos.end());
		TR << select_clashing_pos << ") and " << protocols::jd2::JobDistributor::get_instance()->current_output_name() << std::endl;
	}
	if ( w ) { 
		write_pymol_string_to_pdb( select_clashing_pos ); 
	}
  return( clashing_pos.size() );
} // compute

void ClashCheckFilter::write_to_pdb( core::pose::Pose const & pose, std::string const residue_name, core::Size const residue, std::string const atom_name ) const
{

	protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
	std::string filter_name = this->name();
	std::string user_name = this->get_user_defined_name();
	core::Size output_resi = residue;
	if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
		output_resi = pose.pdb_info()->number( residue );
	}
	std::string unsat_pols_string = filter_name + " " + user_name + ": " + residue_name + ObjexxFCL::string_of(output_resi) + " " + atom_name ;
	job->add_string(unsat_pols_string);
}

void ClashCheckFilter::write_pymol_string_to_pdb( std::string const pymol_selection ) const
{

	protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
	std::string filter_name = this->name();
	std::string user_name = this->get_user_defined_name();
	std::string pymol_string = filter_name + " " + user_name + ": " + pymol_selection + ") and " + protocols::jd2::JobDistributor::get_instance()->current_output_name() ;
	job->add_string(pymol_string);
}

// @brief returns true if the set of residues defined by the TaskOperations have a no clashes. False otherwise.
bool ClashCheckFilter::apply( Pose const & pose ) const
{
	// Get the number of clashes from the compute function and filter
  core::Size const clashes( compute( pose, verbose(), write() ) );

	TR<<"# clashing residues: "<< clashes <<". ";
	if( clashes > threshold_ ){
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
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
  task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
  clash_dist( tag->getOption< core::Real >( "clash_dist", 3.5 ) );
	sym_dof_names_ = tag->getOption< std::string >( "sym_dof_names", "" );
  nsub_bblock_ = tag->getOption<core::Size>("nsub_bblock", 1);
	threshold_ = tag->getOption<core::Size>( "cutoff", 0 );
	verbose_ = tag->getOption< bool >( "verbose", 0 );
	write_ = tag->getOption< bool >("write2pdb", 0);
}

core::Real
ClashCheckFilter::report_sm( core::pose::Pose const & pose ) const
{
  return( compute( pose, false, false ) );
} 

void
ClashCheckFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
  out << "ClashCheckFilter returns " << compute( pose, false, false ) << std::endl;
}

protocols::filters::FilterOP
ClashCheckFilterCreator::create_filter() const { return new ClashCheckFilter; }

std::string
ClashCheckFilterCreator::keyname() const { return "ClashCheck"; }


} // matdes
} // protocols
