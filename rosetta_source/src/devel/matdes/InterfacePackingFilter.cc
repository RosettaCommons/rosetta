// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Jacob Bale (balej@uw.edu)
#include <devel/matdes/InterfacePackingFilter.hh>
#include <devel/matdes/InterfacePackingFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <basic/database/open.hh>
#include <core/id/AtomID_Map.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/symmetry/util.hh>
#include <devel/matdes/util.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/option.hh>

namespace devel {
namespace matdes {

static basic::Tracer TR( "devel.matdes.InterfacePackingFilter" );

///@brief default ctor
InterfacePackingFilter::InterfacePackingFilter() :
	parent( "InterfacePacking" ),
	distance_cutoff_( 9.0 ),
	lower_threshold_( -5 ),
	upper_threshold_( 5 )
{}

core::Real
InterfacePackingFilter::distance_cutoff() const{
	return distance_cutoff_;
}

core::Real
InterfacePackingFilter::lower_threshold() const{
	return lower_threshold_;
}

core::Real
InterfacePackingFilter::upper_threshold() const{
	return upper_threshold_;
}

void
InterfacePackingFilter::distance_cutoff( core::Real const d ){
	distance_cutoff_ = d;
}

void
InterfacePackingFilter::lower_threshold( core::Real const l ){
	lower_threshold_ = l;
}

void
InterfacePackingFilter::upper_threshold( core::Real const u ){
	upper_threshold_ = u;
}

bool
InterfacePackingFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real packing_score(compute( pose ));
	if( (packing_score >= lower_threshold_) && (packing_score <= upper_threshold_) ){
		TR<<"passing."<<std::endl;
		return true;
	} 
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

core::Real
InterfacePackingFilter::compute( core::pose::Pose const & pose ) const{

  utility::vector1<Size> intra_subs; 
	if (!basic::options::option[basic::options::OptionKeys::matdes::num_subs_building_block].user()) {
    utility_exit_with_message("ERROR: You have not set the required option -matdes::num_subs_building_block");
  } else {
		core::Size num_subs = basic::options::option[basic::options::OptionKeys::matdes::num_subs_building_block]();
    for(core::Size intrasub=1; intrasub<=num_subs; intrasub++) {
      intra_subs.push_back(intrasub);
    }
  }
  core::pose::Pose sub_pose = devel::matdes::get_neighbor_subs(pose, intra_subs);
  core::scoring::packing::HolesParams hp(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
  core::scoring::packing::HolesResult hr(core::scoring::packing::compute_holes_score(sub_pose, hp));
  core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
  core::Size nres_monomer = symm_info->num_independent_residues();
  core::Size count = 0; Real if_score = 0;

  core::Real cutoff2 = distance_cutoff_*distance_cutoff_;
  for (core::Size ir=1; ir<=nres_monomer; ir++) {
    for (core::Size ia = 1; ia<=sub_pose.residue(ir).nheavyatoms(); ia++) {
      bool contact = false;
      for (core::Size jr=nres_monomer+1; jr<=sub_pose.n_residue(); jr++) {
        for (core::Size ja = 1; ja<=sub_pose.residue(jr).nheavyatoms(); ja++) {
          if (sub_pose.residue(ir).xyz(ia).distance_squared(sub_pose.residue(jr).xyz(ja)) <= cutoff2)  {
            contact = true;
            break; // ja
          }
        } // ja
        if (contact == true) break;
      } // jr
      if (contact == true) {
        count++;
        if_score += hr.atom_scores[core::id::AtomID(ia, ir)];
      }
    } // ia
  } // ir
	TR << "if_score / count = " << if_score << " / " << count << " = " << (if_score / (Real)count) << std::endl;
  return if_score / (Real)count;
}

core::Real
InterfacePackingFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real packing_score(compute( pose ));
	return( packing_score );
}

void
InterfacePackingFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"InterfacePackingFilter returns "<<compute( pose )<<std::endl;
}

void
InterfacePackingFilter::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
	TR << "InterfacePackingFilter"<<std::endl;
	distance_cutoff( tag->getOption< core::Real >( "distance_cutoff", 9.0 ) );
	lower_threshold( tag->getOption< core::Real >( "lower_cutoff", -5 ) );
	upper_threshold( tag->getOption< core::Real >( "upper_cutoff", 5 ) );
	TR<<"with options lower_threshold: "<<lower_threshold()<<", upper_threshold: "<<upper_threshold()<<", and distance_cutoff: "<<distance_cutoff()<<std::endl;
}

protocols::filters::FilterOP
InterfacePackingFilter::fresh_instance() const{
	return new InterfacePackingFilter();
}

InterfacePackingFilter::~InterfacePackingFilter(){}

protocols::filters::FilterOP
InterfacePackingFilter::clone() const{
	return new InterfacePackingFilter( *this );
}

protocols::filters::FilterOP
InterfacePackingFilterCreator::create_filter() const { return new InterfacePackingFilter; }

std::string
InterfacePackingFilterCreator::keyname() const { return "InterfacePacking"; }

} // matdes
} // devel
