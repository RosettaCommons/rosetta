// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh>
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechainsCreator.hh>

// Project headers
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
//for putting back right variants

//Auto Headers
#include <core/pose/util.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.SaveAndRetrieveSidechains" );

std::string
SaveAndRetrieveSidechainsCreator::keyname() const
{
	return SaveAndRetrieveSidechainsCreator::mover_name();
}

protocols::moves::MoverOP
SaveAndRetrieveSidechainsCreator::create_mover() const {
	return new SaveAndRetrieveSidechains;
}

std::string
SaveAndRetrieveSidechainsCreator::mover_name()
{
	return "SaveAndRetrieveSidechains";
}

SaveAndRetrieveSidechains::SaveAndRetrieveSidechains() :
	simple_moves::DesignRepackMover( SaveAndRetrieveSidechainsCreator::mover_name() )
{
	allsc_ = false; // default
	jumpid_ = 1; //default
	ensure_variant_matching_ = false; //default
	two_step_ = false;
	first_apply_ = new protocols::moves::DataMapObj< bool >;
	first_apply_->obj = true;
	init_pose_ = new core::pose::Pose;
}

SaveAndRetrieveSidechains::SaveAndRetrieveSidechains(
	core::pose::Pose const & pose,
	bool const allsc /*=false*/,
	bool const ensure_variant_matching /*=false*/,
	core::Size const jumpid /*=1*/
) :
	simple_moves::DesignRepackMover( SaveAndRetrieveSidechainsCreator::mover_name() ),
	allsc_( allsc ),
	ensure_variant_matching_(ensure_variant_matching),
	jumpid_( jumpid )
{
	init_pose_ = new core::pose::Pose( pose );
	two_step_ = false;
	first_apply_ = new protocols::moves::DataMapObj< bool >;
	first_apply_->obj = true;
}

SaveAndRetrieveSidechains::~SaveAndRetrieveSidechains() {}

void
SaveAndRetrieveSidechains::apply( Pose & pose )
{
	typedef conformation::Residue Residue;
	if( two_step() && first_apply_->obj ){
		TR<<"Saving sidechains."<<std::endl;
		*init_pose_ = pose;
		first_apply_->obj = false;
		return;
	}
	TR << "Retrieving sidechains..."<<std::endl;
	Size nres = pose.total_residue();
	if (nres != init_pose_->total_residue() && core::pose::symmetry::is_symmetric(pose)) {
		conformation::symmetry::SymmetricConformation & symm_conf (
				dynamic_cast<conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		nres = symm_conf.Symmetry_Info()->num_independent_residues();
	}
	runtime_assert( nres == init_pose_->total_residue() );
	kinematics::Jump new_jump;
	core::Size const rb_jump( jumpid_ );
	if ( jumpid_ > 0 ) {
	new_jump = pose.jump( rb_jump );
	}

	for( core::Size res=1; res<=nres; ++res ) {
		if( allsc_ ) { // replace all sidechains
			pose.replace_residue( res, init_pose_->residue( res ), true/*orient_backbone*/ );
			continue;
		}
		else {
			if( pose.residue( res ).name3() == "ALA" ) // only replace Ala positions
			pose.replace_residue( res, init_pose_->residue( res ), true/*orient_backbone*/ );
		}
	}
	 if (ensure_variant_matching_){
		//make sure variants match, if not put back the initial variants
    using namespace core;
    for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {

      if( !(pose.residue_type( i ).variants_match( init_pose_->residue_type( i ) ) ) ){

        utility::vector1< std::string > const new_var_types( pose.residue_type( i ).variant_types() );
        utility::vector1< std::string > const old_var_types( init_pose_->residue_type( i ).variant_types() );
        for( utility::vector1< std::string >::const_iterator newvars = new_var_types.begin(); newvars  != new_var_types.end(); ++newvars ){
          if( ! (init_pose_->residue_type( i ).has_variant_type( *newvars ) ) ) core::pose::remove_variant_type_from_pose_residue( pose, *newvars, i );
        }

        for( utility::vector1< std::string >::const_iterator oldvars = old_var_types.begin(); oldvars  != old_var_types.end(); ++oldvars ){
          if( !pose.residue_type( i ).has_variant_type( *oldvars ) ) core::pose::add_variant_type_to_pose_residue( pose, *oldvars, i );
        }
      } //if variants don't match
    }
	}
	if ( jumpid_ > 0 ) {
		pose.set_jump( rb_jump, new_jump );
	}
	TR.flush();
}

std::string
SaveAndRetrieveSidechains::get_name() const {
	return SaveAndRetrieveSidechainsCreator::mover_name();
}

void
SaveAndRetrieveSidechains::parse_my_tag( TagPtr const tag, DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	first_apply_->obj = true;
	allsc_ = tag->getOption<bool>( "allsc", 0 );
	two_step( tag->getOption< bool >( "two_step", false ) );
	if( !two_step() )
		init_pose_ = new core::pose::Pose( pose );
}

protocols::moves::MoverOP
SaveAndRetrieveSidechains::clone() const {
  return( protocols::moves::MoverOP( new SaveAndRetrieveSidechains( *this )));
}

} //movers
} //protein_interface_design
} //protocols
