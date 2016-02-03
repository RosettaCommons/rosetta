// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/SavePoseMover.cc
/// @author Florian Richter (floric@u.washington.edu)

// Unit Headers
#include <protocols/rosetta_scripts/SavePoseMover.hh>
#include <protocols/rosetta_scripts/SavePoseMoverCreator.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.SavePoseMover" );

namespace protocols {
namespace rosetta_scripts {

std::string
SavePoseMoverCreator::keyname() const
{
	return SavePoseMoverCreator::mover_name();
}

protocols::moves::MoverOP
SavePoseMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SavePoseMover );
}

std::string
SavePoseMoverCreator::mover_name()
{
	return "SavePoseMover";
}

SavePoseMover::SavePoseMover() :
	Mover( "SavePoseMover" ),
	reference_pose_(/* NULL */),
	restore_pose_(false)
{
}

SavePoseMover::~SavePoseMover() {}

void
SavePoseMover::apply( core::pose::Pose & pose )
{
	if ( !restore_pose_ ) *reference_pose_ = pose;
	else pose = *reference_pose_;

	//make sure this always counts as success
	this->set_last_move_status( protocols::moves::MS_SUCCESS );
}

void
SavePoseMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = saved_reference_pose(tag,data_map );
	} else throw utility::excn::EXCN_RosettaScriptsOption("Need to specify name under which to save pose.");

	if ( tag->hasOption( "pdb_file" ) ) {
		std::string const template_pdb_fname( tag->getOption< std::string >( "pdb_file" ));
		core::import_pose::pose_from_file( *reference_pose_, template_pdb_fname , core::import_pose::PDB_file);
		TR <<"reading in " << template_pdb_fname << " pdb with " << reference_pose_->total_residue() <<" residues"<<std::endl;
	}


	if ( tag->hasOption("restore_pose") ) {
		restore_pose_ = tag->getOption<bool>("restore_pose",1);
	}
}

std::string
SavePoseMover::get_name() const {
	return "SavePoseMover";
}

protocols::moves::MoverOP
SavePoseMover::clone() const{
	return protocols::moves::MoverOP( new SavePoseMover( *this ) );
}

protocols::moves::MoverOP
SavePoseMover::fresh_instance() const{
	return protocols::moves::MoverOP( new SavePoseMover );
}


} // moves
} // protocols

