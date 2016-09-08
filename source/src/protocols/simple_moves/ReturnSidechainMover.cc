// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ReturnSidechainMover.cc
/// @brief ReturnSidechainMover methods implemented
/// @author Steven Lewis


// Unit Headers
#include <protocols/simple_moves/ReturnSidechainMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/util.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string> //making sure that residue_type.name() comparison works


using basic::T;
using basic::Error;
using basic::Warning;


static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.ReturnSidechainMover" );


namespace protocols {
namespace simple_moves {

/// @details this code was copied from protocols/loops/loops_main.cc:187-204, revision 21282
void
ReturnSidechainMover::apply( core::pose::Pose & pose )
{
	core::pose::Pose saved_input_pose(remembered_pose_);
	core::Size nres( end_res_ - start_res_ + 1 );

	//symmetry
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(remembered_pose_) ) {
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	if ( nres != saved_input_pose.size() ) {
		utility_exit_with_message("ReturnSidechainMover used with poses of different length; aborting");
	}

	for ( Size i=start_res_, j=1; i<= end_res_; ++i, ++j ) {
		bool copy_this_residue( false );

		if ( copy_all_chi_ ) {
			copy_this_residue = true;
		} else {
			if ( allow_chi_copy_[i] ) {
				copy_this_residue = true;
			}
		}

		if ( symm_info && !symm_info->bb_is_independent( i ) ) continue;
		if ( pose.residue_type(i).aa() == core::chemical::aa_vrt ) continue;

		if ( copy_this_residue ) {
			core::chemical::ResidueType const & rsd_type ( pose.residue(i).type() );
			core::chemical::ResidueType const & saved_rsd_type( saved_input_pose.residue(j).type() );

			//ensure that there is no sequence change
			if ( rsd_type.name3() != saved_rsd_type.name3() ) {
				utility_exit_with_message("ReturnSidechainMover used with poses of different sequence; aborting");
			}

			//we need to check variant types in case there are cutpoints for loop modeling or whatever
			if ( ! variants_match( rsd_type, saved_rsd_type ) ) {
				utility::vector1< std::string > const & variant_types ( rsd_type.properties().get_list_of_variants() );
				utility::vector1< std::string > missing_variant_types;

				 for ( auto const & variant_type : variant_types ) {
					if ( !saved_rsd_type.has_variant_type( variant_type ) ) {
						missing_variant_types.push_back( variant_type );
					}
				}

				for ( utility::vector1< std::string >::const_iterator it = missing_variant_types.begin(),
						it_end=missing_variant_types.end(); it != it_end; ++it ) {
					core::pose::add_variant_type_to_pose_residue( saved_input_pose,
						core::chemical::ResidueProperties::get_variant_from_string( *it ), i);
				}
			} //checking variants

			pose.replace_residue(i, saved_input_pose.residue(j), true );
		}
	} // residue i

	pose.energies().clear();
	return;
}//apply

std::string
ReturnSidechainMover::get_name() const {
	return "ReturnSidechainMover";
}

// default constructor
ReturnSidechainMover::ReturnSidechainMover() : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "ReturnSidechainMover" );
	copy_all_chi_ = true;
}

// constructor with pose
ReturnSidechainMover::ReturnSidechainMover(
	core::pose::Pose const & pose_in,
	core::Size start_res,
	core::Size end_res ) :
	protocols::moves::Mover(), remembered_pose_(pose_in)
{
	protocols::moves::Mover::type( "ReturnSidechainMover" );
	copy_all_chi_ = true;
	if ( start_res == 0 ) start_res_ = 1;
	else start_res_ = start_res;
	if ( end_res == 0 ) end_res_ = pose_in.size();
	else end_res_ = end_res;
}

// constructor with pose
ReturnSidechainMover::ReturnSidechainMover(
	core::pose::Pose const & pose_in,
	utility::vector1<bool> allow_chi_in,
	core::Size start_res,
	core::Size end_res ) :
	protocols::moves::Mover(),
	allow_chi_copy_( allow_chi_in ),
	remembered_pose_(pose_in)
{
	protocols::moves::Mover::type( "ReturnSidechainMover" );
	copy_all_chi_ = false;
	if ( start_res == 0 ) start_res_ = 1;
	else start_res_ = start_res;
	if ( end_res == 0 ) end_res_ = pose_in.size();
	else end_res_ = end_res;
}

// copy constructor
ReturnSidechainMover::ReturnSidechainMover(ReturnSidechainMover const & ) = default;

ReturnSidechainMover::~ReturnSidechainMover() = default;

void
ReturnSidechainMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "                    start\tend" << std::endl;
	output << "Range of residues:  " << get_start_res() << "\t\t" << get_end_res() << std::endl;
}

protocols::moves::MoverOP
ReturnSidechainMover::clone() const
{
	return protocols::moves::MoverOP( new ReturnSidechainMover(*this) );
}

protocols::moves::MoverOP
ReturnSidechainMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ReturnSidechainMover() );
}

core::Size ReturnSidechainMover::get_start_res() const {
	return start_res_;
}

core::Size ReturnSidechainMover::get_end_res() const {
	return end_res_;
}

std::ostream &operator<< (std::ostream &os, ReturnSidechainMover const &mover)
{
	mover.show(os);
	return os;
}

}//moves
}//protocols
