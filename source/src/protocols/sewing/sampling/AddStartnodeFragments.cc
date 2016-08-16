// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AddStartnodeFragments.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/sampling/AddStartnodeFragments.hh>
#include <protocols/sewing/sampling/AddStartnodeFragmentsCreator.hh>
#include <protocols/sewing/conformation/DisembodiedAssembly.hh>

//Package headers
#include <protocols/sewing/util/io.hh>
#include <protocols/sewing/util/util.hh>

#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/rosetta_scripts/util.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

namespace protocols {
namespace sewing  {

static basic::Tracer TR( "protocols.sewing.AddStartnodeFragments" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
protocols::moves::MoverOP
AddStartnodeFragmentsCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AddStartnodeFragments );
}

std::string
AddStartnodeFragmentsCreator::keyname() const
{
	return AddStartnodeFragmentsCreator::mover_name();
}

std::string
AddStartnodeFragmentsCreator::mover_name()
{
	return "AddStartnodeFragments";
}

protocols::moves::MoverOP
AddStartnodeFragments::clone() const {
	return( protocols::moves::MoverOP( new AddStartnodeFragments( *this ) ) );
}
protocols::moves::MoverOP
AddStartnodeFragments::fresh_instance() const {
	return protocols::moves::MoverOP( new AddStartnodeFragments );
}

std::string
AddStartnodeFragments::get_name() const {
	return "AddStartnodeFragments";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  AddStartnodeFragments function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


AddStartnodeFragments::AddStartnodeFragments():
	Mover()
{}


void
AddStartnodeFragments::apply(
	core::pose::Pose & pose
) {
	using namespace core;
	using namespace core::fragment;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !option[ in::file::native ].user() || !option[in::file::frag9].user() ) {
		utility_exit_with_message("Must specify native and frag9 files for AddStartnodeFragments mover");
	}

	core::pose::PoseOP native_pose( new pose::Pose );
	core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);

	core::Size start_window = protocols::rosetta_scripts::find_nearest_res( pose, *native_pose, start_res_);
	core::Size end_window = protocols::rosetta_scripts::find_nearest_res( pose, *native_pose, end_res_);

	std::string frag9_file = option[ in::file::frag9 ]();
	core::fragment::FragmentIO frag_io(
		//option[ OptionKeys::abinitio::number_9mer_frags ](),
		200,
		option[ OptionKeys::frags::nr_large_copies ](),
		option[ OptionKeys::frags::annotate ]());
	core::fragment::FragSetOP fragset = frag_io.read_data( frag9_file );
	TR << "Stealing native fragments for windows " << start_window << " to " << end_window << std::endl;
	frag_io.write_data(frag9_file+".nonat", *fragset);
	core::fragment::steal_frag_set_from_pose( *native_pose, start_window, end_window, *fragset, core::fragment::FragDataOP( new FragData( SingleResidueFragDataOP( new BBTorsionSRFD ), fragset->max_frag_length() ) ) );
	frag_io.write_data(frag9_file+".nat", *fragset);
}

void
AddStartnodeFragments::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	using namespace basic::options;

	if ( !tag->hasOption("start_res") || !tag->hasOption("end_res") ) {
		utility_exit_with_message("Must specify start_res and end_res");
	}

	start_res_ = tag->getOption<core::Size>("start_res");
	end_res_ = tag->getOption<core::Size>("end_res");

}


} //sewing
} //protocols
