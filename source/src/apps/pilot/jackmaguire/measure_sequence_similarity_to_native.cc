// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jackmaguire/measure_sequence_similarity_to_native.cc
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>

#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/memory.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/Jobs.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/InterGroupInterfaceByVectorSelector.hh>
#include <core/select/residue_selector/LayerSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>

#include <core/import_pose/import_pose.hh>
#include <core/simple_metrics/metrics/SequenceSimilarityMetric.hh>

OPT_1GRP_KEY( String, select, region )
/*
region can be:
all
core
boundary
surface
interface12
interface1
interface2
*/

static basic::Tracer TR( "apps.pilot.jackmaguire.MeasureSequenceSimilarityToNative" );
static basic::Tracer TR_results( "apps.pilot.jackmaguire.MeasureSequenceSimilarityToNative.Results" );

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::select::residue_selector;

ResidueSelectorOP region_for_string( std::string const & region );

int main( int argc, char* argv[] ){

	try {

		NEW_OPT( select::region, "Select which region of the native to use. Options are 'all', 'core', 'boundary', 'surface', 'interface12' (interface between chains 1 and 2), 'interface1' (chain 1 of the interface), and 'interface2' (chain 2 of the interface).", "" );

		devel::init( argc, argv );

		if ( ! option[ select::region ].user() ) {
			TR << "Select which region of the native to use by passing -select:region [region name]. Options are 'all', 'core', 'boundary', 'surface', 'interface', 'interface1' (chain 1 of the interface), and 'interface2' (chain 2 of the interface).";
			return 1;
		}

		if ( ! option[ OptionKeys::in::file::native ].user() ) {
			utility_exit_with_message( "Please pass a reference pose by using the -in:file:native flag." );
		}

		std::string const region = option[ select::region ]();
		ResidueSelectorOP selector = region_for_string( region );

		core::pose::PoseOP native = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ]() );
		core::simple_metrics::metrics::SequenceSimilarityMetric metric( selector );
		metric.set_native_pose( native, true );

		utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();
		for ( auto const & job : input_jobs ) {
			std::string const filename = job->input_tag();
			core::pose::PoseOP pose = core::import_pose::pose_from_file( filename );
			TR_results << filename << " " << metric.calculate( * pose ) << std::endl;
		}

		return 0;

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}



ResidueSelectorOP region_for_string( std::string const & region ){
	using namespace utility::pointer;
	if ( region == "all" ) return make_shared< TrueResidueSelector >();

	//Layers
	if ( region == "core" ) {
		LayerSelectorOP core = make_shared< LayerSelector >();
		core->set_layers( true, false, false );
		return core;
	}
	if ( region == "boundary" ) {
		LayerSelectorOP boundary = make_shared< LayerSelector >();
		boundary->set_layers( false, true, false );
		return boundary;
	}
	if ( region == "surface" ) {
		LayerSelectorOP surface = make_shared< LayerSelector >();
		surface->set_layers( false, false, true );
		return surface;
	}

	//Interfaces
	ChainSelectorOP chain1 = make_shared< ChainSelector >( "1" );
	ChainSelectorOP chain2 = make_shared< ChainSelector >( "2" );
	ResidueSelectorOP interface = make_shared< InterGroupInterfaceByVectorSelector >( chain1, chain2 );
	if ( region == "interface12" ) return interface;
	if ( region == "interface1" ) return make_shared< AndResidueSelector >( chain1, interface );
	if ( region == "interface2" ) return make_shared< AndResidueSelector >( chain2, interface );

	utility_exit_with_message( "No region:select match for " + region );
	return 0;
}
