// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @file   apps/pilot/nobuyasu/fldsgn.cc
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// project headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>

#include <protocols/flxbb/BluePrint.hh>
#include <protocols/fldsgn/BluePrintBDR.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <utility/excn/Exceptions.hh>

// utility headers
#include <utility/file/FileName.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/vector1.hh>

// Key headers
#include <core/options/option.hh>

// C++ headers
#include <string>

// My protocols
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
//#include <basic/options/keys/fldsgn.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( FileVector, blueprint )
OPT_KEY( Boolean, native_bias )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( blueprint, "input blueprint files", "" );
	NEW_OPT( native_bias, "use native bias potentials", false );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief
void * graphics_main( void * ) {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	using namespace basic::options::OptionKeys;
	using basic::options::option;
	using core::Size;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using utility::file::FileName;
 	using protocols::flxbb::BluePrint;
	using protocols::flxbb::BluePrintOP;
 	using protocols::fldsgn::BluePrintBDR;
	using protocols::fldsgn::BluePrintBDR_OP;
	using protocols::moves::SequenceMover;
	using protocols::moves::SequenceMoverOP;
	using protocols::jd2::JobDistributor;

	// read weights
	std::string sfx = option[ score::weights ].user() ? option[ score::weights ].value() : "fldsgn_cen";

	// set SequenceMover of FilterMovers
	SequenceMoverOP seq_mover = new SequenceMover;
	seq_mover->use_mover_status( true );

	utility::vector1< FileName > blueprint_files;
	blueprint_files = option[ blueprint ]().vector();
	for ( utility::vector1< FileName >::iterator it=blueprint_files.begin(), ite=blueprint_files.end(); it != ite; ++it ) {

		// read blueprint
		BluePrintOP blp = new BluePrint( it->name() );

		// setting score function
		ScoreFunctionOP csfxn = ScoreFunctionFactory::create_score_function( sfx );

		// init BDR
		BluePrintBDR_OP bdr = new BluePrintBDR;
		bdr->scorefunction( csfxn );
		bdr->set_blueprint( blp );
		bdr->ss_from_blueprint( true );
		seq_mover->add_mover( bdr );
	} // blueprint_files


#if defined GL_GRAPHICS

	core::pose::Pose pose;
	core::import_pose::pose_from_file( pose, option[ in::file::s ].value().at( 1 ) , core::import_pose::PDB_file);
	protocols::viewer::add_conformation_viewer( pose.conformation(), "fldsgn" );
	return 0;

#else

	JobDistributor::get_instance()->go( seq_mover );

#endif

	return 0;
}


int main( int argc, char * argv [] ) {
	try{
	using namespace basic::options::OptionKeys;
	using basic::options::option;

	// initialize rosetta
	devel::init( argc, argv );

	protocols::viewer::viewer_main( graphics_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


