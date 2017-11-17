// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  Application that reads in dimer (composed of two chains) and outputs the peptide which contibutes most to the interface.

/// @author Nir London
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @date   Nov. 15, 2009

//#define GL_GRAPHICS

// Unit Headers
#include <protocols/moves/FilterReporterMover.hh>
#include <protocols/peptide_deriver/PeptideDeriverFilter.hh>

// Project headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/peptide_deriver.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/NullMover.hh>

// C++ headers
using basic::Error;
using basic::Warning;

// using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR("apps.public.peptide_deriver.PeptideDeriver");


// PeptideDeriver app
//
int
main( int argc, char * argv[] ) {
	try {

		protocols::jd2::register_options();

		option.add_relevant( out::output );
		option.add_relevant( out::nooutput );
		option.add_relevant( out::mute );
		option.add_relevant( out::unmute );
		option.add_relevant( out::user_tag );
		option.add_relevant( out::output_tag );
		option.add_relevant( out::chname );
		option.add_relevant( out::chtimestamp );
		option.add_relevant( out::mpi_tracer_to_file );

		option.add_relevant( in::file::l );
		option.add_relevant( in::file::list );
		option.add_relevant( in::file::s );
		option.add_relevant( in::path::path );

		option.add_relevant( peptide_deriver::pep_lengths );
		option.add_relevant( peptide_deriver::skip_zero_isc );
		option.add_relevant( peptide_deriver::dump_peptide_pose );
		option.add_relevant( peptide_deriver::dump_prepared_pose );
		option.add_relevant( peptide_deriver::dump_cyclic_poses );
		option.add_relevant( peptide_deriver::dump_report_file );
		option.add_relevant( peptide_deriver::restrict_receptors_to_chains );
		option.add_relevant( peptide_deriver::restrict_partners_to_chains );
		option.add_relevant( peptide_deriver::do_minimize );
		option.add_relevant( peptide_deriver::optimize_cyclic_threshold );


		//setup random numbers and options
		devel::init(argc, argv);

		// by default, enable jump generation when missing density
		// this prevents Peptiderive from deriving discontinuous peptides
		if ( !option[ in::missing_density_to_jump ].user() ) {
			option[ in::missing_density_to_jump ].value( true );
		}

		// disable PDB output, method 2
		// SilentFileJobOutputterOP job_out = new SilentFileJobOutputter;
		// job_out->set_write_no_structures();

		// create a pose
		// pose::Pose orig_pose;
		// core::import_pose::pose_from_file(orig_pose, start_file(), core::import_pose::PDB_file);
		protocols::peptide_deriver::PeptideDeriverFilterOP filter =
			protocols::peptide_deriver::PeptideDeriverFilterOP(
			new protocols::peptide_deriver::PeptideDeriverFilter );

		protocols::moves::MoverOP null_mover = protocols::moves::MoverOP(
			new protocols::moves::NullMover );

		protocols::moves::MoverOP mover = protocols::moves::MoverOP(
			new protocols::moves::FilterReporterMover(null_mover, filter,
			1 /*max_tries*/, TR /*out*/) );

		// GO!
		protocols::jd2::JobDistributor::get_instance()->go(mover);

	} catch (utility::excn::Exception const & e) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
