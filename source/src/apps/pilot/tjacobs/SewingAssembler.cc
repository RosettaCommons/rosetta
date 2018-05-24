// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SewAssembler.cc
///
/// @brief Application wrapper for the SEWING protocol. Takes a pre-generated score file and creates backbone assemblies from it
///
/// @author Tim Jacobs

//Package headers
#include <devel/init.hh>
#include <protocols/sewing/sampling/SewGraph.hh>

//Protocol headers
#include <core/types.hh>
#include <core/id/AtomID.hh>

//Utility headers
#include <utility/io/izstream.hh>
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

//C++ headers
#include <map>
#include <set>
#include <iterator>

static basic::Tracer TR("SewAssembler");

namespace SewAssembler {
basic::options::FileOptionKey const model_file( "model_file" );
basic::options::FileOptionKey const score_file( "score_file" );
basic::options::FileOptionKey const starting_pdb( "starting_pdb" );
basic::options::IntegerOptionKey const num_assemblies( "num_assemblies" );
basic::options::IntegerOptionKey const path_size( "path_size" );
}

int
main( int argc, char * argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::sewing;

	option.add( SewAssembler::model_file, "Model file name");
	option.add( SewAssembler::score_file, "Score file to use for assembly");
	option.add( SewAssembler::num_assemblies, "Number of assemblies to generate and print");
	option.add( SewAssembler::path_size, "How many secondary structure segments to add");
	option.add( SewAssembler::starting_pdb, "starting_pdb");

	// initialize core and read options
	devel::init(argc, argv);
	std::string score_filename = option[SewAssembler::score_file];
	std::string model_filename = option[SewAssembler::model_file];
	core::Size num_assemblies = option[SewAssembler::num_assemblies].def(10);
	core::Size path_size = option[SewAssembler::path_size].def(4);

	std::map< int, Model > models = read_model_file(model_filename).second;

	Hasher hash;
	ScoreResults scores = hash.read_scores_from_disk(score_filename);

	SewGraph sew_graph(models);
	sew_graph.generate_from_scores(scores);

	for ( core::Size index=1; index <= num_assemblies; ++index ) {
		utility::vector1<Edge> random_path = sew_graph.random_path(path_size);
		TR << "Found path!" << std::endl;
		for ( core::Size i=1; i<=random_path.size(); ++i ) {
			TR << random_path[i].source_model() << "-->" << random_path[i].target_model() << std::endl;
		}

		sew_graph.print_path(index, random_path);
	}
}
