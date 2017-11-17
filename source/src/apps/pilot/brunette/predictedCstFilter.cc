// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author TJ Brunette <tjbrunette@gmail.edu>

// libRosetta headers

#include <core/types.hh>
#include <devel/init.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

//options
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>


#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


//tj headers
#include <apps/pilot/brunette/tj_util.hh>

using std::string;
using core::pose::Pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring::constraints;
using namespace core::chemical;

static basic::Tracer tr( "brunette.predictedCstFilter" );


int main( int argc, char * argv [] ) {
	try {

		using core::sequence::read_fasta_file;
		devel::init( argc, argv );
		//step 1: read in fasta,alignment,template,csts
		map< string, Pose > templateData = poses_from_cmd_line(
			option[ in::file::template_pdb ]());
		map< string, SequenceAlignment> alnDataMapped = input_alignmentsMapped(true);
		string query_sequence (
			read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence() );
		core::pose::Pose query_pose;
		core::pose::make_pose_from_sequence(
			query_pose, query_sequence, *(rsd_set_from_cmd_line())
		);
		ConstraintSetOP sigmoid_cstset
			= core::scoring::constraints::ConstraintIO::read_constraints(
			core::scoring::constraints::get_cst_file_option(),
			new core::scoring::constraints::ConstraintSet,query_pose);
		map< string, Pose> partialThreadsMapped = generate_partial_threads(alnDataMapped,templateData,query_sequence,true);
		std::cout << "name" << partialThreadsMapped.begin()->first << std::endl;
		partialThreadsMapped.begin()->second;
		vector1<core::Real> burial =  calculate_surface_exposure(partialThreadsMapped.begin()->second);
		//step 2: reduce sigmoid constraints in following ways
		// option 1: gap, top N csts
		// option 2: gap, all csts
		// option 3: gap, top N csts, surface exposed
		// option 4: gap, all csts, surface exposed
		// option 5: gap, all csts, surface exposed, seperated into cohesive contact collections.
		//step 3: output sigmoid csts.

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


//core::scoring::constraints::ConstraintIO::write_constraints(std::cout,*sigmoid_cstset,query_pose);
