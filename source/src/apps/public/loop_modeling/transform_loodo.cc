// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilot/kblacklock/transform_loodo.cc
///
/// @brief Pilot Application for reading CapHit information from LooDo CapHit_RT.txt-style output file
/// and generating CapHits in PDB form.
///
/// @author Kristin Blacklock (kristin.blacklock@rutgers.edu)
/// @note Last Modified: 13 October 2022 by Vikram K. Mulligan. (Promoted to public app.)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <basic/options/keys/loodo.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/Stub.hh>

#include <protocols/sic_dock/util.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using Vector = numeric::xyzVector<core::Real>;
using Matrix = numeric::xyzMatrix<core::Real>;

basic::Tracer TR("apps.public.loop_modeling.transform_loodo");

core::kinematics::Stub
parse_into_stub( utility::vector1<core::Real> const & in_vec )
{
	if ( in_vec.size() < 12 ) {
		utility_exit_with_message("Need 12 values to parse into stub.");
	}
	core::kinematics::Stub to_apply_stub;
	Matrix M;
	Vector v;
	int count (0), i(1), j(1), k(0);
	for ( core::Real val: in_vec ) {
		count++;
		if ( count <= 9 ) {
			if ( j==4 ) {
				i++;
				j=1;
			}
			M(i,j) = val;
			j++;
		} else {
			v[k] = val;
			k++;
		}
	}
	to_apply_stub.M = M;
	to_apply_stub.v = v;
	return to_apply_stub;
}

/// @brief Main method
int main( int argc, char * argv [] ) {

	using namespace basic::options;

	try {

		// Devel init factories
		devel::init(argc, argv);

		////////// Register with Citation Manager //////////

		{
			basic::citation_manager::CitationManager * cm ( basic::citation_manager::CitationManager::get_instance() );
			basic::citation_manager::CitationCollectionOP collection( utility::pointer::make_shared< basic::citation_manager::CitationCollection >( "Transform_LooDo", basic::citation_manager::CitedModuleType::Application ) );
			collection->add_citation( cm->get_citation_by_doi("10.1002/prot.25445") );
			cm->add_citation( collection );
		}

		// Import native insert domain PDB & center.
		std::string pdbfile = option[ OptionKeys::loodo::cap ]();
		core::pose::Pose native_pose;
		core::import_pose::pose_from_file( native_pose, pdbfile , core::import_pose::PDB_file);
		native_pose.center();

		core::pose::Pose moving_pose = native_pose;

		// Get CapHit_RT (or equivalent) file name.
		std::string stubfile = option[ OptionKeys::loodo::caphit_rt_file ]();

		// Open CapHit_RT (or equivalent) File
		utility::io::izstream stream(stubfile);
		if ( !stream ) {
			utility_exit_with_message("File not found");
		}

		// Read CapHit_RT (or equivalent) file lines and separate into CapHit name & stub
		std::string line;
		while ( getline(stream,line) ) {
			//std::cout << line << std::endl;

			core::Size name_index = line.find_first_of(" || ", 0);
			core::Size stub_index = line.rfind("STUB ");

			if ( stub_index == std::string::npos ) {
				utility_exit_with_message("While parsing " + stubfile + " STUB designation not found in line.");
			}
			if ( name_index == std::string::npos ) {
				utility_exit_with_message("While parsing " + stubfile + " name designation not found.");
			}

			std::string name = line.substr(0, name_index);
			//TR << "Name: " << name << std::endl;

			std::string stub = line.substr(stub_index+5);
			//TR << "Stub: " << stub << std::endl;

			std::vector< std::string > stubvec;
			boost::split( stubvec, stub, boost::is_any_of(" "));

			utility::vector1<core::Real> SV;
			for ( auto & i : stubvec ) {
				auto f = boost::lexical_cast<float>(i);
				SV.push_back(f);
			}

			if ( SV.size() < 12 ) {
				utility_exit_with_message("While parsing " + stubfile + " STUB line does not contain sufficient entries.");
			}

			// Parse file information into stub and apply to the centered native pose.
			moving_pose = native_pose;
			core::kinematics::Stub to_apply_stub = parse_into_stub( SV );
			//TR << "Stub after parse_into_stub: " << to_apply_stub << std::endl;
			protocols::sic_dock::xform_pose_rev( moving_pose, to_apply_stub );
			moving_pose.dump_pdb("XFORM_"+name);

		}
		stream.close();

		// Final citation manager output:
		basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info();

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
