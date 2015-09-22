// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/protein_interface_design/ReportPSSMDifference.cc
/// @brief calculation of the difference in PSSM score between mutated and native pose
/// @author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)


#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>


// Utility Headers
#include <utility/vector1.hh>

// Unit Headers
// C++ headers
#include <map>
#include <fstream>
#include <protocols/protein_interface_design/ReportPSSMDifference.hh>

#include <core/pack/task/PackerTask.hh>

//Auto Headers
#include <core/conformation/Residue.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.InterfaceRecapitulationMover" );

bool
protocols::protein_interface_design::ReportPSSMDifferences::load_pssm_data(
	std::string const & native_filename
)
{
	std::string native_substr = native_filename.substr( 0, native_filename.size() - 4 );
	std::string pssm_file_name = native_substr + ".fasta";
	std::cerr << "Openning PSSM File " << pssm_file_name << " " << std::endl;
	std::ifstream pssm_file( pssm_file_name.c_str() );

	utility::vector1< Real > pssm_prob_dist( core::chemical::num_canonical_aas, 0.0 );
	Size linenum( 0 );
	pssm_data_.clear();
	while ( pssm_file ) {
		++linenum;
		char line_aa;
		pssm_file >> line_aa;
		core::chemical::AA aa( core::chemical::aa_from_oneletter_code( line_aa ));
		Real sum( 0.0 );
		for ( Size ii = 1; ii <= core::chemical::num_canonical_aas; ++ii ) {
			pssm_file >> pssm_prob_dist[ ii ];
			sum += pssm_prob_dist[ ii ];
		}
		if ( std::abs( sum - 1 ) > 0.001 ) {
			TR << "Warning: pssm probability distribution does not sum to 1.0: " << sum << std::endl;
			TR << "Problem on line " << linenum << " of " << pssm_file_name << std::endl;
		}
		pssm_data_.push_back( std::make_pair( aa, pssm_prob_dist ));
	}

	if ( pssm_data_.size() == 0 ) { std::cerr << "Did not read file -- possibly not found" << std::endl; return false; }

	return true;
}


core::Real
protocols::protein_interface_design::ReportPSSMDifferences::calculate(
	core::pose::Pose const & pose1_in, core::pose::Pose const & pose2_in, core::pack::task::PackerTaskCOP const & task
)
{
	using namespace core::scoring;

	core::pose::Pose pose1( pose1_in );
	core::pose::Pose pose2( pose2_in );
	core::Real pssm = 0.;

	for ( core::Size i = 1; i <= pose1.total_residue(); ++i ) {
		if ( !pose1.residue(i).is_protein() ) continue;
		core::chemical::AA const restype( pose2.residue(i).aa() );

		if ( task->being_designed( i ) ) {
			if ( pssm_data_[i].first == restype ) {
				pssm += pssm_data_[i].second[ restype ];
			} else {
				TR << "Warning: No pssm data found. Falling back on Sequence comparison." << std::endl;
				if ( pose1.residue(i).aa() ==  pose2.residue(i).aa() ) {
					pssm += 1.;
				}
			}
		}
	}

	return pssm;
}
