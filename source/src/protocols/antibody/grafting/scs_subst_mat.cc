// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/scs_subs_mat.cc
/// @brief Structural Component Selector (SCS) implementation using substition matrices
/// @author Brian D. Weitzner

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/scs_subst_mat.hh>

// #include <protocols/antibody/grafting/chothia_numberer.hh>
#include <protocols/antibody/grafting/scs_functor.hh>

#include <core/sequence/MatrixScoringScheme.hh>

#include <basic/Tracer.hh>
#include <basic/execute.hh>

#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/string_util.hh>
//#include <utility/stream_util.hh>

#include <algorithm>
#include <fstream>
#include <locale>

namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting");



void SCS_SubstitutionMatrix::init_from_options()
{
	SCS_LoopOverSCs::init_from_options();

	// using namespace basic::options;
	// using namespace basic::options::OptionKeys;

	// Probably a good idea to enable different substitution matrices to be set explicitly from an option
}


SCS_ResultsVector populate_results_vector( std::map< string, std::map<string, core::Real> > const & result_map,
                                           std::map< string, std::map<string, string> > const & db )
{
	SCS_ResultsVector results;

	for(auto result : result_map ) {
		// Converting text based results filed to SCS_BlastMetric
		SCS_SubstitutionMatrixResultOP r(new SCS_SubstitutionMatrixResult);
		r->pdb = result.first; // this comes from the antibody.info file; PDB codes are written in the form "2adf"

		r->sid  = result.second.at( "sid" );
		r->score = result.second.at( "score" );

		populate_results_from_db( r, db );

		results.push_back(r);
	}

	return results;
}


void SCS_SubstitutionMatrix::select_template(
	Result & j,
	string const & /*db_to_query*/,
	std::map< string, std::map<string, string> > const & ab_db ) const
{

	// Query sequence is in j.sequence
	string subst_mat_name = (j.name.find("fr") < string::npos or j.name.find("orientation") < string::npos or \
									j.name.find("heavy") < string::npos or j.name.find("light") < string::npos) ? "BLOSUM62" : "PAM30";


	core::sequence::MatrixScoringScheme substitution_matrix_reader;
	substitution_matrix_reader.read_from_database( subst_mat_name ); // probably makes more sense to lazily load and keep these around

	typedef utility::vector1< utility::vector1< core::Real > > SubstitionMatrix;
	SubstitionMatrix subst_mat = substitution_matrix_reader.scoring_matrix();

	std::map< string, std::map<string, core::Real> > search_results;

	core::Size sc_length = j.sequence.length();

	for ( auto & db_entry : ab_db ) {

		string sc_seq = db_entry.second.at( j.name );

		if ( sc_length != sc_seq.length() ) { continue; }

		string::const_iterator res1 = j.sequence.begin();
		string::const_iterator res2 = sc_seq.begin();

		core::Real seq_score = 0.;
		core::Size mismatches = 0;

		// iterate over sequences
		for (; res1 != j.sequence.end(); ++res1, ++res2) {
			core::chemical::AA
			aa1( core::chemical::aa_from_oneletter_code( * res1 ) ),
			aa2( core::chemical::aa_from_oneletter_code( * res2 ) );

			if ( aa1 != aa2 ) { ++mismatches; }
			// check for non-canonicals in sequence
			if ( aa1 != core::chemical::aa_unk &&  aa2 != core::chemical::aa_unk ) { seq_score += subst_mat[ aa1 ][ aa2 ]; }
		}

		core::Real sid = ( core::Real )( sc_length - mismatches ) / ( core::Real ) sc_length;

		// populate the result map
		search_results[ db_entry.first ] = {{"sid", sid}, {"score", seq_score}};
	}

	j.results = populate_results_vector( search_results, ab_db );

}


} // namespace grafting
} // namespace antibody
} // namespace protocols


#endif // __ANTIBODY_GRAFTING__
