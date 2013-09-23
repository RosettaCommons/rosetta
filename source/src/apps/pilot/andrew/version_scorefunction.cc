// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/andrew/apl_msd.cc
/// @brief  Application for determining the fingerprint for a score function and for comparing against
///         an existing fingerprint.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Phil Bradley ()
/// @author Brian Weitzner



/// Core headers
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/LREnergyContainer.hh>

// Protocols headers
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>

/// basic headers
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// Utility headers
#include <utility/sort_predicates.hh>
#include <utility/string_util.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>

#include <protocols/jd2/Job.hh>
#include <utility/vector1.hh>



OPT_1GRP_KEY( String, sfxnfprnt, output_fingerprint_file )
//OPT_1GRP_KEY( String, sfxnfprnt, input_fingerprint_file  )

using namespace core;


class ScoreFunctionFingerprintMover : public protocols::moves::Mover 
{
public:
	ScoreFunctionFingerprintMover();
	virtual ~ScoreFunctionFingerprintMover();
	virtual std::string get_name() const;

	void sfxn( scoring::ScoreFunctionOP sfxn );
	virtual void apply( pose::Pose & pose );

	utility::vector1< std::list< std::string > > lines_for_jobs() const;

private:
	scoring::ScoreFunctionOP sfxn_;
	utility::vector1< std::list< std::string > > lines_for_jobs_;

};

////// Implementation ////////

ScoreFunctionFingerprintMover::ScoreFunctionFingerprintMover() {}
ScoreFunctionFingerprintMover::~ScoreFunctionFingerprintMover() {}

std::string
ScoreFunctionFingerprintMover::get_name() const
{
	return "ScoreFunctionFingerprintMover";
}

void
ScoreFunctionFingerprintMover::sfxn( scoring::ScoreFunctionOP sfxn )
{
	using namespace core;
	using namespace core::scoring;
	using namespace core::scoring::methods;

	sfxn_ = sfxn;
	utility::vector1< std::pair< std::string, std::string > > version_term_lines;
	for ( Size ii = 1; ii <= core::scoring::n_score_types; ++ii ) {
		if ( sfxn_->weights()[ ScoreType(ii) ] == 0.0 ) continue; 
		for ( ScoreFunction::AllMethodsIterator iter = sfxn_->all_energies_begin(),
				iter_end = sfxn_->all_energies_end(); iter != iter_end; ++iter ) {
			ScoreTypes const & iter_types = (*iter)->score_types();
			bool found( false );
			for ( Size jj = 1; jj <= iter_types.size(); ++jj ) {
				if ( iter_types[ jj ] == ScoreType( ii ) ) {
					found = true;
					break;
				}
			}
			if ( ! found ) continue;
			std::string name = ScoreTypeManager::name_from_score_type( ScoreType(ii) );
			std::string version = "Version for term " + name + " " + utility::to_string( (*iter)->version() ) + "\n";
			version_term_lines.push_back( std::make_pair( name, version ));
		}
	}

	// sort output by the name of the score types so that the version stays
	// the same even if the order changes in the ScoreType enumeration
	std::sort( version_term_lines.begin(), version_term_lines.end(), utility::SortFirst< std::string, std::string >() );
	std::list< std::string > lines;
	for ( Size ii = 1; ii <= version_term_lines.size(); ++ii ) {
		lines.push_back( version_term_lines[ ii ].second );
	}
	lines_for_jobs_.push_back( lines );
}

void
ScoreFunctionFingerprintMover::apply( pose::Pose & pose )
{
	using namespace core::graph;
	using namespace core::scoring;

	(*sfxn_)( pose );
	std::list< std::string > lines;
	assert( protocols::jd2::JobDistributor::get_instance()->current_job() );
	lines.push_back( "Begin Fingerprint for " + protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag() + "\n" );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Size jj = 1; jj <= scoring::n_score_types; ++jj ) {
			if ( sfxn_->has_nonzero_weight( ScoreType(jj) ) && pose.energies().onebody_energies( ii )[ ScoreType(jj) ] != 0.0 ) {
				std::string newline = "1b " + utility::to_string( ii ) +
					" " + ScoreTypeManager::name_from_score_type( ScoreType(jj) ) +
					" " + ObjexxFCL::format::F( 9, 3, pose.energies().onebody_energies( ii )[ ScoreType(jj) ] ) + "\n";
				lines.push_back( newline );
			}
		} 
	}
	EnergyGraph const &  eg = pose.energies().energy_graph();
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Node::EdgeListConstIter
				iter = eg.get_node(ii)->const_upper_edge_list_begin(),
				iter_end = eg.get_node(ii)->const_upper_edge_list_end();
				iter != iter_end; ++iter ) {
			Size jj( (*iter)->get_second_node_ind() );

			EnergyEdge const * eedge = static_cast< EnergyEdge const * > ( *iter );
			EnergyMap emap = eedge->fill_energy_map();

			for ( Size kk = 1; kk <= scoring::n_score_types; ++kk ) { // this can be more efficient...
				if ( sfxn_->has_nonzero_weight( ScoreType(kk) ) && emap[ ScoreType(kk) ] != 0.0 ) {
					std::string newline = "2b_sr " + utility::to_string( ii ) +
						" " + utility::to_string( jj ) +
						" " + ScoreTypeManager::name_from_score_type( ScoreType(kk) ) +
						" " + ObjexxFCL::format::F( 9, 3,  emap[ ScoreType(kk) ] ) + "\n";
					lines.push_back( newline );
				}
			}
		}
	}

	for ( ScoreFunction::LR_2B_MethodIterator
			lr_iter = sfxn_->long_range_energies_begin(),
			lr_end  = sfxn_->long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {
		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.
		// Potentially O(N^2) operation...

		for ( Size ii = 1; ii <= 1; ++ ii ) {

			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const jj = rni->upper_neighbor_id();
				EnergyMap emap;

				rni->retrieve_energy( emap );
				for ( Size kk = 1; kk <= (*lr_iter)->score_types().size(); ++kk ) {
					ScoreType kkst = (*lr_iter)->score_types()[ kk ];
					if ( emap[ kkst ] != 0 ) {
						std::string newline = "2b_lr " + utility::to_string( ii ) +
							" " + utility::to_string( jj ) +
							" " + ScoreTypeManager::name_from_score_type( kkst ) +
							" " + ObjexxFCL::format::F( 9, 3,  emap[ kkst ] ) + "\n";
						lines.push_back( newline );
					}
				}
			}
		}
	}
				
	for ( Size ii = 1; ii <= scoring::n_score_types; ++ii ) {
		if ( sfxn_->has_nonzero_weight( ScoreType(ii) ) && pose.energies().total_energies()[ ScoreType(ii) ] != 0.0 ) {
			std::string newline = "total " + ScoreTypeManager::name_from_score_type( ScoreType(ii) ) +
				" " + ObjexxFCL::format::F( 9, 3,  pose.energies().total_energies()[ ScoreType(ii) ] ) + "\n";
			lines.push_back( newline );
		}
	}
	lines_for_jobs_.push_back( lines );

}

utility::vector1< std::list< std::string > >
ScoreFunctionFingerprintMover::lines_for_jobs() const
{
	return lines_for_jobs_;
}

typedef utility::pointer::owning_ptr< ScoreFunctionFingerprintMover > ScoreFunctionFingerprintMoverOP;

int main( int argc, char ** argv )
{
	try {

	using namespace utility;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( sfxnfprnt::output_fingerprint_file, "Fingerprint file destination", "" );
	//NEW_OPT( sfxnfprnt::input_fingerprint_file,  "Input fingerprint to verify against", "" );

	devel::init( argc, argv );
	if ( ! option[ sfxnfprnt::output_fingerprint_file ].user() ) {
		utility_exit_with_message( "Must specify either an output fingerprint file" );
	}


	ScoreFunctionFingerprintMoverOP sffm = new ScoreFunctionFingerprintMover;
	sffm->sfxn( core::scoring::getScoreFunction() );

	protocols::jd2::JobDistributor::get_instance()->go(sffm);

	std::ofstream output_file( option[ sfxnfprnt::output_fingerprint_file ]().c_str() ); 

	utility::vector1< std::list< std::string > > lines = sffm->lines_for_jobs();
	for ( Size ii = 1; ii <= lines.size(); ++ii ) {
		for ( std::list< std::string >::const_iterator iter = lines[ii].begin(), iter_end = lines[ii].end();
				iter != iter_end; ++iter ) {
			output_file << *iter;
		}
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}


