// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <protocols/mpi_refinement/Scheduler.hh>
#include <protocols/mpi_refinement/util.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <fstream> // for ifstream
#include <algorithm> // for sort
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace mpi_refinement {

static THREAD_LOCAL basic::Tracer TR("MPI.LHR.S");

Scheduler::Scheduler()
{
	set_default();
	clear(); // just for initializing
}

Scheduler::~Scheduler()= default;

void
Scheduler::set_default(){
	is_random_ = false;
	niter_ = 1;
	iter_ = 0;
	n_rerelaxed_ = 0;
	n_to_rerelax_ = 0;

	nmethods_enrich_max_ = 3;
	enrich_Zscore_cut_ = -1.0;
	Zdiff_outstand_ = 1.5;
	exclusions_.resize( 0 );
}

void
Scheduler::prepare_search_stage( core::Size const mpi_rank )
{
	using namespace basic::options;
	TR << "prepare for search starge: " << std::endl;

	// Currently only via cmdfile
	if ( option[ OptionKeys::lh::mpi_master_schfile ].user() ) {
		std::string cmdfile = option[ OptionKeys::lh::mpi_master_schfile ]();
		read_cmd( cmdfile, mpi_rank );
	}
}

// Add new schedules if
void
Scheduler::prepare_enrich_stage( protocols::wum::SilentStructStore const &decoys,
	std::string const & scorename )
{
	//std::string roundtype( "enrich" );

	TR << "==================================================================" << std::endl;
	TR << "Pick enrich methods" << std::endl;
	TR << "==================================================================" << std::endl;

	// 1. Pick
	methods_picked_.resize( 0 );
	methods_picked_ = pick_enrich_methods( decoys, scorename );

	TR << "Picked " << methods_picked_.size() << " methods: ";
	for ( core::Size i = 1; i <= methods_picked_.size(); ++i ) TR << " " << methods_picked_[i];
	TR << std::endl;

	// 2. Add enriching methods
	utility::vector1< MethodParams > params_copy( params_ );

	// fill in prv
	params_.resize( 0 );
	for ( core::Size i = 1; i < isch_; ++i ) params_.push_back( params_copy[i] );

	for ( core::Size i = 1; i <= methods_picked_.size(); ++i ) {
		core::Size imethod( methods_picked_[i] );

		bool found( false );
		for ( core::Size j = 1; j <= params_copy.size(); ++j ) {
			if ( params_copy[j].index == imethod ) {
				found = true;
				params_.push_back( params_copy[j] );
				TR << "Adding " << imethod << " == " << params_copy[j].name << " for enriching." << std::endl;
				break;
			}
		}

		if ( !found ) TR << "Warning! " << imethod << " not found!" << std::endl;
	}

	// fill in rest
	for ( core::Size i = isch_+1; i <= params_copy.size(); ++i ) params_.push_back( params_copy[i] );

	// Debug
	// Simulate
	TR << "curr isch_: " << isch_ << std::endl;
	for ( core::Size i = 1; i <= params_copy.size(); ++i ) {
		if ( i == isch_ ) {
			TR << "--------- insertion begin (run from here ) ---------" << std::endl;
		}

		TR << i << " " << params_[i].name << " " << std::endl;
		if ( i == isch_+methods_picked_.size() ) {
			TR << "------------------ insertion end ------------------" << std::endl;
		}

	}

	TR << "==================================================================" << std::endl;
	TR << "End enrich setup" << std::endl;
	TR << "==================================================================" << std::endl;

}

// Clear schedules - I think this should be private...
void
Scheduler::clear(){
	isch_ = 1;
	params_.resize( 0 );
	methods_picked_.resize( 0 );
	methodname_.clear();
}

// Proceed to next round
void
Scheduler::proceed(){
	//TR.Debug << "call proceed" << std::endl;

	isch_++;

	// this may never happen, but make it safer
	if ( isch_ > params_.size() ) {
		isch_ = params_.size();
		return;
	}

	// End of schedule; terminate if no more iteration left
	if ( params_[ isch_ ].roundtype.compare( "done" ) == 0 ) {
		// Check if using iteration; if then revert
		iter_++;
		// also renew rerelax stuffs
		n_to_rerelax_ = 0;
		n_rerelaxed_ = 0;

		if ( iter_ < niter_ ) {
			isch_ = 1;
			TR << "Iter/Niter=" << iter_ << " " << niter_;
			TR << ", Reverting schedule to the initial schedule." << std::endl;
		}
	}
}

void
Scheduler::read_cmd( std::string const & cmdfile,
	core::Size const mpi_rank,
	core::Size const stage_to_run )
{
	// clean before reading
	clear();

	std::ifstream infile( cmdfile.c_str() );
	TR.Debug << "================== Reading script file: ==================" << std::endl;
	if ( !infile.good() ) {
		utility_exit_with_message( "[ERROR] Error opening script file '" + cmdfile + "'" );
	}

	std::string line;
	core::Size nline( 0 );
	std::string roundtype;
	core::Size rank_specific( 999 );
	core::Size nmaster_max( 1 );
	core::Size imethod( 0 );
	//bool run_continuous( false );
	core::Size stage( 1 ); //first round if not specified

	while ( getline(infile,line) ) {
		utility::vector1< std::string > tokens ( utility::split( line ) );
		MethodParams params;
		nline++;

		if ( tokens.size() == 0 ) continue;

		// Skip comments
		if ( tokens[1].compare( 0, 1, "#" ) == 0 ) continue;

		if ( tokens[1] == "stage" ) stage = atoi( tokens[2].c_str() );
		if ( stage != stage_to_run ) continue;

		// Special headers: rank, round
		if ( tokens[1] == "rank" ) {
			rank_specific = atoi( tokens[2].c_str() );
			if ( rank_specific == mpi_rank ) {
				clear();
			}
			if ( rank_specific > nmaster_max ) nmaster_max = rank_specific;
			continue;

		} else if ( tokens[1] == "end_rank" ) {
			rank_specific = 999;
			continue;
		}

		// cas when specified but does not correspond to
		if ( rank_specific != 999 && mpi_rank != rank_specific ) continue;

		// Round control headers
		if ( tokens[1] == "exclude" ) {
			if ( tokens.size() != 3 ) {
				TR << line << ": exclude only works on a pair of names. skip." << std::endl;
			}
			exclusions_.push_back( std::make_pair( tokens[2], tokens[3] ) );
			continue;

		} else if ( tokens[1] == "round" ) {
			roundtype = tokens[2];
			// Just add on a fresh schedule only containing name
			if ( roundtype.compare( "enrich" ) == 0 || roundtype.compare( "wait" ) == 0||
					roundtype.compare( "average") == 0 || roundtype.compare( "calcdev" ) == 0 ||
					roundtype.compare( "nextstage" ) == 0 || roundtype.compare( "nextgen" ) == 0 ||
					roundtype.compare( "cluster" ) == 0 ) {
				add_fresh_param( roundtype );

			} else {
				TR.Warning << "unknown roundtype " << roundtype << " pass." << std::endl;
			}
			continue;

		} else if ( tokens[1] == "done" ) {
			add_fresh_param( tokens[1] );
			continue;

		} else if ( tokens[1] == "iter" ) {
			niter_ = atoi( tokens[2].c_str() );
			continue;

		} else if ( tokens[1] == "pick" ) {
			npick_per_iter_ = atoi( tokens[2].c_str() );
			if ( tokens.size() >= 3 ) pick_strategy_ = tokens[2]; // default is random
			if ( tokens.size() >= 4 ) pick_objfunction_ = tokens[3]; // default is esum
			continue;

		}

		// Will pass here only for "search"
		if ( tokens.size() < 7 ) {
			TR << "Line " << nline << ": ";
			TR << "tokens too few: [name/mover/nrun/scoretype/cstw/relaxtype/rerelaxtype](+/shave1/shave2); ";
			TR << "Skip line." << std::endl;
			continue;
		}

		// Common params:
		// 1:name 2:MoverType 3:nrun 4:nperrun 5:score_type 6:cstw 7.relax_type
		imethod++;
		methodname_[imethod] = tokens[1];
		methods_picked_.push_back( imethod );

		params.index = imethod;
		params.roundtype = roundtype;
		params.name = tokens[1];
		params.movertype = tokens[2];
		params.nrun = atoi( tokens[3].c_str() );
		params.nperrun = atoi( tokens[4].c_str() );
		params.istart = 0;
		params.irun = 0;

		// Token 5
		// Score type is defined by Size, not by its name
		// because WUM is communicating with Master via Size type only
		if ( tokens[5].compare("talaris2013") == 0 ) {
			params.score_type = 0;
		} else if ( tokens[5].compare("scorefacts") == 0 ) {
			params.score_type = 1;
		} else if ( tokens[5].compare("goap") == 0 ) {
			params.score_type = 2;
		} else {
			TR << "Unknown score: " << tokens[5] << " at ";
			TR << nline << "th line! Skip."  << std::endl;
		}

		// Token 6
		params.cstw = atof( tokens[6].c_str() );

		// Token 7
		core::Size relaxtype = atoi( tokens[7].c_str() );
		if ( relaxtype >= 1 && relaxtype <= 20 ) {
			params.relax_type = relaxtype;
		} else {
			TR << "Unknown relaxtype: " << relaxtype << " at ";
			TR << nline << "th line! Skip."  << std::endl;
		}

		// Token 8
		core::Size rerelaxtype = atoi( tokens[8].c_str() );
		if ( rerelaxtype == 0 || (rerelaxtype >= 10 && rerelaxtype <= 20) ) {
			params.rerelax_type = rerelaxtype;
		} else {
			TR << "Unknown rerelaxtype: " << rerelaxtype << " at ";
			TR << nline << "th line! Skip."  << std::endl;
		}

		// Token 9,10,11
		params.fshave1 = 0.0;
		params.fshave2 = 0.0;
		params.fshave3 = 0.5;
		if ( tokens.size() >= 9 ) {
			params.fshave1 = atof( tokens[9].c_str() );
			if ( tokens.size() >= 10 ) {
				params.fshave2 = atof( tokens[10].c_str() );
				if ( tokens.size() >= 11 ) {
					params.fshave3 = atof( tokens[11].c_str() );
				}
			}
		}
		params_.push_back( params );
	}

	if ( params_.size() == 0 ) {
		utility_exit_with_message( "[ERROR] Could not find any scheduler for the Master!" );
	}

	// Final line should be always finished by "done"
	if ( !params_[ params_.size() ].roundtype.compare("done") ) {
		add_fresh_param( "done" );
	}

	// Simulate
	TR << "==================================================================" << std::endl;
	TR << "Simulate scheduler from cmd" << std::endl;
	TR << "==================================================================" << std::endl;

	for ( core::Size imaster = 1; imaster <= nmaster_max; ++imaster ) {
		//if( mpi_rank != imaster ) continue;

		TR << "Schedule to be run for master " << imaster << ":" << std::endl;
		TR << "Total " << params_.size() << " stages read, ";
		TR << " with num. iteration " << niter_ << std::endl;

		for ( core::Size isch = 1; isch <= params_.size(); ++isch ) {
			roundtype = params_[isch].roundtype;
			TR << " " << params_[isch].name << std::endl;
			if ( roundtype.compare( "done" ) == 0 ) break;
		}
	}
	TR << "==================================================================" << std::endl;

	infile.close();
}

core::Size
Scheduler::n_to_gen() const
{
	core::Size n_to_gen( 0 );
	for ( core::Size isch = 1; isch <= params_.size(); ++isch ) {
		MethodParams const &param( params_[isch] );
		core::Real fkeep( 1.0 - param.fshave1 );

		if ( param.movertype.compare( "cartnm" ) == 0 || param.movertype.compare( "cartnmcen" ) == 0 ||
				param.movertype.compare( "torsnm" ) == 0 || param.movertype.compare( "torsnmcen" ) == 0 ) {
			//n_to_gen += fkeep * param.nrun * (param.nperrun + 5)*2; // can we make this more pretty??
			n_to_gen += fkeep * param.nrun * (param.nperrun)*2; // can we make this more pretty??
		} else if ( param.movertype.compare( "fraginsert" ) == 0 ||
				param.movertype.compare( "fraginsertcen" ) == 0 ||
				param.movertype.compare( "kiccloser" ) == 0 ||
				param.movertype.compare( "cartcloser" ) == 0 ||
				param.movertype.compare( "partialabinitio" ) == 0 ||
				param.movertype.compare( "partialrefine" ) == 0
				//param.movertype.compare( "ramapert" ) == 0
				) {
			n_to_gen += fkeep * param.nrun;
		} else {
			n_to_gen += fkeep * param.nrun * param.nperrun;
		}
	}
	return n_to_gen;
}

bool
Scheduler::is_excluded( utility::vector1< std::string > picked, std::string name2 ) const
{

	for ( core::Size ipick = 1; ipick <= picked.size(); ++ipick ) {
		std::string name1 = picked[ipick];
		std::pair< std::string, std::string > pair1 = std::make_pair( name1, name2 );
		std::pair< std::string, std::string > pair2 = std::make_pair( name2, name1 );

		for ( core::Size i = 1; i <= exclusions_.size(); ++i ) {
			if ( exclusions_[i] == pair1 || exclusions_[i] == pair2 ) return true;
		}
	}
	return false;
}

MethodParams const
Scheduler::get_params( core::Size const imethod ) const
{
	for ( core::Size i = 1; i <= params_.size(); ++i ) {
		if ( params_[i].index == imethod ) return params_[i];
	}

	// if failed
	utility_exit_with_message( "[ERROR] cannot find the method index from scheduler" );
}

void
Scheduler::add_fresh_param( std::string const & name ){
	MethodParams param;
	param.name = name;
	param.movertype = "";
	param.roundtype = name;
	param.nrun = 0;
	param.istart = 0;
	param.irun = 0;
	param.cstw = 0.0;
	param.nperrun = 0;
	param.score_type = 0;
	param.relax_type = 0;
	param.rerelax_type = 5;
	param.index = 0;
	param.fshave1 = 0.0;
	param.fshave2 = 0.0;
	params_.push_back( param );
}

utility::vector1< core::Size >
Scheduler::pick_enrich_methods( protocols::wum::SilentStructStore const &decoys,
	std::string const & scorename ) const
{
	// First, get scores based on method_col info
	std::map< core::Size, utility::vector1< core::Real > > scores;
	utility::vector1< core::Real > scores_total;

	core::Size i( 0 );
	for ( auto it = decoys.begin();
			it != decoys.end(); ++it, ++i ) {
		core::Size imethod = (*it)->get_energy( "samplemethod" );
		if ( imethod == 0 ) continue;

		core::Real score = (*it)->get_energy( scorename );
		scores[imethod].push_back( score );
		scores_total.push_back( score );
		TR.Debug << "i, imethod, score: " << i << " " << imethod << " " << score << std::endl;
	}

	core::Size npick = nmethods_enrich_max_ > scores.size() ? scores.size() : nmethods_enrich_max_;

	// Overall pool values
	core::Real mean_tot, stdev_tot, cut, median_tot;
	mean_and_stdev( scores_total, 0.5, median_tot, mean_tot, stdev_tot );
	TR.Debug << "mean/stdev: " << mean_tot << " " << stdev_tot << std::endl;

	std::vector< core::Real > Zscores;
	utility::vector1< std::pair< core::Size, core::Real > > Zscorepair;

	for ( std::map< core::Size, utility::vector1< core::Real > >::const_iterator it = scores.begin();
			it != scores.end(); ++it ) {
		core::Size imethod = it->first;
		utility::vector1< core::Real > const &method_scores = it->second;
		core::Real mean, stdev, Zscore, median;

		// get score cut for Zscore calc
		mean_and_stdev( method_scores, params_[imethod].fshave3, cut, mean, stdev );

		// Filter upper cut
		utility::vector1< core::Real > filtered_scores;
		for ( core::Size i = 1; i<= method_scores.size(); ++i ) {
			if ( method_scores[i] <= cut ) {
				filtered_scores.push_back( method_scores[i] );
			}
		}
		mean_and_stdev( filtered_scores, 0.5, median, mean, stdev );

		Zscore = ( median - median_tot ) / stdev_tot;
		TR << "n/median/Zscore/shave3/ntop: " << imethod << " " << methodname(imethod);
		TR << " " << filtered_scores.size() << " " << median << " " << Zscore;
		TR << " " << params_[imethod].fshave3 << " " << std::endl;

		Zscores.push_back( Zscore );
		Zscorepair.push_back( std::make_pair( imethod, Zscore ) );
	}

	std::sort( Zscores.begin(), Zscores.end() );
	core::Real Zscorecut = (Zscores[npick-1] > enrich_Zscore_cut_)
		? enrich_Zscore_cut_ : Zscores[npick-1];
	core::Real Zscoremin = Zscores[npick-1];

	TR << "max pick/ Zscorecut: " << npick << " " << Zscorecut << std::endl;

	// check exclusion
	bool exclusion_found( false );
	utility::vector1< std::string > names_picked;

	for ( utility::vector1< std::pair< core::Size, core::Real > >::const_iterator it
			= Zscorepair.begin(); it != Zscorepair.end(); ++it ) {

		core::Size const imethod = it->first;
		std::string name = methodname_.at( imethod );
		core::Real const Zscore = it->second;

		if ( is_excluded( names_picked, name ) ) {
			TR << "Exclude " << name << " based on exclusion rule!" << std::endl;
			exclusion_found = true;
			continue;
		}

		if ( Zscore <= Zscorecut && Zscore - Zscoremin < Zdiff_outstand_ ) {
			names_picked.push_back( name );
		}
	}

	if ( exclusion_found && Zscores.size() > npick+1 ) {
		Zscorecut = (Zscores[npick-2] > enrich_Zscore_cut_)
			? enrich_Zscore_cut_ : Zscores[npick-2];
	}

	// Pick methods based on methods' Zscores
	utility::vector1< core::Size > methods_picked;
	Zscoremin = 10000.0; // just to make sure
	core::Size imin( 1 );
	for ( utility::vector1< std::pair< core::Size, core::Real > >::const_iterator it
			= Zscorepair.begin(); it != Zscorepair.end(); ++it ) {

		core::Size const imethod = it->first;
		std::string name = methodname_.at( imethod );
		core::Real const Zscore = it->second;

		if ( is_excluded( names_picked, name ) ) continue;

		if ( Zscore <= Zscorecut ) {
			TR << "Pick method based on Zscore, method/Zscore: " << name << " " <<  Zscore << std::endl;
			methods_picked.push_back( imethod );
		}

		if ( Zscore < Zscoremin ) {
			imin = imethod;
			Zscoremin = Zscore;
		}

		if ( methods_picked.size() >= nmethods_enrich_max_ ) break;
	}

	// Return lowest Zscore if nothing picked
	if ( methods_picked.size() == 0 ) {
		TR << "No method with Zscore lower than " << Zscorecut << ", pick method the lowest..." << std::endl;
		methods_picked.push_back( imin );

	}

	return methods_picked;
}

}
}
