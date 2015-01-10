// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/mpi_refinement/MultiObjective.cc
/// @brief
/// @author Hahnbeom Park

#include <protocols/mpi_refinement/util.hh>
#include <protocols/mpi_refinement/MultiObjective.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pose/Pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// Utility headers
#include <math.h>
#include <numeric/random/random.hh>

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace mpi_refinement {

static basic::Tracer TR("MPI.LHR.O");

// Methods should be REGISTERED here in order to use as multi-objective
MultiObjective::MultiObjective()
{
	set_defaults();
}

MultiObjective::~MultiObjective(){}

void
MultiObjective::set_defaults()
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	objsfxnOPs_.resize( 0 );
	fobjnames_.resize( 0 );

	fobjnames_.push_back( "score" );

	// Add extra multi objective functions 
	if( !option[ OptionKeys::lh::multi_objective_functions ].user() ){
		TR << "Add default extra objective function." << std::endl;
		fobjnames_.push_back( "goap" );
		fobjnames_.push_back( "similarity" );

	} else {
		TR << "Setup multi-objective function defined by user." << std::endl;
		utility::vector1< std::string > fobjnames_add
			= option[ OptionKeys::lh::multi_objective_functions ]();

		for( core::Size i_obj = 1; i_obj <= fobjnames_add.size(); ++i_obj )
			fobjnames_.push_back( fobjnames_add[i_obj] );
	}

	// use only if specified by user
	//if( option[ constraints::cst_file ].user() ) fobjnames_.push_back( "cst_cen" );
	//if( option[ constraints::cst_fa_file ].user() ) fobjnames_.push_back( "cst_fa" );

	for( core::Size i_obj = 1; i_obj <= fobjnames_.size(); ++i_obj ){
		std::string const score_name( fobjnames_[i_obj] );

		if( score_name.compare("goap") == 0 ){
			TR << "- " << i_obj << ". Added Goap potential as 'goap'" << std::endl;
			core::scoring::ScoreFunctionCOP sfxn =
				core::scoring::ScoreFunctionFactory::create_score_function( "goap" );
			objsfxnOPs_.push_back( sfxn );

		} else if( score_name.compare("score") == 0){
			TR << "- " << i_obj << ". Added input score as 'score'" << std::endl;
			core::scoring::ScoreFunctionOP sfxn = 
				core::scoring::ScoreFunctionFactory::create_score_function( option[ score::weights ]() );
			// modify cart_bonded term to be less sensitive to small non-ideality
			// note that this is just for "evaluation", not for sampling
			//sfxn->set_weight( core::scoring::cart_bonded, 0.1 );
			objsfxnOPs_.push_back( sfxn );
		} else if( score_name.compare("similarity") == 0){
			TR << "- " << i_obj << ". Added similarity as 'similarity'" << std::endl;
			// Just to ensure the length, but never been used
			core::scoring::ScoreFunctionOP sfxn = 
				core::scoring::ScoreFunctionFactory::create_score_function( option[ score::weights ]() );
			objsfxnOPs_.push_back( sfxn );

		} else if( score_name.compare("cst_fa") == 0 ||
							 score_name.compare("cst_cen") == 0 ){
			TR << "- " << i_obj << ". Added " << score_name << " as '" << score_name << "'" << std::endl;
			core::scoring::ScoreFunctionOP sfxn_cst( new core::scoring::ScoreFunction );
			sfxn_cst->set_weight( core::scoring::atom_pair_constraint, 1.0 );
			sfxn_cst->set_weight( core::scoring::coordinate_constraint, 1.0 );
			objsfxnOPs_.push_back( sfxn_cst );

		} else {
			TR << "WARNING: Skip unrecognized multi-objective function: " << score_name << std::endl;
		}
  }

	// dominant cut initialize
	utility::vector1< core::Real > objective_cut_from_cmd =
		option[ lh::objective_dominate_cut ]();

	utility::vector1< core::Real > objective_inc_from_cmd =
		option[ lh::objective_cut_increment ]();

	if( nobjs() != objective_cut_from_cmd.size() ){
		utility_exit_with_message( "lh::objective_dominate_cut num. args are different from actual objective function definition!");
	}
	if( nobjs() != objective_inc_from_cmd.size() ){
		utility_exit_with_message( "lh::objective_dominate_cut num. args are different from actual objective function definition!");
	}

	obj_dominant_cut_ = objective_cut_from_cmd;
	obj_cut_increment_ = objective_inc_from_cmd;
}

// Check if ss1 is better than ss2 by any obj function
bool
MultiObjective::is_dominant( core::io::silent::SilentStructCOP ss1,
														 core::io::silent::SilentStructCOP ss2 )
{
	utility::vector1< std::string > fobj_loc( fobjnames_ );

  for( core::Size iobj = 1; iobj <= fobj_loc.size(); ++iobj ){
    std::string const score_name( fobj_loc[iobj] );
    core::Real dobj = ss1->get_energy( score_name ) - ss2->get_energy( score_name );

    // when ss1 no more dominates ss2
		if( dobj > obj_dominant_cut_[iobj] ) return false;
  }

  return true;
}

	/*
core::Real
MultiObjective::fobj_density( protocols::wum::SilentStructStore const &ref_structs,
															protocols::wum::SilentStructStore &structs,
															bool const is_symmetric )
{

	core::Size const nstruct_ref( ref_structs.size() );
	core::Size const nstruct    ( structs.size() );

	// First, fill in Similarity matrix
	utility::vector1< utility::vector1< core::Real > > sim_matrix( nstruct );
	for( core::Size i = 1; i <= nstruct; ++i ) sim_matrix[i].resize( nstruct_ref );

  for( core::Size i = 1; i <= nstruct; ++i ){
    SilentStructCOP p( totalpool.get_struct( i_p ) );
		for( core::Size j = 1; i <= nstruct; ++i ){
		
		}
	}

}
	*/

// Takes care of objective function "similarity"
void
MultiObjective::calculate_pool_diversity(
    protocols::wum::SilentStructStore &structs ) const
{
	core::Size const n( structs.size() );

	for( core::Size i = 0; i < n; ++i ){
		core::io::silent::SilentStructOP ss1 = structs.get_struct( i );
		calculate_structure_diversity( ss1, structs );
	}
}

void
MultiObjective::calculate_pool_diversity(
				 protocols::wum::SilentStructStore &structs1, 
				 protocols::wum::SilentStructStore &structs2 
																				 ) const
{
	core::Size const n( structs1.size() );

	for( core::Size i = 0; i < n; ++i ){
		core::io::silent::SilentStructOP ss1 = structs1.get_struct( i );
		calculate_structure_diversity( ss1, structs2 );
	}
}

void
MultiObjective::calculate_structure_diversity( 
		core::io::silent::SilentStructOP ss1,
    protocols::wum::SilentStructStore &structs ) const
{

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	//core::Real simlimit = option[ lh::similarity_reference ]();
	core::Real simlimit = option[ lh::rms_limit ]();
	std::string similarity_method = option[ lh::similarity_method ](); 
	std::string similarity_measure = option[ lh::similarity_measure ](); 
	core::Size const n( structs.size() );
	core::Real const simtol = option[ lh::similarity_tolerance ]();

	//if measured by similarity to Emin
	utility::vector0< core::Real > dvector( n, 0.0 );
	core::Size nrel( dvector.size() );
	for( core::Size j = 0; j < n; ++j ){
		core::io::silent::SilentStructOP ss2 = structs.get_struct( j );
		std::string ssname = ss2->decoy_tag();

		// avoid itself
		if( ssname.compare( ss1->decoy_tag() ) == 0 ){
			nrel--;
			dvector[j] = 0.0;
		} else {
			core::Real dumm, dist( 0.0 );
			if( similarity_measure.compare( "Sscore" ) == 0 ){
				dist = CA_Sscore( ss1, ss2, dumm, 2.0 );
				//dist = CA_Sscore( ss1, ss2, dumm, 1.0 );

			} else if (similarity_measure.compare( "rmsd" ) == 0 ){
				dumm = CA_Sscore( ss1, ss2, dist, 2.0 );
				//dumm = CA_Sscore( ss1, ss2, dist, 1.0 );

			} else if (similarity_measure.compare( "looprmsd" ) == 0 ){
				std::string loopstr = option[ lh::loop_string ]();
				utility::vector1< core::Size > loopres = loopstring_to_loopvector( loopstr );
				dumm = CA_Sscore( ss1, ss2, dist, loopres, 2.0 );
				//dumm = CA_Sscore( ss1, ss2, dist, loopres, 1.0 );
			} else {
				
			}
			dvector[j] = dist;
			//TR << "d: " << j << " " << dist << std::endl; 
		}
	}

	core::Real similarity( 0.0 );
	// Diversity based on similarity sum
	if( similarity_method.compare( "sum" ) == 0 ){
		for( core::Size j = 0; j < n; ++j ){
			similarity += dvector[j];
		}
		similarity *= 100.0/nrel;

	} else if( similarity_method.compare( "sigsum" ) == 0 ){
		for( core::Size j = 0; j < n; ++j ){
			if( similarity_measure.compare( "Sscore" ) == 0 ){
				if( dvector[j] >= simlimit ){
					similarity += 1.0;
				} else if( dvector[j] > simlimit-simtol ){
					similarity += 0.5 + 0.5*(dvector[j]-simlimit)/simtol;
				}
			} else {
				if( dvector[j] <= simlimit ){
					similarity += 1.0;
				} else if( dvector[j] <= simlimit+simtol ){
					similarity += 0.5 - 0.5*(dvector[j]-simlimit)/simtol;
				}
			}
		}
		similarity *= 100.0/nrel;

		// Diversity based on closest one
	} else if( similarity_method.compare( "min" ) == 0 ||
						 similarity_method.compare( "frontiermin" ) == 0 ) {
		for( core::Size j = 0; j < n; ++j ){
			if( dvector[j] != 1.0 && dvector[j] > similarity ) similarity = dvector[j]*100.0;
		}
	}

	ss1->add_energy( "similarity", similarity );

}

// Update library pool based on NSGAII rule, the multiobj GA 
bool
MultiObjective::update_library_NSGAII(protocols::wum::SilentStructStore &structs, 
																			protocols::wum::SilentStructStore &new_structs,
																			core::Size const nmax,
																			bool const update_obj_cut
																			)
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::io::silent;

	core::Real const simlimit = option[ OptionKeys::lh::rms_limit ]();
	std::string const measure = option[ lh::similarity_measure ]();

  // Total pool
  protocols::wum::SilentStructStore totalpool;
  totalpool.clear();
	// Parent pool
	protocols::wum::SilentStructStore prvpool( structs ); 
	std::string const sim_replace_obj = option[ OptionKeys::lh::sim_replace_obj ]();
	std::string const similarity_method = option[ lh::similarity_method ](); 

	// increment on obj_dominant cut
	if( update_obj_cut ){
		for( core::Size i_obj = 1; i_obj <= nobjs(); ++i_obj )
			obj_dominant_cut_[i_obj] += obj_cut_increment_[i_obj];
	}
	// report current obj_cut
	TR << "Current obj_dominate_cut: ";
	for( core::Size i_obj = 1; i_obj <= nobjs(); ++i_obj )
		TR << "/ " << fobjnames_[i_obj] << " " << F(6,1,obj_dominant_cut_[i_obj]);
	TR << std::endl;

	core::Size nbefore;
	nbefore = new_structs.size();
	filter_similar( new_structs, measure, simlimit, sim_replace_obj );

	//calculate_pool_diversity( new_structs );
	//filter_similar( new_structs, measure, simlimit, "similarity" );
	TR << "New structures, initial/filtered_by_similiarity:" << nbefore << " " << new_structs.size() << std::endl;
  totalpool.add( new_structs );

	// then filter out within the whole pool
	// but only for unexpired parents
	core::Size const expire_cut = option[ OptionKeys::lh::expire_after_rounds ]();
  for( core::Size i_p = 0; i_p < structs.size(); ++i_p ){
		core::io::silent::SilentStructOP ss = structs.get_struct( i_p );
		if( ss->get_energy( "nuse" ) <= expire_cut ) {
			totalpool.add( ss );
		} else {
			TR << "Limit parent " << ss->decoy_tag() << " based on nuse cut." << std::endl;
		}
	}

	nbefore = totalpool.size();

	core::Size nmax_filt( nmax );
	//if( option[ lh::filter_up_to_maxlib ]() ) nmax_filt = nmax;
	filter_similar( totalpool, measure, simlimit, sim_replace_obj, nmax_filt );

	TR << "Total pool, initial/filtered_by_similarity: " << nbefore << " " << totalpool.size() << std::endl;

  if( totalpool.size() <= nmax_filt ){
    structs = totalpool;
		new_structs.clear();
		TR << "NSGAII call, but Direct addition due to empty spaces" << std::endl;
    return true;
  }

	// get list of unused pools

	// entering diversity

	// org:
	totalpool.all_add_energy( "similarity", 0.0 ); //just initialize
	if( similarity_method.compare( "frontiermin" ) != 0 ){
		calculate_pool_diversity( totalpool, totalpool );
		//calculate_pool_diversity( totalpool, prvpool );
	}

	TR << "update_library based on NSGAII, nprv/nnew/nfilt/nmax_filt = ";
	TR << structs.size() << "/" << new_structs.size() << "/";
	TR << totalpool.size() << "/" << nmax_filt << std::endl;

	// Be careful - all the indices start from 0
  utility::vector0< protocols::wum::SilentStructStore > frontier( 1 );
  utility::vector0< utility::vector0< core::Size > > ifront( 1 );
  core::Size nstruct( totalpool.size() );
  utility::vector0< utility::vector0< core::Size > > Sp( nstruct );
  utility::vector0< core::Size > np( nstruct );
	utility::vector0< bool > is_left( nstruct, true );

  // Fast non-dominated-sort
  for( core::Size i_p = 0; i_p < nstruct; ++i_p ){
    SilentStructCOP p( totalpool.get_struct( i_p ) );
		Sp[i_p].resize( 0 );
		np[i_p] = 0;

    for( core::Size i_q = 0; i_q < nstruct; ++i_q ){
      if( i_p == i_q ) continue;
      SilentStructCOP q( totalpool.get_struct( i_q ) );

			/*
			TR << "i_p/i_q/score/goap/sim: " << i_p << " " << i_q;
			for( core::Size iobj = 1; iobj <= fobjnames_.size(); ++iobj ) 
				TR << " " << p->get_energy( fobjnames_[iobj] ) - q->get_energy( fobjnames_[iobj] );
			TR << " " << is_dominant( p, q) << " " << is_dominant( q, p ) << std::endl;
			*/

      if( is_dominant( p, q ) ){ // p dominates q
				Sp[i_p].push_back( i_q );
      } else if ( is_dominant( q, p ) ) { // q dominates p
				np[i_p]++; // N non-dominating 
      }
    }

    // if not dominiated by any, set as first-frontier
    if( np[i_p] == 0 ){
			frontier[0].add( totalpool.get_struct( i_p ) );
			ifront[0].push_back( i_p );
			for( core::Size j = 0; j < ifront[0].size(); ++j ) is_left[ ifront[0][j] ] = false;
		}
  }
	frontier[0].all_add_energy( "frontier", 0 );

	TR.Debug << "frontier1 (size " << ifront[0].size() << ": ";
	for( core::Size j = 0; j < ifront[0].size(); ++j ) TR << " " << ifront[0][j];
	TR.Debug << ")" << std::endl;

	// Debug
	for( core::Size j = 0; j < nstruct; ++j ){
		TR.Debug << "i_p/np/Sp: " << j << " " << np[j] << " | ";
		for( core::Size k = 0; k < Sp[j].size(); ++k ) TR.Debug << " " << Sp[j][k];
		TR.Debug << std::endl;
	}

  // Assign frontiers
  core::Size i( 0 );
  core::Size n_in_frontier( frontier[0].size() );
	protocols::wum::SilentStructStore sorted_pool;

  while( true ){
    // Stop if assigned enough 
    if( n_in_frontier >= nmax )	break;

	  sorted_pool.add( frontier[i] );

		// update diversity on currently added so far if sorted_pool is large enough
		// experimental
		if( similarity_method.compare("frontiermin") == 0 )
			calculate_pool_diversity( totalpool, sorted_pool );

    protocols::wum::SilentStructStore Q; // Next frontier
		utility::vector0< core::Size > iQ;
		protocols::wum::SilentStructStore const &P = frontier[i];
		utility::vector0< core::Size > const &iP = ifront[i];

    for( core::Size istr = 0; istr < P.size(); ++istr ){
      core::Size i_p = iP[istr];

      // When p pops out, p's next frontier members' dominant number reduces by 1
      for( core::Size jstr = 0; jstr < Sp[i_p].size(); ++jstr ){
				core::Size i_q = Sp[i_p][jstr];
				np[i_q]--; // mark that 

				if( np[i_q] == 0 ){ // If all the dominating members are already poped out
					Q.add( totalpool.get_struct( i_q ) );
					iQ.push_back( i_q );
				}
      }
    }

		i++;
		Q.all_add_energy( "frontier", i );
    frontier.push_back( Q );
		ifront.push_back( iQ );
		for( core::Size j = 0; j < iQ.size(); ++j ) is_left[ iQ[j] ] = false;
    n_in_frontier += Q.size();

		TR.Debug << "frontier" << i+1 << "(size " << ifront[i].size() << ": ";
		for( core::Size j = 0; j < ifront[i].size(); ++j ) TR.Debug << " " << ifront[i][j];
		TR.Debug << "), left: " << int(nmax) - int(n_in_frontier) << std::endl;
  }

  // Add into sorted_pool
  // up to n-1 frontier, store all
  //for( core::Size i_front = 0; i_front < frontier.size()-1; ++i_front ){
  //  sorted_pool.add( frontier[i_front] );
	//}

  // for the last frontier, sort by certain objective (e.g. the first objective)
  // Make sure your objective is reasonable for doing this; 
  protocols::wum::SilentStructStore &last_frontier = frontier[ frontier.size()-1 ];
	//utility::vector0< core::Size > &ilast = ifront[ frontier.size()-1 ];

  //last_frontier.sort_by( fobjnames(1) );
  //last_frontier.sort_by( "goap" );
	// try using "similarity to the already selected pool" for picking among the last frontier
	for( core::Size j = 0; j < last_frontier.size(); ++j ){
		if( sorted_pool.size() > 0 ){
			calculate_structure_diversity( last_frontier.get_struct(j), sorted_pool );
		} else {
			calculate_structure_diversity( last_frontier.get_struct(j), last_frontier );
		}
	}
  last_frontier.sort_by( "similarity" );

  core::Size const nleft( nmax - sorted_pool.size() );

	if( frontier.size() == 1 ){
		TR << "NSGAII: Adding on pool, " << sorted_pool.size() << " from first frontier" << std::endl;
	} else {
		TR << "NSGAII: Adding on pool, " << sorted_pool.size() << " from first " << frontier.size() - 2 << " frontiers";
		TR << " and " << nleft << " from the last frontier" << std::endl;
	}

	// new_structs will return the removed ones
	new_structs.clear();
  for( core::Size i_left = 0; i_left < last_frontier.size(); ++i_left ){
		if( i_left < nleft ){
			sorted_pool.add( last_frontier.get_struct( i_left ) );
		} else {
			new_structs.add( last_frontier.get_struct( i_left ) );
		}
	}

	// remaining frontiers
	core::Size const final_frontier( frontier.size() );
	for( core::Size j = 0; j < nstruct; ++j ){
		if( is_left[j] ){
			totalpool.get_struct( j )->add_energy( "frontier", final_frontier );
			new_structs.add( totalpool.get_struct( j ) );
		}
	}

  // Update library
  structs.clear();
  structs.add( sorted_pool );
	
	// finally add diversity info : different from number used for frontier decision,
	// but will give estimation to the next iteration
	//calculate_pool_diversity( structs );

	return true;
}

void
MultiObjective::filter_similar( protocols::wum::SilentStructStore &structs,
																std::string const measure,
																core::Real const criteria,
																std::string const score_for_priority,
																core::Size const nmax
																) 
{

  using namespace core::io::silent;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  core::Size const nstruct( structs.size() );
  utility::vector0< bool > filtered( nstruct, false );
  utility::vector0< core::Size > nuse( nstruct, 0 );

	core::Size nfilt( 0 );
	core::Size nfilt_max = (nstruct >= nmax) ? nstruct - nmax : 0;

  for( core::Size i = 1; i < nstruct; ++i ){
    SilentStructOP ss1 = structs.get_struct( i );
    core::Real score1 = ss1->get_energy( score_for_priority );
		std::string const iname = ss1->decoy_tag();

    for( core::Size j = 0; j < i; ++j ){
      SilentStructOP ss2 = structs.get_struct( j );
      core::Real score2 = ss2->get_energy( score_for_priority );
			std::string const jname = ss2->decoy_tag();

 			if( nfilt >= nfilt_max ) break;

      if( filtered[j] ) continue;

      // only Sscore for now
      bool is_similar( false );
      core::Real dist( 0.0 );
			core::Real dumm( 0.0 );

      if( measure.compare( "Sscore" ) == 0 ){
				dist = CA_Sscore( ss1, ss2, dumm, 2.0 );
				if( dist > criteria ) is_similar = true;

      } else if (measure.compare( "rmsd" ) == 0 ){
				dumm = CA_Sscore( ss1, ss2, dist, 2.0 );
				if( dist < criteria ) is_similar = true;

			} else if (measure.compare( "looprmsd" ) == 0 ){
				std::string loopstr = option[ lh::loop_string ]();
				utility::vector1< core::Size > loopres = loopstring_to_loopvector( loopstr );
				dumm = CA_Sscore( ss1, ss2, dist, loopres, 2.0 );
				if( dist < criteria ) is_similar = true;

			}

			//TR << "check similarity, " << i << " " << j << " with measure " << measure << ": " << dist << " " << criteria << std::endl;

      if( !is_similar ) continue;

			// get larger nuse if similar
      if( score1 < score2 ){
				filtered[j] = true;
				nfilt++;
				nuse[i] = nuse[j] > nuse[i] ? nuse[j] : nuse[i];
				TR << "remove " << jname << " on " << iname << ", dist/dE " << dist << "/" << score2 - score1 << std::endl;
				break;
      } else {
				filtered[i] = true;
				nfilt++;
				nuse[j] = nuse[j] > nuse[i] ? nuse[j] : nuse[i];
				TR << "remove " << iname << " on " << jname << ", dist/dE " << dist << "/" << score1 - score2 << std::endl;
      }
    } // j

		if( nfilt >= nfilt_max ) break;
  } // i

  core::Size i( 0 );
	protocols::wum::SilentStructStore outstructs;
  for( protocols::wum::SilentStructStore::const_iterator it = structs.begin();
       it != structs.end(); ++it, ++i ){
    if( !filtered[i] ){
			core::io::silent::SilentStructOP ss( *it );
			ss->add_energy( "nuse", nuse[i] ); //updated based on max within similar
			outstructs.add( ss );
		}
  }

	structs = outstructs;
}

void
MultiObjective::add_objective_function_info( 
																						core::io::silent::SilentStructOP ss,
																						protocols::wum::SilentStructStore & ) const
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

	// Get a copy
	core::io::silent::SilentStructCOP ss_copy = ss->clone();
	
	// First cleanup useless information so as to make it "readable"
	ss->clear_energies();

	// Equally add information from ss no matter how they are generated
	ss->add_energy( "state", ss_copy->get_energy("state") );
	ss->add_energy( "samplemethod", ss_copy->get_energy("samplemethod") );
	ss->add_energy( "round", ss_copy->get_energy("round") );
	ss->add_energy( "nuse", ss_copy->get_energy("nuse") );
	ss->add_energy( "iter", ss_copy->get_energy("iter") );
	ss->add_energy( "NMmode", ss_copy->get_energy("NMmode") );
	ss->add_energy( "NMscale", ss_copy->get_energy("NMscale") );

	// copy loopinfo
	core::Size nloop = (core::Size)( ss->get_energy( "nloop" ) );
	if( nloop > 0 ){
		//ss->add_string_value( "loops_nsampled", ss_copy->get_string_value( "loops_nsampled" ) );
		for( core::Size iloop = 1; iloop <= nloop; ++iloop ){
			std::stringstream sstream( "" );
			sstream << "nsampled_loop" << iloop;
			ss->add_energy( sstream.str(), ss_copy->get_energy( sstream.str() ) );
		}
	}

	// Get native info
	core::chemical::ResidueTypeSetCOP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::chemical::ResidueTypeSetCOP rsd_set_cen
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	if( option[ in::file::native ].user() ){
		core::pose::Pose native_pose;
    core::import_pose::pose_from_pdb( native_pose, *rsd_set, option[ in::file::native ]() );
		add_nativeinfo_to_ss( *ss, native_pose );
	}

	protocols::moves::MoverOP tocen 
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

  // Iter over multi objective functions
  for( core::Size i_obj = 1; i_obj <= nobjs(); ++ i_obj ){
    std::string score_name( fobjnames( i_obj ) );
    core::Real score_i = ss->get_energy( score_name );
		//TR << "Fobj: " << score_name << std::endl;
    // Calculate energy if empty yet
		if( score_name.compare("score") == 0 ){
			core::pose::Pose pose_tmp;
			ss->fill_pose( pose_tmp, *rsd_set );
			objsfxnOPs_[i_obj]->score( pose_tmp );
			ss->energies_from_pose( pose_tmp );

		} else if( score_name.compare("similarity") == 0 ){
			continue;

		} else if( score_name.compare("cst_cen") == 0 ){
			core::pose::Pose pose_tmp;
			ss->fill_pose( pose_tmp, *rsd_set_cen );
			core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose_tmp );
      score_i = objsfxnOPs_[i_obj]->score( pose_tmp );
			ss->add_energy( "cst_cen", score_i );

		} else if( score_name.compare("cst_fa") == 0 ){
			core::pose::Pose pose_tmp;
			ss->fill_pose( pose_tmp, *rsd_set );
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose_tmp );
      score_i = objsfxnOPs_[i_obj]->score( pose_tmp );
			ss->add_energy( "cst_fa", score_i );

		} else if( score_i == 0.0 && score_name.compare("score") != 0 ){
      core::pose::Pose pose_tmp;
      ss->fill_pose( pose_tmp );
      score_i = objsfxnOPs_[i_obj]->score( pose_tmp );
      ss->add_energy( score_name, score_i );
    }
  }

	// add ediff if both score and goap exists
	if( fobjnames_.contains("score") && fobjnames_.contains("goap") ){
		ss->add_energy( "ediff", ss->get_energy("goap") - ss->get_energy("score") );
		ss->add_energy( "esum", ss->get_energy("goap") + ss->get_energy("score") );
	}

	// constraints

	// finally add diversity info
	// skip?
	//calculate_structure_diversity( ss, sstore ); 
}

void
MultiObjective::add_objective_function_info( protocols::wum::SilentStructStore & sstore ) const 
{
  core::Size const nstruct( sstore.size() );
  for( core::Size i = 0; i < nstruct; ++i ){
		add_objective_function_info( sstore.get_struct( i ), sstore );
	}
}


std::string
MultiObjective::formatted_objs_values( core::io::silent::SilentStruct const &ss ) const
{
  std::stringstream sstream;

  for( core::Size iobj = 1; iobj <= nobjs(); ++iobj ){
		std::string const score_name( fobjnames(iobj) );
    sstream << " | " << F(6,1, ss.get_energy( score_name ) );
  }
	return sstream.str();
}

std::string
MultiObjective::formatted_objs_names() const
{
  std::stringstream sstream;

  for( core::Size iobj = 1; iobj <= nobjs(); ++iobj ){
		std::string const scorename( fobjnames(iobj) );
    sstream << " " << std::setw(7) << scorename;
  }
	return sstream.str();
}

}
}
