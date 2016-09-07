// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: LoopModeler
///
///
/// @author Vatsan Raman

//devel headers
#include <devel/integrated_loop/LoopManager.hh>

//core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

//numeric headers
#include <numeric/random/random.hh>

//protocols headers
#include <protocols/loops/looprelax_protocols.hh>

//C++ headers

using basic::T;
using basic::Error;
using basic::Warning;


namespace protocols {
namespace moves {


static THREAD_LOCAL basic::Tracer TR( "devel.IntegratedLoop.LoopManager" );
using namespace core;

//v typedef utility::vector1< protocols::Loop > Loops;
//v typedef utility::vector1< protocols::Loop >::iterator LoopsIt;
//v typedef utility::vector1< protocols::Loop >::const_iterator LoopsConsIt;

//////////////////////////////////////////////////////////////////////////
/*Size LoopManager::NumSingleLoops()
{
	return LoopList_.size();
}
*/
//////////////////////////////////////////////////////////////////////////
loops::Loops LoopManager::LoopsToPerturb()
{

	loops tmpLoops;

	std::cout << "LoopManager " << pose_.total_residue() << std::endl;

	for ( LoopsIt it = LoopList_.begin(), it_end = LoopList_.end(); it != it_end; ++it; )
		{
			float skip_rate_random = numeric::random::rg().uniform();
			std::cout << "skip_rate_random " << skip_rate_random << std::endl;
			if( skip_rate_random > it->skip_rate() )
				//				tmpLoops.push_back( VaryStems( VaryCutpoint( *it ) ) );
				if( skip_rate_random > it->skip_rate() ) tmpLoops.push_back( VaryStems( VaryCutpoint( *it ) ) );
		}

	return tmpLoops;

}
///////////////////////////////////////////////////////////////////////////
bool LoopManager::IsNtermLoop(
	Loop const & ThisLoop
);
{

	if( ThisLoop.loop_begin() == 1 ) return true;
	return false;

}
///////////////////////////////////////////////////////////////////////////
bool LoopManager::IsCtermLoop(
	Loop const & ThisLoop
)
{

	if( ThisLoop.loop_end() == pose_.total_residue() ) return true;
	return false;

}
///////////////////////////////////////////////////////////////////////////
void LoopManager::ReorderLoops()
{

	if( LoopList_.size() <= 1 ) return;
	loops new_loop_list;
	LoopsIt it_begin = LoopList_.begin();
	new_loop_list.push_back( *it_begin );

	for ( LoopsIt it = ++it_begin, it_end = LoopList_.end();
				it != it_end; ++it ) {
		bool inserted = false;
		for( LoopsIt it2 = new_loop_list.begin(), it2_end = new_loop_list.end();
				 it2 != it2_end; ++it2 ) {
			if ( it->loop_begin() < it2->loop_begin() ) {
				new_loop_list.insert( it2, *it );
				inserted = true;
				break;
			}
		}
		if ( ! inserted ) new_loop_list.push_back( *it );
	}

	assert( LoopList_.size() == new_loop_list.size() );
	LoopList_ = new_loop_list;

	return;
}
///////////////////////////////////////////////////////////////////////////
Loop LoopManager::PreviousLoop(
	Loop const & ThisLoop
)
{

	Loop tmp_loop;
	for( LoopsIt it = LoopList_.begin(), it_end = LoopList_.end();
			 it != it_end; ++it )
		{
			if( it->loop_begin() == ThisLoop.loop_begin() &&
				it->loop_end() == ThisLoop.loop_end() )
				tmp_loop = *( it - 1 );
		}

	return tmp_loop;
}
///////////////////////////////////////////////////////////////////////////
Loop LoopManager::NextLoop(
	Loop const & ThisLoop
)
{

	Loop tmp_loop;
	for( LoopsIt it = LoopList_.begin(), it_end = LoopList_.end();
			 it != it_end; ++it )
		{
			if( it->loop_begin() == ThisLoop.loop_begin() &&
				it->loop_end() == ThisLoop.loop_end() )
				tmp_loop = *( it + 1 );
		}
	return tmp_loop;
}

///////////////////////////////////////////////////////////////////////////
bool LoopManager::IsFirstLoop(
	Loop const & ThisLoop
)
{

	if( LoopList_.front().loop_begin() == ThisLoop.loop_begin() &&
		LoopList_.front().loop_end() == ThisLoop.loop_end() )
		return true;
	return false;

}


///////////////////////////////////////////////////////////////////////////
bool LoopManager::IsLastLoop(
	Loop const & ThisLoop
)
{

	if( LoopList_.back().loop_begin() == ThisLoop.loop_begin() &&
		LoopList_.back().loop_end() == ThisLoop.loop_end() )
		return true;
	return false;


}
///////////////////////////////////////////////////////////////////////////
Loop LoopManager::VaryStems(
	Loop const & ThisLoop
)
{

	int window_size = 5;//Not sure if this is the right size.
	Size stem_vary_window_n = 0;
	Size stem_vary_window_c = 0;
	if( !IsFirstLoop( ThisLoop ) && !IsLastLoop( ThisLoop ) ) {

		stem_vary_window_n = std::min( window_size, int( ThisLoop.loop_begin() -
				PreviousLoop( ThisLoop ).loop_end() ) );

		stem_vary_window_c = std::min( window_size, int( NextLoop( ThisLoop ).loop_begin() -
				ThisLoop.loop_end() ) );
	} else if ( !IsFirstLoop( ThisLoop )  && IsLastLoop( ThisLoop ) ) {
		stem_vary_window_n = std::min( window_size, int( ThisLoop.loop_begin() -
				PreviousLoop( ThisLoop ).loop_end() ) );
		if( !IsCtermLoop( ThisLoop ) ) {
			stem_vary_window_c = std::min( window_size, int( pose_.total_residue()  -
					ThisLoop.loop_end() ) );
		} else {
			stem_vary_window_c = 0;
		}
	} else if ( !IsLastLoop( ThisLoop )  && IsFirstLoop( ThisLoop ) ) {

		stem_vary_window_c = std::min( window_size, int( NextLoop( ThisLoop ).loop_begin() -
				ThisLoop.loop_end() ) );
		if( !IsNtermLoop( ThisLoop ) ) {
			stem_vary_window_n = std::min( window_size, int( ThisLoop.loop_begin() ) );
		} else {
			stem_vary_window_n = 0;
		}
	} else if ( IsFirstLoop( ThisLoop)  && IsLastLoop( ThisLoop ) ) {

		stem_vary_window_n = std::min( window_size, int( ThisLoop.loop_begin() ) );
		stem_vary_window_c = std::min( window_size, int( pose_.total_residue() -
				ThisLoop.loop_end() ) );
	}


	Size loop_begin = ThisLoop.loop_begin() - Size( stem_vary_window_n*loopmanager_RG.uniform() );
	Size loop_end = ThisLoop.loop_end() + Size( stem_vary_window_c*loopmanager_RG.uniform() );

	TR << "Loop in LoopManager " << loop_begin << " " << loop_end << " " << ThisLoop.cutpoint() << " " << ThisLoop.skip_rate() << "\n";
	std::cout << "Loop in LoopManager " << loop_begin << " " << loop_end << " " << ThisLoop.cutpoint() << " " << ThisLoop.skip_rate() << std::endl;
	return Loop( loop_begin, loop_end, ThisLoop.cutpoint(), ThisLoop.skip_rate() );

}
///////////////////////////////////////////////////////////////////////////
Loop LoopManager::VaryCutpoint(
	Loop const & ThisLoop
)
{
	if( ThisLoop.cutpoint() != 0 ) return ThisLoop;

	Size const loop_size( ThisLoop.loop_end() - ThisLoop.loop_begin() );
	Size const n_cutpoints( loop_size - 1 );
	std::vector< Real > cut_weight;
	Real total_cut_weight( 0.0 );

	for( Size i = 1; i <= n_cutpoints; ++i ) {
		Real const weight( std::max( i, n_cutpoints - i + 1 ) );
		total_cut_weight += weight;
		cut_weight.push_back( total_cut_weight );
	}


	Size cutpoint( 0 );
	if( ThisLoop.loop_begin() > 1 && ThisLoop.loop_end() < pose_.total_residue() ) {

		Size nfail( 0 );
		do {
			nfail++;
			Real const weight( total_cut_weight*loopmanager_RG.uniform() );
			for( Size i = 1; i <= n_cutpoints; ++i ) {
				if( weight <= cut_weight.at( i - 1 ) ) {
					cutpoint = ThisLoop.loop_begin() + i - 1;
					break;
				}

			}

		} while( nfail < 20 ); //Need to add secondary structure logic here
		if( cutpoint == 0 ) {//not sure if it is actually needed though
			cutpoint = ThisLoop.loop_begin() + n_cutpoints/2;
			Warning() << "cutpoint choice problem; setting cutpoint = "
								<< cutpoint << "\n";
		}
	} else if ( ThisLoop.loop_begin() == 1 ) {
		cutpoint = 1;
	} else if ( ThisLoop.loop_end() == pose_.total_residue() ) {
		cutpoint = pose_.total_residue();
	}


	return Loop( ThisLoop.loop_begin(), ThisLoop.loop_end(), cutpoint, ThisLoop.skip_rate() );

}

///////////////////////////////////////////////////////////////////////////

}//moves
}//protocols
