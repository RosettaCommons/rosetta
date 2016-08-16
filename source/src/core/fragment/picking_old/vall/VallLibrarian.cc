// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/VallLibrarian.cc
/// @brief  Librarian that picks fragments from the Vall
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/VallLibrarian.hh>

#include <core/fragment/picking_old/vall/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


// static initialization
static THREAD_LOCAL basic::Tracer TR( "core.fragment.picking_old.vall.VallLibrarian" );


/// @brief default constructor
VallLibrarian::VallLibrarian() :
	Super(),
	preallocate_( true )
{}


/// @brief default destructor
VallLibrarian::~VallLibrarian()
{}


/// @brief get top 'N' fragments from prior catalog()
VallLibrarian::FragDataOPs VallLibrarian::top_fragments(
	Size const n,
	BBTorsionSRFD const & srfd_type
) const {
	return fragments( 1, n, srfd_type );
}


/// @brief get fragments from prior catalog() [from, to]
/// @param from index of the first fragment in the list, indexing starts from '1'
/// @param to index of the last fragment in the list (inclusive)
/// @return filled FragDataOPs if sort() was called successfully, otherwise empty FragDataOPs
VallLibrarian::FragDataOPs VallLibrarian::fragments(
	Size from,
	Size to,
	BBTorsionSRFD const & srfd_type
) const
{
	FragDataOPs fdl;

	if ( scores().size() == 0 ) {
		return fdl;
	}

	// auto correct range in case of error
	if ( from < 1 ) {
		from = 1;
	}
	if ( to > scores().size() ) {
		to = scores().size();
	}

	Scores::const_iterator begin = scores().begin() + ( from - 1 );
	Scores::const_iterator end = scores().begin() + to;

	TR << "best fragment:   " << ( *begin ) << std::endl;
	TR << "worst fragment:  " << ( *( end - 1) ) << std::endl;

	core::Real max_allowed_score( basic::options::option[ basic::options::OptionKeys::frags::picking_old_max_score ].value() );
	//flo debug
	//TR << "allscores: ";
	core::Size num_frags_picked(1);
	//for( Scores::const_iterator tmp_it( scores().begin() + ( from - 1 ) ); tmp_it != scores().begin() + to -1; ++tmp_it, ++hack ){
	// TR << "(" << hack << "/" << *tmp_it << "),";
	//}
	//TR << std::endl;
	//debug over

	for ( Scores::const_iterator i = begin; i != end; ++i, ++num_frags_picked ) {
		core::Real this_score(i->score);
		if ( this_score < max_allowed_score ) {
			fdl.push_back( extent_to_fragdata( i->extent_begin, i->extent_end, this_score, srfd_type ) );
		} else {
			if ( fdl.size() == 0 ) { //in case every fragment was too bad we only take the first to prevent crash
				TR << "no fragments had good score, only picking 1. this is probably a bad sign..." << std::endl;
				fdl.push_back( extent_to_fragdata( begin->extent_begin, begin->extent_end,this_score, srfd_type ) );
			} else TR << "fragments with bad score (" << this_score << ") observed after fragment " << num_frags_picked << ", skipping everything afterwards." << std::endl;
			break;
		}
	}

	return fdl;
}


/// @brief this function runs before main routine in catalog() starts
void VallLibrarian::pre_catalog_ops( VallLibrary const & library ) {
	// do any Evaluator pre-catalog operations
	for ( VallFragmentEvalOPs::iterator i = Super::extent_eval().begin(), ie = Super::extent_eval().end(); i != ie; ++i ) {
		(**i).pre_catalog_op( library );
	}

	// preallocate heuristic
	if ( preallocate_ ) {
		scores().reserve( library.size() );
	}
}


/// @brief this function runs after main routine catalog() finishes
void VallLibrarian::post_catalog_ops( VallLibrary const & library ) {
	// do any Evaluator post-catalog operations
	for ( VallFragmentEvalOPs::iterator i = Super::extent_eval().begin(), ie = Super::extent_eval().end(); i != ie; ++i ) {
		(**i).post_catalog_op( library );
	}
}


} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core
