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
/// @author gsmurphy
// Unit Headers
#include <devel/denovo_protein_design/SSClass.hh>
// Project Headers
// C++ Headers
#include <string>
#include <iostream>
#include <sstream>
// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
#include <numeric/random/random.fwd.hh>
#include <basic/Tracer.hh>

// ObjexxFCL Headers

namespace devel {
namespace denovo_protein_design {

using namespace core;

static THREAD_LOCAL basic::Tracer tr( "SS" );


std::ostream & operator<< ( std::ostream & os, const SSs & sss_ ) {
	os << "SS  begin  end  type" << std::endl;
	for ( auto const & ss : sss_ ) {
		os << ss << std::endl;
	}
	return os;
}


/////////////////////////////////////////////////////////////////////
void
SSs::write_ss_to_file(
	std::string const & filename
)
{

	utility::io::ozstream data;
	data.open( filename );
	if ( !data ) {
		utility_exit_with_message( "Couldn't write check point ss file " );
	}

	for ( auto const & it : *this ) {
		data << "SS " << it.start() << " " << it.stop() << " " << it.sstype() << std::endl;
	}

	data.close();
	data.clear();
}


void
SSs::add_ss( SS ss_ ) {
	ss_.start();
	ss_.stop();
	ss_.sstype();
	sss_.push_back( ss_ );
}


//////////////////////////////////////////////////////////////////////
void
SSs::add_ss(
	core::Size const start,
	core::Size const stop,
	char const sstype
)
{
	add_ss( SS( start, stop, sstype ));
}
/////////////////////////////////////////////////////////////////////////////
void
SSs::add_ss( const SSs::const_iterator & it ) {
	add_ss( it->start(), it->stop(), it->sstype() );
}
/////////////////////////////////////////////////////////////////////////////
void
SSs::add_ss( const SSs::iterator & it ) {
	add_ss( it->start(), it->stop(), it->sstype() );
}

//////////////////////////////////////////////////////////////////////////////
void
SSs::delete_ss(
	Size const start,
	Size const stop
)
{
	runtime_assert( start < stop );

	for ( auto it=sss_.begin(), it_end=sss_.end();
			it != it_end; ++it ) {
		if ( start == it->start() && stop == it->stop() ) {
			sss_.erase( it );
			break;
		}
	}
}
/////////////////////////////////////////////////////////////////////////////
SSs::const_iterator
SSs::one_random_ss_element() const {
	Size const size = sss_.size();
	runtime_assert( size > 0 );
	Size index =0;
	Size const end = static_cast< Size >( numeric::random::uniform()*size );
	auto it = sss_.begin();
	while ( index != end ) { ++index; ++it; }
	return it;

}
/////////////////////////////////////////////////////////////////////////////
/*
Size
SSs::ss_size(
Size const num
) const {
runtime_assert( num > 0 && num <= sss_.size() );
return sss_[num-1].size();
}
/////////////////////////////////////////////////////////////////////////////
Size
SSs::ss_size() const {
Size sssize = 0;
for( const_iterator it=sss_.begin(), it_end=sss_.end();
it != it_end; ++it ) {
sssize += it->size();
}
return sssize;
}
*/

void
SSs::clear(){
	sss_.clear();
}

/*
void SSs::read_ss_file(
std::string filename
) {
clear();
std::ifstream infile( filename.c_str() );

if (!infile.good()) {
utility_exit_with_message( "[ERROR] Error opening SS file '" + filename + "'" );
}

}
*/

} // namespace denovo_protein_design
} // namespace devel

