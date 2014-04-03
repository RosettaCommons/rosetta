// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/frag_picker/CSTalosIO.hh
/// @brief Class that converts Talos object into secondary shifts list
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_CS2ndShift_hh
#define INCLUDED_protocols_frag_picker_CS2ndShift_hh

// utility headers
#include <core/types.hh>
// AUTO-REMOVED #include <utility/vector0.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <protocols/frag_picker/CSTalosIO.hh>



#include <string>
#include <map>

// boost headers
#include <boost/tuple/tuple.hpp>

namespace protocols {
namespace frag_picker {

class CS2ndShift {
public:

	CS2ndShift() {}

	//CS2ndShift( CSTalosIO & input_data );

	CS2ndShift( CSTalosIO & input_data, bool use_sslimit );

	//

	inline utility::vector1< utility::vector1< std::pair< core::Size, core::Real > > > shifts() {
		return secondary_shifts_;
	}

private:

	std::map< char, std::map< std::string, core::Real > >
	read_adjust_table( std::string const & file_name );

	std::map<char,std::map<std::string,std::pair< core::Real, core::Real > > >
	read_sslimit_table(std::string const & file_name);

	std::string also_check_fix_disulf( std::string instring );

	utility::vector1< utility::vector1< std::pair< core::Size, core::Real > > > secondary_shifts_;

};

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_CSTalosIO_HH */
