// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/segment_file_generator.hh
/// @brief a generator of sewing segment files
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_apps_pilot_frankdt_segment_file_generator_hh
#define INCLUDED_apps_pilot_frankdt_segment_file_generator_hh

#include <apps/public/sewing/segment_file_generator.fwd.hh>
#include <string>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

namespace apps {
namespace pilot {
namespace frankdt {

///@brief a generator of smart sewing segment files
class segment_file_generator : public utility::pointer::ReferenceCount {

public:

	segment_file_generator();
	segment_file_generator(segment_file_generator const & src);

	virtual ~segment_file_generator();

	segment_file_generatorOP
	clone() const;

private:

};



struct SecondaryStruct{
	SecondaryStruct( std::string dssp, core::Size min, core::Size max):
		dssp_(dssp),
		min_(min),
		max_(max)
	{}
	SecondaryStruct(){}
	std::string dssp_;
	core::Size min_;
	core::Size max_;
};

struct Motif{
	Motif( std::string abv_motif_string, std::string motif_string, utility::vector1< SecondaryStruct > secondary_structs ):
		abv_motif_string_( abv_motif_string ),
		motif_string_( motif_string ),
		secondary_structs_( secondary_structs )
	{}
	Motif(){}
	std::string abv_motif_string_;
	std::string motif_string_;
	utility::vector1< SecondaryStruct > secondary_structs_;
};

} //apps
} //pilot
} //frankdt

#endif //INCLUDED_apps_pilot_frankdt_segment_file_generator_hh





