// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rotamers/PeptoidRotamerLibrarySpecification.cc
/// @brief  The PeptoidRotamerLibrarySpecification class specifies building PeptoidRotamers.
/// @author Rocco Moretti (rmorettase@gmail.com)

// Unit headers
#include <core/chemical/rotamers/PeptoidRotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/PeptoidRotamerLibrarySpecificationCreator.hh>

#include <core/chemical/ResidueType.fwd.hh>

// Utility headers
#include <utility>
#include <utility/exit.hh>
#include <basic/Tracer.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {
namespace rotamers {

static basic::Tracer TR("core.chemical.rotamers.PeptoidRotamerLibrarySpecification");

// Creator Functions

RotamerLibrarySpecificationOP
PeptoidRotamerLibrarySpecificationCreator::create() const {
	return RotamerLibrarySpecificationOP( new PeptoidRotamerLibrarySpecification );
}

RotamerLibrarySpecificationOP
PeptoidRotamerLibrarySpecificationCreator::create( std::istream & input ) const {
	return RotamerLibrarySpecificationOP( new PeptoidRotamerLibrarySpecification( input ) );
}

std::string
PeptoidRotamerLibrarySpecificationCreator::keyname() const {
	return PeptoidRotamerLibrarySpecification::library_name();
}

// Specification Functions

PeptoidRotamerLibrarySpecification::PeptoidRotamerLibrarySpecification() = default;

PeptoidRotamerLibrarySpecification::PeptoidRotamerLibrarySpecification(std::string const & peptoid_rotlib_path ) :
	peptoid_rotlib_path_( peptoid_rotlib_path )
{}

PeptoidRotamerLibrarySpecification::PeptoidRotamerLibrarySpecification( std::istream & input )
{
	input >> peptoid_rotlib_path_;
	if ( ! input ) {
		utility_exit_with_message("Must provide rotamer library path for Peptoid rotamer library specification.");
	}
	core::Size nbins;
	input >> nbins;
	while ( input ) {
		peptoid_rotlib_n_bins_per_rot_.push_back( nbins );
		input >> nbins;
	}
	TR.Debug << "Read " << peptoid_rotlib_n_bins_per_rot_.size() << " entries for Peptoid rotamer library specification with path " << peptoid_rotlib_path_ << std::endl;
}

PeptoidRotamerLibrarySpecification::~PeptoidRotamerLibrarySpecification() = default;

std::string
PeptoidRotamerLibrarySpecification::keyname() const {
	return library_name();
}

std::string
PeptoidRotamerLibrarySpecification::library_name() {
	return "PEPTOID";
}

std::string
PeptoidRotamerLibrarySpecification::cache_tag(core::chemical::ResidueType const &) const {
	std::ostringstream ss;
	debug_assert( peptoid_rotlib_path_.size() );
	ss << peptoid_rotlib_path_;
	for ( core::Size ii(1); ii <= peptoid_rotlib_n_bins_per_rot_.size(); ++ii ) {
		ss << "|" << peptoid_rotlib_n_bins_per_rot_[ii];
	}
	return ss.str();
}

} //namespace rotamers
} //namespace chemical
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::rotamers::PeptoidRotamerLibrarySpecification::save( Archive & arc ) const {
	arc( cereal::base_class< core::chemical::rotamers::RotamerLibrarySpecification >( this ) );
	arc( CEREAL_NVP( peptoid_rotlib_path_ ) ); // std::string
	arc( CEREAL_NVP( peptoid_rotlib_n_bins_per_rot_ ) ); // utility::vector1<Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::rotamers::PeptoidRotamerLibrarySpecification::load( Archive & arc ) {
	arc( cereal::base_class< core::chemical::rotamers::RotamerLibrarySpecification >( this ) );
	arc( peptoid_rotlib_path_ ); // std::string
	arc( peptoid_rotlib_n_bins_per_rot_ ); // utility::vector1<Size>
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::rotamers::PeptoidRotamerLibrarySpecification );
CEREAL_REGISTER_TYPE( core::chemical::rotamers::PeptoidRotamerLibrarySpecification )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_rotamers_PeptoidRotamerLibrarySpecification )
#endif // SERIALIZATION
