// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/dna/DNAParameters.cc
/// @brief A class to query for base-paired partners as well as base pair and base step parameters
/// @author Jim Havranek

#include <protocols/dna/DNAParameters.hh>

#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/basic.hh>

#include <core/chemical/ResidueType.hh>

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>

//#include <ObjexxFCL/format/formatted.o.hh> // I()

#include <utility/vector1.hh>

#include <numeric/random/random.hh>


using utility::vector1;

namespace protocols {
namespace dna {

using namespace core;

// Note that the DNABasepair and DNABasestep classes are meant to managed by the DNAParameters
// class, and not used directly.  The main reason for this is that they have little in the way
// of error-checking - most of this is done by the DNAParameters class, which creates them as
// needed.  Even the DNAParameters class is not meant to be created by user code.  The intent
// is that it will be created and updated as necessary by the DNAParameterCalculator PoseMetric,
// which does all the instantiation and recalculation work and logic.

// This is the case where we pass in the residue
// and call out to a function in core/scoring/dna/base_geometry.hh
DNABase::DNABase( conformation::Residue const & rsd )
{
	assert( rsd.is_DNA() );

	alpha_ = rsd.mainchain_torsion(1);
	beta_ = basic::unsigned_periodic_range( rsd.mainchain_torsion(2), 360.0 );
	gamma_ = rsd.mainchain_torsion(3);
	delta_ = rsd.mainchain_torsion(4);
	epsilon_ = basic::unsigned_periodic_range( rsd.mainchain_torsion(5), 360.0 );
	zeta_ = basic::unsigned_periodic_range( rsd.mainchain_torsion(6), 360.0 );
	chi_ = rsd.chi(1);


	Real pseudorotation( 0.0 );
	Real amplitude( 0.0 );
	std::pair< std::string, int > pucker;

	//scoring::dna::get_base_pucker( rsd, pucker, pseudorotation, amplitude );

	pseudorotation_ = pseudorotation;
	amplitude_ = amplitude;
}

// This is the case where we pass in the parameters in a vector,
// with the ordering as defined in core/scoring/dna/base_geometry.hh
DNABasepair::DNABasepair( utility::vector1< core::Real > const & values ) :
	stretch_( values[ 5 ] ), stagger_( values[ 6 ] ), shear_( values[ 4 ] ),
	propeller_( values[ 1 ] ), opening_( values[ 3 ] ), buckle_( values[ 2 ] ) {}

// This is the case where we pass in the residues that form the base pair
// and call out to a function in core/scoring/dna/base_geometry.hh
DNABasepair::DNABasepair( conformation::Residue const & rsd1, conformation::Residue const & rsd2 )
{
	assert( rsd1.is_DNA() );
	assert( rsd2.is_DNA() );

	utility::vector1< core::Real > values( 6 );

	scoring::dna::get_base_pair_params( rsd1, rsd2, values );

	stretch_ = values[ 5 ];
	stagger_ = values[ 6 ];
	shear_ = values[ 4 ];
	propeller_ = values[ 1 ];
	opening_ = values[ 3 ];
	buckle_ = values[ 2 ];
}

// This is the case where we pass in the parameters in a vector,
// with the ordering as defined in core/scoring/dna/base_geometry.hh
DNABasestep::DNABasestep( utility::vector1< core::Real > const & values ) :
	slide_( values[ 4 ] ), shift_( values[ 6 ] ), rise_( values[ 5 ] ),
	roll_( values[ 2 ] ), twist_( values[ 1 ] ), tilt_( values[ 3 ] ) {}

// This is the case where we pass in the residues that form the base pair
// and call out to a function in core/scoring/dna/base_geometry.hh
DNABasestep::DNABasestep( conformation::Residue const & rsd1, conformation::Residue const & rsd2,
													conformation::Residue const & rsd1_next, conformation::Residue const & rsd2_prev )
{
	assert( rsd1.is_DNA() );
	assert( rsd2.is_DNA() );
	assert( rsd1_next.is_DNA() );
	assert( rsd2_prev.is_DNA() );

	utility::vector1< Real > values( 6, 0.0 );

	scoring::dna::get_base_step_params( rsd1, rsd2, rsd1_next, rsd2_prev, values );

	slide_ = values[ 4 ];
	shift_ = values[ 6 ];
	rise_ = values[ 5 ];
	roll_ = values[ 2 ];
	twist_ = values[ 1 ];
	tilt_ = values[ 3 ];
}


// End of the bs/bp data holders
//
//
// Begin DNAParameters class methods

// Return the torsions for this residue, if it's a base
DNABase const &
DNAParameters::base( core::Size resid ) const
{
	std::map< core::Size, DNABase >::const_iterator find_itr( bases_.find( resid ) );
	assert( find_itr != bases_.end() );

	return find_itr->second;
}

// Return the parameters for a base pair including this residue
DNABasepair const &
DNAParameters::basepair( core::Size resid ) const
{
	std::map< core::Size, DNABasepair >::const_iterator find_itr( basepairs_.find( resid ) );
	assert( find_itr != basepairs_.end() );

	return find_itr->second;
}

// Return the parameters for the base step consisting of this residue, the next residue,
// and their two base paired partners
DNABasestep const &
DNAParameters::basestep( core::Size resid ) const
{
	std::map< core::Size, DNABasestep >::const_iterator find_itr( basesteps_.find( resid ) );
	assert( find_itr != basesteps_.end() );

	return find_itr->second;
}

// Check to see if this residue is in base pair where "in a base pair" is defined as
// "in my map of base pairs"
bool
DNAParameters::is_base_paired( core::Size resid ) const
{
	return( basepairs_.find( resid) != basepairs_.end() );
}

// Check to see if this residue is a valid starting point for a base step - that is, do
// this residue, its base paired partner, the residue at resid+1, and its base paired partner
// form a base step
bool
DNAParameters::valid_basestep_start( core::Size resid ) const
{
	return( basesteps_.find( resid) != basesteps_.end() );
}

// Find the base paired partner for a residue.  If none exists, you get the constructor
// default value in the vector1, which is at the time of this typing 0.
core::Size
DNAParameters::find_partner( core::Size resid ) const
{
	return partners_[ resid ];
}

core::Size
DNAParameters::random_basepair() const
{
	return unique_basepairs_[ numeric::random::rg().random_range(1, unique_basepairs_.size() ) ];
}

core::Size
DNAParameters::random_basestep() const
{
	return unique_basestep_starts_[ numeric::random::rg().random_range(1, unique_basestep_starts_.size() ) ];
}

// More or less trying to mimic the functionality of Phil's code in scoring::dna with as
// little code duplication as necessary.  One difference is that I am double storing
// base pair information, once for each of the bases.  I think certain parameters change
// sign upon this kind of inversion, but I didn't want to leave any missing pieces in case
// they get queried by someone later.
void
DNAParameters::calculate( core::pose::Pose const & pose )
{
	basepairs_.clear();
	basesteps_.clear();
	bases_.clear();
	partners_.clear();
	dna_base_positions_.clear();
	unique_basepairs_.clear();
	unique_basestep_starts_.clear();

	// This function gets just what we want, including the desired value of zero
	// for unpartnered bases/non-DNA
	scoring::dna::find_basepairs( pose, partners_ );

	// Fill in the single base torsions
	for( core::Size resid( 1 ), endid( pose.total_residue() ) ; resid <= endid ; ++resid ) {
		if( pose.residue( resid ).is_DNA() ) {
			dna_base_positions_.push_back( resid );
			bases_[ resid ] =  DNABase( pose.residue( resid ) ) ;
		}
	}

	// Now get the base pairs
	for( core::Size resid( 1 ), endid( pose.total_residue() ) ; resid <= endid ; ++resid ) {
		core::Size partner_check( partners_[ resid ] );
		if( partner_check > resid ) {
			basepairs_[ resid ] =  DNABasepair( pose.residue( resid ), pose.residue( partner_check ) ) ;
			basepairs_[ partner_check] =  DNABasepair( pose.residue( partner_check ), pose.residue( resid ) ) ;
			unique_basepairs_.push_back( resid );
		}
	}

	// And the base steps - note the range of the loop is foreshortened, since a base step
	// can't start on the last residue
	for( core::Size resid( 1 ), endid( pose.total_residue() - 1 ) ; resid <= endid ; ++resid ) {

		core::Size partner_check( partners_[ resid ] );
		core::Size next_residue( resid + 1 );
		core::Size next_residue_partner( partners_[ next_residue ] );

		// To form a base step, both resid and next_residue must be base paired.
		// Note that I do one insertion max in the inner loop.  I let the same base step be
		// found twice, once in each direction.  This is different from the base pairs above.
		if( partner_check != 0 &&
				next_residue_partner != 0 &&
				!pose.residue( resid ).is_upper_terminus() ) {
			basesteps_[ resid ] =  DNABasestep( pose.residue( resid ), pose.residue( partner_check ),
																							pose.residue( next_residue ), pose.residue( next_residue_partner ) ) ;
			if( resid < partner_check ) {
				unique_basestep_starts_.push_back( resid );
			}
		}
	}


	return;
}

} // namespace dna
} // namespace protocols
