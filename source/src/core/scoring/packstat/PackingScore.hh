// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/PackingScore.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_PackingScore_hh
#define INCLUDED_core_scoring_packstat_PackingScore_hh


#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1_bool.hh>


namespace core {
namespace scoring {
namespace packstat {


struct PackingScoreResData : public utility::pointer::ReferenceCount {
	PackingScoreResData( Size nrad, Size npr ) : nrad_(nrad), npr_(npr), msa_(nrad*npr,0.0) {}
	Size  npr() const { return npr_; }
	Size nrad() const { return nrad_; }
	core::Real       & msa( Size rad, Size pr )       { return msa_[ (rad-1)*npr_ + pr ]; }
	core::Real const & msa( Size rad, Size pr ) const { return msa_[ (rad-1)*npr_ + pr ]; }
private:
	Size const nrad_,npr_;
	utility::vector1<core::Real> msa_;
};
std::ostream & operator<< ( std::ostream & out, PackingScoreResData const & dat );
typedef utility::pointer::shared_ptr< PackingScoreResData >       PackingScoreResDataOP;
typedef utility::pointer::shared_ptr< PackingScoreResData const > PackingScoreResDataCOP;

struct PackingScore : public utility::pointer::ReferenceCount {
	PackingScore( Size nrad, Size npr, bool /*compprob = false*/ ) :
		nrad_(nrad), npr_(npr), weights_(nrad*npr,0.0), centers_(nrad*npr,0.0), compprob_( false /*compprob_*/ ) {}
	Size  npr() const { return npr_; }
	Size nrad() const { return nrad_; }
	core::Real       & weight( Size rad, Size pr )       { return weights_[ (rad-1)*npr_ + pr ]; }
	core::Real const & weight( Size rad, Size pr ) const { return weights_[ (rad-1)*npr_ + pr ]; }
	core::Real       & center( Size rad, Size pr )       { return centers_[ (rad-1)*npr_ + pr ]; }
	core::Real const & center( Size rad, Size pr ) const { return centers_[ (rad-1)*npr_ + pr ]; }
	core::Real       &   rho( )       { return rho_; }
	core::Real const &   rho( ) const { return rho_; }
	core::Real       & probA( )       { return probA_; }
	core::Real const & probA( ) const { return probA_; }
	core::Real       & probB( )       { return probB_; }
	core::Real const & probB( ) const { return probB_; }
	bool       & compprob( )       { return compprob_; }
	bool const & compprob( ) const { return compprob_; }
	core::Real score( PackingScoreResDataCOP dat ) const;
	core::Real score( utility::vector1<PackingScoreResDataCOP> dats ) const;
private:
	Size const nrad_,npr_;
	utility::vector1<core::Real> weights_, centers_;
	core::Real rho_,probA_,probB_;
	bool compprob_;
};

typedef utility::pointer::shared_ptr< PackingScore >       PackingScoreOP;
typedef utility::pointer::shared_ptr< PackingScore const > PackingScoreCOP;


} // namespace packstat
} // namespace scoring
} // namespace core


#endif
