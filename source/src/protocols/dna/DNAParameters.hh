// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/dna/DNAParameters.hh
/// @brief A class to query for base-paired partners as well as base pair and base step parameters
/// @author Jim Havranek

#ifndef INCLUDED_protocols_dna_DNAParameters
#define INCLUDED_protocols_dna_DNAParameters

#include <protocols/dna/DNAParameters.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <iosfwd>
#include <map>

namespace protocols {
namespace dna {

class DNABase {
public: // constructor
	DNABase() :  alpha_( 0.0 ), beta_( 0.0 ), gamma_( 0.0 ), delta_( 0.0 ), epsilon_( 0.0 ), zeta_( 0.0 ),
		chi_( 0.0 ), pseudorotation_( 0.0 ), amplitude_( 0.0 ) {}
	// These next two use the formats from core::scoring::dna, which are ultimately used to generate the data
	DNABase( core::conformation::Residue const & rsd );

public:
	// const getter methods
	core::Real alpha() const { return alpha_; }
	core::Real beta() const { return beta_; }
	core::Real gamma() const { return gamma_; }
	core::Real delta() const { return delta_; }
	core::Real epsilon() const { return epsilon_; }
	core::Real zeta() const { return zeta_; }
	core::Real chi() const { return chi_; }
	core::Real pseudorotation() const { return pseudorotation_; }
	core::Real amplitude() const { return amplitude_; }
	// setter methods
	void alpha( core::Real in_value ) { alpha_ = in_value; }
	void beta( core::Real in_value ) { beta_ = in_value; }
	void gamma( core::Real in_value ) { gamma_ = in_value; }
	void delta( core::Real in_value ) { delta_ = in_value; }
	void epsilon( core::Real in_value ) { epsilon_ = in_value; }
	void zeta( core::Real in_value ) { zeta_ = in_value; }
	void chi( core::Real in_value ) { chi_ = in_value; }
	void pseudorotation( core::Real in_value ) { pseudorotation_ = in_value; }
	void amplitude( core::Real in_value ) { amplitude_ = in_value; }

private:
	core::Real alpha_;
	core::Real beta_;
	core::Real gamma_;
	core::Real delta_;
	core::Real epsilon_;
	core::Real zeta_;
	core::Real chi_;
	core::Real pseudorotation_;
	core::Real amplitude_;
};

class DNABasepair {
public: // constructor
	DNABasepair() :  stretch_( 0.0 ), stagger_( 0.0 ), shear_( 0.0 ), propeller_( 0.0 ), opening_( 0.0 ), buckle_( 0.0 ) {}
	// These next two use the formats from core::scoring::dna, which are ultimately used to generate the data
	DNABasepair( utility::vector1<core::Real> const & init_values );
	DNABasepair( core::conformation::Residue const & rsd1, core::conformation::Residue const & rsd2 );

public:
	// const getter methods
	core::Real stretch() const { return stretch_; }
	core::Real stagger() const { return stagger_; }
	core::Real shear() const { return shear_; }
	core::Real propeller() const { return propeller_; }
	core::Real opening() const { return opening_; }
	core::Real buckle() const { return buckle_; }
	// setter methods
	void stretch( core::Real in_value ) { stretch_ = in_value; }
	void stagger( core::Real in_value ) { stagger_ = in_value; }
	void shear( core::Real in_value ) { shear_ = in_value; }
	void propeller( core::Real in_value ) { propeller_ = in_value; }
	void opening( core::Real in_value ) { opening_ = in_value; }
	void buckle( core::Real in_value ) { buckle_ = in_value; }

private:
	core::Real stretch_;
	core::Real stagger_;
	core::Real shear_;
	core::Real propeller_;
	core::Real opening_;
	core::Real buckle_;
};


class DNABasestep {
public: // constructor
	DNABasestep() :  slide_( 0.0 ), shift_( 0.0 ), rise_( 0.0 ), roll_( 0.0 ), twist_( 0.0 ), tilt_( 0.0 ) {}
	// These next two use the formats from core::scoring::dna, which are ultimately used to generate the data
	DNABasestep( utility::vector1<core::Real> const & init_values );
	DNABasestep( core::conformation::Residue const & rsd1, core::conformation::Residue const & rsd2,
		core::conformation::Residue const & rsd1_next, core::conformation::Residue const & rsd2_prev );

public:
	// const getter methods
	core::Real slide() const { return slide_; }
	core::Real shift() const { return shift_; }
	core::Real rise() const { return rise_; }
	core::Real roll() const { return roll_; }
	core::Real twist() const { return twist_; }
	core::Real tilt() const { return tilt_; }
	// setter methods
	void slide( core::Real in_value ) { slide_ = in_value; }
	void shift( core::Real in_value ) { shift_ = in_value; }
	void rise( core::Real in_value ) { rise_ = in_value; }
	void roll( core::Real in_value ) { roll_ = in_value; }
	void twist( core::Real in_value ) { twist_ = in_value; }
	void tilt( core::Real in_value ) { tilt_ = in_value; }

private:
	core::Real slide_;
	core::Real shift_;
	core::Real rise_;
	core::Real roll_;
	core::Real twist_;
	core::Real tilt_;
};


class DNAParameters : public utility::pointer::ReferenceCount {
public: // constructors
	DNAParameters(){};
	DNAParameters( core::pose::Pose const & pose ) { calculate( pose ); }
	~DNAParameters(){};

	core::Size number_of_bases() { return bases_.size(); }
	core::Size number_of_basepairs() { return basepairs_.size(); }
	core::Size number_of_basesteps() { return basesteps_.size(); }
	core::Size random_basepair() const;
	core::Size random_basestep() const;

	core::Size  const & nth_dna_base( core::Size index ) const { return dna_base_positions_[ index ]; }
	DNABase     const & base( core::Size resid ) const;
	DNABasepair const & basepair( core::Size resid ) const;
	DNABasestep const & basestep( core::Size resid ) const;

	std::map< core::Size, DNABase >::const_iterator bases_begin() const { return bases_.begin(); }
	std::map< core::Size, DNABase >::const_iterator bases_end() const { return bases_.end(); }
	std::map< core::Size, DNABasepair >::const_iterator basepairs_begin() const { return basepairs_.begin(); }
	std::map< core::Size, DNABasepair >::const_iterator basepairs_end() const { return basepairs_.end(); }
	std::map< core::Size, DNABasestep >::const_iterator basesteps_begin() const { return basesteps_.begin(); }
	std::map< core::Size, DNABasestep >::const_iterator basesteps_end() const { return basesteps_.end(); }

	core::Size find_partner( core::Size resid ) const;
	bool is_base_paired( core::Size resid ) const;
	bool valid_basestep_start( core::Size resid ) const;

	void calculate( core::pose::Pose const & pose );

private:
	// vectors to hold the actual data
	std::map< core::Size, DNABasepair > basepairs_;
	std::map< core::Size, DNABasestep > basesteps_;
	std::map< core::Size, DNABase > bases_;
	utility::vector1< core::Size > partners_;
	utility::vector1< core::Size > unique_basepairs_;
	utility::vector1< core::Size > unique_basestep_starts_;
	utility::vector1< core::Size > dna_base_positions_;

};

} // namespace dna
} // namespace protocols

#endif
