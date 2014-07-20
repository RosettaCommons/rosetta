// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSamplerOneValueComb.hh
/// @brief Aggregate of multiple rotamer samplers for sampling combinatorially.
/// @author  Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerSamplerOneValueComb_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerSamplerOneValueComb_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerSamplerOneValueComb.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerSamplerOneValue.hh>
#include <protocols/rotamer_sampler/RotamerSamplerSizedComb.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerSamplerOneValueComb : public RotamerSamplerSizedComb {

public:

	RotamerSamplerOneValueComb();

	virtual ~RotamerSamplerOneValueComb();

	ValueList const &
	get_value_list( Size const id );

	ValueList const &
	get_value_list( utility::vector1< Size > const & id_list );

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_external_loop_rotamer( RotamerSamplerOneValueOP const & rotamer ) {
		RotamerSamplerSizedComb::add_external_loop_rotamer( rotamer );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerSamplerOneValueComb"; }

	/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
	virtual RotamerSamplerType type() const { return ONE_VALUE_COMB; }

private:

	/// @brief Add one more rotamer sampler to this sampler -- note that this is now made 'private', so that users can only add RotamerSamplerOneValueOPs
	using	RotamerSamplerSizedComb::add_external_loop_rotamer;

private:

	ValueList value_list_cached_, id_list_cached_;

};

} //rotamer_sampler
} //protocols

#endif

