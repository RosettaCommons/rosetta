// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief very simple--pulls apart non-DNA (i.e. protein and friends) from DNA to facilitate DNA binding score calculation

/// @details sets an *appropriate fold tree* for a non-DNA/DNA interface (including multichain non-DNA and/or DNA) and then simply separates non-DNA chain(s) from the DNA chain(s) with a single translation. More basic than RigidBodyMover (which maybe needs a new base class?)

/// @author ashworth

#ifndef INCLUDED_protocols_dna_SeparateDnaFromNonDna_hh
#define INCLUDED_protocols_dna_SeparateDnaFromNonDna_hh

#include <protocols/dna/SeparateDnaFromNonDna.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>

#ifdef WIN32
#include <utility/tag/Tag.hh>
#endif


namespace protocols {
namespace dna {

class SeparateDnaFromNonDna : public moves::Mover {
public:
	SeparateDnaFromNonDna();
	SeparateDnaFromNonDna( core::Real x, core::Real y, core::Real z );
	SeparateDnaFromNonDna( numeric::xyzVector< core::Real > const & xyz );
	SeparateDnaFromNonDna( SeparateDnaFromNonDna const & other );

	~SeparateDnaFromNonDna() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
	numeric::xyzVector< core::Real > translation() const { return translation_; }

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & ) override;
	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP fresh_instance() const override;
	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP clone() const override;

private:
	numeric::xyzVector< core::Real > translation_;

};

} // dna
} // protocols

#endif
