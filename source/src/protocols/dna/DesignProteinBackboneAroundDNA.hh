// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DesignProteinBackboneAroundDNA.hh
/// @brief performs a round flexible backbone sampling/design.  (Interim(?) encapsulation of some loose code in a mover class.)
/// @author ashworth

#ifndef INCLUDED_protocols_dna_DesignProteinBackboneAroundDNA_hh
#define INCLUDED_protocols_dna_DesignProteinBackboneAroundDNA_hh

#include <protocols/dna/DesignProteinBackboneAroundDNA.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/dna/DnaDesignDef.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.fwd.hh>

#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace dna {

// conceptually this not really a PackRotamersMover, but deriving this way makes sense on a practical level
class DesignProteinBackboneAroundDNA : public protocols::simple_moves::PackRotamersMover {
public:
	typedef loops::Loops Loops;

public:
	DesignProteinBackboneAroundDNA();
	~DesignProteinBackboneAroundDNA() override;

	DesignProteinBackboneAroundDNA( std::string , ScoreFunctionCOP );

	void targeted_dna( DnaDesignDefOPs const & );
	DnaDesignDefOPs const & targeted_dna() const;

	void apply( Pose & ) override;

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

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // methods
	void set_loop_info( Pose const &, Loops const & );
	void ccd( Pose &, loops::LoopsOP const, TaskFactoryCOP );
	void backrub( Pose &, loops::LoopsOP const, TaskFactoryCOP );

private: // data
	std::string type_;
	DnaDesignDefOPs targeted_dna_;
	core::Size gapspan_, spread_, cycles_outer_, cycles_inner_, repack_rate_;
	core::Real temp_initial_, temp_final_;
	bool designable_second_shell_;
};

} // namespace dna
} // namespace protocols

#endif
