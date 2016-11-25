// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DnaInterfaceMinMover.hh
/// @brief A mover that minimizes protein residues near DNA
/// @author ashworth

#ifndef INCLUDED_protocols_dna_DnaInterfaceMinMover_hh
#define INCLUDED_protocols_dna_DnaInterfaceMinMover_hh

#include <protocols/dna/DnaInterfaceMinMover.fwd.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <protocols/dna/DnaInterfaceFinder.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace dna {

class DnaInterfaceMinMover : public protocols::simple_moves::MinMover {
public:
	typedef utility::tag::TagCOP TagCOP;
public:
	DnaInterfaceMinMover();
	DnaInterfaceMinMover( DnaInterfaceMinMover const & );
	DnaInterfaceMinMover & operator = ( DnaInterfaceMinMover const & );
	~DnaInterfaceMinMover() override;
	DnaInterfaceMinMover( DnaInterfaceFinderOP );

	void apply( core::pose::Pose & ) override;
	// XRW TEMP  std::string get_name() const override;

	void use_interface( DnaInterfaceFinderOP );
	void chi( bool value ) { chi_ = value; }
	void bb( bool value ) { bb_ = value; }
	void reset_from_interface();

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & ) override;

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


private:
	DnaInterfaceFinderOP interface_;
	bool chi_;
	bool bb_;
};

} // namespace dna
} // namespace protocols

#endif
