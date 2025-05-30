// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ligand_docking/HeavyAtomFilter.hh
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_ligand_docking_HeavyAtomFilter_hh
#define INCLUDED_protocols_ligand_docking_HeavyAtomFilter_hh


#include <core/types.hh>
#include <protocols/filters/Filter.hh>

#include <utility/tag/Tag.fwd.hh>



namespace protocols {
namespace ligand_docking {

class HeavyAtomFilter : public protocols::filters::Filter
{
public:
	HeavyAtomFilter() :
		//utility::VirtualBase(),
		protocols::filters::Filter( "HeavyAtom" )
	{}

	HeavyAtomFilter(std::string const & chain, core::Size heavy_atom_limit ) :
		//utility::VirtualBase(),
		protocols::filters::Filter( "HeavyAtom" ),
		chain_(chain),
		heavy_atom_limit_(heavy_atom_limit)
	{}

	HeavyAtomFilter( HeavyAtomFilter const & init ) :
		//utility::VirtualBase(),
		protocols::filters::Filter( init ),
		chain_(init.chain_),
		heavy_atom_limit_(init.heavy_atom_limit_)

	{};

	~HeavyAtomFilter() override= default;

	bool apply( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override {
		return utility::pointer::make_shared< HeavyAtomFilter >( *this );
	}

	core::Real report_sm( core::pose::Pose const & ) const override;

	protocols::filters::FilterOP fresh_instance() const override{
		return utility::pointer::make_shared< HeavyAtomFilter >();
	}

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string chain_;
	core::Size heavy_atom_limit_;
};
} // ligand_docking
} // protocols

#endif //INCLUDED_protocols_ProteinInterfaceDesign_ligand_docking_HeavyAtomFilter_HH_
