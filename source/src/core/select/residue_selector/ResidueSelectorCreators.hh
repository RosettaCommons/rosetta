// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueSelectorCreators.hh
/// @brief  Class declarations for the ResidueSelectorCreators for a set of simple ResidueSelectors
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_ResidueSelectorCreators_HH
#define INCLUDED_core_select_residue_selector_ResidueSelectorCreators_HH

// Package headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace select {
namespace residue_selector {

class AndResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class BinSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class BondedResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class ChainSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class DensityFitResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class InterGroupInterfaceByVectorSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class GlycanPositionSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class GlycanResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class GlycanSequonsSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class LayerSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class NotResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class ResidueIndexSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class ResidueSpanSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class ResidueNameSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class ResidueInMembraneSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};


class NeighborhoodResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class NumNeighborsSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class OrResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class PhiSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class JumpUpstreamSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class JumpDownstreamSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class SecondaryStructureSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class SSElementSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class TrueResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class FalseResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class ResiduePDBInfoHasLabelSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class RandomResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class RandomGlycanFoliageSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class ScoreTermValueBasedSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class ResidueInSequenceMotifSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

class SliceResidueSelectorCreator : public ResidueSelectorCreator {
public:
	ResidueSelectorOP create_residue_selector() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const override;
};

} //namespace residue_selector
} //namespace select
} //namespace core


#endif
