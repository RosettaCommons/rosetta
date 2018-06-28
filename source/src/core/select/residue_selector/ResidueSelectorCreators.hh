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
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class BinSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class BondedResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class ChainSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class DensityFitResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class InterGroupInterfaceByVectorSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class GlycanPositionSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class GlycanResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class GlycanSequonsSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class LayerSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class NotResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class ResidueIndexSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class ResidueSpanSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class ResidueNameSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class ResidueInMembraneSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};


class NeighborhoodResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class NumNeighborsSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class OrResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class PhiSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class JumpUpstreamSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class JumpDownstreamSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class SecondaryStructureSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class SSElementSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class TrueResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class ResiduePDBInfoHasLabelSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class RandomResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class RandomGlycanFoliageSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class ScoreTermValueBasedSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

class ResidueInSequenceMotifSelectorCreator : public ResidueSelectorCreator {
public:
	virtual ResidueSelectorOP create_residue_selector() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};


} //namespace residue_selector
} //namespace select
} //namespace core


#endif
