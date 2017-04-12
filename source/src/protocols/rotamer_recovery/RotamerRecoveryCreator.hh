// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/RotamerRecoveryCreator.hh
/// @brief  Base class for RotamerRecoveryCreators for the RotamerRecovery load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_rotamer_recovery_RotamerRecoveryCreator_hh
#define INCLUDED_protocols_rotamer_recovery_RotamerRecoveryCreator_hh

// Unit Headers
#include <protocols/rotamer_recovery/RotamerRecoveryCreator.fwd.hh>

// Package Headers
#include <protocols/rotamer_recovery/RRProtocol.fwd.hh>
#include <protocols/rotamer_recovery/RRComparer.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>

namespace protocols {
namespace rotamer_recovery {

/// @brief The Creator class is responsible for creating a particular
/// mover class.
class RRProtocolCreator : public utility::pointer::ReferenceCount {
public:
	RRProtocolCreator() {}
	~RRProtocolCreator() override = default;

	virtual RRProtocolOP create_protocol() const = 0;
	virtual std::string type_name() const = 0;
	virtual void append_attributes( utility::tag::AttributeList & attlist ) const = 0;
};

class RRProtocolReferenceStructureCreator : public RRProtocolCreator {
public:
	RRProtocolReferenceStructureCreator() {}
	~RRProtocolReferenceStructureCreator() override = default;

	RRProtocolOP create_protocol() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRProtocolRTMinCreator : public RRProtocolCreator {
public:
	RRProtocolRTMinCreator() {}
	~RRProtocolRTMinCreator() override = default;

	RRProtocolOP create_protocol() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRProtocolRotamerTrialsCreator : public RRProtocolCreator {
public:
	RRProtocolRotamerTrialsCreator() {}
	~RRProtocolRotamerTrialsCreator() override = default;

	RRProtocolOP create_protocol() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRProtocolMinPackCreator : public RRProtocolCreator {
public:
	RRProtocolMinPackCreator() {}
	~RRProtocolMinPackCreator() override = default;

	RRProtocolOP create_protocol() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRProtocolPackRotamersCreator : public RRProtocolCreator {
public:
	RRProtocolPackRotamersCreator() {}
	~RRProtocolPackRotamersCreator() override = default;

	RRProtocolOP create_protocol() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRProtocolRelaxCreator : public RRProtocolCreator {
public:
	RRProtocolRelaxCreator() {}
	~RRProtocolRelaxCreator() override = default;

	RRProtocolOP create_protocol() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRProtocolMoverCreator : public RRProtocolCreator {
public:
	RRProtocolMoverCreator() {}
	~RRProtocolMoverCreator() override = default;

	RRProtocolOP create_protocol() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};


/// @brief The Creator class is responsible for creating a particular
/// mover class.
class RRComparerCreator : public utility::pointer::ReferenceCount {
public:
	RRComparerCreator() {}
	~RRComparerCreator() override = default;

	virtual RRComparerOP create_comparer() const = 0;
	virtual std::string type_name() const = 0;
	virtual void append_attributes( utility::tag::AttributeList & attlist ) const = 0;
};

class RRComparerAutomorphicRMSDCreator : public RRComparerCreator {
public:
	RRComparerAutomorphicRMSDCreator() {}
	~RRComparerAutomorphicRMSDCreator() override = default;

	RRComparerOP create_comparer() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRComparerRotBinsCreator : public RRComparerCreator {
public:
	RRComparerRotBinsCreator() {}
	~RRComparerRotBinsCreator() override = default;

	RRComparerOP create_comparer() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRComparerChiDiffCreator : public RRComparerCreator {
public:
	RRComparerChiDiffCreator() {}
	~RRComparerChiDiffCreator() override = default;

	RRComparerOP create_comparer() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRComparerElecDensDiffCreator : public RRComparerCreator {
public:
	RRComparerElecDensDiffCreator() {}
	~RRComparerElecDensDiffCreator() override = default;

	RRComparerOP create_comparer() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

/// @brief The Creator class is responsible for creating a particular
/// mover class.
class RRReporterCreator : public utility::pointer::ReferenceCount {
public:
	RRReporterCreator() {}
	~RRReporterCreator() override = default;

	virtual RRReporterOP create_reporter() const = 0;
	virtual std::string type_name() const = 0;
	virtual void append_attributes( utility::tag::AttributeList & attlist ) const = 0;
};

class RRReporterSimpleCreator : public RRReporterCreator {
public:
	RRReporterSimpleCreator() {}
	~RRReporterSimpleCreator() override = default;

	RRReporterOP create_reporter() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRReporterHumanCreator : public RRReporterCreator {
public:
	RRReporterHumanCreator() {}
	~RRReporterHumanCreator() override = default;

	RRReporterOP create_reporter() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};

class RRReporterSQLiteCreator : public RRReporterCreator {
public:
	RRReporterSQLiteCreator() {}
	~RRReporterSQLiteCreator() override = default;

	RRReporterOP create_reporter() const override;
	std::string type_name() const override;
	void append_attributes( utility::tag::AttributeList & attlist ) const override;
};


} //namespace
} //namespace

#endif
