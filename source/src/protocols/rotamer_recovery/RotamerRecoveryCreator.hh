// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#include <string>

namespace protocols {
namespace rotamer_recovery {

/// @brief The Creator class is responsible for creating a particular
/// mover class.
class RRProtocolCreator : public utility::pointer::ReferenceCount {
public:
	RRProtocolCreator() {}
	virtual ~RRProtocolCreator() {}

	virtual RRProtocolOP create_protocol() const = 0;
	virtual std::string type_name() const = 0;
};

class RRProtocolReferenceStructureCreator : public RRProtocolCreator {
public:
	RRProtocolReferenceStructureCreator() {}
	~RRProtocolReferenceStructureCreator() {}

	RRProtocolOP create_protocol() const;
	std::string type_name() const;
};

class RRProtocolRTMinCreator : public RRProtocolCreator {
public:
	RRProtocolRTMinCreator() {}
	~RRProtocolRTMinCreator() {}

	RRProtocolOP create_protocol() const;
	std::string type_name() const;
};

class RRProtocolRotamerTrialsCreator : public RRProtocolCreator {
public:
	RRProtocolRotamerTrialsCreator() {}
	~RRProtocolRotamerTrialsCreator() {}

	RRProtocolOP create_protocol() const;
	std::string type_name() const;
};

class RRProtocolMinPackCreator : public RRProtocolCreator {
public:
	RRProtocolMinPackCreator() {}
	~RRProtocolMinPackCreator() {}

	RRProtocolOP create_protocol() const;
	std::string type_name() const;
};

class RRProtocolPackRotamersCreator : public RRProtocolCreator {
public:
	RRProtocolPackRotamersCreator() {}
	~RRProtocolPackRotamersCreator() {}

	RRProtocolOP create_protocol() const;
	std::string type_name() const;
};

class RRProtocolRelaxCreator : public RRProtocolCreator {
public:
	RRProtocolRelaxCreator() {}
	~RRProtocolRelaxCreator() {}

	RRProtocolOP create_protocol() const;
	std::string type_name() const;
};

class RRProtocolMoverCreator : public RRProtocolCreator {
public:
	RRProtocolMoverCreator() {}
	~RRProtocolMoverCreator() {}

	RRProtocolOP create_protocol() const;
	std::string type_name() const;
};


/// @brief The Creator class is responsible for creating a particular
/// mover class.
class RRComparerCreator : public utility::pointer::ReferenceCount {
public:
	RRComparerCreator() {}
	virtual ~RRComparerCreator() {}

	virtual RRComparerOP create_comparer() const = 0;
	virtual std::string type_name() const = 0;
};

class RRComparerAutomorphicRMSDCreator : public RRComparerCreator {
public:
	RRComparerAutomorphicRMSDCreator() {}
	~RRComparerAutomorphicRMSDCreator() {}

	RRComparerOP create_comparer() const;
	std::string type_name() const;
};

class RRComparerRotBinsCreator : public RRComparerCreator {
public:
	RRComparerRotBinsCreator() {}
	~RRComparerRotBinsCreator() {}

	RRComparerOP create_comparer() const;
	std::string type_name() const;
};

class RRComparerChiDiffCreator : public RRComparerCreator {
public:
	RRComparerChiDiffCreator() {}
	~RRComparerChiDiffCreator() {}

	RRComparerOP create_comparer() const;
	std::string type_name() const;
};

class RRComparerElecDensDiffCreator : public RRComparerCreator {
public:
	RRComparerElecDensDiffCreator() {}
	~RRComparerElecDensDiffCreator() {}

	RRComparerOP create_comparer() const;
	std::string type_name() const;
};

/// @brief The Creator class is responsible for creating a particular
/// mover class.
class RRReporterCreator : public utility::pointer::ReferenceCount {
public:
	RRReporterCreator() {}
	virtual ~RRReporterCreator() {}

	virtual RRReporterOP create_reporter() const = 0;
	virtual std::string type_name() const = 0;
};

class RRReporterSimpleCreator : public RRReporterCreator {
public:
	RRReporterSimpleCreator() {}
	~RRReporterSimpleCreator() {}

	RRReporterOP create_reporter() const;
	std::string type_name() const;
};

class RRReporterHumanCreator : public RRReporterCreator {
public:
	RRReporterHumanCreator() {}
	~RRReporterHumanCreator() {}

	RRReporterOP create_reporter() const;
	std::string type_name() const;
};

class RRReporterSQLiteCreator : public RRReporterCreator {
public:
	RRReporterSQLiteCreator() {}
	~RRReporterSQLiteCreator() {}

	RRReporterOP create_reporter() const;
	std::string type_name() const;
};


} //namespace
} //namespace

#endif
