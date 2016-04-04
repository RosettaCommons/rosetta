// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	protocols/toolbox/PyReturnValuePolicyTest.hh
/// @brief	A few functions test how PyRosetta handle boost ReturnValuePolicy
/// @author	Sergey Lyskov

#ifndef INCLUDED_protocols_toolbox_PyReturnValuePolicyTest_hh
#define INCLUDED_protocols_toolbox_PyReturnValuePolicyTest_hh

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <utility/pointer/owning_ptr.fwd.hh>

#include <core/types.hh>

namespace protocols {
namespace toolbox {

class DummyClass;

typedef utility::pointer::shared_ptr< DummyClass > DummyClassOP;
typedef utility::pointer::shared_ptr< DummyClass const > DummyClassCOP;

//typedef utility::pointer::access_ptr< DummyClass > DummyClassAP;
//typedef utility::pointer::access_ptr< DummyClass const > DummyClassCAP;


class DummyClass : public utility::pointer::ReferenceCount {
public:
	DummyClass() {};
	DummyClass(int) {};
	~DummyClass() {};

	DummyClass(DummyClass const &) {};

	DummyClass & operator=( DummyClass const & ) { return *this; };

	static DummyClass * create() { return new DummyClass(); };

	virtual DummyClassOP clone() const { return DummyClassOP(new DummyClass()); };

	void test() {};

private:
};

inline DummyClassOP PyReturnValuePolicyTest_DummyClassOP(void) { return DummyClassOP( DummyClass::create() ); };
inline DummyClassCOP PyReturnValuePolicyTest_DummyClassCOP(void) { return DummyClassCOP( DummyClass::create() ); };
//inline DummyClassAP PyReturnValuePolicyTest_DummyClassAP(void)  { return DummyClass::create(); };
//inline DummyClassCAP PyReturnValuePolicyTest_DummyClassCAP(void)  { return DummyClass::create(); };


class SF_Replica;

typedef utility::pointer::shared_ptr< SF_Replica > SF_ReplicaOP;
typedef utility::pointer::shared_ptr< SF_Replica const > SF_ReplicaCOP;

//typedef utility::pointer::access_ptr< SF_Replica > SF_ReplicaAP;
//typedef utility::pointer::access_ptr< SF_Replica const > SF_ReplicaCAP;


// SF_Replica
class SF_Replica : public utility::pointer::ReferenceCount
{
public:

	SF_Replica() {};

	~SF_Replica() {};

	SF_Replica &
	operator=( SF_Replica const & ) { return *this; };

	SF_Replica( SF_Replica const & ) {};

	virtual SF_ReplicaOP clone() const { return 0; };

	virtual core::Real operator ()( core::pose::Pose & pose ) const { return 0.0; };

private:
	int some_private_int_;
};

inline SF_ReplicaOP PyReturnValuePolicyTest_SF_ReplicaOP(void) { return SF_ReplicaOP( new SF_Replica ); };
inline SF_ReplicaCOP PyReturnValuePolicyTest_SF_ReplicaCOP(void) { return SF_ReplicaCOP( new SF_Replica ); };
//inline SF_ReplicaAP PyReturnValuePolicyTest_SF_ReplicaAP(void)  { return new SF_Replica; };
//inline SF_ReplicaCAP PyReturnValuePolicyTest_SF_ReplicaCAP(void)  { return new SF_Replica; };


inline core::pose::PoseOP PyReturnValuePolicyTest_PoseOP(void) { return core::pose::PoseOP(new core::pose::Pose() ); };
inline core::pose::PoseCOP PyReturnValuePolicyTest_PoseCOP(void) { return core::pose::PoseCOP( new core::pose::Pose() ); };
//inline core::pose::PoseAP PyReturnValuePolicyTest_PoseAP(void)  { return new core::pose::Pose(); };
//inline core::pose::PoseCAP PyReturnValuePolicyTest_PoseCAP(void)  { return new core::pose::Pose(); };

inline core::scoring::ScoreFunctionOP PyReturnValuePolicyTest_ScoreFunctionOP(void) { return core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ); };
inline core::scoring::ScoreFunctionCOP PyReturnValuePolicyTest_ScoreFunctionCOP(void) { return core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ); };
inline core::scoring::ScoreFunctionCOP PyReturnValuePolicyTest_ScoreFunctionCOP2(void) { return core::scoring::ScoreFunctionCOP( new core::scoring::ScoreFunction() ); };


} //toolbox
} //protocols

#endif  // INCLUDED_protocols_toolbox_PyReturnValuePolicyTest_hh
