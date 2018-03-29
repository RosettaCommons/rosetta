// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @author Andy Watkins

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/UTracer.hh>

// Unit headers
#include <protocols/parser/BluePrint.hh>
#include <protocols/task_operations/CrystalContactsOperation.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <map>

static basic::Tracer TR("protocols.task_operations.CrystalContactsOperationTests.cxxtest");

using namespace protocols::task_operations;
using namespace core::pack::task;

class CrystalContactsOperationTests : public CxxTest::TestSuite {

public:

public:

	void setUp(){
		protocols_init();
	}

	// This guy requires a very particularly set up pose and so I'd like someone else
	// to write this one
	void test_CrystalContactsOperation() {
	}

};
