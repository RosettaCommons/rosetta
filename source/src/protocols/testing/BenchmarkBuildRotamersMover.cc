// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/testing/BenchmarkBuildRotamersMover.cc
/// @brief For benchmarking purposes only. Builds data required to run pack rotamers but does not do any sampling. Input pose is untouched.
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/testing/BenchmarkBuildRotamersMover.hh>
#include <protocols/testing/BenchmarkBuildRotamersMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace testing {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;

using basic::Warning;
using basic::t_warning;
static basic::Tracer TR( "protocols.testing.BenchmarkBuildRotamersMover" );

// BenchmarkBuildRotamersMover

std::string
BenchmarkBuildRotamersMoverCreator::keyname() const
{
	return BenchmarkBuildRotamersMover::mover_name();
}

protocols::moves::MoverOP
BenchmarkBuildRotamersMoverCreator::create_mover() const {
	return utility::pointer::make_shared< BenchmarkBuildRotamersMover >();
}

void BenchmarkBuildRotamersMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	BenchmarkBuildRotamersMover::provide_xml_schema( xsd );
}


BenchmarkBuildRotamersMover::BenchmarkBuildRotamersMover() :
	minimization_packing::PackRotamersMover()
{}

BenchmarkBuildRotamersMover::BenchmarkBuildRotamersMover(
	utility::options::OptionCollection const & options
) :
	minimization_packing::PackRotamersMover( options )
{}

BenchmarkBuildRotamersMover::BenchmarkBuildRotamersMover( std::string const & type_name ) :
	minimization_packing::PackRotamersMover( type_name )
{}

// constructors with arguments
BenchmarkBuildRotamersMover::BenchmarkBuildRotamersMover(
	ScoreFunctionCOP scorefxn,
	PackerTaskCOP task,
	Size nloop
) :
	minimization_packing::PackRotamersMover( scorefxn, task, nloop )
{}

BenchmarkBuildRotamersMover::~BenchmarkBuildRotamersMover() = default;

BenchmarkBuildRotamersMover::BenchmarkBuildRotamersMover( BenchmarkBuildRotamersMover const & other ) :
	minimization_packing::PackRotamersMover( other )
{}

void
BenchmarkBuildRotamersMover::apply( Pose & pose )
{
	//ONLY initializes data, no sampling
	this->setup( pose );//call PackRotamersMover::setup()
}

std::string
BenchmarkBuildRotamersMoverCreator::mover_name() {
	return "BenchmarkBuildRotamersMover";
}

std::string
BenchmarkBuildRotamersMover::mover_name() {
	return BenchmarkBuildRotamersMoverCreator::mover_name();
}


/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BenchmarkBuildRotamersMover::fresh_instance() const
{
	return utility::pointer::make_shared< BenchmarkBuildRotamersMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BenchmarkBuildRotamersMover::clone() const
{
	return utility::pointer::make_shared< protocols::testing::BenchmarkBuildRotamersMover >( *this );
}

void
BenchmarkBuildRotamersMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = minimization_packing::PackRotamersMover::complex_type_generator_for_pack_rotamers_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "For benchmarking purposes only. Builds data required to run pack rotamers but does not do any sampling. Input pose is untouched." )
		.write_complex_type_to_schema( xsd );
}


} // moves
} // protocols

