// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/testing/BenchmarkBuildRotamersMover.hh
/// @brief For benchmarking purposes only. Builds data required to run pack rotamers but does not do any sampling. Input pose is untouched.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_testing_BenchmarkBuildRotamersMover_hh
#define INCLUDED_protocols_testing_BenchmarkBuildRotamersMover_hh

// Unit headers
#include <protocols/testing/BenchmarkBuildRotamersMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

// Project headers
#include <core/types.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace testing {

/// @brief A protocols::moves::Mover that packs the side-chains using a rotamer library
/// It uses a ScoreFunction for packing and a PackerTask,
/// or a TaskFactory that generates a PackerTask, for instructions on
/// what rotamer sets are allowed at each residue position during packing
///
/// Common Methods:
///     BenchmarkBuildRotamersMover.apply
/// @note please derive from BenchmarkBuildRotamersMover instead of attempting to add protocol-specific stuff here!
class BenchmarkBuildRotamersMover : public minimization_packing::PackRotamersMover {
public:
	/// @brief default constructor; reads the nloop_ value from the global options system.
	BenchmarkBuildRotamersMover();

	/// @brief constructor that reads from the (possibly local) options collection object
	BenchmarkBuildRotamersMover( utility::options::OptionCollection const & options );

	/// @brief constructor with typename; reads the nloop_ value from the global options system.
	BenchmarkBuildRotamersMover( std::string const & );

	/// @brief Constructs a BenchmarkBuildRotamersMover with PackerTask  <task>
	/// evaluated using  <scorefxn>
	///
	/// ScoreFunction  scorefxn   /function to minimize while changine rotamers
	/// PackerTask     task       /object specifying what to design/pack
	/// Size (int)     nloop      /number of rounds to run packing
	BenchmarkBuildRotamersMover(
		core::scoring::ScoreFunctionCOP scorefxn,
		core::pack::task::PackerTaskCOP task = nullptr,
		core::Size nloop = 1
	);

	// destructor (important for properly forward-declaring smart-pointer members)
	~BenchmarkBuildRotamersMover() override;

	// copy constructor
	BenchmarkBuildRotamersMover( BenchmarkBuildRotamersMover const & other );

	//just calls PackRotamersMover::setup
	void apply( Pose & pose ) override;

	std::string get_name() const override {
		return mover_name();
	}

	static std::string mover_name();

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP clone() const override;

	//using minimization_packing::PackRotamersMover::provide_xml_schema;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	//static utility::tag::XMLSchemaComplexTypeGeneratorOP complex_type_generator_for_benchmark_build_rotamers_mover( utility::tag::XMLSchemaDefinition & xsd );

};

} // moves
} // protocols

#endif
