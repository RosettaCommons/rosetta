// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/evaluation/EvaluatorFactory.hh
/// @brief Factory for creating Evaluator objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_evaluation_EvaluatorFactory_hh
#define INCLUDED_protocols_evaluation_EvaluatorFactory_hh

// Unit Headers
#include <protocols/evaluation/EvaluatorFactory.fwd.hh>

// Project Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/EvaluatorCreator.fwd.hh>

// Platform Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/factory/WidgetRegistrator.hh>

// C++ Headers
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace evaluation {

/// Create Evaluator Reporters
class EvaluatorFactory {

	// Private constructor to make it singleton managed
	EvaluatorFactory();
	EvaluatorFactory(const EvaluatorFactory & src); // unimplemented

	EvaluatorFactory const &
	operator=( EvaluatorFactory const & ); // unimplemented

public:

	// Warning this is not called because of the singleton pattern
	virtual ~EvaluatorFactory();

	static EvaluatorFactory * get_instance();

	void factory_register( EvaluatorCreatorCOP creator );
	void add_evaluators(
		std::string const & type_name,
		MetaPoseEvaluator & eval);

	void add_all_evaluators(
		MetaPoseEvaluator & eval);

	/// @brief Replace the load-time EvaluatorCreator with another creator.
	// undefined, commenting out to fix PyRosetta build  void replace_creator( EvaluatorCreatorCOP creator );

	// undefined, commenting out to fix PyRosetta build  EvaluatorCreatorCOP
	// get_creator( std::string const & type_name );

private:

	static EvaluatorFactory * instance_;

	// This is a vector rather than a map so the canonical order of the
	// evaluators is determined by how they are inserted.
	typedef std::vector< std::pair< std::string, EvaluatorCreatorCOP > > EvaluatorCreatorMap;
	EvaluatorCreatorMap types_;
};


/// @brief This templated class will register an instance of an
/// EvaluatorCreator (class T) with the
/// EvaluatorFactory.  It will ensure that no
/// EvaluatorCreator is registered twice, and, centralizes this
/// registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class EvaluatorRegistrator : public utility::factory::WidgetRegistrator< EvaluatorFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< EvaluatorFactory, T > parent;
public:
	EvaluatorRegistrator() : parent() {}
};



} // namespace
} // namespace

#endif
