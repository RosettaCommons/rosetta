// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_evaluation_PoseEvaluator_hh
#define INCLUDED_protocols_evaluation_PoseEvaluator_hh


// Unit Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>

// Package Headers

// Project Headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
//// C++ headers

// due to template function
#include <core/io/silent/SilentStruct.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace evaluation {

class PoseEvaluator : public utility::pointer::ReferenceCount {
public:
	PoseEvaluator() {};
	virtual ~PoseEvaluator() {};
	// we may put a cache for user values into the pose, then it will just return a pose
	// or return just the Map from strings to values...

	/// @brief evaluate pose and store values in Silent_Struct
	virtual void apply( core::pose::Pose&, std::string tag, core::io::silent::SilentStruct &pss) const = 0;

	/// @brief direct application to SilentStruct...
	/// default implementation makes pose and calls "apply", you can overload if you don't need the pose-step
	virtual void apply( core::io::silent::SilentStruct &pss) const;
	virtual bool applicable( core::pose::Pose const& ) const { return true; }

	virtual core::Size size() const = 0;
	virtual std::string name( core::Size ) const = 0;

};

template <class T >
class SingleValuePoseEvaluator : public PoseEvaluator {
public:
	SingleValuePoseEvaluator( std::string name = "UNSPECIFIED_SingleValuePoseEvaluator" ) : name_( name ) {};

	/// @brief evaluate pose and store values in Silent_Struct
	/// why is this specific to a specific type of SilentStruct? that seems needlessly pointless and overly constraining.
	virtual void apply( core::pose::Pose&, std::string tag, core::io::silent::SilentStruct &pss) const;
	//	void apply( core::pose::Pose&, std::string tag, core::io::silent::SilentStruct &pss) const;

	using PoseEvaluator::apply;

	/// @brief evaluate pose
	virtual T apply( core::pose::Pose& ) const = 0;
	virtual bool applicable( core::pose::Pose const& ) const { return true; }

	//	void apply( core::pose::Pose&, std::string tag, core::io::silent::SilentStruct &pss) const;

	virtual core::Size size() const { return 1; }
	virtual std::string name( core::Size ) const {
		return name_;
	}

private:
	std::string name_;

};

//now also named PoseEvaluators ( see typedef below )
//typedef MetaPoseEvaluator PoseEvaluators; //moved to .fwd.hh, SML 09/09/2009

class MetaPoseEvaluator : public PoseEvaluator {
public:
	typedef PoseEvaluator Parent;
	typedef utility::vector1< PoseEvaluatorOP > EvaluatorList;

	MetaPoseEvaluator() {};

	virtual void
	apply( core::pose::Pose& pose, std::string tag, core::io::silent::SilentStruct &pss) const;

	void
	add_evaluation( PoseEvaluatorOP pe) {
		evaluators_.push_back( pe );
	}

	/// @brief clear the list of evaluators
	void
	clear() {
		evaluators_.clear();
	}

	//removes last element
	void
	pop_back() {
		evaluators_.pop_back();
	}

	MetaPoseEvaluator& operator<< ( PoseEvaluatorOP pe) {
		add_evaluation( pe );
		return *this;
	}

	Size size() const {
		Size s( 0 );
		for ( EvaluatorList::const_iterator it = evaluators_.begin(); it != evaluators_.end(); ++it ) {
			s += (*it)->size();
		}
		return s;
	}

	virtual std::string name( core::Size ind ) const;

	EvaluatorList const& evaluators() { return evaluators_; }

private:
	EvaluatorList evaluators_;
};


//@brief evaluate pose and store values in Silent_Struct
template< class T >
void SingleValuePoseEvaluator< T >::apply( core::pose::Pose &pose, std::string, core::io::silent::SilentStruct &pss) const {
	if ( applicable( pose ) )	pss.add_energy ( name_, apply( pose ) );
}


// PyRosetta concreate classes
class SingleValuePoseEvaluator_Size : public SingleValuePoseEvaluator<core::Size> {};
class SingleValuePoseEvaluator_SSize: public SingleValuePoseEvaluator<core::SSize> {};


}
}
#endif
