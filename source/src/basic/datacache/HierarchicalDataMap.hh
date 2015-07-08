// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_basic_datacache_HierarchicalDataMap_hh
#define INCLUDED_basic_datacache_HierarchicalDataMap_hh

// Project headers
#include <basic/datacache/HierarchicalDataMap.fwd.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace basic {
namespace datacache {

/// @brief A data map designed to help create hierarchies with complicated 
/// shared data requirements.
///
/// @details To explain what exactly this class is supposed to do, consider the 
/// following mover hierarchy (which comes from the loop modeling framework and 
/// was the actual motivation for this code):
///
/// <pre>
///                  A                         ~
///                  |                         ~
///   +-------+------+------+------+           ~
///   |       |      |      |      |           ~
///   B       C      D      E      F           ~
///          / \    / \           / \          ~
///         G   H  I   J         K   L         ~
///                |   |        /   /|\        ~
///                M   N       O   P Q R       ~
/// </pre>
/// 
/// To explain this diagram a little, each letter corresponds to a mover.  The 
/// movers are mostly of different classes, although they all inherit from the 
/// same abstract base class.  This whole hierarchy is constructed by A.  When 
/// A.apply() is called, it calls B.apply(), C.apply(), D.apply(), E.apply(), 
/// and F.apply() in that order.  When B.apply() is called, it calls G.apply() 
/// and H.apply(), and so on.  If you are interested in seeing exactly how this 
/// hierarchy is implemented, here are the classes the letters correspond to:
///
/// A. protocols::loop_modeling::LoopModeler
/// B. protocols::loop_modeling::utilities::PrepareForCentroid
/// C. protocols::loop_modeling::LoopBuilder
/// D. protocols::loop_modeling::LoopProtocol
/// E. protocols::loop_modeling::utilities::PrepareForFullatom
/// F. protocols::loop_modeling::LoopProtocol
/// G. protocols::kinematic_closure::KicMover
/// H. protocols::loop_modeling::refiners::MinimizationRefiner
/// I. protocols::loop_modeling::utilities::LoopMoverGroup
/// J. protocols::loop_modeling::utilities::LoopMoverGroup
/// K. protocols::loop_modeling::utilities::LoopMoverGroup
/// L. protocols::loop_modeling::utilities::LoopMoverGroup
/// M. protocols::kinematic_closure::KicMover
/// N. protocols::loop_modeling::refiners::MinimizationRefiner
/// O. protocols::kinematic_closure::KicMover
/// P. protocols::loop_modeling::refiners::RepackingRefiner
/// Q. protocols::loop_modeling::refiners::RotamerTrialsRefiner
/// R. protocols::loop_modeling::refiners::MinimizationRefiner
/// 
/// This hierarchy had the following shared data requirements:
/// 
/// 1. Everything in the hierarchy has to share a loops object.
/// 2. A has to know about centroid and fullatom score functions.
/// 3. B, D, and all their children have to use the centroid score function.
/// 4. F and all its children have to use the fullatom score function.
/// 5. A has to provide a default task factory to all its children.
/// 6. Only P and Q have any use for a task factory, and they may want to 
///    override the default provided by A.
/// 
/// Many of these requirements have a similar form: movers higher in the 
/// hierarchy must be able to provide default attributes to all their children, 
/// and movers lower in the hierarchy must be able to override those defaults.  
/// This class supports that pattern by being a data map that asks its parent 
/// for missing values.
///
/// More specifically, this class stores two pieces of information: a DataMap 
/// and a pointer to a parent HierarchicalDataMap.  The DataMap provides the 
/// basic attribute lookup capability: it's just a string to owning pointer 
/// map.  The parent pointer allows the HierarchicalDataMap to recursively 
/// access default values from parent maps for key that aren't in it's own 
/// DataMap.
///
/// Going back to the example hierarchy above: A needs to provide a default 
/// task factory that P and Q can override.  To accomplish this, we start by 
/// giving all the objects in the hierarchy HierarchicalDataMap objects that 
/// have been properly connected to their parents.  Then A sets the default 
/// "task_ops" attribute in it's own HierarchicalDataMap as soon as it's 
/// constructed.  P and Q provide convenience methods like get_task_ops() and 
/// set_task_ops() to allow the "task_ops" field of their HierarchicalDataMaps 
/// to be accessed publicly.  Since the data maps in P and Q are connected to A 
/// via L and F, the attributes stored in A will serve as defaults for P and Q.
/// 
/// It's worth noting that in this example, it is also possible to get the 
/// right behavior by putting a task operations data member in the shared base 
/// class.  However, this is a poor solution for two reasons.  The first is 
/// that it requires many classes that have no need for task operations to have 
/// a task operations attribute.  Second, it requires that you write the 
/// default logic for each attribute that needs it, which is a duplication of 
/// effort.  Using HierarchicalDataMap addresses both of these problems.
///
/// When using this class to help build a hierarchy like the one illustrated 
/// above, you might find it useful to move public accessor methods (like 
/// get_task_ops() and set_task_ops() in P and Q) into a reusable mixin class.  
/// This takes a little bit of C++ template magic, because you have to use the 
/// "curiously recurring template pattern" (CRTP).  But you can see how this 
/// works by looking at protocols::loop_modeling::RepackingLoopMover (P in the 
/// example) and protocols::loop_modeling::utilites::TaskFactoryMixin (the 
/// mixin class).

class HierarchicalDataMap : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor.
	HierarchicalDataMap();

	/// @brief Default destructor.
	virtual ~HierarchicalDataMap();

public:

	/// @brief Set a parent for this data map.
	/// @details If a key is requested and not found in this data map, the search 
	/// will continue (recursively) in the parent.  Having a parent is optional.
	void set_parent(HierarchicalDataMapCAP parent);

	/// @brief Unset the parent for this data map.
	/// @details When a data map doesn't have a parent, it will only return keys 
	/// it is holding itself.  Having a parent is optional.
	void unset_parent();

public:

	/// @brief Return the value associated with the given key.  If the key is not 
	/// found, throw an exception.
	///
	/// @details The key is first searched for in this data map, and then 
	/// subsequently in parent data maps.  The first key that is found is 
	/// returned.  The key is composed of a type and a name.  Furthermore, the 
	/// name "" is handled specially: it matches any name that's not already in 
	/// the data map.  For example, consider a data map with the following keys: 
	/// ("spam", "") and ("spam", "eggs").  If you request ("spam", "bacon"), the 
	/// value associated with ("spam", "") will be returned.  If ("spam", "") 
	/// hadn't been present, ("spam", "bacon") wouldn't have been found and the 
	/// search would have continued.
	template< typename ValueOP >
	ValueOP get(std::string const & type, std::string const & name) const {
		ValueOP value = get<ValueOP>(type, name, ValueOP());
		if (! value) {
			std::stringstream message;
			message << "HierarchicalDataMap has no attribute with type '" << type << "'";
			if (! name.empty()) { message << "' and name '" << name << "'"; }
			message << "." << std::endl;
			throw utility::excn::EXCN_Msg_Exception(message.str());
		}
		return value;
	}

	/// @brief Return the value associated with the given key.  If the key is not 
	/// found, return the given fallback value.
	///
	/// @details The key is first searched for in this data map, and then 
	/// subsequently in parent data maps.  The first key that is found is 
	/// returned.  The key is composed of a type and a name.  Furthermore, the 
	/// name "" is handled specially: it matches any name that's not already in 
	/// the data map.  For example, consider a data map with the following keys: 
	/// ("spam", "") and ("spam", "eggs").  If you request ("spam", "bacon"), the 
	/// value associated with ("spam", "") will be returned.  If ("spam", "") 
	/// hadn't been present, ("spam", "bacon") wouldn't have been found and the 
	/// search would have continued.
	template< typename ValueOP >
	ValueOP get(std::string const & type, std::string const & name, ValueOP fallback) const {

		// If the requested key is present in this data map, return it.

		if (data_map_.has(type, name)) {
			return data_map_.get_ptr<typename ValueOP::element_type>(type, name);
		}

		else if (data_map_.has(type, "")) {
			return data_map_.get_ptr<typename ValueOP::element_type>(type, "");
		}

		// If not, recursively see if any parent data maps have it.

		else if (HierarchicalDataMapCOP parent = parent_.lock()) {
			return parent->template get<ValueOP>(type, name, fallback);
		}

		// Return the fallback value if the requested key was not found.
		
		else {
			return fallback;
		}
	}

	/// @brief Return the value associated with the given key.  If the key is not 
	/// found, return a null owning pointer.
	///
	/// @details The key is first searched for in this data map, and then 
	/// subsequently in parent data maps.  The first key that is found is 
	/// returned.  The key is composed of a type and a name.  Furthermore, the 
	/// name "" is handled specially: it matches any name that's not already in 
	/// the data map.  For example, consider a data map with the following keys: 
	/// ("spam", "") and ("spam", "eggs").  If you request ("spam", "bacon"), the 
	/// value associated with ("spam", "") will be returned.  If ("spam", "") 
	/// hadn't been present, ("spam", "bacon") wouldn't have been found and the 
	/// search would have continued.
	template< typename ValueOP >
	ValueOP get_or_null(std::string const & type, std::string const & name) const {
		return get<ValueOP>(type, name, ValueOP());
	}

	/// @brief Set the value of the given key in this map.  This value may also 
	/// be accessed by children maps.
	template< typename ValueOP >
	ValueOP set( std::string const & type, std::string const & name, ValueOP value ) {
		data_map_[type][name] = value;
		return value;
	}

private:
	HierarchicalDataMapCAP parent_;
	DataMap data_map_;

};

}
}

#endif
