#ifndef INCLUDED_ObjexxFCL_ObserverGraph_hh
#define INCLUDED_ObjexxFCL_ObserverGraph_hh


// ObserverGraph: Observer Graph Representation
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/Observer.fwd.hh>

// C++ Headers
#include <cstddef>
#include <map>
#include <vector>


namespace ObjexxFCL {
namespace internal {


/// @brief ObserverGraph: Observer Graph Representation
class ObserverGraph
{


public: // Types


	// STL style
	typedef  std::size_t  size_type;

	// C++ style
	typedef  std::size_t  Size;

	typedef  std::map< Observer *, size_type >  Graph; // Maps Observers to their in-degree counts
	typedef  std::vector< Graph::iterator >  Sources; // Iterators to zero in-degree graph nodes


public: // Creation


	/// @brief Subject Constructor
	ObserverGraph( Subject const & s );


	/// @brief Destructor
	inline
	~ObserverGraph()
	{}


public: // Inspector


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return graph_.empty();
	}


public: // Modifier


	/// @brief Push a Subject's Transitive Observers onto Graph and Return Acyclicity
	bool
	push( Subject const & s_root, Subject const & s );


	/// @brief Pop a Source Observer from Graph
	Observer *
	pop();


private: // Data


	/// @brief Graph representation
	Graph graph_;

	/// @brief Source Observers with in-degree == zero
	Sources sources_;


}; // ObserverGraph


} // namespace internal
} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_ObserverGraph_HH
