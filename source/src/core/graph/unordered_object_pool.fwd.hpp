// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef BOOST_unordered_object_pool_fwd_hpp
#define BOOST_unordered_object_pool_fwd_hpp


// std::size_t

// boost::details::pool::default_mutex

#include <boost/pool/poolfwd.hpp>

namespace boost {


// Location: <boost/pool/unordered_object_pool.hpp>
//
template <typename T, typename UserAllocator = default_user_allocator_new_delete>
class unordered_object_pool;

} // namespace boost

#endif
