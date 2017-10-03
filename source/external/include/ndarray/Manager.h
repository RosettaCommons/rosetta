// -*- c++ -*-
/*
 * Copyright (c) 2010-2012, Jim Bosch
 * All rights reserved.
 *
 * ndarray is distributed under a simple BSD-like license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/ndarray/ndarray
 */
#ifndef NDARRAY_Manager_h_INCLUDED
#define NDARRAY_Manager_h_INCLUDED

/** 
 *  @file ndarray/Manager.h
 *
 *  @brief Definition of Manager, which manages the ownership of array data.
 */

#include "ndarray_fwd.h"
#include <boost/noncopyable.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/scoped_array.hpp>

namespace ndarray {

class Manager : private boost::noncopyable {
public:

    typedef boost::intrusive_ptr<Manager> Ptr;

    friend inline void intrusive_ptr_add_ref(Manager const * manager) {
        ++manager->_rc;
    }
 
    friend inline void intrusive_ptr_release(Manager const * manager) {
        if ((--manager->_rc)==0) delete manager;
    }

    int getRC() const { return _rc; }

    virtual bool isUnique() const { return false; }

protected:

    virtual ~Manager() {}

    explicit Manager() : _rc(0) {}

private:
    mutable int _rc;
};

template <typename T>
class SimpleManager : public Manager {
    typedef typename boost::remove_const<T>::type U;
public:
    
    static std::pair<Manager::Ptr,T*> allocate(Size size) {
        boost::intrusive_ptr<SimpleManager> r(new SimpleManager(size));
        return std::pair<Manager::Ptr,T*>(r, r->_p.get());
    }

    virtual bool isUnique() const { return true; }

private:
    explicit SimpleManager(Size size) : _p() {
        if (size > 0) _p.reset(new U[size]);
    }
    boost::scoped_array<U> _p;
};

template <typename T> Manager::Ptr makeManager(T const & owner);

template <typename U>
class ExternalManager : public Manager, private U {
public:
    typedef U Owner;

    template <typename T> friend Manager::Ptr makeManager(T const & owner);

    Owner const & getOwner() const { return *static_cast<Owner const *>(this); }

private:
    explicit ExternalManager(Owner const & owner) : Owner(owner) {}
};

template <typename T>
inline Manager::Ptr makeManager(T const & owner) {
    return Manager::Ptr(new ExternalManager<T>(owner));
}

// A no-op overload for makeManager to avoid unnecessary levels of indirection.
inline Manager::Ptr makeManager(Manager::Ptr const & owner) {
    return owner;
}

} // namespace ndarray

#endif // !NDARRAY_Manager_h_INCLUDED
