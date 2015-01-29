// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/Filter.hh
/// @brief  True/False query class which handles basic boolean logic
/// @author Will Sheffler (willsheffler@gmail.com)
/// @date   Thu Aug  9 19:49:23 2007
///

#ifndef INCLUDED_utility_query_Filter_HH
#define INCLUDED_utility_query_Filter_HH

#include <utility/query/types.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

namespace utility {
namespace query {

	using utility::pointer::ReferenceCount;


template<class T>
class FilterBase : public utility::pointer::ReferenceCount {
public:

	// virtual ~FilterBase() {}

	virtual bool operator()( T arg ) = 0;

	virtual std::string description(std::string indent="") = 0;

	virtual owning_ptr<FilterBase<T> > clone()  = 0;

};

template<class T>
class Filter : public FilterBase<T> {
public:
	Filter() : wrapped_(NULL) {}
	Filter(FilterBase<T> & tf) { *this = tf; }
	Filter(FilterBase<T> * tf) { *this = tf; }
	Filter(owning_ptr<FilterBase<T> > tf) { *this = tf; }
	// virtual ~Filter() {}
	void operator=(owning_ptr<FilterBase<T> > p) { wrapped_ = p; }
	void operator=(FilterBase<T> & tf) { wrapped_ = tf.clone(); }
	void operator=(FilterBase<T> * tf) { wrapped_ = tf->clone(); }
	inline bool operator()( T arg ) {
		if( wrapped_ ) return (*wrapped_)(arg);
		std::cerr << "WARNING: calling Filter w/ no held query!!!" << std::endl;
		return false;
		}
	std::string description(std::string indent="") {
		if( wrapped_ ) return indent + wrapped_->description();
		std::cerr << "WARNING: calling Filter w/ no held query!!!" << std::endl;
		return indent + "You Must Define string description() in Your FilterBase<T> Subclass!";
	}
	owning_ptr<FilterBase<T> > clone() {
	debug_assert( wrapped_ ); // foo!
		return wrapped_;
	}
protected:
	owning_ptr<FilterBase<T> > wrapped_;
};

template<class T>
class Filter0Param : public FilterBase<T> {
public:
	// virtual ~Filter0Param() {};
	Filter0Param( owning_ptr<Converter<T,bool> > conv) : converter_(conv) {}
	inline bool operator() (T arg) { return converter_->convert(arg); }
	std::string description(std::string indent="") { return indent + converter_.description(); }
	owning_ptr<FilterBase<T> > clone() { return owning_ptr<FilterBase<T> >(new Filter0Param<T>(converter_)); }
private:
	owning_ptr<Converter<T,bool> > converter_;
};

template<class T, typename A>
class Filter1Param : public FilterBase<T> {
public:
	Filter1Param( C1P_TAB_OP conv, A param) : converter_(conv),param_(param) {}
	// virtual ~Filter1Param() {};
	inline bool operator() (T arg) {
		return converter_->convert(arg,param_);
	}
	std::string description(std::string indent="") {
		std::ostringstream o; o << param_;
		return indent + converter_->description() + " " + o.str();
	}
	owning_ptr<FilterBase<T> > clone() {
		return owning_ptr<FilterBase<T> >(new Filter1Param<T,A>(converter_,param_));
	}
private:
	C1P_TAB_OP converter_;
	A param_;
};

template<class T, typename A>
class Filter1ParamGenerator {
public:
	Filter1ParamGenerator( C1P_TAB_OP conv ) : converter_(conv) {}
	owning_ptr<FilterBase<T> > operator()(A param) { return owning_ptr<FilterBase<T> >( new Filter1Param<T,A>(converter_,param) ); }
private:
	C1P_TAB_OP converter_;
};

template<class T>
class Filter_And : public FilterBase<T> {
public:
	Filter_And( owning_ptr<FilterBase<T> > lhs          , owning_ptr<FilterBase<T> > rhs          ) : lhs_(lhs)         , rhs_(rhs)         { }
	Filter_And( FilterBase<T> & lhs , owning_ptr<FilterBase<T> > rhs          ) : lhs_(lhs.clone()) , rhs_(rhs)         { }
	Filter_And( owning_ptr<FilterBase<T> > lhs          , FilterBase<T> & rhs ) : lhs_(lhs)         , rhs_(rhs.clone()) { }
	Filter_And( FilterBase<T> & lhs , FilterBase<T> & rhs ) : lhs_(lhs.clone()) , rhs_(rhs.clone()) { }
	// virtual ~Filter_And() {}
	inline bool operator()( T arg ) {
		return ( (*lhs_)(arg) && (*rhs_)(arg) );
	}
	std::string description(std::string indent="") {
		std::string s = "";
		s += lhs_->description(indent+"    ")+"\n";
		s += indent + "AND\n";
		s += rhs_->description(indent+"    ");
		return s;
	}
	owning_ptr<FilterBase<T> > clone() {
		return new Filter_And( lhs_, rhs_ );
	}
protected:
	owning_ptr<FilterBase<T> > lhs_, rhs_;
};

template<class T>
class Filter_Or  : public FilterBase<T> {
public:
	Filter_Or( owning_ptr<FilterBase<T> > lhs          , owning_ptr<FilterBase<T> > rhs          ) : lhs_(lhs)         , rhs_(rhs)         { }
	Filter_Or( FilterBase<T> & lhs , owning_ptr<FilterBase<T> > rhs          ) : lhs_(lhs.clone()) , rhs_(rhs)         { }
	Filter_Or( owning_ptr<FilterBase<T> > lhs          , FilterBase<T> & rhs ) : lhs_(lhs)         , rhs_(rhs.clone()) { }
	Filter_Or( FilterBase<T> & lhs , FilterBase<T> & rhs ) : lhs_(lhs.clone()) , rhs_(rhs.clone()) { }
	// virtual ~Filter_Or() {}
	inline bool operator()( T arg ) {
		return ( (*lhs_)(arg) || (*rhs_)(arg) );
	}
	std::string description(std::string indent="") {
		std::string s = "";
		s += lhs_->description(indent+"    ")+"\n";
		s += indent + "OR\n";
		s += rhs_->description(indent+"    ");
		return s;
	}
	owning_ptr<FilterBase<T> > clone() { return new Filter_Or( lhs_, rhs_ ); }
protected:
	owning_ptr<FilterBase<T> > lhs_, rhs_;
};

template<class T>
class Filter_Not : public FilterBase<T> {
public:
	Filter_Not(    owning_ptr<FilterBase<T> >     rhs ) : rhs_(rhs) { }
	Filter_Not( FilterBase<T> & rhs ) : rhs_(rhs.clone()) { }
	// virtual ~Filter_Not() {}
	inline bool operator()( T arg ) { return !( (*rhs_)(arg) ); }
	std::string description(std::string indent="") {
		std::string s = indent + "NOT\n";
		s += rhs_->description(indent+"    ");
		return s;
	}
	owning_ptr<FilterBase<T> > clone() { return new Filter_Not( rhs_ ); }
protected:
	owning_ptr<FilterBase<T> > rhs_;
};

namespace filter_operators {
#define FILT_OP utility::pointer::owning_ptr< FilterBase< T > > // dunno what else to do..

template<class T> FILT_OP create_ander(FILT_OP lhs, FILT_OP rhs) { return (FILT_OP)( new Filter_And<T>( lhs, rhs ) ); }
template<class T> FILT_OP create_orer (FILT_OP lhs, FILT_OP rhs) { return (FILT_OP)( new Filter_Or <T>( lhs, rhs ) ); }
template<class T> FILT_OP create_noter(FILT_OP rhs           ) { return (FILT_OP)( new Filter_Not<T>( rhs ) );}

template<class T> FILT_OP operator & ( FILT_OP lhs, FILT_OP rhs ) { return create_ander<T>(lhs,rhs); }
template<class T> FILT_OP operator | ( FILT_OP lhs, FILT_OP rhs ) { return create_orer <T>(lhs,rhs); }
template<class T> FILT_OP operator ~ ( FILT_OP rhs )            { return create_noter <T>(rhs); }

template<class T> FILT_OP operator & ( FilterBase<T> & lhs, FILT_OP rhs ) { return create_ander<T>(lhs.clone(),rhs); }
template<class T> FILT_OP operator | ( FilterBase<T> & lhs, FILT_OP rhs ) { return create_orer <T>(lhs.clone(),rhs); }

template<class T> FILT_OP operator & ( FILT_OP lhs, FilterBase<T> & rhs ) { return create_ander<T>(lhs,rhs.clone()); }
template<class T> FILT_OP operator | ( FILT_OP lhs, FilterBase<T> rhs ) { return create_orer <T>(lhs,rhs.clone()); }

template<class T> FILT_OP operator ~ ( FilterBase<T> & rhs ) { return create_noter<T>(rhs.clone()); }
template<class T> FILT_OP operator & ( FilterBase<T> & lhs, FilterBase<T> & rhs ) { return create_ander<T>(lhs.clone(),rhs.clone()); }
template<class T> FILT_OP operator | ( FilterBase<T> & lhs, FilterBase<T> & rhs ) { return create_orer <T>(lhs.clone(),rhs.clone()); }
}
using namespace filter_operators;

} // end namespace query
} // end namespace utility




#endif
