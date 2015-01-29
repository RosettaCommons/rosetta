// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (C) 199x-2007 Vadderbilt University

/// @file   utility/Metric.hh
/// @brief  scalar metric class which allows adding / subtracting mult etc...
/// @author Will Sheffler (willsheffler@gmail.com)
/// @date   Thu Aug  9 19:49:23 2007
///

#ifndef INCLUDED_utility_query_Metric_HH
#define INCLUDED_utility_query_Metric_HH

#include <core/types.hh>

#include <utility/query/Filter.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/string.functions.hh>

#include <string>
#include <sstream>
#include <math.h>


namespace utility {
namespace query {

	using core::Real;
	using utility::pointer::owning_ptr;
	using utility::pointer::ReferenceCount;

template<class T>
class MetricBase : public utility::pointer::ReferenceCount {
public:

	// virtual ~MetricBase() {}

	virtual Real operator()( T arg ) = 0;

	virtual std::string description(std::string indent="") = 0;

	virtual owning_ptr<MetricBase<T> > clone()  = 0;

};

template<class T>
class Metric : public MetricBase<T> {
public:
	Metric() : wrapped_(NULL) {}
	Metric(MetricBase<T> & tf) { *this = tf; }
	Metric(MetricBase<T> * tf) { *this = tf; }
	Metric(owning_ptr<MetricBase<T> > tf) { *this = tf; }
	// virtual ~Metric() {}
	void operator=(owning_ptr<MetricBase<T> > p) { wrapped_ = p; }
	void operator=(MetricBase<T> & tf) { wrapped_ = tf.clone(); }
	void operator=(MetricBase<T> * tf) { wrapped_ = tf->clone(); }
	inline Real operator()( T arg ) {
		if( wrapped_ ) return (*wrapped_)(arg);
		std::cerr << "WARNING: calling Metric w/ no held query!!!" << std::endl;
		return false;
		}
	std::string description(std::string indent="") {
		if( wrapped_ ) return indent + wrapped_->description();
		std::cerr << "WARNING: calling Metric w/ no held query!!!" << std::endl;
		return indent + "You Must Define string description() in Your MetricBase<T> Subclass!";
	}
	owning_ptr<MetricBase<T> > clone() {
	debug_assert( wrapped_ ); // foo!
		return wrapped_;
	}
protected:
	owning_ptr<MetricBase<T> > wrapped_;
};

template<class T>
class MetricConstant : public MetricBase<T> {
public:
	MetricConstant(Real const value) : value_(value) {}
	// virtual ~MetricConstant() {}
	inline Real operator()( T ) { return value_; }
	std::string description(std::string indent="") { return indent+ObjexxFCL::string_of(value_); };
	owning_ptr<MetricBase<T> > clone() { return owning_ptr<MetricBase<T> >(new MetricConstant(value_)); }
private:
	Real const value_;
};

template<class T>
class Metric0Param : public MetricBase<T> {
public:
	// virtual ~Metric0Param() {};
	Metric0Param() : converter_(new ImplicitConverter<T,Real>) {}
	Metric0Param( owning_ptr<Converter<T,Real> > conv) : converter_(conv) {}
	inline Real operator() (T arg) { return converter_->convert(arg); }
	std::string description(std::string indent="") { return indent + converter_->description(); }
	owning_ptr<MetricBase<T> > clone() { return owning_ptr<MetricBase<T> >(new Metric0Param<T>(converter_)); }
private:
	owning_ptr<Converter<T,Real> > converter_;
};

template<class T, typename A>
class Metric1Param : public MetricBase<T> {
public:
	Metric1Param( C1P_TAR_OP conv, A param) : converter_(conv),param_(param) {}
	// virtual ~Metric1Param() {};
	inline Real operator() (T arg) {
		return converter_->convert(arg,param_);
	}
	std::string description(std::string indent="") {
		std::ostringstream o; o << param_;
		return indent + converter_->description() + " " + o.str();
	}
	owning_ptr<MetricBase<T> > clone() {
		return owning_ptr<MetricBase<T> >(new Metric1Param<T,A>(converter_,param_));
	}
private:
	C1P_TAR_OP converter_;
	A param_;
};

template<class T, typename A>
class Metric1ParamGenerator {
public:
	Metric1ParamGenerator( C1P_TAR_OP conv ) : converter_(conv) {}
	owning_ptr<MetricBase<T> > operator()(A param) { return owning_ptr<MetricBase<T> >( new Metric1Param<T,A>(converter_,param) ); }
private:
	C1P_TAR_OP converter_;
};

///////////////////////////////////////////
// operator nodes (fromed by expressions)
///////////////////////////////////////////

template<class T>
class Metric_Add : public MetricBase<T> {
public:
	Metric_Add( owning_ptr<MetricBase<T> > lhs          , owning_ptr<MetricBase<T> > rhs          ) : lhs_(lhs)         , rhs_(rhs)         { }
	Metric_Add( MetricBase<T> & lhs , owning_ptr<MetricBase<T> > rhs          ) : lhs_(lhs.clone()) , rhs_(rhs)         { }
	Metric_Add( owning_ptr<MetricBase<T> > lhs          , MetricBase<T> & rhs ) : lhs_(lhs)         , rhs_(rhs.clone()) { }
	Metric_Add( MetricBase<T> & lhs , MetricBase<T> & rhs ) : lhs_(lhs.clone()) , rhs_(rhs.clone()) { }
	// virtual ~Metric_Add() {}
	inline Real operator()( T arg ) {
		return ( (*lhs_)(arg) + (*rhs_)(arg) );
	}
	std::string description(std::string indent="") {
		return lhs_->description(indent+"    ")+"\n" + indent + "+\n" + rhs_->description(indent+"    ");
	}
	owning_ptr<MetricBase<T> > clone() {
		return new Metric_Add( lhs_, rhs_ );
	}
protected:
	owning_ptr<MetricBase<T> > lhs_, rhs_;
};

template<class T>
class Metric_Sub : public MetricBase<T> {
public:
	Metric_Sub( owning_ptr<MetricBase<T> > lhs, owning_ptr<MetricBase<T> > rhs ) : lhs_(lhs)         , rhs_(rhs)         { }
	Metric_Sub( MetricBase<T> &            lhs, owning_ptr<MetricBase<T> > rhs ) : lhs_(lhs.clone()) , rhs_(rhs)         { }
	Metric_Sub( owning_ptr<MetricBase<T> > lhs, MetricBase<T> &            rhs ) : lhs_(lhs)         , rhs_(rhs.clone()) { }
	Metric_Sub( MetricBase<T> &            lhs, MetricBase<T> &            rhs ) : lhs_(lhs.clone()) , rhs_(rhs.clone()) { }
	// virtual ~Metric_Add() {}
	inline Real operator()( T arg ) {
		return ( (*lhs_)(arg) - (*rhs_)(arg) );
	}
	std::string description(std::string indent="") {
		return lhs_->description(indent+"    ")+"\n" + indent + "-\n" + rhs_->description(indent+"    ");
	}
	owning_ptr<MetricBase<T> > clone() {
		return new Metric_Sub( lhs_, rhs_ );
	}
protected:
	owning_ptr<MetricBase<T> > lhs_, rhs_;
};

template<class T>
class Metric_Mult  : public MetricBase<T> {
public:
	Metric_Mult( owning_ptr<MetricBase<T> > lhs, owning_ptr<MetricBase<T> > rhs ) : lhs_(lhs)         , rhs_(rhs)         { }
	Metric_Mult( MetricBase<T> &            lhs, owning_ptr<MetricBase<T> > rhs ) : lhs_(lhs.clone()) , rhs_(rhs)         { }
	Metric_Mult( owning_ptr<MetricBase<T> > lhs, MetricBase<T> & rhs            ) : lhs_(lhs)         , rhs_(rhs.clone()) { }
	Metric_Mult( MetricBase<T> &            lhs, MetricBase<T> & rhs            ) : lhs_(lhs.clone()) , rhs_(rhs.clone()) { }
	// virtual ~Metric_Mult() {}
	inline Real operator()( T arg ) {
		return ( (*lhs_)(arg) * (*rhs_)(arg) );
	}
	std::string description(std::string indent="") {
		return lhs_->description(indent+"    ")+"\n" + indent + "*\n" + rhs_->description(indent+"    ");
	}
	owning_ptr<MetricBase<T> > clone() { return new Metric_Mult( lhs_, rhs_ ); }
protected:
	owning_ptr<MetricBase<T> > lhs_, rhs_;
};

template<class T>
class Metric_Div  : public MetricBase<T> {
public:
	Metric_Div( owning_ptr<MetricBase<T> > lhs, owning_ptr<MetricBase<T> > rhs ) : lhs_(lhs)         , rhs_(rhs)         { }
	Metric_Div( MetricBase<T> &            lhs, owning_ptr<MetricBase<T> > rhs ) : lhs_(lhs.clone()) , rhs_(rhs)         { }
	Metric_Div( owning_ptr<MetricBase<T> > lhs, MetricBase<T> &            rhs ) : lhs_(lhs)         , rhs_(rhs.clone()) { }
	Metric_Div( MetricBase<T> &            lhs, MetricBase<T> &            rhs ) : lhs_(lhs.clone()) , rhs_(rhs.clone()) { }
	// virtual ~Metric_Mult() {}
	inline Real operator()( T arg ) {
		return ( (*lhs_)(arg) / (*rhs_)(arg) );
	}
	std::string description(std::string indent="") {
		return lhs_->description(indent+"    ")+"\n" + indent + "/\n" + rhs_->description(indent+"    ");
	}
	owning_ptr<MetricBase<T> > clone() { return new Metric_Div( lhs_, rhs_ ); }
protected:
	owning_ptr<MetricBase<T> > lhs_, rhs_;
};

template<class T>
class Metric_Pow  : public MetricBase<T> {
public:
	Metric_Pow( owning_ptr<MetricBase<T> > lhs, owning_ptr<MetricBase<T> > rhs ) : lhs_(lhs)         , rhs_(rhs)         { }
	Metric_Pow( MetricBase<T> &            lhs, owning_ptr<MetricBase<T> > rhs ) : lhs_(lhs.clone()) , rhs_(rhs)         { }
	Metric_Pow( owning_ptr<MetricBase<T> > lhs, MetricBase<T> &            rhs ) : lhs_(lhs)         , rhs_(rhs.clone()) { }
	Metric_Pow( MetricBase<T> &            lhs, MetricBase<T> &            rhs ) : lhs_(lhs.clone()) , rhs_(rhs.clone()) { }
	// virtual ~Metric_Mult() {}
	inline Real operator()( T arg ) {
		return pow( (*lhs_)(arg) , (*rhs_)(arg) );
	}
	std::string description(std::string indent="") {
		return lhs_->description(indent+"    ")+"\n" + indent + "^ (pow)\n" + rhs_->description(indent+"    ");
	}
	owning_ptr<MetricBase<T> > clone() { return new Metric_Pow( lhs_, rhs_ ); }
protected:
	owning_ptr<MetricBase<T> > lhs_, rhs_;
};

template<class T>
class Metric_Neg : public MetricBase<T> {
public:
	Metric_Neg( owning_ptr<MetricBase<T> > rhs ) : rhs_(rhs) { }
	Metric_Neg( MetricBase<T> & rhs ) : rhs_(rhs.clone()) { }
	// virtual ~Metric_Neg() {}
	inline Real operator()( T arg ) { return !( (*rhs_)(arg) ); }
	std::string description(std::string indent="") {
		std::string s = indent + "-\n";
		s += rhs_->description(indent+"    ");
		return s;
	}
	owning_ptr<MetricBase<T> > clone() { return new Metric_Neg( rhs_ ); }
protected:
	owning_ptr<MetricBase<T> > rhs_;
};

enum MetricInequalityType {
	EQ, NEQ, GT, GTE, LT, LTE
};

template<class T>
class MetricInequality : public FilterBase<T> {
public:
	MetricInequality( MetricInequalityType op, owning_ptr<MetricBase<T> > lhs , owning_ptr<MetricBase<T> > rhs ) : op_(op), lhs_(lhs)         , rhs_(rhs)         { }
	MetricInequality( MetricInequalityType op, MetricBase<T> &            lhs , owning_ptr<MetricBase<T> > rhs ) : op_(op), lhs_(lhs.clone()) , rhs_(rhs)         { }
	MetricInequality( MetricInequalityType op, owning_ptr<MetricBase<T> > lhs , MetricBase<T> &            rhs ) : op_(op), lhs_(lhs)         , rhs_(rhs.clone()) { }
	MetricInequality( MetricInequalityType op, MetricBase<T> &            lhs , MetricBase<T> &            rhs ) : op_(op), lhs_(lhs.clone()) , rhs_(rhs.clone()) { }
	// virtual ~MetricInequality() {}
	bool operator()( T arg ) {
		switch(op_) {
			case EQ : return (*lhs_)(arg) == (*rhs_)(arg); break;
			case NEQ: return (*lhs_)(arg) != (*rhs_)(arg); break;
			case GT : return (*lhs_)(arg) >  (*rhs_)(arg); break;
			case GTE: return (*lhs_)(arg) >= (*rhs_)(arg); break;
			case LT : return (*lhs_)(arg) <  (*rhs_)(arg); break;
			case LTE: return (*lhs_)(arg) <= (*rhs_)(arg); break;
		}
	debug_assert(false);
		return false;
	}
	std::string description(std::string indent="") {
		std::string op = "";
		switch(op_) {
			case EQ : op += " == "; break;
			case NEQ: op += " != "; break;
			case GT : op += " > " ; break;
			case GTE: op += " >= "; break;
			case LT : op += " < " ; break;
			case LTE: op += " <= "; break;
		}
		return lhs_->description(indent+"    ")+"\n" + indent + op + "\n" + rhs_->description(indent+"    ");
	}
	FILT_OP clone() {
		return new MetricInequality( op_, lhs_, rhs_ );
	}
protected:
	MetricInequalityType op_;
	owning_ptr<MetricBase<T> > lhs_, rhs_;
};


namespace metric_operators { // operators
#define MET_OP utility::pointer::owning_ptr< MetricBase< T > > // dunno what else to do..

template<class T> MET_OP create_sub (MET_OP lhs, MET_OP rhs) { return (MET_OP)( new Metric_Sub<T>( lhs, rhs ) ); }
template<class T> MET_OP operator - ( MET_OP          lhs, MET_OP          rhs ) { return create_sub<T>( lhs                               , rhs                                ); }
template<class T> MET_OP operator - ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_sub<T>( lhs.clone()                       , rhs                                ); }
template<class T> MET_OP operator - ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_sub<T>( lhs                               , rhs.clone()                        ); }
template<class T> MET_OP operator - ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_sub<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> MET_OP operator - ( double          lhs, MET_OP          rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator - ( double          lhs, MetricBase<T> & rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator - ( MET_OP          lhs, double          rhs ) { return create_sub<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( MetricBase<T> & lhs, double          rhs ) { return create_sub<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( float           lhs, MET_OP          rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator - ( float           lhs, MetricBase<T> & rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator - ( MET_OP          lhs, float           rhs ) { return create_sub<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( MetricBase<T> & lhs, float           rhs ) { return create_sub<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( int             lhs, MET_OP          rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator - ( int             lhs, MetricBase<T> & rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator - ( MET_OP          lhs, int             rhs ) { return create_sub<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( MetricBase<T> & lhs, int             rhs ) { return create_sub<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( long            lhs, MET_OP          rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator - ( long            lhs, MetricBase<T> & rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator - ( MET_OP          lhs, long            rhs ) { return create_sub<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( MetricBase<T> & lhs, long            rhs ) { return create_sub<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( unsigned int    lhs, MET_OP          rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator - ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator - ( MET_OP          lhs, unsigned int    rhs ) { return create_sub<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_sub<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( unsigned long   lhs, MET_OP          rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator - ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_sub<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator - ( MET_OP          lhs, unsigned long   rhs ) { return create_sub<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator - ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_sub<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> MET_OP create_mult(MET_OP lhs, MET_OP rhs) { return (MET_OP)( new Metric_Mult<T>( lhs, rhs ) ); }
template<class T> MET_OP operator * ( MET_OP          lhs, MET_OP          rhs ) { return create_mult<T>( lhs                               , rhs                                ); }
template<class T> MET_OP operator * ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_mult<T>( lhs.clone()                       , rhs                                ); }
template<class T> MET_OP operator * ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_mult<T>( lhs                               , rhs.clone()                        ); }
template<class T> MET_OP operator * ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_mult<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> MET_OP operator * ( double          lhs, MET_OP          rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator * ( double          lhs, MetricBase<T> & rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator * ( MET_OP          lhs, double          rhs ) { return create_mult<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( MetricBase<T> & lhs, double          rhs ) { return create_mult<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( float           lhs, MET_OP          rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator * ( float           lhs, MetricBase<T> & rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator * ( MET_OP          lhs, float           rhs ) { return create_mult<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( MetricBase<T> & lhs, float           rhs ) { return create_mult<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( int             lhs, MET_OP          rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator * ( int             lhs, MetricBase<T> & rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator * ( MET_OP          lhs, int             rhs ) { return create_mult<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( MetricBase<T> & lhs, int             rhs ) { return create_mult<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( long            lhs, MET_OP          rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator * ( long            lhs, MetricBase<T> & rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator * ( MET_OP          lhs, long            rhs ) { return create_mult<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( MetricBase<T> & lhs, long            rhs ) { return create_mult<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( unsigned int    lhs, MET_OP          rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator * ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator * ( MET_OP          lhs, unsigned int    rhs ) { return create_mult<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_mult<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( unsigned long   lhs, MET_OP          rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator * ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_mult<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator * ( MET_OP          lhs, unsigned long   rhs ) { return create_mult<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator * ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_mult<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> MET_OP create_div (MET_OP lhs, MET_OP rhs) { return (MET_OP)( new Metric_Div<T>( lhs, rhs ) ); }
template<class T> MET_OP operator / ( MET_OP          lhs, MET_OP          rhs ) { return create_div<T>( lhs                               , rhs                                ); }
template<class T> MET_OP operator / ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_div<T>( lhs.clone()                       , rhs                                ); }
template<class T> MET_OP operator / ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_div<T>( lhs                               , rhs.clone()                        ); }
template<class T> MET_OP operator / ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_div<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> MET_OP operator / ( double          lhs, MET_OP          rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator / ( double          lhs, MetricBase<T> & rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator / ( MET_OP          lhs, double          rhs ) { return create_div<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( MetricBase<T> & lhs, double          rhs ) { return create_div<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( float           lhs, MET_OP          rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator / ( float           lhs, MetricBase<T> & rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator / ( MET_OP          lhs, float           rhs ) { return create_div<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( MetricBase<T> & lhs, float           rhs ) { return create_div<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( int             lhs, MET_OP          rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator / ( int             lhs, MetricBase<T> & rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator / ( MET_OP          lhs, int             rhs ) { return create_div<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( MetricBase<T> & lhs, int             rhs ) { return create_div<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( long            lhs, MET_OP          rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator / ( long            lhs, MetricBase<T> & rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator / ( MET_OP          lhs, long            rhs ) { return create_div<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( MetricBase<T> & lhs, long            rhs ) { return create_div<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( unsigned int    lhs, MET_OP          rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator / ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator / ( MET_OP          lhs, unsigned int    rhs ) { return create_div<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_div<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( unsigned long   lhs, MET_OP          rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator / ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_div<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator / ( MET_OP          lhs, unsigned long   rhs ) { return create_div<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator / ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_div<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> MET_OP create_add (MET_OP lhs, MET_OP rhs) { return (MET_OP)( new Metric_Add<T>( lhs, rhs ) ); }
template<class T> MET_OP operator + ( MET_OP          lhs, MET_OP          rhs ) { return create_add<T>( lhs                               , rhs                                ); }
template<class T> MET_OP operator + ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_add<T>( lhs.clone()                       , rhs                                ); }
template<class T> MET_OP operator + ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_add<T>( lhs                               , rhs.clone()                        ); }
template<class T> MET_OP operator + ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_add<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> MET_OP operator + ( double          lhs, MET_OP          rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator + ( double          lhs, MetricBase<T> & rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator + ( MET_OP          lhs, double          rhs ) { return create_add<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( MetricBase<T> & lhs, double          rhs ) { return create_add<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( float           lhs, MET_OP          rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator + ( float           lhs, MetricBase<T> & rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator + ( MET_OP          lhs, float           rhs ) { return create_add<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( MetricBase<T> & lhs, float           rhs ) { return create_add<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( int             lhs, MET_OP          rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator + ( int             lhs, MetricBase<T> & rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator + ( MET_OP          lhs, int             rhs ) { return create_add<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( MetricBase<T> & lhs, int             rhs ) { return create_add<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( long            lhs, MET_OP          rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator + ( long            lhs, MetricBase<T> & rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator + ( MET_OP          lhs, long            rhs ) { return create_add<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( MetricBase<T> & lhs, long            rhs ) { return create_add<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( unsigned int    lhs, MET_OP          rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator + ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator + ( MET_OP          lhs, unsigned int    rhs ) { return create_add<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_add<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( unsigned long   lhs, MET_OP          rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator + ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_add<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator + ( MET_OP          lhs, unsigned long   rhs ) { return create_add<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator + ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_add<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> MET_OP create_pow (MET_OP lhs, MET_OP rhs) { return (MET_OP)( new Metric_Pow<T>( lhs, rhs ) ); }
template<class T> MET_OP operator ^ ( MET_OP          lhs, MET_OP          rhs ) { return create_pow<T>( lhs                               , rhs                                ); }
template<class T> MET_OP operator ^ ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_pow<T>( lhs.clone()                       , rhs                                ); }
template<class T> MET_OP operator ^ ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_pow<T>( lhs                               , rhs.clone()                        ); }
template<class T> MET_OP operator ^ ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_pow<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> MET_OP operator ^ ( double          lhs, MET_OP          rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator ^ ( double          lhs, MetricBase<T> & rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator ^ ( MET_OP          lhs, double          rhs ) { return create_pow<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( MetricBase<T> & lhs, double          rhs ) { return create_pow<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( float           lhs, MET_OP          rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator ^ ( float           lhs, MetricBase<T> & rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator ^ ( MET_OP          lhs, float           rhs ) { return create_pow<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( MetricBase<T> & lhs, float           rhs ) { return create_pow<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( int             lhs, MET_OP          rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator ^ ( int             lhs, MetricBase<T> & rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator ^ ( MET_OP          lhs, int             rhs ) { return create_pow<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( MetricBase<T> & lhs, int             rhs ) { return create_pow<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( long            lhs, MET_OP          rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator ^ ( long            lhs, MetricBase<T> & rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator ^ ( MET_OP          lhs, long            rhs ) { return create_pow<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( MetricBase<T> & lhs, long            rhs ) { return create_pow<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( unsigned int    lhs, MET_OP          rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator ^ ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator ^ ( MET_OP          lhs, unsigned int    rhs ) { return create_pow<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_pow<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( unsigned long   lhs, MET_OP          rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> MET_OP operator ^ ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_pow<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> MET_OP operator ^ ( MET_OP          lhs, unsigned long   rhs ) { return create_pow<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> MET_OP operator ^ ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_pow<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> MET_OP create_neg (MET_OP rhs) { return (MET_OP)( new Metric_Neg <T>( rhs ) );}
template<class T> MET_OP operator - ( MET_OP          rhs ) { return create_neg <T>(rhs); }
template<class T> MET_OP operator - ( MetricBase<T> & rhs ) { return create_neg<T>(rhs.clone()); }

// boolean conversions
template<class T> FILT_OP create_eq (MET_OP lhs, MET_OP rhs) { return (FILT_OP)( new MetricInequality<T>( EQ, lhs, rhs ) ); }
template<class T> FILT_OP operator == ( MET_OP          lhs, MET_OP          rhs ) { return create_eq<T>( lhs                               , rhs                                ); }
template<class T> FILT_OP operator == ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_eq<T>( lhs.clone()                       , rhs                                ); }
template<class T> FILT_OP operator == ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_eq<T>( lhs                               , rhs.clone()                        ); }
template<class T> FILT_OP operator == ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_eq<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> FILT_OP operator == ( double          lhs, MET_OP          rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator == ( double          lhs, MetricBase<T> & rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator == ( MET_OP          lhs, double          rhs ) { return create_eq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( MetricBase<T> & lhs, double          rhs ) { return create_eq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( float           lhs, MET_OP          rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator == ( float           lhs, MetricBase<T> & rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator == ( MET_OP          lhs, float           rhs ) { return create_eq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( MetricBase<T> & lhs, float           rhs ) { return create_eq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( int             lhs, MET_OP          rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator == ( int             lhs, MetricBase<T> & rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator == ( MET_OP          lhs, int             rhs ) { return create_eq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( MetricBase<T> & lhs, int             rhs ) { return create_eq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( long            lhs, MET_OP          rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator == ( long            lhs, MetricBase<T> & rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator == ( MET_OP          lhs, long            rhs ) { return create_eq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( MetricBase<T> & lhs, long            rhs ) { return create_eq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( unsigned int    lhs, MET_OP          rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator == ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator == ( MET_OP          lhs, unsigned int    rhs ) { return create_eq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_eq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( unsigned long   lhs, MET_OP          rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator == ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_eq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator == ( MET_OP          lhs, unsigned long   rhs ) { return create_eq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator == ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_eq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> FILT_OP create_neq(MET_OP lhs, MET_OP rhs) { return (FILT_OP)( new MetricInequality<T>( NEQ, lhs, rhs ) ); }
template<class T> FILT_OP operator != ( MET_OP          lhs, MET_OP          rhs ) { return create_neq<T>( lhs                               , rhs                                ); }
template<class T> FILT_OP operator != ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_neq<T>( lhs.clone()                       , rhs                                ); }
template<class T> FILT_OP operator != ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_neq<T>( lhs                               , rhs.clone()                        ); }
template<class T> FILT_OP operator != ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_neq<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> FILT_OP operator != ( double          lhs, MET_OP          rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator != ( double          lhs, MetricBase<T> & rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator != ( MET_OP          lhs, double          rhs ) { return create_neq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( MetricBase<T> & lhs, double          rhs ) { return create_neq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( float           lhs, MET_OP          rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator != ( float           lhs, MetricBase<T> & rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator != ( MET_OP          lhs, float           rhs ) { return create_neq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( MetricBase<T> & lhs, float           rhs ) { return create_neq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( int             lhs, MET_OP          rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator != ( int             lhs, MetricBase<T> & rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator != ( MET_OP          lhs, int             rhs ) { return create_neq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( MetricBase<T> & lhs, int             rhs ) { return create_neq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( long            lhs, MET_OP          rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator != ( long            lhs, MetricBase<T> & rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator != ( MET_OP          lhs, long            rhs ) { return create_neq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( MetricBase<T> & lhs, long            rhs ) { return create_neq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( unsigned int    lhs, MET_OP          rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator != ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator != ( MET_OP          lhs, unsigned int    rhs ) { return create_neq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_neq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( unsigned long   lhs, MET_OP          rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator != ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_neq<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator != ( MET_OP          lhs, unsigned long   rhs ) { return create_neq<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator != ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_neq<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> FILT_OP create_gt(MET_OP lhs, MET_OP rhs) { return (FILT_OP)( new MetricInequality<T>( GT, lhs, rhs ) ); }
template<class T> FILT_OP operator > ( MET_OP          lhs, MET_OP          rhs ) { return create_gt<T>( lhs                               , rhs                                ); }
template<class T> FILT_OP operator > ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_gt<T>( lhs.clone()                       , rhs                                ); }
template<class T> FILT_OP operator > ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_gt<T>( lhs                               , rhs.clone()                        ); }
template<class T> FILT_OP operator > ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_gt<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> FILT_OP operator > ( double          lhs, MET_OP          rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator > ( double          lhs, MetricBase<T> & rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator > ( MET_OP          lhs, double          rhs ) { return create_gt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( MetricBase<T> & lhs, double          rhs ) { return create_gt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( float           lhs, MET_OP          rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator > ( float           lhs, MetricBase<T> & rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator > ( MET_OP          lhs, float           rhs ) { return create_gt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( MetricBase<T> & lhs, float           rhs ) { return create_gt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( int             lhs, MET_OP          rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator > ( int             lhs, MetricBase<T> & rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator > ( MET_OP          lhs, int             rhs ) { return create_gt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( MetricBase<T> & lhs, int             rhs ) { return create_gt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( long            lhs, MET_OP          rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator > ( long            lhs, MetricBase<T> & rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator > ( MET_OP          lhs, long            rhs ) { return create_gt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( MetricBase<T> & lhs, long            rhs ) { return create_gt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( unsigned int    lhs, MET_OP          rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator > ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator > ( MET_OP          lhs, unsigned int    rhs ) { return create_gt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_gt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( unsigned long   lhs, MET_OP          rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator > ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_gt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator > ( MET_OP          lhs, unsigned long   rhs ) { return create_gt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator > ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_gt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> FILT_OP create_gte(MET_OP lhs, MET_OP rhs) { return (FILT_OP)( new MetricInequality<T>( GTE, lhs, rhs ) ); }
template<class T> FILT_OP operator >= ( MET_OP          lhs, MET_OP          rhs ) { return create_gte<T>( lhs                               , rhs                                ); }
template<class T> FILT_OP operator >= ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_gte<T>( lhs.clone()                       , rhs                                ); }
template<class T> FILT_OP operator >= ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_gte<T>( lhs                               , rhs.clone()                        ); }
template<class T> FILT_OP operator >= ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_gte<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> FILT_OP operator >= ( double          lhs, MET_OP          rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator >= ( double          lhs, MetricBase<T> & rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator >= ( MET_OP          lhs, double          rhs ) { return create_gte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( MetricBase<T> & lhs, double          rhs ) { return create_gte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( float           lhs, MET_OP          rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator >= ( float           lhs, MetricBase<T> & rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator >= ( MET_OP          lhs, float           rhs ) { return create_gte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( MetricBase<T> & lhs, float           rhs ) { return create_gte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( int             lhs, MET_OP          rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator >= ( int             lhs, MetricBase<T> & rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator >= ( MET_OP          lhs, int             rhs ) { return create_gte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( MetricBase<T> & lhs, int             rhs ) { return create_gte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( long            lhs, MET_OP          rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator >= ( long            lhs, MetricBase<T> & rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator >= ( MET_OP          lhs, long            rhs ) { return create_gte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( MetricBase<T> & lhs, long            rhs ) { return create_gte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( unsigned int    lhs, MET_OP          rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator >= ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator >= ( MET_OP          lhs, unsigned int    rhs ) { return create_gte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_gte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( unsigned long   lhs, MET_OP          rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator >= ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_gte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator >= ( MET_OP          lhs, unsigned long   rhs ) { return create_gte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator >= ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_gte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> FILT_OP create_lt(MET_OP lhs, MET_OP rhs) { return (FILT_OP)( new MetricInequality<T>( LT, lhs, rhs ) ); }
template<class T> FILT_OP operator < ( MET_OP          lhs, MET_OP          rhs ) { return create_lt<T>( lhs                               , rhs                                ); }
template<class T> FILT_OP operator < ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_lt<T>( lhs.clone()                       , rhs                                ); }
template<class T> FILT_OP operator < ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_lt<T>( lhs                               , rhs.clone()                        ); }
template<class T> FILT_OP operator < ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_lt<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> FILT_OP operator < ( double          lhs, MET_OP          rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator < ( double          lhs, MetricBase<T> & rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator < ( MET_OP          lhs, double          rhs ) { return create_lt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( MetricBase<T> & lhs, double          rhs ) { return create_lt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( float           lhs, MET_OP          rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator < ( float           lhs, MetricBase<T> & rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator < ( MET_OP          lhs, float           rhs ) { return create_lt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( MetricBase<T> & lhs, float           rhs ) { return create_lt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( int             lhs, MET_OP          rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator < ( int             lhs, MetricBase<T> & rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator < ( MET_OP          lhs, int             rhs ) { return create_lt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( MetricBase<T> & lhs, int             rhs ) { return create_lt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( long            lhs, MET_OP          rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator < ( long            lhs, MetricBase<T> & rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator < ( MET_OP          lhs, long            rhs ) { return create_lt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( MetricBase<T> & lhs, long            rhs ) { return create_lt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( unsigned int    lhs, MET_OP          rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator < ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator < ( MET_OP          lhs, unsigned int    rhs ) { return create_lt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_lt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( unsigned long   lhs, MET_OP          rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator < ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_lt<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator < ( MET_OP          lhs, unsigned long   rhs ) { return create_lt<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator < ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_lt<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }

template<class T> FILT_OP create_lte(MET_OP lhs, MET_OP rhs) { return (FILT_OP)( new MetricInequality<T>( LTE, lhs, rhs ) ); }
template<class T> FILT_OP operator <= ( MET_OP          lhs, MET_OP          rhs ) { return create_lte<T>( lhs                               , rhs                                ); }
template<class T> FILT_OP operator <= ( MetricBase<T> & lhs, MET_OP          rhs ) { return create_lte<T>( lhs.clone()                       , rhs                                ); }
template<class T> FILT_OP operator <= ( MET_OP          lhs, MetricBase<T> & rhs ) { return create_lte<T>( lhs                               , rhs.clone()                        ); }
template<class T> FILT_OP operator <= ( MetricBase<T> & lhs, MetricBase<T> & rhs ) { return create_lte<T>( lhs.clone()                       , rhs.clone()                        ); }
template<class T> FILT_OP operator <= ( double          lhs, MET_OP          rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator <= ( double          lhs, MetricBase<T> & rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator <= ( MET_OP          lhs, double          rhs ) { return create_lte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( MetricBase<T> & lhs, double          rhs ) { return create_lte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( float           lhs, MET_OP          rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator <= ( float           lhs, MetricBase<T> & rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator <= ( MET_OP          lhs, float           rhs ) { return create_lte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( MetricBase<T> & lhs, float           rhs ) { return create_lte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( int             lhs, MET_OP          rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator <= ( int             lhs, MetricBase<T> & rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator <= ( MET_OP          lhs, int             rhs ) { return create_lte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( MetricBase<T> & lhs, int             rhs ) { return create_lte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( long            lhs, MET_OP          rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator <= ( long            lhs, MetricBase<T> & rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator <= ( MET_OP          lhs, long            rhs ) { return create_lte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( MetricBase<T> & lhs, long            rhs ) { return create_lte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( unsigned int    lhs, MET_OP          rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator <= ( unsigned int    lhs, MetricBase<T> & rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator <= ( MET_OP          lhs, unsigned int    rhs ) { return create_lte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( MetricBase<T> & lhs, unsigned int    rhs ) { return create_lte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( unsigned long   lhs, MET_OP          rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs                                ); }
template<class T> FILT_OP operator <= ( unsigned long   lhs, MetricBase<T> & rhs ) { return create_lte<T>( MET_OP(new MetricConstant<T>(lhs)), rhs.clone()                        ); }
template<class T> FILT_OP operator <= ( MET_OP          lhs, unsigned long   rhs ) { return create_lte<T>( lhs                               , MET_OP(new MetricConstant<T>(rhs)) ); }
template<class T> FILT_OP operator <= ( MetricBase<T> & lhs, unsigned long   rhs ) { return create_lte<T>( lhs.clone()                       , MET_OP(new MetricConstant<T>(rhs)) ); }
}
using namespace metric_operators;




} // end namespace query
} // end namespace utility




#endif

