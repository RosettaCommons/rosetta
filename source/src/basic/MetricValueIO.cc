// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Colin A. Smith

#include <basic/MetricValueIO.hh>

#include <basic/MetricValue.hh>

#include <string>
#include <sstream>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <ostream>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/Tracer.hh>


using utility::vector1;

namespace basic {

template <class T>
bool
write_metric_value_scalar(
	std::ostream & os,
	MetricValue<T> const & metric_value,
	std::string const & type_name
)
{
	os << type_name << " " << metric_value.value();

	return true;
}

template <class T>
bool
write_metric_value_vector1(
	std::ostream & os,
	MetricValue<vector1<T> > const & metric_value,
	std::string const & type_name
)
{
	vector1<T> const & vec(metric_value.value());

	os << type_name << "[";
	for ( typename vector1<T>::const_iterator iter(vec.begin()), iter_end(vec.end()); iter != iter_end; ++iter ) {
		os << ' ' << *iter;
	}
	os << " ]";

	return true;
}

/// @brief central write introspection and delegation function
bool
write_metric_value(
	std::ostream & os,
	MetricValueBase const & metric_value
)
{
	MetricValue<double> const * const metric_value_real( dynamic_cast<MetricValue<double> const * > (&metric_value) );
	if ( metric_value_real ) return write_metric_value_scalar(os, *metric_value_real, "Real");

	MetricValue<int> const * const metric_value_int( dynamic_cast<MetricValue<int> const * > (&metric_value) );
	if ( metric_value_int ) return write_metric_value_scalar(os, *metric_value_int, "Int");

	MetricValue<size_t> const * const metric_value_size( dynamic_cast<MetricValue<size_t> const * > (&metric_value) );
	if ( metric_value_size ) return write_metric_value_scalar(os, *metric_value_size, "Size");

	MetricValue<bool> const * const metric_value_bool( dynamic_cast<MetricValue<bool> const * > (&metric_value) );
	if ( metric_value_bool ) return write_metric_value_scalar(os, *metric_value_bool, "Bool");

	MetricValue<vector1<double> > const * const metric_value_vector_real(
		dynamic_cast<MetricValue<vector1<double> > const * > (&metric_value)
	);
	if ( metric_value_vector_real ) return write_metric_value_vector1(os, *metric_value_vector_real, "Real");

	MetricValue<vector1<int> > const * const metric_value_vector_int(
		dynamic_cast<MetricValue<vector1<int> > const * > (&metric_value)
	);
	if ( metric_value_vector_int ) return write_metric_value_vector1(os, *metric_value_vector_int, "Int");

	MetricValue<vector1<size_t> > const * const metric_value_vector_size(
		dynamic_cast<MetricValue<vector1<size_t> > const * > (&metric_value)
	);
	if ( metric_value_vector_size ) return write_metric_value_vector1(os, *metric_value_vector_size, "Size");

	MetricValue<vector1<bool> > const * const metric_value_vector_bool(
		dynamic_cast<MetricValue<vector1<bool> > const * > (&metric_value)
	);
	if ( metric_value_vector_bool ) return write_metric_value_vector1(os, *metric_value_vector_bool, "Bool");

	return false;
}

bool
handles_metric_value(
	MetricValueBase const & metric_value
)
{
	// this isn't the fastest but avoids duplication
	std::ostringstream oss;
	return write_metric_value(oss, metric_value);
}

template <class T>
bool
read_metric_value_scalar(
	std::istream & is,
	MetricValueBase & metric_value
)
{
	MetricValue<T> * const metric_value_typed( dynamic_cast<MetricValue<T> * > (&metric_value) );
	if ( !metric_value_typed ) return false;

	T value;
	if ( !(is >> value) ) return false;

	metric_value_typed->set(value);

	return true;
}

template <class T>
bool
read_metric_value_scalar(
	std::istream & is,
	MetricValueBaseOP & metric_value
)
{
	metric_value = MetricValueBaseOP( new MetricValue<T> );

	return read_metric_value_scalar<T>(is, *metric_value);
}

template <class T>
bool
read_metric_value_vector1(
	std::istream & is,
	MetricValueBase & metric_value
)
{
	MetricValue<vector1<T> > * const metric_value_typed( dynamic_cast<MetricValue<vector1<T> > * > (&metric_value) );
	if ( !metric_value_typed ) return false;

	vector1<T> vec;

	std::string word;
	while ( is >> word ) {
		if ( word == "]" ) {
			metric_value_typed->set(vec);
			return true;
		}
		std::istringstream iss(word);
		T value;
		if ( !(iss >> value) ) return false;
		vec.push_back(value);
	}

	return false;
}

template <class T>
bool
read_metric_value_vector1(
	std::istream & is,
	MetricValueBaseOP & metric_value
)
{
	metric_value = MetricValueBaseOP( new MetricValue<vector1<T> > );

	return read_metric_value_vector1<T>(is, *metric_value);
}

/// @brief central read introspection and delegation function
template <class MetricValueBase_or_MetricValueBaseOP>
bool
read_metric_value_template(
	std::istream & is,
	MetricValueBase_or_MetricValueBaseOP & metric_value
)
{
	std::string word;
	if ( !(is >> word) ) return false;

	if ( word == "Real" ) return read_metric_value_scalar<double>(is, metric_value);
	if ( word == "Int" ) return read_metric_value_scalar<int>(is, metric_value);
	if ( word == "Size" ) return read_metric_value_scalar<size_t>(is, metric_value);
	if ( word == "Bool" ) return read_metric_value_scalar<bool>(is, metric_value);

	if ( word == "Real[" ) return read_metric_value_vector1<double>(is, metric_value);
	if ( word == "Int[" ) return read_metric_value_vector1<int>(is, metric_value);
	if ( word == "Size[" ) return read_metric_value_vector1<size_t>(is, metric_value);
	if ( word == "Bool[" ) return read_metric_value_vector1<bool>(is, metric_value);

	return false;
}

bool
read_metric_value(
	std::istream & is,
	MetricValueBase & metric_value
)
{
	return read_metric_value_template(is, metric_value);
}

MetricValueBaseOP
read_metric_value(
	std::istream & is
)
{
	MetricValueBaseOP metric_value;

	if ( read_metric_value_template(is, metric_value) ) return metric_value;

	return NULL;
}


} // namespace basic

