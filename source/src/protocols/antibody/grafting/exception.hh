// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/exception.hh
/// @brief Exceptions classes for antibody grafting code
/// @author Sergey Lyskov


#ifndef INCLUDED_protocols_antibody_grafting_exception_hh
#define INCLUDED_protocols_antibody_grafting_exception_hh


#include <utility/excn/Exceptions.hh>

namespace protocols {
namespace antibody {
namespace grafting {


class Grafting_Base_Exception : public utility::excn::Exception
{
public:
	using utility::excn::Exception::Exception;
};


class _AE_grafting_failed_ : public Grafting_Base_Exception
{
public:
	using Grafting_Base_Exception::Grafting_Base_Exception;
};



class _AE_cdr_detection_failed_ : public Grafting_Base_Exception
{
public:
	using Grafting_Base_Exception::Grafting_Base_Exception;
};

class _AE_json_cdr_detection_failed_ : public Grafting_Base_Exception
{
public:
	using Grafting_Base_Exception::Grafting_Base_Exception;
};


class _AE_cdr_undefined_ : public Grafting_Base_Exception
{
public:
	using Grafting_Base_Exception::Grafting_Base_Exception;
};


class _AE_invalid_cdr_region_ : public Grafting_Base_Exception
{
public:
	using Grafting_Base_Exception::Grafting_Base_Exception;
};



class _AE_numbering_failed_ : public Grafting_Base_Exception
{
public:
	using Grafting_Base_Exception::Grafting_Base_Exception;
};


class _AE_unexpected_region_length_ : public _AE_numbering_failed_
{
public:
	using _AE_numbering_failed_::_AE_numbering_failed_;
};

class _AE_scs_failed_ : public Grafting_Base_Exception
{
public:
	using Grafting_Base_Exception::Grafting_Base_Exception;
};



} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // INCLUDED_protocols_antibody_grafting_exception_hh
