// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/options/PyRosettaFunctions.hh
/// @brief  Additional functions to set/get Options from PyRosetta
/// @author Sergey Lyskov

#ifndef INCLUDED_basic_options_PyRosettaFunctions_hh
#define INCLUDED_basic_options_PyRosettaFunctions_hh

#include <iosfwd>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>


namespace basic {
namespace options {

bool        get_boolean_option(std::string const & id);
int         get_integer_option(std::string const & id);
double      get_real_option(std::string const & id);
std::string get_string_option(std::string const & id);
std::string get_file_option(std::string const & id);

void set_boolean_option(std::string const & id, bool );
void set_integer_option(std::string const & id, int         );
void set_real_option(std::string const & id, double      );
void set_string_option(std::string const & id, std::string const &);
void set_file_option(std::string const & id, std::string const &);


utility::vector1<bool>        get_boolean_vector_option(std::string const & id);
utility::vector1<int>         get_integer_vector_option(std::string const & id);
utility::vector1<double>      get_real_vector_option(std::string const & id);
utility::vector1<std::string> get_string_vector_option(std::string const & id);
utility::vector1<std::string> get_file_vector_option(std::string const & id);

void set_boolean_vector_option(std::string const & id, utility::vector1<bool> const &);
void set_integer_vector_option(std::string const & id, utility::vector1<int> const &);
void set_real_vector_option(std::string const & id, utility::vector1<double> const &);
void set_string_vector_option(std::string const & id, utility::vector1<std::string> const &);
void set_file_vector_option(std::string const & id, utility::vector1<std::string> const &);

} // namespace options
} // namespace basic


#endif  // INCLUDED_basic_options_PyRosettaFunctions_hh
