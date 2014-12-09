// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/UTracer.hh
/// @brief  Unit test Tracer/Diff system
/// @author Sergey Lyskov
///

#ifndef INCLUDED_UTracer_HH
#define INCLUDED_UTracer_HH

#include <test/UTools.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
//#include <utility/stream_util.hh>
//#include <utility/vector1.hh>
//#include <utility/exit.hh>

#include <map>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace test {

/// Helper macro for "standard" usage.
#define UTRACE UT.file(__FILE__).line(__LINE__)

/// @brief special verison of Tracer that record and compare io data to the given file.
class UTracer : public basic::otstream
{
public:
	UTracer(std::string const & file_name) :
	  u_file_name_(file_name),
	  file_original_(file_name.c_str(), std::fstream::in),
	  file_new_((file_name+"._tmp_").c_str(), std::fstream::out),
	  abs_tolerance_(1e-200), rel_tolerance_(1e-200),
	  error_(false) //error_signal(false), ufile_line_(-1),
	{
		static char const * nfi = "NoFileInfo";
		line_ = 0;
		file_ = nfi;

		if( !file_original_.good() ) {
			std::string m( "Unable to open file: " + file_name + " for reading!\nPlease verify that it has been added to the 'testinputfiles' block in appropriate rosetta_source/test/<LIBRARY>.test.settings file\n" );
			TS_FAIL( m );
			error_ = true;
		}

		if( !file_new_.good() ) {
			std::string m("Unable to open file: " + file_name + ".new for writing!\n" );
			TS_FAIL( m );
			error_ = true;
		}
	};

	virtual ~UTracer() {
		// We want to insure that last symbol in output is a new line,
		// and that inner buffer is flushed.
		(*this) << std::endl;
		file_new_.close();

		if( ! error_ ) {

			file_new_.open((u_file_name_+"._tmp_").c_str(), std::fstream::in);

			int ufile_line = 1;  // keeping track of line number for generating error message.

			std::string original_line, new_line;
			for(; getline(file_new_, new_line); ufile_line++) {
				if( !getline(file_original_, original_line) ) { ///< no lines remain in the buffer
					std::ostringstream buf;
					buf << "(UTrace: " << u_file_name_ << ":" << ufile_line <<
							 ") Not enough lines in the file to comapre with: " << new_line << std::endl;
					_TS_FAIL(file_, line_, buf.str() );
					error_ = true;
				}

				std::string error_message;
				if( !test::utools::isEq(new_line, original_line, abs_tolerance_, rel_tolerance_, error_message)
				   && (!error_) ) {
					std::ostringstream msg;

					std::cerr << "UTracer(" << u_file_name_ << ") line " << ufile_line << " not equal:" << std::endl
					    << error_message << std::endl
						<< "old: " << original_line << std::endl
						<< "new: " << new_line;

					//_TS_FAIL(file_, line_, msg.str() );
					//TS_FAIL(msg.str() );
					TS_FAIL( error_message );
					std::cerr << std::endl;

					error_ = true;
				}
				/*if( (new_line != original_line) && (!error_) ) {
					std::ostringstream buf;  buf << ufile_line;
					std::string err = "UTracer("+u_file_name_+") line " + buf.str() + " not equal:\n" +
					"  " + original_line + "\n" +
					"  " + new_line + "\n";

					_TS_FAIL(file_, line_, err.c_str() );

					//error_signal = false;
					error_ = true;
				}
				*/
			}
		}

		if( ! error_ ) {
			std::string line;
			if( getline(file_original_, line) ) { ///< Some lines remain in the buffer

				line = "(UTrace: "+ u_file_name_ + ") Extra line in the end of a file to comapre with: " + line + "\n";
				_TS_FAIL(file_, line_, line.c_str() );
			}
		}

		file_original_.close();
	};

	/// Set line file/line number for the next output
	UTracer & file(const char * file) { file_= (char *)file; return *this; };

	UTracer & line(int line) {
		line_ = line;
		return *this;
	};

	/// Set absolute delta for real types asserts
	//UTracer & abs_tolerance(double d) {	abs_tolerance_ = d; return *this; };
	UTracer & abs_tolerance(double d) {	(*this) << "set_abs_tolerance(" <<  d << ")"; return *this; };

	/// Set relative delta for real types asserts
	//UTracer & rel_tolerance(double r) { rel_tolerance_ = r; return *this; };
	UTracer & rel_tolerance(double r) {	(*this) << "set_rel_tolerance(" <<  r << ")"; return *this; };

protected:
	/// @brief overload member function.
	virtual void t_flush(std::string const & s) {
		file_new_ << s;
	}

private:
	std::string u_file_name_;
	std::fstream file_original_;
	std::fstream file_new_;

	/// Source file name/line number of current output. We need this to generate informative error messages.
	char const * file_;
	int line_;

	/// @brief delta for real asserts, absolute precision
	double abs_tolerance_;

	/// @brief relative for real asserts, relative precision
	double rel_tolerance_;

	/// @brief flag indicating that error already happened and reported,
	///        and futher errors should not be reported.
	bool error_;

	/// @brief flag indicating that error happend at current line.
	//bool error_signal_;

	/// @brief buffers that holds original and new lines, we use it for error reporting.
	//std::string original_line_;
	//std::ostringstream new_line_;

	/// @brief number of line in utracer file.
	//int ufile_line_;

	/// Link to file tracer
	//basic::Tracer * tr_;
}; // UTracer


} // namespace test

#endif // INCLUDED_test_utracer_hh
