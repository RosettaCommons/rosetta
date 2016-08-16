// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/Emitter.hh
///
/// @brief  Lightweight class to ease writting YAML documents
/// @author Ian W. Davis

#ifndef INCLUDED_basic_Emitter_hh
#define INCLUDED_basic_Emitter_hh

#include <basic/Tracer.hh>

#include <cstddef>
#include <iosfwd>
#include <ostream>


namespace basic {

class Emitter; // fwd declaration
typedef utility::pointer::shared_ptr< Emitter > EmitterOP;
typedef utility::pointer::shared_ptr< Emitter const > EmitterCOP;

class YamlEmitter; // fwd declaration
typedef utility::pointer::shared_ptr< YamlEmitter > YamlEmitterOP;
typedef utility::pointer::shared_ptr< YamlEmitter const > YamlEmitterCOP;

class JsonEmitter; // fwd declaration
typedef utility::pointer::shared_ptr< JsonEmitter > JsonEmitterOP;
typedef utility::pointer::shared_ptr< JsonEmitter const > JsonEmitterCOP;

/// @brief Lightweight class to ease writting YAML documents (http://yaml.org)
/// @details YAML is a structured data format, a more human-readable
///  alternative to XML.
///
///  I use the terms "map" and "list" here, but you may also see "mapping" and "sequence".
///  The former is a series of key : value pairs, while the latter is a simple series of items.
///  This class is not entirely able to enforce the generation of valid YAML --
///  you can still e.g. write a key : value pair in the middle of a list.
///  It will print a warning message but will otherwise try to continue on.
///
///  YAML documents have optional explicit start/end markers; if the emitter
///  supports them, they will be auto-generated when the class is initialized
///  and when you're finished and you call end().
///
///  Whitespace is YAML documents is semi-significant.  By default, this class
///  will try to pretty-print, showing the depth of nesting through indentation.
///  When starting a new list or map, you can request that it not be indented
///  (i.e., be printed all on one line), but all lists and maps nested inside it
///  will also be printed without linebreaks, regardless of the requested indentation.
///  YAML refers to these two styles as "block" and "flow", respectively.
class Emitter : public utility::pointer::ReferenceCount
{
public:

	Emitter(std::ostream & out):
		out_(out),
		indent_str_(" "),
		in_map_(),
		first_(),
		indent_(),
		indent_depth_()
	{}

	virtual ~Emitter() {}

	/// @brief Flush the underlying output stream.
	void flush() { out_.flush(); }

	/// @brief Write method for use inside lists / arrays
	template <typename T>
	void write(T const & data);
	/// @brief Write method for use inside lists / arrays
	void start_list(bool indent=true);
	/// @brief Write method for use inside lists / arrays
	void start_map(bool indent=true);

	/// @brief Write method for use inside maps / dicts / objects
	template <typename T>
	void write(std::string const & label, T const & data);

	// The const char* variants are required here so that string literals
	// won't be cast to booleans (resulting in calls to the list methods above!).

	/// @brief Write method for use inside maps / dicts / objects
	void start_list(std::string const & label, bool indent=true)
	{ start_list(label.c_str(), indent); }
	/// @brief Write method for use inside maps / dicts / objects
	void start_list(const char * label, bool indent=true);
	/// @brief Write method for use inside maps / dicts / objects
	void start_map(std::string const & label, bool indent=true)
	{ start_map(label.c_str(), indent); }
	/// @brief Write method for use inside maps / dicts / objects
	void start_map(const char * label, bool indent=true);

	/// @brief Counterpart to start_list() -- writes closing bracket.
	void end_list();
	/// @brief Counterpart to start_map() -- writes closing brace.
	void end_map();
	/// @brief By default, closes all open maps/lists, ending the document.
	void end(size_t desired_depth=0);
	/// @brief Number of closing brackets/braces required to end document (non-negative)
	int depth() const { return in_map_.size(); }
	/// @brief Number of spaces used for indenting.  Default is one.
	void set_indent(int num_spaces)
	{
		std::ostringstream s;
		for ( int i = 0; i < num_spaces; ++i ) s << ' ';
		indent_str_ = s.str();
	}

	/// @brief Start a new document, ending the previous one first if necessary.
	/// A no-op for some types of output (e.g. JSON) that don't allow multiple documents.
	virtual void start_doc() = 0;

protected:

	Emitter(); // no null ctor
	Emitter(Emitter const &); // no copy ctor

	/// @brief Check that we're in the expected context (either map or list)
	bool assert_in(bool in_map, std::string const & msg)
	{
		if ( in_map_.empty() || in_map_.back() != in_map ) {
			basic::Warning() << "Bad YAML: " << msg << "\n";
			return false;
		} else return true;
	}

	/// @brief Write the key part of a key-value pair.
	void write_label(std::string const & label)
	{
		// Space is required after the colon for it to also be legal YAML
		out_ << quote_string(label) << ": ";
	}
	void write_raw(bool b)
	{ out_ << (b ? "true" : "false"); }
	void write_raw(int i)
	{ out_ << i; }
	void write_raw(long l)
	{ out_ << l; }
	void write_raw(float f)
	{ out_ << f; }
	void write_raw(double d)
	{ out_ << d; }
	void write_raw(std::string const & s)
	{ out_ << quote_string(s); }

	/// @brief Converts special characters (newlines, etc) to escape codes (\n, etc).
	std::string escape_string(std::string const & s, bool & needs_quotes_out);

	/// @brief Quotes strings as much as needed for this format (e.g. always for JSON).
	virtual
	std::string quote_string(std::string const & s) = 0;

	/// @brief Handle pretty-printing indentation.
	/// Don't want to use commas for opening/closing brace/bracket.
	virtual
	void do_indent(bool write_comma=true) = 0;

	/// @brief Used for traditional YAML output, writes the "-" marker.
	virtual
	void write_list_marker() = 0;

	/// @brief Actual implementation of start_map() and start_list().
	virtual
	void start_raw(bool is_map, bool indent) = 0;

	/// @brief Actual implementation of end_map() and end_list().
	virtual
	void end_raw() = 0;

protected:
	std::ostream & out_;
	std::string indent_str_;

	// stacks (pushed down by start_list/map)
	std::vector<bool> in_map_;      //< if false, we're in a list, not a map
	std::vector<bool> first_;       //< first element in this list/map?
	std::vector<bool> indent_;      //< if false, print elements on one line
	std::vector<int> indent_depth_; //< depth of indent (if being used)
}; // Emitter


// Templated functions must go in header file:
template <typename T>
void Emitter::write(T const & data)
{
	assert_in(false, "Tried to write list data inside a map");
	do_indent();
	write_list_marker();
	write_raw(data);
}

template <typename T>
void Emitter::write(std::string const & label, T const & data)
{
	assert_in(true, "Tried to write map data inside a list");
	do_indent();
	write_label(label);
	write_raw(data);
}


//////////////////////////////////////////////////////////////////////////////


/// @brief Emitter for more classically formatted YAML
class YamlEmitter : public Emitter
{
public:

	YamlEmitter(std::ostream & out):
		Emitter(out)
	{ start_raw(true, true); } // can't call virtual func from base class ctor

	virtual ~YamlEmitter() {}

	/// @brief Start a new YAML document, ending the previous one first if necessary
	virtual void start_doc()
	{ end(); out_ << "\n---\n"; start_raw(true, true); }

private:

	YamlEmitter(); // no null ctor
	YamlEmitter(YamlEmitter const &); // no copy ctor

protected:

	/// @brief YAML only quotes strings when they contain special characters.
	virtual
	std::string quote_string(std::string const & s)
	{
		bool needs_quotes(false); // dummy val; will be overwritten below
		std::string t = escape_string(s, needs_quotes);
		if ( needs_quotes ) return "\""+t+"\"";
		else return t;
	}

	virtual
	void do_indent(bool write_comma=true)
	{
		bool indent = indent_.back();
		if ( indent ) {
			out_ << "\n";
			int depth = (indent_depth_.empty() ? 0 : indent_depth_.back());
			for ( int i = 0; i < depth; ++i ) out_ << indent_str_;
		} else {
			if ( write_comma ) {
				bool first = first_.back();
				if ( first ) first_[ first_.size()-1 ] = false;
				else out_ << ",";
			}
			out_ << " ";
		}
	}

	/// @brief YAML uses "-" for list items when in block (indented) mode.
	virtual
	void write_list_marker()
	{
		bool indent = indent_.back();
		if ( indent ) out_ << "- ";
	}

	/// @details YAML only uses brackets and braces if data is not being indented.
	virtual
	void start_raw(bool is_map, bool indent)
	{
		if ( !indent ) {
			if ( is_map ) out_ << "{ ";
			else out_ << "[ ";
		}
		in_map_.push_back(is_map);
		first_.push_back(true);
		// Once you start not indenting, no children can be indented
		if ( indent_.empty() ) indent_.push_back(indent);
		else indent_.push_back( indent_.back() & indent );
		if ( indent_depth_.empty() ) indent_depth_.push_back(1);
		else indent_depth_.push_back(indent_depth_.back() + 1);
	}

	/// @details YAML only uses brackets and braces if data is not being indented.
	virtual
	void end_raw()
	{
		if ( in_map_.empty() ) return; // bad op
		bool indent = indent_.back();
		bool is_map = in_map_.back();
		// This is popped first to get the closing brace/bracket to line up right.
		// Other pops go after do_indent() or weird stuff could happen.
		indent_depth_.pop_back();
		if ( !indent ) {
			//do_indent(false /* no trailing comma */);
			if ( is_map ) out_ << " }";
			else out_ << " ]";
		}
		in_map_.pop_back();
		first_.pop_back();
		indent_.pop_back();
	}
}; // YamlEmitter


//////////////////////////////////////////////////////////////////////////////


/// @brief Lightweight class to ease writting JSON documents, a subset of YAML (http://json.org)
/// @details Using this class is not recommended -- use YamlEmitter instead.
///
///  JSON is JavaScript Object Notation, the format for
///  object and array literals in the JavaScript language.
///  It also looks a great deal like nested Python dicts and lists,
///  except for the special JSON values true, false, and null.
///  (Python uses True, False, and None instead.
///  Otherwise, JSON can be embedded directly into Python scripts.)
///
///  JSON is legal subset of YAML, but leaves out much of YAML's complexity and advanced features.
///  At the present time, there appears to be more library support
///  for JSON than for YAML, especially for reading and parsing.
///  (Both are fairly easy to write, as demonstrated here.)
///
///  The topmost level of every valid JSON document consists of exactly one map,
///  which is automatically started for you when the class is initialized.
///  When you're finished you should call end(), and then probably flush().
///
///  Whitespace is JSON documents is not significant.  By default, this class
///  will try to pretty-print, showing the depth of nesting through indentation.
///  When starting a new list or map, you can request that it not be indented
///  (i.e., be printed all on one line), but all lists and maps nested inside it
///  will also be printed without linebreaks, regardless of the requested indentation.
class JsonEmitter : public Emitter
{
public:

	JsonEmitter(std::ostream & out):
		Emitter(out)
	{ start_raw(true, true); } // can't call virtual func from base class ctor

	virtual ~JsonEmitter() {}

	/// @brief JSON doesn't support multiple documents, so this just calls end(1)
	/// to return you to the top (least-nested) level.
	virtual void start_doc()
	{ end(1); }

private:

	JsonEmitter(); // no null ctor
	JsonEmitter(JsonEmitter const &); // no copy ctor

protected:

	/// @brief JSON always quotes strings, regardless.
	virtual
	std::string quote_string(std::string const & s)
	{
		bool needs_quotes(false); // dummy val; will be overwritten below
		std::string t = escape_string(s, needs_quotes);
		return "\""+t+"\""; // JSON always uses quotes, whether needed or not
	}

	virtual
	void do_indent(bool write_comma=true)
	{
		if ( write_comma ) {
			bool first = first_.back();
			if ( first ) first_[ first_.size()-1 ] = false;
			else out_ << ",";
		}
		bool indent = indent_.back();
		if ( indent ) {
			out_ << "\n";
			int depth = (indent_depth_.empty() ? 0 : indent_depth_.back());
			for ( int i = 0; i < depth; ++i ) out_ << indent_str_;
		} else out_ << " ";
	}

	/// @brief JSON has no special marker for list items
	virtual
	void write_list_marker() {} // not used in JSON

	/// @details JSON always uses brackets and braces, regardless of whether we're indenting.
	virtual
	void start_raw(bool is_map, bool indent)
	{
		if ( is_map ) out_ << "{";
		else out_ << "[";
		//if( !indent ) out_ << " "; // not needed
		in_map_.push_back(is_map);
		first_.push_back(true);
		// Once you start not indenting, no children can be indented
		if ( indent_.empty() ) indent_.push_back(indent);
		else indent_.push_back( indent_.back() & indent );
		if ( indent_depth_.empty() ) indent_depth_.push_back(1);
		else indent_depth_.push_back(indent_depth_.back() + 1);
	}

	/// @details JSON always uses brackets and braces, regardless of whether we're indenting.
	virtual
	void end_raw()
	{
		if ( in_map_.empty() ) return; // bad op
		bool is_map = in_map_.back();
		// This is popped first to get the closing brace/bracket to line up right.
		// Other pops go after do_indent() or weird stuff could happen.
		indent_depth_.pop_back();
		do_indent(false /* no trailing comma */);
		if ( is_map ) out_ << "}";
		else out_ << "]";
		in_map_.pop_back();
		first_.pop_back();
		indent_.pop_back();
	}
}; // JsonEmitter


} // basic

#endif // INCLUDED_basic_Emitter_HH
