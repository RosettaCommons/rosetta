// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief An object for reading/writing a simple xml-like language.
/// @details
///  Reads simple xml-like files, such as:
///
///  <loop_design>
///     <params max_cycles=100 />
///     <loop begin=A:10 end=A:15 length=5-6/>
///  </loop_design>
///
///  When read into a tag, the following code would be valid:
///    std::ifstream fin;
///    fin.open("loop_design.tag")
///    TagCOP tag = Tag::create( fin );
///    runtime_assert( tag->getName() == "loop_design" );
///    int max_cycles = tag->getTag("params")->getOption<int>("max_cycles");
///
///  Files must obey the following grammar:
///
///  EBNF for a simple xml-like language:
///  -----------------------------------
///  top              := !xml_schema_tag misc* tag misc*
///  xml_schema_tag   := <?xml (char - '?')* '?>'
///  misc             := comment_tag | (char - '<')
///  comment_tag      := '<!--' ((char - '-') | ('-' (char - '-')))* '-->'
///  tag              := leaf_tag | branch_tag
///  leaf_tag         := '<' name options* '/>'
///  option           := name '=' ( name | quote | squote )
///  name             := (alnum | '_' | '-' | '.' | '*' | ',' )+
///  quote            := '"' (alnum - '"')* '"'
///  squote           := '\'' (alnum - '\'')* '\''
///  branch_tag       := '<' name options* '>' ( tag | misc )* '</' name '>'
///  *whitespace allowed between all tokens
///
///  Less complex than XML, also: no need to quote options,
///  and text outside of tags is ignored. Implemented with
///  the boost spirit library.
///
///
/// To debug the schema:
/// 1) uncomment this:
///
/// #define BOOST_SPIRIT_DEBUG
///
/// 2) add in the definition function below
///
/// BOOST_SPIRIT_DEBUG_NODE(rule_name);
///
/// for each rule defined.
///
/// NOTE: As of r53105, Rosetta doesn't compile with line 162 in boost_1_46_1/boost/spirit/home/debug/debug_node.hpp when BOOST_SPIRIT_DEBUG is enabled.
///
/// @author Paul Murphy
/// @author (schema modified by Matthew O'Meara)

#include <utility/tag/Tag.hh>

#include <iostream>
#include <cstdio>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>


#include <string>
#include <sstream>

#include <boost/spirit/home/classic/attribute.hpp> // AUTO IWYU For closure<>::context_t, closure
#include <boost/spirit/home/classic/core/composite/actions.hpp> // AUTO IWYU For action
#include <boost/spirit/home/classic/core/composite/alternative.hpp> // AUTO IWYU For alternative, operator|
#include <boost/spirit/home/classic/core/composite/kleene_star.hpp> // AUTO IWYU For operator*, kleene_star
#include <boost/spirit/home/classic/core/composite/optional.hpp> // AUTO IWYU For operator!
#include <boost/spirit/home/classic/core/composite/positive.hpp> // AUTO IWYU For operator+
#include <boost/spirit/home/classic/core/composite/sequence.hpp> // AUTO IWYU For operator>>, sequence
#include <boost/spirit/home/classic/core/non_terminal/grammar.hpp> // AUTO IWYU For grammar
#include <boost/spirit/home/classic/core/non_terminal/rule.hpp> // AUTO IWYU For abstract_parser, rule
#include <boost/spirit/home/classic/iterator/position_iterator.hpp> // AUTO IWYU For file_position_base, position_iterator
#include <boost/spirit/home/classic/phoenix/binders.hpp> // AUTO IWYU For bind
#include <boost/spirit/home/classic/utility/functor_parser.hpp> // AUTO IWYU For functor_parser

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector0.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace boost {
void throw_exception(std::exception const&) {
	std::cerr << "some kind of exception was thrown" << std::endl;
}
}

namespace utility {
namespace tag {

/// @details Auto-generated virtual destructor
Tag::~Tag() = default;

using namespace std;

utility::vector0<TagCOP> const Tag::vEmpty_; // need to return this from getTags

Tag::Tag() :
	parentTag_(/* NULL */),
	quote_options_( true )
{}

void Tag::clear() {
	name_.clear();
	mOptions_.clear();
	vTags_.clear();
	mvTags_.clear();
	parentTag_.reset();
}

Tag::options_t const&
Tag::getOptions() const { return mOptions_; }

#ifdef OLDER_CLANG
/// @brief Variant that will only be defined for Sizes.
/// @note Older clang compilers have trouble with the general case being
/// deleted and then specialized cases being defined, even though this is supposed to be supported
/// by the cxx11 standard.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template< class T >
T
Tag::getOption( std::string const &, int const ) const {
	utility_exit_with_message( "Program error: an integer default was specified for a non-integer type.  This is effectively a compilation error, but it can only appear at runtime with older standard libraries." );
	return T(0);
}
#endif

#ifdef PYROSETTA
bool Tag::get_option_bool(std::string const& key) const {
	return Tag::getOption<bool>(key);
}

bool Tag::get_option_bool(std::string const& key, char const * default_as_string_literal) const {
	return Tag::getOption<bool>(key, default_as_string_literal);
}

bool Tag::get_option_bool(std::string const& key, bool const& t_default) const {
	return Tag::getOption<bool>(key, t_default);
}
int Tag::get_option_int(std::string const& key) const {
	return Tag::getOption<int>(key);
}
int Tag::get_option_int(std::string const& key, int const& t_default) const {
	return Tag::getOption<int>(key, t_default);
}
std::string Tag::get_option_string(std::string const& key) const {
	return Tag::getOption<std::string>(key);
}
std::string Tag::get_option_string(std::string const& key, string const& t_default) const {
	return Tag::getOption<std::string>(key, t_default);
}
platform::Real Tag::get_option_real(std::string const& key) const {
	return Tag::getOption<platform::Real>(key);
}
platform::Real Tag::get_option_real(std::string const& key, platform::Real const& t_default) const {
	return Tag::getOption<platform::Real>(key, t_default);
}
platform::Size Tag::get_option_size(std::string const& key) const {
	return Tag::getOption<platform::Size>(key);
}
platform::Size Tag::get_option_size(std::string const& key, platform::Size const& t_default) const {
	return Tag::getOption<platform::Size>(key, t_default);
}
#endif


void Tag::setOptions( options_t const& options ) {
	mOptions_ = options;
	accessed_options_.clear();
}

TagCOP const &
Tag::operator[](std::string const& key) const {
	return getTag(key);
} // operator[]

Tag &
Tag::operator=(Tag const &other) {

	clear();
	name_ = other.name_;
	mOptions_ = other.mOptions_;

	for ( size_t i = 0; i < other.vTags_.size(); ++i ) {
		TagOP tag( new Tag() );
		*tag = *other.vTags_[i];
		addTag( tag );
	}

	return *this;
}

void Tag::reset_accessed_options() const
{
	accessed_options_.clear();
	for ( auto child : vTags_ ) {
		child->reset_accessed_options();
	}
}

void Tag::die_for_unaccessed_options() const {
	auto option= mOptions_.begin();
	for ( ; option != mOptions_.end(); ++option ) {
		std::string const & option_key= option->first;
		options_t::const_iterator const & found= accessed_options_.find(option_key);
		if ( found == accessed_options_.end() ) {
			std::cerr<<"In tag\n"<<*this<<"\n available options are: ";
			for ( options_t::const_iterator available_option = accessed_options_.begin(); available_option != accessed_options_.end(); ++available_option ) {
				std::cerr<<available_option->first<<", ";
			}
			std::cerr<<std::endl;
			utility_exit_with_message("'"+option->first+"' is not a valid option for "+name_);
		}
	}
}

void Tag::die_for_unaccessed_options_recursively() const {
	die_for_unaccessed_options();
	auto begin= vTags_.begin();
	for ( ; begin != vTags_.end(); ++begin ) {
		(*begin)->die_for_unaccessed_options_recursively();
	}
}


void Tag::setName(string const& name) {
	name_ = name;
}

bool
Tag::hasOption( string const& key ) const {
	auto i = mOptions_.find(key);
	return i != mOptions_.end() && i->second.size() != 0 ;
}

void Tag::addTag( TagOP tag ) {
	tag->parentTag_ = get_self_weak_ptr();
	vTags_.push_back(tag);
	mvTags_[ tag->getName() ].push_back( tag );
}

utility::vector0< TagCOP > const &
Tag::getTags() const {
	return vTags_;
}

utility::vector0< TagCOP > const &
Tag::getTags( string const& name ) const {
	auto i = mvTags_.find(name);
	if ( i == mvTags_.end() ) {
		return vEmpty_;
	}
	return i->second;
}

TagCOP const &
Tag::getTag( string const& name ) const
{
	utility::vector0< TagCOP > const& v = getTags(name);
	if ( v.size() != 1 ) {
		std::stringstream error_msg;
		error_msg << "Tag::getTag - name=" << name << ", appears " << v.size() << " times in xml file. Only allowed once.\n";
		error_msg << *this << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, error_msg.str() );
	}
	return v[0];
} // getTag

bool
Tag::hasTag(std::string const& name ) const
{
	return getTags(name).size() != 0;
}

void Tag::write(std::ostream& out, int num_tabs ) const
{
	string tabs(num_tabs,'\t');

	out << tabs << "<" << name_;
	for ( auto const & mOption : mOptions_ ) {
		out << " " << mOption.first << "=";
		if ( quote_options_ ) out << "\"";
		out << mOption.second;
		if ( quote_options_ ) out << "\"";
	}

	if ( vTags_.size() == 0 ) {
		out << "/>\n";
	} else {
		out << ">\n";
		for ( auto tag : vTags_ ) {
			tag->write(out,num_tabs+1);
		}
		out << tabs << "</" << name_ << ">\n";
	}

} // Tag::write

std::string
Tag::to_string( int num_tabs ) const {
	std::stringstream ss;
	write( ss, num_tabs );
	return ss.str();
}


size_t Tag::size() const {
	size_t rval = 1;
	for ( auto const & vTag : vTags_ ) {
		rval += (*vTag).size();
	}
	return rval;
} // Tag::size

/// @details This is explicit specialization for boolean values
/// to allow for use of "true" "false" etc. in addition to 1 and 0
/// @throws Throws a utility::excn::EXCN_Msg_Exception if the option
/// is proviced, but the string that's given is not a valid true/false
/// string (either "true", "false", "1" or "0").
template<>
bool
Tag::getOption<bool>(std::string const& key, bool const& t_default) const {
	auto i = mOptions_.find(key);
	if ( i == mOptions_.end() ) {
		accessed_options_[key]= key;
		return t_default;
	}
	accessed_options_[key]= i->second;
	if ( utility::is_true_string( i->second ) ) { return true; }
	if ( utility::is_false_string( i->second ) ) { return false; }
	std::stringstream error_message;
	error_message << "getOption: key= " << key << " stream extraction for boolean value failed! Tried to parse '" << i->second << "\n";
	throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
	return t_default;
}

/// @brief This is for the variant in which someone has specified a default
/// using "true" instead of true, or "false" instead of false.
template<>
bool
Tag::getOption<bool>(std::string const& key, char const * default_as_string_literal) const {
	auto i = mOptions_.find(key);
	if ( i != mOptions_.end() ) {
		return getOption<bool>(key);
	}

	runtime_assert_string_msg( default_as_string_literal != nullptr, "Program error: the developer has erroneously provided 0 or nullptr as the default for boolean option key \"" + key + "\".  This is effectively a compilation error, but it is only detectable at runtime.  The temporary workaround is to provide a value for the " + key + " option.  Please also inform a developer." );
	std::string const default_string( default_as_string_literal );

	if ( utility::is_true_string( default_string ) ) { return true; }
	if ( utility::is_false_string( default_string ) ) { return false; }
	std::stringstream error_message;
	error_message << "getOption: key= " << key << " stream extraction for default boolean value failed! Tried to parse \"" << default_string << "\", but the developer has provided a default string that can't be interpreted as a boolean.  This is really coding problem that would ideally be a compilation error, but it can only show up at runtime.\n";
	throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
	return true; //Keep the compiler satsified.
}

/// @throws Throws a utility::excn::EXCN_Msg_Exception if none of the
/// accepted boolean strings is provided.
template<>
bool
Tag::getOption<bool>(std::string const& key) const {
	auto i = mOptions_.find(key);
	if ( i == mOptions_.end() ) {
		std::stringstream error_message;
		error_message << "Option '" << key << "' not found in Tag named '" << this->getName() << "'.";
		throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
	}
	accessed_options_[key]= i->second;
	if ( utility::is_true_string( i->second ) ) { return true; }
	if ( utility::is_false_string( i->second ) ) { return false; }
	std::stringstream error_message;
	error_message << "getOption: key= " << key << " stream extraction for boolean value failed! Tried to parse '" << i->second << "\n";
	throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
	return false; // appease compiler
}

template<>
AutoBool
Tag::getOption< AutoBool >(std::string const& key, AutoBool const& t_default) const {
	auto i = mOptions_.find(key);
	if ( i == mOptions_.end() ) {
		accessed_options_[key]= key;
		return t_default;
	}
	accessed_options_[key]= i->second;
	if ( utility::is_true_string( i->second ) ) { return AutoBool::True; }
	if ( utility::is_false_string( i->second ) ) { return AutoBool::False; }
	if ( utility::upper( i->second ) == "AUTO" ) { return AutoBool::Auto; }
	std::stringstream error_message;
	error_message << "getOption: key= " << key << " stream extraction for autoboolean value failed! Tried to parse '" << i->second << "\n";
	throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
	return t_default;
}

/// @brief This is for the variant in which someone has specified a default
/// using "true" instead of true, or "false" instead of false.
template<>
AutoBool
Tag::getOption<AutoBool>(std::string const& key, char const * default_as_string_literal) const {
	auto i = mOptions_.find(key);
	if ( i != mOptions_.end() ) {
		return getOption<AutoBool>(key);
	}
	runtime_assert_string_msg( default_as_string_literal != nullptr, "Program error: the developer has erroneously provided 0 or nullptr as the default for AutoBool option key \"" + key + "\".  This is effectively a compilation error, but it is only detectable at runtime.  The temporary workaround is to provide a value for the " + key + " option.  Please also inform a developer." );
	std::string const default_string( default_as_string_literal );

	if ( utility::is_true_string( default_string ) ) { return AutoBool::True; }
	if ( utility::is_false_string( default_string ) ) { return AutoBool::False; }
	if ( utility::upper( default_string ) == "AUTO" ) { return AutoBool::Auto; }
	std::stringstream error_message;
	error_message << "getOption: key= " << key << " stream extraction for default AutoBool value failed! Tried to parse \"" << default_string << "\", but the developer has provided a default string that can't be interpreted as an AutoBool.  This is really coding problem that would ideally be a compilation error, but it can only show up at runtime.\n";
	throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
	return AutoBool::True; //keep the compiler happy
}

/// @throws Throws a utility::excn::EXCN_Msg_Exception if none of the
/// accepted boolean strings is provided.
template<>
AutoBool
Tag::getOption< AutoBool >(std::string const& key) const {
	auto i = mOptions_.find(key);
	if ( i == mOptions_.end() ) {
		std::stringstream error_message;
		error_message << "Option '" << key << "' not found in Tag named '" << this->getName() << "'.";
		throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
	}
	accessed_options_[key]= i->second;
	if ( utility::is_true_string( i->second ) ) { return AutoBool::True; }
	if ( utility::is_false_string( i->second ) ) { return AutoBool::False; }
	if ( utility::upper( i->second ) == "AUTO" ) { return AutoBool::Auto; }
	std::stringstream error_message;
	error_message << "getOption: key= " << key << " stream extraction for autoboolean value failed! Tried to parse '" << i->second << "\n";
	throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str() );
	return AutoBool::False; // appease compiler
}

// @brief If this were uncommented, this would add a special-case treatment to ensure that integer defaults
// get interpreted as Reals when setting a Real option.  So, for example, a developer could write
// tag->getOption<core::Real>("myoption", 1) instead of tag->getOption<core::Real>("myoption", 1.0).
// @details This has been deliberately REMOVED because there's value in not allowing a Real's default to be set with a Size
// -- removing this revealed a number of actual errors.  I'm leaving this here, commented out, in case we someday
// decide to allow a Real option to be given a Size default.
// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
// template<>
// platform::Real
// Tag::getOption<platform::Real>(
//  std::string const & key,
//  int const default_int
// ) const {
//  return getOption< platform::Real >(key, static_cast<platform::Real>(default_int));
// }

/// @brief Special-casing the string literal version for string options.  In this case,
/// there shouldn't be an error thrown.  A string literal should be allowed to set the
/// default value for a string.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
std::string
Tag::getOption<std::string>( std::string const & key, char const * default_as_string_literal ) const {
	runtime_assert_string_msg( default_as_string_literal != nullptr, "Program error: the developer has erroneously provided 0, nullptr, or false as the default for string option key \"" + key + "\".  This is effectively a compilation error, but it is only detectable at runtime.  Please inform a developer." );
	return getOption< std::string >( key, std::string( default_as_string_literal ) );
}

/// @brief Special-casing to ensure that 0 gets interpreted as Size(0) rather than nullptr.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
platform::Size
Tag::getOption<platform::Size>( std::string const & key, int const default_int ) const {
	runtime_assert_string_msg( default_int >= 0, "Error in Tag::getOption(): A platform::Size cannot have a negative default.  The platform::Size type is unsigned." );
	return getOption< platform::Size >( key, static_cast<platform::Size>(default_int) );
}

/// @brief Special-casing for signed ints.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
int
Tag::getOption<int>( std::string const & key, int const default_int ) const {
	auto i = mOptions_.find(key);
	if ( i == mOptions_.end() ) {
		accessed_options_[key]= key;
		return default_int;
	}
	accessed_options_[key]= i->second;

	try{
		return boost::lexical_cast<int>(i->second);
	} catch(boost::bad_lexical_cast &) {
		std::stringstream error_message;
		error_message << "getOption: key= " << key << " stream extraction failed! Tried to parse '" << i->second << "' as a signed integer, but could not!\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
	}
	return default_int;
}

/// @brief Special-casing for int64_t.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
template<>
int64_t
Tag::getOption<int64_t>(
	std::string const & key,
	int const default_int
) const {
	return getOption<int64_t>( key, static_cast<int64_t>(default_int) );
}

// ____________________ <boost::spirit parser definition> ____________________

namespace {

using namespace boost::spirit::classic;
using namespace phoenix;

struct string_closure : public boost::spirit::classic::closure<string_closure,string> {
	member1 value;
}; // string_closure

typedef pair<string,string> option_value_type;
struct option_closure : public boost::spirit::classic::closure<option_closure,option_value_type > {
	member1 value;
}; // option_closure

typedef map<string,string> options_value_type;
struct options_closure : public boost::spirit::classic::closure<options_closure,options_value_type >
{
	member1 value;
}; // options_closure

typedef std::pair<string,options_value_type> name_and_options_value_type;
struct name_and_options_closure : public boost::spirit::classic::closure<name_and_options_closure,name_and_options_value_type>
{
	member1 value;
}; // name

struct tag_closure : public boost::spirit::classic::closure< tag_closure,TagOP >
{
	member1 value;
}; // tag_closure

void set_name_and_options( TagOP & tag, name_and_options_value_type const & v )
{
	tag = utility::pointer::make_shared< Tag >();
	tag->setName( v.first );
	for ( auto const & i : v.second ) {
		tag->setOption( i.first, i.second );
	}
} // set_name_and_options

void add_tag( TagOP tag, TagOP other ) {
	tag->addTag( other );
} // set_name_and_options

void assert_matching_tag_names( TagCOP tag, string const& closing_tag_name )
{
	if ( tag->getName() != closing_tag_name ) {
		std::cerr << "error - tag names do not match! (Possibly due to unclosed tag.) " << tag->getName() << " != " << closing_tag_name << endl;
	}
} // set_name_and_options

struct error_report_parser {
	//   char const* eol_msg;
	//   char const* msg;
	vector< file_position >& errors;
	explicit error_report_parser( vector<file_position>& error_vec ) : errors( error_vec ) {}

	using result_t = boost::spirit::classic::nil_t;

	template <typename ScannerT>
	int operator()(ScannerT const& scan, result_t&) const {
		errors.push_back( scan.first.get_position() );
		return -1;
	}

};

//typedef functor_parser<error_report_parser> error_report_p;

struct tag_grammar : public grammar<tag_grammar,tag_closure::context_t> {

	vector< file_position > errors;
	functor_parser<error_report_parser> error_report_p;

	tag_grammar() : error_report_p( error_report_parser(errors) ) {}

	template< class Scanner >
	struct definition {

		rule<Scanner> top;
		rule<Scanner> xml_schema_tag;
		rule<Scanner> misc;
		rule<Scanner> comment_tag;
		rule<Scanner,tag_closure::context_t> tag;
		rule<Scanner,name_and_options_closure::context_t> leaf_tag;
		rule<Scanner,name_and_options_closure::context_t> name_and_options;
		rule<Scanner,string_closure::context_t> name;
		rule<Scanner,options_closure::context_t> options;
		rule<Scanner,option_closure::context_t> option;
		rule<Scanner,string_closure::context_t> name_or_quote;
		rule<Scanner,string_closure::context_t> quote;
		rule<Scanner,string_closure::context_t> squote;
		rule<Scanner,tag_closure::context_t> branch_tag;
		rule<Scanner,name_and_options_closure::context_t> open_tag;
		rule<Scanner,string_closure::context_t> close_tag;

		explicit definition( tag_grammar const& self) {

			top
				= !xml_schema_tag >> *misc >> tag[ self.value = arg1 ] >> *misc
				| self.error_report_p;

			xml_schema_tag = str_p("<?xml") >> *~ch_p('?') >> str_p("?>");

			misc = ( ~ch_p('<') | comment_tag );

			comment_tag =
				str_p("<!--") >>
				*( ~ch_p('-') | ( ch_p('-') >> ~ch_p('-') ) ) >>
				str_p("-->");

			tag
				= leaf_tag[ phoenix::bind(set_name_and_options)(tag.value,arg1) ]
				| branch_tag[ tag.value = arg1 ];


			leaf_tag =
				ch_p('<') >> *space_p >>
				name_and_options[ leaf_tag.value = arg1 ] >>
				ch_p('/') >> *space_p >> ch_p('>');

			name_and_options = name[ phoenix::bind( &name_and_options_value_type::first )( name_and_options.value ) = arg1 ] >> *space_p >> options[ phoenix::bind( &name_and_options_value_type::second)(name_and_options.value) = arg1 ] >> *space_p;

			name = (+(alnum_p | ch_p('_') | ch_p(':') | ch_p('-') | ch_p('.') | ch_p('*') | ch_p(',') ))[ name.value = construct_<string>(arg1,arg2) ];

			typedef pair< options_value_type::iterator,bool> (options_value_type::*insert_t)( pair<const string,string> const&);
			insert_t ins = &options_value_type::insert;
			options = *option[ phoenix::bind( ins )( options.value, arg1 ) ]; // without the typedef, C++ can't figure out the type of &::insert since it is templatized

			option  = name[ phoenix::bind( &option_value_type::first)(option.value) = arg1 ] >> *space_p >> ch_p('=') >> *space_p >> name_or_quote[ phoenix::bind( &option_value_type::second)(option.value) = arg1 ] >> *space_p;

			name_or_quote = name[name_or_quote.value = arg1 ] | quote[name_or_quote.value = arg1 ] | squote[name_or_quote.value = arg1 ];

			quote  = ch_p('"' ) >> (*~ch_p('"' ))[  quote.value = construct_<string>(arg1,arg2) ] >> ch_p('"' ) >> *space_p;
			squote = ch_p('\'') >> (*~ch_p('\''))[ squote.value = construct_<string>(arg1,arg2) ] >> ch_p('\'') >> *space_p;


			branch_tag =
				open_tag[ phoenix::bind(set_name_and_options)(branch_tag.value,arg1) ] >>
				*( tag[ phoenix::bind(add_tag)(branch_tag.value,arg1) ] | misc ) >>
				close_tag[ phoenix::bind(assert_matching_tag_names)(branch_tag.value,arg1) ];

			open_tag =
				ch_p('<') >> *space_p >>
				name_and_options[ open_tag.value = arg1 ] >>
				ch_p('>');

			close_tag =
				ch_p('<') >> *space_p >>
				ch_p('/') >> *space_p >>
				name[ close_tag.value = arg1 ] >> *space_p >>
				ch_p('>');

		} // definition

		rule<Scanner> const& start() const { return top; }

	}; // definition

}; // tag_grammar

// ____________________ </boost::spirit parser definition> ____________________

void print_error( ostream& out, string const& str, file_position const& lc ) {

	size_t line_begin = 0;
	size_t line_end = str.find('\n',line_begin);

	for ( int line = 1; line < lc.line; ++line ) {
		line_begin = line_end+1;
		if ( line_begin < str.size() ) {
			line_end = str.find('\n',line_begin);
		}
	}

	if ( line_end == std::string::npos ) {
		line_end = str.size();
	}


	if ( line_begin < str.size() ) {
		out << "Tag::read - parse error - file:" << lc.file << " line:" << lc.line << " column:" << lc.column << " - ";
		for ( size_t column = line_begin; column < line_end; ++column ) {
			out << str[column];
		}
		out << "\n";

		out << "Tag::read - parse error - file:" << lc.file << " line:" << lc.line << " column:" << lc.column << " - ";
		for ( int k = 0; k < (lc.column-1); ++k ) {
			size_t column = line_begin + k;
			if ( column < str.size() && str[ column ] == '\t' ) {
				out << "\t";
			} else {
				out << " ";
			}
		}
		out << "^\n" << endl;
	} else {
		out << "Tag::read - parse error - file:" << lc.file << " line:" << lc.line << " column:" << lc.column << " - somewhere in file\n" << endl;
	}
}

} // namespace

void Tag::read( std::istream& in ) {
	string str;
	while ( in.peek() != EOF ) str.push_back( in.get() );
	read( str );
}

void Tag::read( std::string const & str ) {

	using iterator_t = position_iterator<const char *>;
	iterator_t begin(str.c_str(), str.c_str() + str.size(), "istream"); // is this remotely correct?
	iterator_t end;
	//begin.set_tabchars(8); // ??? what does this do

	TagOP tag; // don't need to initialize this at all
	tag_grammar g;
	auto parsed = parse(begin, end, g[ var(tag) = arg1 ] ); // on a separate line to get around GCC 3.12 peculiarities
	bool full = parsed.full;

	if ( tag && full ) {
		*this = *tag;
	} else {
		stringstream err_msg;
		err_msg << "Tag::read - parse error, printing backtrace.\n" << endl;
		for ( auto & error : g.errors ) {
			print_error(err_msg,str,error);
		}

		throw CREATE_EXCEPTION(utility::excn::BadInput, err_msg.str() );
	}

} // NOLINT -- Needed as the spirit parse() function somehow stashes stuff into a global variable

/*
void Tag::read(std::istream& in ) {
string str;
while( in.peek() != EOF ) str.push_back(in.get());

tag_grammar g;

TagCOP tag; // don't need to initialize this at all

if( parse( str.c_str(), g[ var(tag) = arg1 ] ).full ) {
*this = *tag;
}
else {
runtime_assert( false );
}

} // read
*/

TagOP
Tag::create(std::istream& in) {
	TagOP tag( new Tag() );
	tag->read(in);
	return tag;
} // creates a new tag and reads into it

TagOP
Tag::create(std::string const &instring) {
	std::stringstream in(instring);

	return create(in);
} // creates a new tag and reads into it

TagOP
Tag::clone() const {

	TagOP rval( new Tag() );
	*rval = *this;
	return rval;

} // Tag::clone

/// @brief if true, options will be quoted when the tag is outputted
///        if false, options will be left as-is (default)
/// @param[in] quote_options_val Whether or not option values should be quoted.
///                              Default=true
void
Tag::set_quote_options( bool const quote_options_val ) {
	quote_options_ = quote_options_val;
} // Tag::set_quote_options

std::ostream& operator<<(std::ostream& out, Tag const& tag) {
	tag.write(out);
	return out;
}

std::ostream& operator<<(std::ostream& out, TagCOP const& tag_ptr) {
	tag_ptr->write(out);
	return out;
}


} // namespace tag
} // namespace utility

#ifdef    SERIALIZATION

//template< class Archive >
//void deserialize_tags_t( Archive & arc, utility::tag::Tag::tags_t & tags )
//{
//
//}

/// @brief Automatically generated serialization method
template< class Archive >
void
utility::tag::Tag::save( Archive & arc ) const {
	arc( CEREAL_NVP( name_ ) ); // std::string
	arc( CEREAL_NVP( mOptions_ ) ); // options_t
	arc( CEREAL_NVP( accessed_options_ ) ); // options_t
	arc( CEREAL_NVP( vTags_ ) ); // tags_t
	arc( CEREAL_NVP( mvTags_ ) ); // std::map<std::string, tags_t>
	arc( CEREAL_NVP( parentTag_ ) ); // TagCAP
	arc( CEREAL_NVP( quote_options_ ) ); // bool

}

/// @brief Automatically generated deserialization method
template< class Archive >
void
utility::tag::Tag::load( Archive & arc ) {
	arc( name_ ); // std::string
	arc( mOptions_ ); // options_t
	arc( accessed_options_ ); // options_t

	utility::vector0< std::shared_ptr< utility::tag::Tag > > local_vTags;
	arc( local_vTags ); // tags_t
	vTags_ = local_vTags; // copy the non-const pointer(s) into the const pointer(s)

	std::map< std::string, utility::vector0< TagOP > > local_mvTags;
	arc( local_mvTags ); // std::map<std::string, tags_t>
	for ( std::map< std::string, utility::vector0< TagOP > >::const_iterator
			iter = local_mvTags.begin(), iter_end = local_mvTags.end();
			iter != iter_end; ++iter ) {
		mvTags_[ iter->first ] = iter->second;
	}

	TagAP local_parentTag;
	arc( local_parentTag ); // TagCAP
	parentTag_ = local_parentTag;

	arc( quote_options_ );
}

SAVE_AND_LOAD_SERIALIZABLE( utility::tag::Tag );
CEREAL_REGISTER_TYPE( utility::tag::Tag )

CEREAL_REGISTER_DYNAMIC_INIT( utility_tag_Tag )
#endif // SERIALIZATION
