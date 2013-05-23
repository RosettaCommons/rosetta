// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief An object for reading/writing a simple xml-like language.
/// @detailed
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
///    TagPtr tag = Tag::create( fin );
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
///  option           := name '=' ( name | quote )
///  name             := (alnum | '_' | '-' | '.' | '*' | ',' )+
///  quote            := '"' (alnum - '"')* '"'
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
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_attribute.hpp>
#include <boost/spirit/include/classic_closure.hpp>
#include <boost/spirit/include/phoenix1_primitives.hpp>
#include <boost/spirit/include/phoenix1_operators.hpp>
#include <boost/spirit/include/phoenix1_functions.hpp>
#include <boost/spirit/include/phoenix1_binders.hpp>
#include <boost/spirit/include/phoenix1_casts.hpp>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/spirit/include/classic_position_iterator.hpp> // for reporting where errors occur
#include <boost/spirit/include/classic_functor_parser.hpp>
#include <platform/types.hh>

#include <string>
#include <iostream>
#include <sstream>

namespace boost {
	void throw_exception(std::exception const&) {
		std::cerr << "some kind of exception was thrown" << std::endl;
	}
}

namespace utility {
namespace tag {

/// @details Auto-generated virtual destructor
Tag::~Tag() {}

using namespace std;

size_t Tag::num_tags = 0;

utility::vector0<TagPtr> const Tag::vEmpty_; // need to return this from getTags

Tag::Tag() {
	++num_tags;
}

void Tag::clear() {
	name_.clear();
	mOptions_.clear();
	vTags_.clear();
	mvTags_.clear();
}

void Tag::die_for_unaccessed_options(){
	options_t::const_iterator option= mOptions_.begin();
	for(; option != mOptions_.end(); ++option){
		std::string const & option_key= option->first;
		options_t::const_iterator const & found= accessed_options_.find(option_key);
		if( found == accessed_options_.end()){
			std::cerr<<"In tag\n"<<*this<<"\n available options are: ";
			for( options_t::const_iterator available_option = accessed_options_.begin(); available_option != accessed_options_.end(); ++available_option )
				std::cerr<<available_option->first<<", ";
			std::cerr<<std::endl;
			utility_exit_with_message("'"+option->first+"' is not a valid option for "+name_);
		}
	}
}


void Tag::setName(string const& name) {
	name_ = name;
}

bool
Tag::hasOption( string const& key ) const {
	options_t::const_iterator i = mOptions_.find(key);
	return i != mOptions_.end() && i->second.size() != 0 ;
}

void Tag::addTag( TagPtr const& tag ) {
	vTags_.push_back(tag);
	mvTags_[ tag->getName() ].push_back( tag );
}

utility::vector0< TagPtr > const &
Tag::getTags() const {
	return vTags_;
}

utility::vector0< TagPtr > const &
Tag::getTags( string const& name ) const {
	map< string, utility::vector0< TagPtr > >::const_iterator i = mvTags_.find(name);
	if( i == mvTags_.end() ) {
		return vEmpty_;
	}
	return i->second;
}

TagPtr const &
Tag::getTag( string const& name ) const
{
	utility::vector0< TagPtr > const& v = getTags(name);
	if( v.size() != 1 ) {
		std::cerr << "Tag::getTag - name=" << name << ", appears " << v.size() << " times in xml file. Only allowed once." << std::endl;
		std::cerr << *this << std::endl;
		runtime_assert( false );
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
	for( options_t::const_iterator i = mOptions_.begin(), end_i = mOptions_.end(); i != end_i; ++i ) {
		out << " " << i->first << "=" << i->second;
	}

	if( vTags_.size() == 0 ) {
		out << "/>\n";
	} else {
		out << ">\n";
		for( tags_t::const_iterator i = vTags_.begin(), end_i = vTags_.end(); i != end_i; ++i ) {
			TagPtr tag = *i;
			tag->write(out,num_tabs+1);
		}
		out << tabs << "</" << name_ << ">\n";
	}

} // Tag::write

size_t Tag::size() const {
	size_t rval = 1;
	for( tags_t::const_iterator i = vTags_.begin(), end_i = vTags_.end(); i != end_i; ++i ) {
		rval += (**i).size();
	}
	return rval;
} // Tag::size

//This is explicit specialization for boolean values
//to allow for use of "true" "false" etc. in addition to 1 and 0
template<>
bool
Tag::getOption<bool>(std::string const& key, bool const& t_default) const {
	options_t::const_iterator i = mOptions_.find(key);
	if( i == mOptions_.end() ) {
		accessed_options_[key]= key;
		return t_default;
	}
	accessed_options_[key]= i->second;
	if ( utility::is_true_string( i->second ) ) { return true; }
	if ( utility::is_false_string( i->second ) ) { return false; }
	std::cerr << "getOption: key= " << key << " stream extraction for boolean value failed! Tried to parse '" << i->second <<
				"' returning default value: '" << t_default << "'." << std::endl;
	return t_default;
}

template<>
bool
Tag::getOption<bool>(std::string const& key) const {
	options_t::const_iterator i = mOptions_.find(key);
	if( i == mOptions_.end() ) {
		std::cerr << "Option " << key << " not found." << std::endl;
		runtime_assert( false );
	}
	accessed_options_[key]= i->second;
	if ( utility::is_true_string( i->second ) ) { return true; }
	if ( utility::is_false_string( i->second ) ) { return false; }
	std::cerr << "getOption: key= " << key << " stream extraction for boolean value failed! Tried to parse '" << i->second <<
				"' returning false." << std::endl;
	return false;
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

struct tag_closure : public boost::spirit::classic::closure<tag_closure,TagPtr>
{
	member1 value;
}; // tag_closure

void set_name_and_options( TagPtr& tag, name_and_options_value_type const & v )
{
	tag = TagPtr( new Tag() );
	tag->setName( v.first );
	for(  map<string,string>::const_iterator i = v.second.begin(), end_i = v.second.end(); i != end_i ; ++i ) {
		tag->setOption( i->first, i->second );
	}
} // set_name_and_options

void add_tag( TagPtr& tag, TagPtr const& other ) {
	tag->addTag( other );
} // set_name_and_options

void assert_matching_tag_names( TagPtr& tag, string const& closing_tag_name )
{
	if( tag->getName() != closing_tag_name ) {
		std::cerr << "error - tag names do not match! (Possibly due to unclosed tag.) " << tag->getName() << " != " << closing_tag_name << endl;
	}
} // set_name_and_options

struct error_report_parser {
	//   char const* eol_msg;
	//   char const* msg;
	vector< file_position >& errors;
	error_report_parser( vector<file_position>& errors ) : errors(errors) {}

	typedef boost::spirit::classic::nil_t result_t;

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
		rule<Scanner,tag_closure::context_t> branch_tag;
		rule<Scanner,name_and_options_closure::context_t> open_tag;
		rule<Scanner,string_closure::context_t> close_tag;

		definition( tag_grammar const& self) {

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
				= leaf_tag[ bind(set_name_and_options)(tag.value,arg1) ]
				| branch_tag[ tag.value = arg1 ];


			leaf_tag =
				ch_p('<') >> *space_p >>
				name_and_options[ leaf_tag.value = arg1 ] >>
				ch_p('/') >> *space_p >> ch_p('>');

			name_and_options = name[ bind( &name_and_options_value_type::first )( name_and_options.value ) = arg1 ] >> *space_p >> options[ bind( &name_and_options_value_type::second)(name_and_options.value) = arg1 ] >> *space_p;

			name = (+(alnum_p | ch_p('_') | ch_p(':') | ch_p('-') | ch_p('.') | ch_p('*') | ch_p(',') ))[ name.value = construct_<string>(arg1,arg2) ];

			typedef pair< options_value_type::iterator,bool> (options_value_type::*insert_t)( pair<const string,string> const&);
			insert_t ins = &options_value_type::insert;
			options = *option[bind( ins )( options.value, arg1 ) ]; // without the typedef, C++ can't figure out the type of &::insert since it is templatized

			option  = name[ bind( &option_value_type::first)(option.value) = arg1 ] >> *space_p >> ch_p('=') >> *space_p >> name_or_quote[ bind( &option_value_type::second)(option.value) = arg1 ] >> *space_p;

			name_or_quote = name[name_or_quote.value = arg1 ] | quote[name_or_quote.value = arg1 ];

			quote = ch_p('"') >> (*~ch_p('"'))[ quote.value = construct_<string>(arg1,arg2) ] >> ch_p('"') >> *space_p;


			branch_tag =
				open_tag[ bind(set_name_and_options)(branch_tag.value,arg1) ] >>
				*( tag[ bind(add_tag)(branch_tag.value,arg1) ] | misc ) >>
				close_tag[ bind(assert_matching_tag_names)(branch_tag.value,arg1) ];

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

	for( int line = 1; line < lc.line; ++line ) {
		line_begin = line_end+1;
		if( line_begin < str.size() ) {
			line_end = str.find('\n',line_begin);
		}
	}

	if( line_end == std::string::npos ){
		line_end = str.size();
	}


	if( line_begin < str.size() ) {
		out << "Tag::read - parse error - file:" << lc.file << " line:" << lc.line << " column:" << lc.column << " - ";
		for( size_t column = line_begin; column < line_end; ++column ) {
			out << str[column];
		}
		out << "\n";

		out << "Tag::read - parse error - file:" << lc.file << " line:" << lc.line << " column:" << lc.column << " - ";
		for( int k = 0; k < (lc.column-1); ++k ) {
			size_t column = line_begin + k;
			if( column < str.size() && str[ column ] == '\t' ) {
				out << "\t";
			}
			else {
				out << " ";
			}
		}
		out << "^\n" << endl;
	}
	else {
		out << "Tag::read - parse error - file:" << lc.file << " line:" << lc.line << " column:" << lc.column << " - somewhere in file\n" << endl;
	}
}

} // namespace

void Tag::read(std::istream& in ) {
	string str;
	while( in.peek() != EOF ) str.push_back(in.get());

	typedef position_iterator<char const*> iterator_t;
	iterator_t begin(str.c_str(), str.c_str() + str.size(), "istream"); // is this remotely correct?
	iterator_t end;
	//begin.set_tabchars(8); // ??? what does this do

	TagPtr tag; // don't need to initialize this at all
	tag_grammar g;
	bool full = parse(begin, end, g[ var(tag) = arg1 ] ).full;

	if( tag && full ) {
		*this = *tag;
	} else {
		stringstream err_msg;
		err_msg << "Tag::read - parse error, printing backtrace.\n" << endl;
		for( size_t k = 0; k < g.errors.size(); ++k ) {
			print_error(err_msg,str,g.errors[k]);
		}

		throw utility::excn::EXCN_BadInput( err_msg.str() );
	}

} // read

/*
void Tag::read(std::istream& in ) {
string str;
while( in.peek() != EOF ) str.push_back(in.get());

tag_grammar g;

TagPtr tag; // don't need to initialize this at all

if( parse( str.c_str(), g[ var(tag) = arg1 ] ).full ) {
		*this = *tag;
}
else {
		runtime_assert( false );
}

} // read
*/

TagPtr
Tag::create(std::istream& in) {
	TagPtr tag( new Tag() );
	tag->read(in);
	return tag;
} // creates a new tag and reads into it

TagPtr
Tag::create(std::string instring) {
	std::stringstream in(instring);

	return create(in);
} // creates a new tag and reads into it

TagPtr
Tag::clone() const {

	TagPtr rval = TagPtr( new Tag() );

	rval->name_ = name_;
	rval->mOptions_ = mOptions_;

	for( size_t i = 0; i < vTags_.size(); ++i ) {
		rval->addTag( vTags_[i]->clone() );
	}

	return rval;

} // Tag::clone

std::ostream& operator<<(std::ostream& out, Tag const& tag) {
	tag.write(out);
	return out;
}

std::ostream& operator<<(std::ostream& out, TagPtr const& tag_ptr) {
	tag_ptr->write(out);
	return out;
}

} // namespace tag
} // namespace utility
