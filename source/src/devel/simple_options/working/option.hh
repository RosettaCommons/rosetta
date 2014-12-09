// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef option_H
#define option_H

#include <map>
#include <string>
#include <list>
#include <iostream>
#include <vector>

#define CONTINUE    '\\'
#define COMMENTSIGN ';'
#define STRLEN 2000

/* #############################################################
	 ################## OPTION HANDLING ##########################
	 user oriented description:
	 to provide options for your module you need to inherit from
	 the virtual class OptionModule. Put any variable you need
	 to set options for in this class and overload the method
	 register_options to connect the variables with names, descriptions
	 and default values, as they will appear in the option file.
	 This task is actually simplified by using the provided macros
	 START_OPT,REG,BEGIN_VARLIST,END_OPT which have to be used in that
	 sequence, as in the example given below:

	 START_OPT(t_fabrelax_opt,"fabrelax","options for ab-initio relaxation")
	 REG(COPT("GENERAL OPTIONS"))
	 REG(SOPT("string_option",title,"default_name","that is how you define string options"));
	 REG(IOPT("cycles",cycles,3,"how many cycles to run"));
	 BEGIN_VAR_LIST
	 // the same variables as used in the REG lines have to be used
	 std::string string_option;
	 int cycles;
	 END_OPT

	 to use the options provided in this 'option_module' anywhere in code
	 do the following:

	 t_fabrelax_opt opt; //use the class you have defined with START_OPT
	 option_file.get_options(opt);
	 ...
	 for (i=0;i<opt.cycles;i++) {...}
	 ...

	 you can even use the same option_module more than once. It will
	 only be read from file once. option_file keeps a copy and gives
	 you a new copy of it:
	 ... { <some local code>
			 t_fabrelax_opt opt_copy;
		 option_file.get_options(opt_copy);
		 if (opt_copy.cycles==3) cerr << "yipee" << endl;
	 };

	 ...


	 the option file is read once and all matching entries are
	 used to fill the registered options. Entries not found in the
	 input file are set with their default values.
	 as soon as option_file.write_options() is called all
	 registerd options are written in an output option-file.
	 ===========================================================

	 how does it work:
	 option_file contains to containers that correspond to two
	 distinct layers:

	 the first layer -- the code-side (OptionModule)
	 the OptionModules that contain the variables used in the program
	 code and that provided virtual functions to connect the variables
	 with strings and to set default values

	 the second layer -- the file-side (t_option_block/t_option_line)
	 t_option_file contains a list of option_blocks, corresponding
	 to the '[ block ]' entries in the file. All options following one
	 '[ block ]' statement are contained in the 'lines' container of
	 a given instance of t_option_block.
	 if blocks or lines not contained in the file are registered from
	 OptionModule they will be put into these lists and are
	 also written to file by the method write_options().

	 when a child of OptionModule is registered all value-str of the
	 lines in the corresponding t_option_block will be used to determine
	 the values for the variables in the OptionModule child.
	 Because we do not know how to transform a string into a variable
	 of the correct type the polymorphy in set_value() of the
	 base_opt class tree is used for that step.
	 Any type used for options must be encapsulated into a child of
	 base_opt and overload set_value(std::string &valstr) to
	 set its variable and to set the default value if the string
	 is empty.

	 how to extend it:
	 * extend the list of types that can be used as options:
	 derive a new class from base_opt and overload set_default()

	 * in the moment options can be only string in one line.

	 * if we have many different types for options it
			might be worthwhile to use the iostream library throughout
		 and to use overloaded <<,>> operators for io of objects.
		 i think there is something like a stringbuf that could be
		 used for that.

*/

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace devel {
namespace option {

class OptionBlock;
class OptionModule;
class base_opt;
typedef utility::pointer::owning_ptr< OptionBlock > OptionBlockOP;
typedef utility::pointer::owning_ptr< OptionModule > OptionModuleOP;


 /*################# ###################################
		stuff used for file parsing
		#################################################### */
	class OptionEntry {
	public:
		OptionEntry() {};
		OptionEntry(std::string _name ) : name(_name), bKnown(false), collect_(false), entry_number(-1) {};
		OptionEntry(std::string _name,std::string _val) :
			name(_name),comment(""),bKnown(false), collect_(false), entry_number(-1) { values.push_back(_val); };
		OptionEntry(std::string _name, std::string _val, std::string _comment) :
			name(_name), collect_(false), comment(_comment) { values.push_back(_val); };
		void set_boolean( bool setting ) { truefalse = setting; };
		int process( std::string const& opt, int pos, int &argc, char* argv[] );
	public:
			bool operator> (const OptionEntry &lin) const {return entry_number>lin.entry_number;};
	private:
		friend class OptionBlock;
		base_opt theOption;

		bool bKnown;
		bool truefalse;
		bool collect_;
		int entry_number; /* this is set by when the lines are registered */
	};

class base_opt;

class OptionBlock : public utility::pointer::ReferenceCount {
	public:
		OptionBlock()
			: number_of_known_entries(0), bKnown(false) {};
		OptionBlock(std::string _name) : block_name(_name), number_of_known_entries(0), bKnown(false) {};
		virtual void process_entries( int &argc, char* argv[] );
	//    virtual bool read_entries( FILE *in_file, std::string &next_block );
		virtual void write_entries( std::ostream& out );
		virtual void write_lines( FILE* out );
		const std::string& name() {return block_name;};
		const std::string& descr() {return block_description;};
		void set_descr(const std::string& descr) {block_description=descr;};
		void reg_opt(const base_opt& opt);
		bool isKnown() const {return bKnown;};
	void register_module( OptionModuleOP mm) { my_module =mm;} ;
	protected:
		std::string block_name;
		std::string block_description;
		typedef std::map<std::string,OptionEntry> t_map_lines;
		t_map_lines lines;
private:
	// OptionBlock(const OptionBlock&) {};
	// OptionBlock& operator= (const OptionBlock&) {};
	int  number_of_known_entries;
	bool bKnown;
	OptionModuleOP my_module;
};

class OptionFileBlock : public OptionBlock {
public:
	OptionFileBlock ( std::string name ) : OptionBlock( name ) {};
	bool read_entries(FILE *in_file, std::string &next_block);
	void write_entries(FILE *out_file);
};




class OptionBackend : public   utility::pointer::ReferenceCount {
public:
	virtual void get_options(OptionModule& module) = 0;
	virtual void register_opt(OptionModule const& module);
	virtual void write_options( std::ostream& );
protected:
	virtual OptionBlockOP new_block( std::string mid ) = 0;

typedef std::map<std::string,OptionBlockOP> t_map_blocks;
	typedef std::map<std::string, OptionModuleOP> t_map_modules;
	t_map_modules modules;
	t_map_blocks blocks;
};

class OptionFile : public OptionBackend {
	public:
		/* provide filenames for input and output(control) option files */
	OptionFile(std::string in_file,std::string out_file);
	~OptionFile(); /* destructor calls write_options if not already done manually */
	void get_options(OptionModule& block);
	void write_options(); //* writes options, if called multiple times file is overwritten */
protected:
	void read_options();
	OptionBlockOP new_block( std::string mid ) { return new OptionFileBlock( mid ); }
private:
	std::string in_file;
	std::string out_file;
	bool bFileWritten;
};

class CommandLineOptions : public OptionBackend {
public:
	CommandLineOptions( ) : processed_( false ) { std::cerr << "CommandLineOptions initialized" << std::endl; };
	void process_options( int& argc, char* argv[] );
	void get_options(OptionModule&);
	void write_options() { OptionBackend::write_options(std::cout); };
protected:
	OptionBlockOP new_block( std::string mid ) { return new OptionBlock( mid ); };
private:
	std::string cmdline;
	bool processed_;
};

extern CommandLineOptions cmd_options;

	/* virtual base class, derive classes for your options using
		 the START_OPT,REG,BEGIN_VAR_LIST,END_OPT macros
	*/
class OptionModule : public  utility::pointer::ReferenceCount {
protected:
	OptionModule(std::string _id,std::string _descr) :
		id(_id),module_description(_descr) { };
	virtual OptionModuleOP clone() const = 0;
	virtual OptionModule& operator=(const OptionModule&) { return *this;};
public:
	const std::string get_id() const {return id;};
	virtual ~OptionModule() {};
	virtual void register_block(OptionBlock&) = 0;
private:
	std::string id;
	std::string module_description;
	friend class OptionBackend;
};

	/* the base_opt class tree is used to define
		 the string-2-yourtype mapping for individual option entries
		 basic types are provided, for your own types create a new
		 class XXOPT  : public base_opt {...}
	 */

	class base_opt {
	public:
		base_opt(std::string _id,std::string comment) : id(_id),_comment(comment) {};
		virtual ~base_opt() {};

		const std::string name() const {return id;};
		virtual const std::string full_name() const { return name(); };
		const std::string comment() const {return _comment;};
		std::vector<std::string>& values() { return values; };
		virtual void set_value(std::string value ) const {};
		virtual void set_name(std::string name ) const {};

	protected:
		std::string id;
		std::string _comment;
		std::vector<std::string> values;
	};

	/* no option - prints just a comment */
	class COPT : public base_opt {
	public:
		COPT(std::string comment) : base_opt(";"+comment,"") {};
		void set_value(std::vector< std::string> &values, bool) const { values.clear();};
	};


	/* string - option */
	class SOPT : public base_opt {
	public:
		SOPT(std::string name, std::string &_var, const std::string& __def, std::string comment) :
			base_opt(name,comment),var(_var),def(__def) {};
		void set_value(std::vector< std::string> &values, bool) const {
	if ( values.empty() ) values.push_back(def); var=values.front();};
	private:
		mutable std::string &var;
		const std::string def;
	};


	/* integer - option */
	class IOPT : public base_opt {
	public:
		IOPT(std::string name, int &_var, const int _def,std::string comment) :
			base_opt(name,comment),var(_var),def(_def) {};
		void set_value(std::vector< std::string> &values, bool) const { if (values.empty() ) {
		char buf[20]; sprintf(buf,"%d",def);	values.push_back( std::string(buf) );}
			var=atoi(values.front().c_str()); };
	private:
		mutable int &var;
		const int def;
	};

	/* float - option */
	class FOPT : public base_opt {
	public:
		FOPT(std::string name, float &_var, const float _def,std::string comment) :
			base_opt(name,comment),var(_var),def(_def) {};
		void set_value(std::vector< std::string> &values, bool) const { if (values.empty() ) {
	char buf[20]; sprintf(buf,"%g",def);  values.push_back( std::string(buf) ); }
			var=atof(values.front().c_str());};
	private:
			mutable float &var;
			const float def;
	};

	/* bool - option */
	class YNOPT : public base_opt {
	public:
		YNOPT(std::string name, bool &_var, const bool _def,std::string comment,std::string _true_str, std::string _false_str) :
			base_opt(name,comment),var(_var),def(_def),true_str(_true_str),false_str(_false_str) {};
		void set_value(std::vector< std::string> &values, bool) const { if ( values.empty() ) {
	char buf[20]; sprintf(buf,"%s",def ? true_str.c_str() : false_str.c_str()); values.push_back( std::string(buf) ); }
		var=(values.front()==true_str); };
	private:
		mutable bool &var;
		const bool def;
		std::string true_str;
		std::string false_str;
	};

	/* bool - option */
	class BOPT : public base_opt {
	public:
		BOPT(std::string name, bool &_var, const bool _def,std::string comment) :
			base_opt(name,comment),var(_var),def(_def) { var = def; };
		void set_value(std::vector< std::string> &, bool setting ) const { var = setting; };
		const std::string full_name() { return (var ? "" : "no") + name(); };
	private:
		mutable bool &var;
		const bool def;
		std::string true_str;
		std::string false_str;
	};

}
}
/* definition of macros to create option block */
#define __NS devel::option
#define START_OPT(classname,blockname,descr)	\
	class classname : public __NS::OptionModule {	\
	public:								\
	classname() : __NS::OptionModule(blockname,descr)        {	\
		__NS::cmd_options.get_options( *this );				\
};									\
	__NS::OptionModule& operator= (const __NS::OptionModule &_data) {			\
		const classname &data=dynamic_cast<const classname &>(_data);	\
		return *this=data;							\
	};									\
	protected:								\
	__NS::OptionModuleOP clone() const {return new classname;};	\
	public:								\
	void register_block( __NS::OptionBlock& block) {				\
		using namespace __NS;					\
		OptionModule::register_block(block);

#define REG(object) block.reg_opt( __NS::object);
#define BEGIN_VAR_LIST } public:
#define END_OPT(classname) };				\
	namespace ns_##classname {					\
		class _the_registrator_ {				\
		public:						\
		_the_registrator_() {					\
			__NS::cmd_options.register_opt( classname() ); };		\
		};								\
		static _the_registrator_ _register_;			\
	}





#endif
