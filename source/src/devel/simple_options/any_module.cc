// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#define CPLUSPLUS
#include "option.hh"

using namespace devel;
using namespace devel::option;
using namespace std;


/* the example class my_test_options shows
	 how the subclasses created by the START_OPT,REG,END_OPT actually look like in full C++,
	 however, for standard applications it is easier to use the provided macros
*/
class my_test_options : public devel::option::OptionModule {
	public:
	 /* intialize the parent-class with "module_id" and a description of the module */
	 my_test_options() : OptionModule("test_options","this is just for testing"),number(0) {
		 cmd_options.register_opt(my_test_options());
	 };

	 /* provide a copy operator */
	 OptionModule& operator= (const OptionModule &_data) {
	const my_test_options &data=dynamic_cast<const my_test_options &>(_data);
	return *this=data;
	 };

 protected:
	 /* provide a polymorhpic generator for instances of this class */
	 OptionModuleOP clone() const{return new my_test_options;};

	 /* the function that registers all the options and connects them with the variables */
	 void register_block( OptionBlock&);
			/* the option fields */
	public:
			std::string title;
			std::string cpp;
			std::string include;
			std::string define;
			int number;
	 /* the options are left public and mutable
			this is done to enable direct and
			easy read out by the module-code
			we trust the code-developer
			not to change the values manually

			it would be nice if the fields could be made constant
			use (const_cast to change them in register ?
					 but gives problem with copy operator=() )

			NOTE: all types used here should
			have a copy-semantic, i.e.,
			operator=() should be overloaded.

	 */

 };


void my_test_basic::options::register_block(OptionBlock& block) {
	OptionModule::register_block(block);
	block.reg_opt(COPT("VARIOUS TEST OPTIONS"));
	block.reg_opt(SOPT("title",title,"ubiquitin","this is the title"));
	block.reg_opt(IOPT("cycles",number,3,"this defines the number of runs"));
};


START_OPT(module2_opt,"block1","another test module")
	REG(COPT("NOCH MEHR KOMMENTAR"))
	REG(SOPT("title",title,"another","this is the title"));
	REG(IOPT("bla",number,3,"this d"));
REG(BOPT("have coffee",bCoffee,true,"this should always be switched on","yes","no"));
REG(FOPT("temp",temp,1.23,"control the monte-carlo sampling"));
BEGIN_VAR_LIST
std::string title;
int number;
bool bCoffee;
float temp;
END_OPT(module2_opt)


using namespace devel::option;

int main( int argc, char** argv) {


	std::cout << "vorher " << std::endl;
	for (int i=0;i<argc;i++)
		std::cout << argv[i] << " " << std::endl;


	std::cout << "nacher " << std::endl;
	for (int i=0;i<argc;i++)
		std::cout << argv[i] << " " << std::endl;

}


/*
int main(int argc, char **argv) {
		OptionFile option_file("in.mdp","out.mdp");

		my_test_options op1;
		my_test_options op2;

		option_file.get_options(op1);
		option_file.get_options(op2);
		option_file.write_options();
		module2_opt op3;
		option_file.get_options(op3);
		cout << "op1 " << op1.number << "    op2 " << op2.number << endl;
		cout << "block1.number " << op3.number << endl;
		ObjexxFCL::FArray2D<double> M;
		ObjexxFCL::FArray2D<double> eigvec;
		ObjexxFCL::FArray1D<double> eigval;
		LapackFunctions<double,Symmetric<double> >::eigsolv(M,eigvec,eigval);
		Seigsolv_<double,Symmetric<double> > (M,eigvec,eigval);

};*/
