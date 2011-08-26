#ifndef __cxxtest__ProcessCommandLine_h__
#define __cxxtest__ProcessCommandLine_h__

//
// The ProcessCommandLine is a function to remove test arguments from a command line represented
// by command_line_argc and command_line_argv variables.
// Test arguments consider to be: 'test_*'
//

#include <string>

#include <iostream>

extern int real_command_line_argc;  extern char ** real_command_line_argv;
extern int command_line_argc;       extern char ** command_line_argv;

enum RunType { _NormalRun_, _OneTest_, _OneSuite_, _ListAllTests_};

//#include <iostream>

namespace CxxTest
{
	extern RunType run_type;

	/// @brief remove first comman line argument from and store results in command_line_arg[c/v].
	inline void __remove_first_command_line_arg(void)
	{
		if( real_command_line_argc > 1 ) {
			command_line_argc = real_command_line_argc - 1;
			command_line_argv = new char * [ command_line_argc ];

			command_line_argv[0] = real_command_line_argv[0];
			for(int i=1; i<command_line_argc; i++) command_line_argv[i] = real_command_line_argv[i+1];
		}
	}


	inline void ProcessCommandLine(void)
	{
		if( real_command_line_argc > 1 ) {
			/*
			// testing if first command line agrment is 'test_*'
			// if( strcmp(real_command_line_argv[1], "all") ) {
			std::string test_("test_");
			std::string s( real_command_line_argv[1] );
			if( s.size() > test_.size() ) s.resize( test_.size() );

			if( s == test_ ) { /// shcecking if run mode should be _OneTest_
				run_type = _OneTest_;

				/// deleting first arg from command line.
				command_line_argc = real_command_line_argc - 1;
				command_line_argv = new char * [ command_line_argc ];

				command_line_argv[0] = real_command_line_argv[0];
				for(int i=1; i<command_line_argc; i++) command_line_argv[i] = real_command_line_argv[i+1];
			}
			*/
			std::string s( real_command_line_argv[1] );

			RealWorldDescription wd;

			if( s == "_ListAllTests_" ) {
				run_type = _ListAllTests_;
				__remove_first_command_line_arg();

				for( SuiteDescription *sd = wd.firstSuite(); sd; sd = sd->next() ) {
					std::cout << sd->suiteName() << " ";
					//if( s == sd->suiteName() ) {
				}
				std::cout << std::endl;

				return;
			}
			else {
				for( SuiteDescription *sd = wd.firstSuite(); sd; sd = sd->next() ) {
					//std::cout << "\n" << sd->suiteName() << "\n";
					if( s == sd->suiteName() ) {
						// Running one suite
						run_type = _OneSuite_;
						__remove_first_command_line_arg();
						return;
					}
					for( TestDescription *td = sd->firstTest(); td; td = td->next() ) {
						if( s == std::string(sd->suiteName()) + ":" + td->testName() ) {
							// Running one test
							run_type = _OneTest_;
							__remove_first_command_line_arg();
							return;
						}
					}
				}
			}
		}
	}
}

#endif // __cxxtest__ProcessCommandLine_h__
