#ifndef __cxxtest__ErrorFormatter_h__
#define __cxxtest__ErrorFormatter_h__

//
// The ErrorFormatter is a TestListener that
// prints reports of the errors to an output
// stream.  Since we cannot rely ou the standard
// iostreams, this header defines a base class
// analogout to std::ostream.
//

#include <cxxtest/TestRunner.h>
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/ValueTraits.h>

#include <fstream>

extern int real_command_line_argc;  extern char ** real_command_line_argv;
extern int command_line_argc;       extern char ** command_line_argv;

namespace CxxTest
{
    class OutputStream
    {
    public:
        virtual ~OutputStream() {}
        virtual void flush() {};
        virtual OutputStream &operator<<( unsigned /*number*/ ) { return *this; }
        virtual OutputStream &operator<<( const char * /*string*/ ) { return *this; }

        typedef void (*Manipulator)( OutputStream & );

        virtual OutputStream &operator<<( Manipulator m ) { m( *this ); return *this; }
        static void endl( OutputStream &o ) { (o << "\n").flush(); }
    };

    class ErrorFormatter : public TestListener
    {
    public:
        ErrorFormatter( OutputStream *o, const char *preLine = ":", const char *postLine = "" ) :
            _dotting( true ),
            _reported( false ),
            _o(o),
            _preLine(preLine),
            _postLine(postLine)
        {
        }

        int run()
        {
	    // Initially we write json report and mark tests as failed and later owerwrite it with real results. This is needed for cases when test terminate with hard exceptions
	    std::string file_name;
	    if( run_type == _OneSuite_ ) {
		file_name = std::string(real_command_line_argv[0]) + '.' + std::string(real_command_line_argv[1]) + ".yaml";
		write_json(file_name, getAllTestsNames(), getAllSuiteNames( std::string(real_command_line_argv[1]) ));
	    }
	    else {
		file_name = std::string(real_command_line_argv[0]) + ".yaml";
		write_json(file_name, getAllTestsNames(), getAllTestsNames());
	    }

	    TestRunner::runAllTests( *this );
	    write_json(file_name, getAllTestsNames(), failed_tests_);

            return tracker().failedTests();
        }

		void write_json(std::string file_name,
					   std::vector< std::string > const & all_tests,
					   std::vector< std::string > const & failed_tests)
		{
			std::fstream f( file_name.c_str(), std::fstream::out);

			f << "{" << std::endl;
			f << "  \"ALL_TESTS\" : [";
			std::string delim = "";
			for(unsigned int i=0; i<all_tests.size(); i++) {
			    f << delim;  delim = ", ";
			    f << all_tests[i];
			}
			f << "]," <<std::endl;

			f << "  \"FAILED_TESTS\" : [";
			delim = "";
			for(unsigned int i=0; i<failed_tests.size(); i++) {
			    f << delim;  delim = ", ";
			    f << failed_tests[i];
			}
			f << "]" << std::endl;
			f << "}" << std::endl;
			f.close();
		}

		std::vector< std::string > getAllTestsNames() {
			std::vector< std::string > res;

			RealWorldDescription wd;
			for( SuiteDescription *sd = wd.firstSuite(); sd; sd = sd->next() ) {
				for( TestDescription *td = sd->firstTest(); td; td = td->next() ) {
					res.push_back( std::string("\"") + sd->suiteName() + ":" + td->testName()+"\"" );
				}
			}
			return res;
		}

		std::vector< std::string > getAllSuiteNames(std::string suite) {
			std::vector< std::string > res;

			RealWorldDescription wd;
			for( SuiteDescription *sd = wd.firstSuite(); sd; sd = sd->next() ) {
				if( sd->suiteName() == suite ) {
				    for( TestDescription *td = sd->firstTest(); td; td = td->next() ) {
				    	res.push_back( std::string("\"") + sd->suiteName() + ":" + td->testName()+"\"" );
				    }
				}
			}
			return res;
		}





        void enterWorld( const WorldDescription & /* desc*/ )
        {
			//if( real_command_line_argc != command_line_argc ) { // check if first command line arg was test_*
			if( run_type == _NormalRun_ ) (*_o) << "Running " << totalTests << endl;
			if( run_type == _OneTest_ )	(*_o) << "Running one test: " << real_command_line_argv[1] << endl;
			if( run_type == _OneSuite_ ) (*_o) << "Running one suite: " << real_command_line_argv[1] << endl;

            _o->flush();
            _dotting = true;
            _reported = false;
        }

        static void totalTests( OutputStream &o )
        {
            char s[WorldDescription::MAX_STRLEN_TOTAL_TESTS];
            const WorldDescription &wd = tracker().world();
            o << wd.strTotalTests( s ) << (wd.numTotalTests() == 1 ? " test" : " tests");
        }

        void enterSuite( const SuiteDescription & desc )
        {
			(*_o) << "Test suite: " << desc.suiteName() << " (" << desc.file() << ")" << endl;
            _reported = false;
        }

        void enterTest( const TestDescription & /* desc */)
        {
			/* (*_o) << "- " << desc.testName() << "..." << endl;
			_o->flush();
			*/
            _reported = false;
        }

        void leaveTest( const TestDescription & desc)
        {
            if ( !tracker().testFailed() && !tracker().warnings() ) {
                /* ((*_o) << "- " << desc.testName() << "... Ok." << endl).flush();
				((*_o) << ".").flush();
				*/
                _dotting = true;
            }
			else {
                (*_o) << "CXXTEST_ERROR: " << desc.testName() << " Failed!" << endl;
				_o->flush();
				failed_tests_.push_back( std::string("\"") + desc.suiteName() +":"+desc.testName() + std::string("\""));
			}
        }

        void leaveWorld( const WorldDescription &desc )
        {
            if ( !tracker().failedTests() ) {
                (*_o) << "All tests passed!" << endl;
                return;
            }
            newLine();
            (*_o) << "Failed " << tracker().failedTests() << " of " << totalTests << endl;
            unsigned numPassed = desc.numTotalTests() - tracker().failedTests();
            (*_o) << "Success rate: " << (numPassed * 100 / desc.numTotalTests()) << "%" << endl;
        }

        void trace( const char *, unsigned, const char *expression )
        {
            (*_o) << "Trace: " << expression << endl;
        }

        void warning( const char *, unsigned, const char *expression )
        {
            (*_o) << "Warning: " << expression << endl;
        }

        void failedTest( const char *file, unsigned line, const char *expression )
        {
            stop( file, line ) << "Error: Test failed: " << expression << endl;
        }

        void failedAssert( const char *file, unsigned line, const char *expression )
        {
            stop( file, line ) << "Error: Assertion failed: " << expression << endl;
        }

        void failedAssertEquals( const char *file, unsigned line,
                                 const char *xStr, const char *yStr,
                                 const char *x, const char *y )
        {
            stop( file, line ) << "Error: Expected (" << xStr << " == " << yStr << "), found (" <<
                x << " != " << y << ")" << endl;
        }

        void failedAssertSameData( const char *file, unsigned line,
                                   const char *xStr, const char *yStr,
                                   const char *sizeStr, const void *x,
                                   const void *y, unsigned size )
        {
            stop( file, line ) << "Error: Expected " << sizeStr << " (" << size << ") bytes to be equal at (" <<
                xStr << ") and (" << yStr << "), found:" << endl;
            dump( x, size );
            (*_o) << "     differs from" << endl;
            dump( y, size );
        }

        void failedAssertDelta( const char *file, unsigned line,
                                const char *xStr, const char *yStr, const char *dStr,
                                const char *x, const char *y, const char *d )
        {
            stop( file, line ) << "Error: Expected (" <<
                xStr << " == " << yStr << ") up to " << dStr << " (" << d << "), found (" <<
                x << " != " << y << ")" << endl;
        }

        void failedAssertDiffers( const char *file, unsigned line,
                                  const char *xStr, const char *yStr,
                                  const char *value )
        {
            stop( file, line ) << "Error: Expected (" <<
                xStr << " != " << yStr << "), found (" <<
                value << ")" << endl;
        }

        void failedAssertLessThan( const char *file, unsigned line,
                                   const char *xStr, const char *yStr,
                                   const char *x, const char *y )
        {
            stop( file, line ) << "Error: Expected (" <<
                xStr << " < " << yStr << "), found (" <<
                x << " >= " << y << ")" << endl;
        }

        void failedAssertLessThanEquals( const char *file, unsigned line,
                                         const char *xStr, const char *yStr,
                                         const char *x, const char *y )
        {
            stop( file, line ) << "Error: Expected (" <<
                xStr << " <= " << yStr << "), found (" <<
                x << " > " << y << ")" << endl;
        }

        void failedAssertRelation( const char *file, unsigned line,
                                   const char *relation, const char *xStr, const char *yStr,
                                   const char *x, const char *y )
        {
            stop( file, line ) << "Error: Expected " << relation << "( " <<
                xStr << ", " << yStr << " ), found !" << relation << "( " << x << ", " << y << " )" << endl;
        }

        void failedAssertPredicate( const char *file, unsigned line,
                                    const char *predicate, const char *xStr, const char *x )
        {
            stop( file, line ) << "Error: Expected " << predicate << "( " <<
                xStr << " ), found !" << predicate << "( " << x << " )" << endl;
        }

        void failedAssertThrows( const char *file, unsigned line,
                                 const char *expression, const char *type,
                                 bool otherThrown )
        {
            stop( file, line ) << "Error: Expected (" << expression << ") to throw (" <<
                type << ") but it " << (otherThrown ? "threw something else" : "didn't throw") <<
                endl;
        }

        void failedAssertThrowsNot( const char *file, unsigned line, const char *expression )
        {
            stop( file, line ) << "Error: Expected (" << expression << ") not to throw, but it did" <<
                endl;
        }

	public:
		std::vector< std::string > failed_tests_;

    protected:
        OutputStream *outputStream() const
        {
            return _o;
        }

    private:
        ErrorFormatter( const ErrorFormatter & );
        ErrorFormatter &operator=( const ErrorFormatter & );

        OutputStream &stop( const char *file, unsigned line )
        {
            newLine();
            reportTest();
            return (*_o) << file << _preLine << line << _postLine << ": ";
        }

        void newLine( void )
        {
            if ( _dotting ) {
                (*_o) << endl;
                _dotting = false;
            }
        }

        void reportTest( void )
        {
            if( _reported )
                return;
            (*_o) << "In " << tracker().suite().suiteName() << "::" << tracker().test().testName() << ":" << endl;
            _reported = true;
        }

        void dump( const void *buffer, unsigned size )
        {
            if ( !buffer )
                dumpNull();
            else
                dumpBuffer( buffer, size );
        }

        void dumpNull()
        {
            (*_o) << "   (null)" << endl;
        }

        void dumpBuffer( const void *buffer, unsigned size )
        {
            unsigned dumpSize = size;
            if ( maxDumpSize() && dumpSize > maxDumpSize() )
                dumpSize = maxDumpSize();

            const unsigned char *p = (const unsigned char *)buffer;
            (*_o) << "   { ";
            for ( unsigned i = 0; i < dumpSize; ++ i )
                (*_o) << byteToHex( *p++ ) << " ";
            if ( dumpSize < size )
                (*_o) << "... ";
            (*_o) << "}" << endl;
        }

        static void endl( OutputStream &o )
        {
            OutputStream::endl( o );
        }

        bool _dotting;
        bool _reported;
        OutputStream *_o;
        const char *_preLine;
        const char *_postLine;
    };
}

#endif // __cxxtest__ErrorFormatter_h__
