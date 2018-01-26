#ifndef __cxxtest__TestSuite_cpp__
#define __cxxtest__TestSuite_cpp__

#include <cxxtest/TestSuite.h>

#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>

namespace CxxTest
{
    //
    // TestSuite members
    //
    TestSuite::~TestSuite() {}
    void TestSuite::setUp() {}
    void TestSuite::tearDown() {}

    //
    // Test-aborting stuff
    //
    static bool currentAbortTestOnFail = false;

    bool abortTestOnFail()
    {
        return currentAbortTestOnFail;
    }

    void setAbortTestOnFail( bool value )
    {
        currentAbortTestOnFail = value;
    }

    void doAbortTest()
    {
#   if defined(_CXXTEST_HAVE_EH)
        if ( currentAbortTestOnFail )
            throw AbortTest();
#   endif // _CXXTEST_HAVE_EH
    }

    //
    // Max dump size
    //
    static unsigned currentMaxDumpSize = CXXTEST_MAX_DUMP_SIZE;

    unsigned maxDumpSize()
    {
        return currentMaxDumpSize;
    }

    void setMaxDumpSize( unsigned value )
    {
        currentMaxDumpSize = value;
    }

    //
    // Some non-template functions
    //
    void doTrace( const char *file, unsigned line, const char *message )
    {
        tracker().trace( file, line, message );
    }

    void doWarn( const char *file, unsigned line, const char *message )
    {
        tracker().warning( file, line, message );
    }

    void doFailTest( const char *file, unsigned line, const char *message )
    {
        tracker().failedTest( file, line, message );
        TS_ABORT();
    }

    void doFailAssert( const char *file, unsigned line,
                       const char *expression, const char *message )
    {
        if ( message )
            tracker().failedTest( file, line, message );
        tracker().failedAssert( file, line, expression );
        TS_ABORT();
    }

    bool sameData( const void *x, const void *y, unsigned size )
    {
        if ( size == 0 )
            return true;

        if ( x == y )
            return true;

        if ( !x || !y )
            return false;

        const char *cx = (const char *)x;
        const char *cy = (const char *)y;
        while ( size -- )
            if ( *cx++ != *cy++ )
                return false;

        return true;
    }

    void doAssertSameData( const char *file, unsigned line,
                           const char *xExpr, const void *x,
                           const char *yExpr, const void *y,
                           const char *sizeExpr, unsigned size,
                           const char *message )
    {
        if ( !sameData( x, y, size ) ) {
            if ( message )
                tracker().failedTest( file, line, message );
            tracker().failedAssertSameData( file, line, xExpr, yExpr, sizeExpr, x, y, size );
            TS_ABORT();
        }
    }

	/// Local helper functions.
	static std::string msg_(const char *message, std::string message2)
	{
		if( message ) return message;
		else return message2;
	}
	/// @details split given std::string using ' ' or '\n' symbols.
	static std::vector< std::string > split(const std::string &s)
	{
		std::vector<std::string> r;
		unsigned int start=0, i=0;
		while( start < s.size() ) {
			if( s[i] == ' ' || s[i]== '\n' ) {
				r.push_back( std::string(s.begin()+start, s.begin()+i) );
				start = i+1;
			}
			i++;
			if( i == s.size() ) {
				r.push_back( std::string(s.begin()+start, s.begin()+i) );
				break;
			}
		}

		for(std::vector< std::string >::iterator it=r.end(); it!=r.begin(); ) {  /// removing empty lines
			it--;
			if( it->size() == 0 ) r.erase( it );
		}
		return r;
	}

	void doAssertStringContains(
		const char *file, unsigned line,
		const char *xExpr, std::string const &x,
		const char *yExpr, std::string const &y,
		const char *message )
	{
		bool contains = x.find(y) != std::string::npos;

        if ( !contains ) {
            if ( message )
                tracker().failedTest( file, line, message );
            tracker().failedAssertContains( file, line, xExpr, yExpr, TS_AS_STRING(x), TS_AS_STRING(y) );
            TS_ABORT();
        }
	}


	/// New, for MiniRosetta unit tests
	/// Assert for comparing files (ie file eq)
    void doAssertFileEQ( const char *file, unsigned line,
                         const char *file1Expr, const char *file1,
                         const char *file2Expr, const char *file2,
                         const char *message )
	{
		using std::string;
		std::ifstream hFile1(file1);
		std::ifstream hFile2(file2);

		if( !hFile1.good() ) {
			doFailTest(file, line, msg_(message, "Unable to open file: " + string(file1Expr) + " for reading!").c_str() );
			return;
		}
		if( !hFile2.good() ) {
			doFailTest(file, line, msg_(message, "Unable to open file: " + string(file2Expr) + " for reading!").c_str() );
			return;
		}

		std::string line1, line2;
		while( getline(hFile1, line1) ) {
			if( getline(hFile2, line2) ) {
				if( line1 != line2 ) {
					string m = "Files are not equal!\nLine1: "+line1+"\nLine2: "+line2+"\n";
					doFailTest(file, line, m.c_str() );
					break;
				}
			} else {
				doFailTest(file, line, msg_(message, "Files have different size!").c_str() );
				return;
			}
		}
		if( getline(hFile2, line2) ) {
			doFailTest(file, line, msg_(message, "Files have different size!").c_str() );
			return;
		}
	}

	/// New, for MiniRosetta unit tests
	/// Load file, treat it as number separated by spaces, and return vector of double.
	void loadFileAsDoubleVector(const char *sourceFile, unsigned line,
								const char *filename, std::vector<double> & result,
								const char *message)
	{
		/// reading file
		std::ifstream file(filename, std::ios::in | std::ios::binary | std::ios::ate);
		if( file.is_open() )  {
			std::string res;

			int fsize = file.tellg();
			res.resize( fsize );
			file.seekg(0, std::ios::beg);
			file.read(&res[0], fsize);
			file.close();

			/// Spliting file in to separate strings using ' ' as separator.
			std::vector< std::string > dv = split(res);

			result.resize(0);
			for(unsigned int i=0; i<dv.size(); i++) {
				result.push_back( atof( dv[i].c_str() ) );
			}
		}
		else {
			doFailTest(sourceFile, line, msg_(message, "Unable to open file: " + std::string(filename) + " for reading!").c_str() );
		}
	}

	/// Assert for comparing files as double vectors/
	/// abs_p - specify absolute presision
	/// rel_p - specify relative precision
	/// resulting presision for each number calcualted as max(abs_p, rel_p)
    void doAssertFileEQ_AsDouble( const char *file, unsigned line,
                         const char *file1Expr, const char *file1,
                         const char *file2Expr, const char *file2,
						 double abs_p, double rel_p,
                         const char *message )
	{
		std::vector<double> f1, f2;
		loadFileAsDoubleVector(file, line, file1, f1, message);
		loadFileAsDoubleVector(file, line, file2, f2, message);

		if( f1.size() != f2.size() ) {
			doFailTest(file, line, msg_(message,
				"Files: "+std::string(file1Expr)+" and "+
				std::string(file2Expr)+" have different size!").c_str() );
			return;
		}
		for(unsigned int i=0; i<f1.size(); i++) {
			double p = (f1[i]+f2[i])/2. * rel_p;
			p = std::max(abs_p, p);
			//_TS_ASSERT_DELTA(file, line, f1[i], f2[i], p, "");
			std::ostringstream m1, m2, mp;
			m1 << file1Expr << "[" << i << "]";
			m2 << file2Expr << "[" << i << "]";
			mp << "Presision: a(" << abs_p << ") or r(" << rel_p << ")";
			doAssertDelta(file, line,
				  m1.str().c_str(), f1[i],
				  m2.str().c_str(), f2[i],
				  mp.str().c_str(), p,
                  message);
		}
	}


    void doFailAssertThrows( const char *file, unsigned line,
                             const char *expr, const char *type,
                             bool otherThrown,
                             const char *message )
    {
        if ( message )
            tracker().failedTest( file, line, message );

        tracker().failedAssertThrows( file, line, expr, type, otherThrown );
        TS_ABORT();
    }

    void doFailAssertThrowsNot( const char *file, unsigned line,
                                const char *expression, const char *message )
    {
        if ( message )
            tracker().failedTest( file, line, message );

        tracker().failedAssertThrowsNot( file, line, expression );
        TS_ABORT();
    }
}

#endif // __cxxtest__TestSuite_cpp__
