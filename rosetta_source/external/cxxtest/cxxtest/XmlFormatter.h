// Licenced under the GLPL, see http://www.gnu.org/licenses/lgpl.html

#ifndef __CXXTEST__XMLFORMATTER_H
#define __CXXTEST__XMLFORMATTER_H

//
// The XmlFormatter is a TestListener that
// prints reports of the errors to an output
// stream in the form of an XML document.
//

// The following definitions are used if stack trace support is enabled,
// to give the traces an easily-parsable XML format.  If stack tracing is
// not enabled, then these definitions will be ignored.
#define CXXTEST_STACK_TRACE_ESCAPE_AS_XML
#define CXXTEST_STACK_TRACE_NO_ESCAPE_FILELINE_AFFIXES

#define CXXTEST_STACK_TRACE_INITIAL_PREFIX "<stack-frame function=\""
#define CXXTEST_STACK_TRACE_INITIAL_SUFFIX "\"/>\n"
#define CXXTEST_STACK_TRACE_OTHER_PREFIX CXXTEST_STACK_TRACE_INITIAL_PREFIX
#define CXXTEST_STACK_TRACE_OTHER_SUFFIX CXXTEST_STACK_TRACE_INITIAL_SUFFIX
#define CXXTEST_STACK_TRACE_ELLIDED_MESSAGE ""
#define CXXTEST_STACK_TRACE_FILELINE_PREFIX "\" location=\""
#define CXXTEST_STACK_TRACE_FILELINE_SUFFIX ""


#include <cxxtest/TestRunner.h>
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/ValueTraits.h>
#include <cxxtest/ErrorFormatter.h>
#include <iostream>
#include <string>
#include <sstream>

namespace CxxTest
{
    class XmlFormatter : public TestListener
    {
        public:
        XmlFormatter( OutputStream *o ) : _o(o) { }

        int run()
        {
                (*_o) << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
                _o->flush();

            TestRunner::runAllTests( *this );
            return tracker().failedTests();
        }

        void enterWorld( const WorldDescription & /*desc*/ )
        {
            (*_o) << "<world>" << endl;
            _o->flush();
        }

        static void totalTests( OutputStream &o )
        {
            char s[WorldDescription::MAX_STRLEN_TOTAL_TESTS];
            const WorldDescription &wd = tracker().world();
            o << wd.strTotalTests( s ) << (wd.numTotalTests() == 1 ? " test" : " tests");
        }

        void enterSuite( const SuiteDescription& desc )
        {
                (*_o) << "    <testsuite name=\"" << desc.suiteName() << "\" ";
                (*_o) << "file=\"" << desc.file() << "\" ";
                (*_o) << "line=\"" << desc.line() << "\"";
                (*_o) << ">"<< endl;
                _o->flush();
        }

        void leaveSuite( const SuiteDescription & )
        {
                (*_o) << "    </testsuite>" << endl;
                _o->flush();
        }

        void enterTest( const TestDescription & desc )
        {
                (*_o) << "        <testcase name=\"" << desc.testName() << "\" ";
                (*_o) << "line=\"" << desc.line() << "\"";
                (*_o) << ">" << endl;
                _o->flush();
        }

        void leaveTest( const TestDescription & )
        {
                (*_o) << "        </testcase>" << endl;
                _o->flush();
        }

        void leaveWorld( const WorldDescription &desc )
        {
                (*_o) << "</world>" << endl;
                _o->flush();
        }

        void trace( const char *file, unsigned line, const char *expression )
        {
            startTag( "trace", file, line );
            attribute( "message", expression );
            endSingletonTag();
        }

        void suiteInitError( const char *file, unsigned line, const char *expression )
        {
            (*_o) << "        <testsuite-error type=\"init\" line=\"" << line << "\">" << endl;
            (*_o) << expression;
            (*_o) << "        </testsuite-error>" << endl;
            _o->flush();
        }

        void warning( const char *file, unsigned line, const char *expression )
        {
            (*_o) << "            <warning line=\"" << line << "\">" << endl;
            (*_o) << expression;
            (*_o) << "            </warning>" << endl;
            _o->flush();
        }

        void failedTest( const char *file, unsigned line, const char *expression )
        {
            testFailure( file, line, "failure", expression );
        }

        void failedAssert( const char *file, unsigned line, const char *expression )
        {
            testFailure( file, line, "failedAssert",
              ( std::string( "Assertion failed: " ) + expression ).c_str() );
        }

        void failedAssertEquals( const char *file, unsigned line,
                                 const char *expectedStr, const char *actualStr,
                                 const char *expected, const char *actual )
        {
            testFailure( file, line, "failedAssertEquals",
                ( std::string( "expected:<" ) + expected + "> but was:<" + actual + ">" ).c_str() );

            // TODO FIXME - use xStr and yStr too, which are the original strings typed into assertequals -
            // Sample output -- Error: Expected (1 == 2 + 3), found (1 != 5)
            // Note: less JUnit-like, but very useful!
        }

        void failedAssertSameData( const char *file, unsigned line,
                                   const char *xStr, const char *yStr, const char *sizeStr,
                                   const void *x, const void *y, unsigned size )
        {
            startTag( "failed-assert-same-data", file, line );
            attribute( "lhs-desc", xStr );
            attributeBinary( "lhs-value", x, size );
            attribute( "rhs-desc", yStr );
            attributeBinary( "rhs-value", y, size );
            attribute( "size-desc", sizeStr );
            attribute( "size-value", size );
            endSingletonTag();
        }

        void failedAssertDelta( const char *file, unsigned line,
                                const char *xStr, const char *yStr, const char *dStr,
                                const char *x, const char *y, const char *d )
        {
            startTag( "failed-assert-delta", file, line );
            attribute( "lhs-desc", xStr );
            attribute( "lhs-value", x );
            attribute( "rhs-desc", yStr );
            attribute( "rhs-value", y );
            attribute( "delta-desc", dStr );
            attribute( "delta-value", d );
            endSingletonTag();
        }

        void failedAssertDiffers( const char *file, unsigned line,
                                  const char *xStr, const char *yStr,
                                  const char *value )
        {
            startTag( "failed-assert-ne", file, line );
            attribute( "lhs-desc", xStr );
            attribute( "rhs-desc", yStr );
            attribute( "value", value );
            endSingletonTag();
        }

        void failedAssertLessThan( const char *file, unsigned line,
                                   const char *xStr, const char *yStr,
                                   const char *x, const char *y )
        {
            startTag( "failed-assert-lt", file, line );
            attribute( "lhs-desc", xStr );
            attribute( "lhs-value", x );
            attribute( "rhs-desc", yStr );
            attribute( "rhs-value", y );
            endSingletonTag();
        }

        void failedAssertLessThanEquals( const char *file, unsigned line,
                                         const char *xStr, const char *yStr,
                                         const char *x, const char *y )
        {
            startTag( "failed-assert-le", file, line );
            attribute( "lhs-desc", xStr );
            attribute( "lhs-value", x );
            attribute( "rhs-desc", yStr );
            attribute( "rhs-value", y );
            endSingletonTag();
        }

        void failedAssertRelation( const char *file, unsigned line,
                                   const char *relation, const char *xStr, const char *yStr,
                                   const char *x, const char *y )
        {
            startTag( "failed-assert-relation", file, line );
            attribute( "relation", relation );
            attribute( "lhs-desc", xStr );
            attribute( "lhs-value", x );
            attribute( "rhs-desc", yStr );
            attribute( "rhs-value", y );
            endSingletonTag();
        }

        void failedAssertPredicate( const char *file, unsigned line,
                                    const char *predicate, const char *xStr, const char *x )
        {
            startTag( "failed-assert-predicate", file, line );
            attribute( "predicate", predicate );
            attribute( "arg-desc", xStr );
            attribute( "arg-value", x );
            endSingletonTag();
        }

        void failedAssertThrows( const char *file, unsigned line,
                                 const char *expression, const char *type,
                                 bool otherThrown )
        {
            startTag( "failed-assert-throws", file, line );
            attribute( "expression", expression );
            attribute( "type", type );
            attribute( "threw", otherThrown ? "other" : "none" );
            endSingletonTag();
        }

        void failedAssertThrowsNot( const char *file, unsigned line, const char *expression )
        {
            startTag( "failed-assert-nothrow", file, line );
            attribute( "expression", expression );
            endSingletonTag();
        }

    protected:
        OutputStream *outputStream() const
        {
            return _o;
        }

    private:
        XmlFormatter( const XmlFormatter & );
        XmlFormatter &operator=( const XmlFormatter & );

        void testFailure( const char *file, unsigned line, const char *failureType, const char *message = NULL )
        {
            std::cerr << "Test " <<  file << " FAILED" << std::endl;
            std::stringstream output;
            const char* tagName = "failure";

            output << startTagText( tagName, file, line );
            if ( message != NULL )
            {
                output << attributeText( "message", message );
            }
            output << attributeText( "type", failureType );
            output << bodyText( file, line, failureType, message );
            output << endTagText( tagName );

            writeOut( output.str() );
        }

        // xxx kill me
        void failedAssertGeneric( const char *file, unsigned line,
            const char *failureType, const char *message = NULL )
        {
            const char* tagName = "failure";

            startTag( tagName, file, line );
            if ( message != NULL )
            {
                attribute( "message", message );
            }
            attribute( "type", failureType );
            body( file, line, failureType, message );
            endTag( tagName );
        }

        static const char * escape(const std::string& str)
        {
            std::string escStr = "";
            for(size_t i = 0; i < str.length(); i++)
            {
                switch(str[i])
                {
                    case '"':  escStr += "&quot;"; break;
                    case '\'': escStr += "&apos;"; break;
                    case '<':  escStr += "&lt;"; break;
                    case '>':  escStr += "&gt;"; break;
                    case '&':  escStr += "&amp;"; break;
                    default:   escStr += str[i]; break;
                }
            }
            return escStr.c_str();
        }

        void startTag( const char* name, const char *file, unsigned line )
        {
            writeOut( startTagText( name, file, line ) );
        }

        static std::string startTagText( const char* name, const char *file, unsigned line )
        {
            std::stringstream output;
            output << "<" << name;
            return output.str();
        }

        void attribute( const char* name, const char *value )
        {
            writeOut( attributeText( name, value ) );
        }

        static std::string attributeText( const char* name, const char *value )
        {
            std::stringstream output;
            output << " " << name;
            output << "=\"" << escape(value) << "\"";
            return output.str();
        }

        void attribute( const char* name, unsigned value )
        {
            (*_o) << name;
            (*_o) << "=\"" << value << "\" ";
        }

        void attributeBinary( const char* name, const void *value, unsigned size )
        {
            (*_o) << name;
            (*_o) << "=\"";
            dump(value, size);
            (*_o) << "\" ";
        }

        void body( const char *file, unsigned line, const char* failureType, const char* message )
        {
            writeOut( bodyText( file, line, failureType, message ) );
        }

        static std::string bodyText( const char *file, unsigned line, const char* failureType, const char* message )
        {
            std::stringstream output;
            output << ">";    //finish the start tag before writing body
            output << escape( failureType );
            if ( message != NULL )
            {
                output << ": " << escape( message );
            }
            output << std::endl;
            output << "at (" << file << ":" << line << ")" << std::endl;

            return output.str();
        }

        void endSingletonTag()
        {
            (*_o) << " />" << endl;
            _o->flush();
        }

        void endTag( const char* name )
        {
            writeOut( endTagText( name ) );
        }

        static std::string endTagText( const char* name )
        {
            std::stringstream output;
            output << "</" << name << ">" << std::endl;
            return output.str();
        }

        void writeOut( std::string output )
        {
            (*_o) << output.c_str();
            _o->flush();
        }

        void dump( const void *buffer, unsigned size )
        {
            unsigned dumpSize = size;
            if ( maxDumpSize() && dumpSize > maxDumpSize() )
                dumpSize = maxDumpSize();

            const unsigned char *p = (const unsigned char *)buffer;
            for ( unsigned i = 0; i < dumpSize; ++ i )
                (*_o) << byteToHex( *p++ ) << " ";
            if ( dumpSize < size )
                (*_o) << "... ";
        }

        static void endl( OutputStream &o )
        {
            OutputStream::endl( o );
        }

        OutputStream *_o;
    };
}

#endif // __CXXTEST__ERRORFORMATTER_H
