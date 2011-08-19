// Licenced under the GLPL, see http://www.gnu.org/licenses/lgpl.html

#ifndef __cxxtest__XmlPrinter_h__
#define __cxxtest__XmlPrinter_h__

//
// The XmlPrinter is a simple TestListener that
// prints JUnit style xml to the output stream
//


#include <cxxtest/Flags.h>

#ifndef _CXXTEST_HAVE_STD
#   define _CXXTEST_HAVE_STD
#endif // _CXXTEST_HAVE_STD

#include <cxxtest/XmlFormatter.h>
#include <cxxtest/StdValueTraits.h>

#ifdef _CXXTEST_OLD_STD
#   include <iostream.h>
#else // !_CXXTEST_OLD_STD
#   include <iostream>
#endif // _CXXTEST_OLD_STD

namespace CxxTest
{
    class XmlPrinter : public XmlFormatter
    {
    public:
        XmlPrinter( CXXTEST_STD(ostream) &o = CXXTEST_STD(cout), const char *preLine = ":", const char *postLine = "" ) :
            XmlFormatter( new Adapter(o) ) {}
        virtual ~XmlPrinter() { delete outputStream(); }

    private:
        class Adapter : public OutputStream
        {
            CXXTEST_STD(ostream) &_o;
        public:
            Adapter( CXXTEST_STD(ostream) &o ) : _o(o) {}
            void flush() { _o.flush(); }
            OutputStream &operator<<( const char *s ) { _o << s; return *this; }
            OutputStream &operator<<( Manipulator m ) { return OutputStream::operator<<( m ); }
            OutputStream &operator<<( unsigned i )
            {
                char s[1 + 3 * sizeof(unsigned)];
                numberToString( i, s );
                _o << s;
                return *this;
            }
        };
    };
}

#endif // __cxxtest__XmlPrinter_h__
