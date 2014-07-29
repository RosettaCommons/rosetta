#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/xpressive/xpressive.hpp>

#include <iostream>
using namespace std;

int main()
{
	// Regex

	static const int N = 3;
	string inputs[N] = {"12", "20x", "x"};

	using namespace boost::xpressive;

	smatch match;
	sregex regex = (s1 = +digit) >> (s2 = !as_xpr('x'));

	for(int i = 0; i < N; ++i) {
		if (regex_match(inputs[i], match, regex)) {
			cout << "'" << match[1] << "' '" << match[2] << "'" << endl;
		} else {
			cout << "Bad input: '" << inputs[i] << "'" << endl;
		}
	}

	// Lexical cast

	using boost::lexical_cast;
	using boost::bad_lexical_cast;

	try { cout << "true: " << lexical_cast<bool>("true") << endl; }
	catch (const bad_lexical_cast &) { cout << "Can't cast 'true'." << endl; }

	try { cout << "false: " << lexical_cast<bool>("false") << endl; }
	catch (const bad_lexical_cast &) { cout << "Can't cast 'false'." << endl; }

	try { cout << "yes: " << lexical_cast<bool>("yes") << endl; }
	catch (const bad_lexical_cast &) { cout << "Can't cast 'yes'." << endl; }

	try { cout << "no: " << lexical_cast<bool>("no") << endl; }
	catch (const bad_lexical_cast &) { cout << "Can't cast 'no'." << endl; }

	try { cout << "1: " << lexical_cast<bool>("1") << endl; }
	catch (const bad_lexical_cast &) { cout << "Can't cast '1'." << endl; }

	try { cout << "0: " << lexical_cast<bool>("0") << endl; }
	catch (const bad_lexical_cast &) { cout << "Can't cast '0'." << endl; }


	return 0;
}
