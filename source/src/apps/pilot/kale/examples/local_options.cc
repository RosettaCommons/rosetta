#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <iostream>

using namespace std;
using namespace basic::options;

OPT_2GRP_KEY(String, kale, in, pdb)

// Application {{{1

int main(int argc, char* argv []) {
	NEW_OPT(kale::in::pdb, "A friendly greeting", "Hello world!");

	devel::init(argc, argv);

	cout << option[OptionKeys::kale::in::pdb]() << endl;
}
// }}}1
