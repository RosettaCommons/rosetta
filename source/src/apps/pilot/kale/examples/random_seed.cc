#include <devel/init.hh>
#include <numeric/random/random.hh>

using namespace std;

int main(int argc, char** argv) {
	devel::init(argc, argv);
	cout << numeric::random::rg().get_seed() << endl;
}
