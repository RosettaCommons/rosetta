// Unit headers
#include <apps/pilot/kale/kic_refactor/KicSandbox.hh>

// C++ headers
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
	apps::pilot::KicSandbox kic(argc, argv);
	kic.protocol->apply(kic.pose);
}

