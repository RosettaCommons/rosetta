/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

int main() {
 return CxxTest::ErrorPrinter().run();
}
#include "TraitsTest.h"

static TestFunky suite_TestFunky;

static CxxTest::List Tests_TestFunky = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestFunky( "TraitsTest.h", 53, "TestFunky", suite_TestFunky, Tests_TestFunky );

static class TestDescription_TestFunky_testPets : public CxxTest::RealTestDescription {
public:
 TestDescription_TestFunky_testPets() : CxxTest::RealTestDescription( Tests_TestFunky, suiteDescription_TestFunky, 56, "testPets" ) {}
 void runTest() { suite_TestFunky.testPets(); }
} testDescription_TestFunky_testPets;

#include <cxxtest/Root.cpp>
