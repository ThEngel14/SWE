/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/TestMain.h>
#include <cxxtest/ErrorPrinter.h>

int main( int argc, char *argv[] ) {
 int status;
    CxxTest::ErrorPrinter tmp;
    CxxTest::RealWorldDescription::_worldName = "cxxtest";
    status = CxxTest::Main< CxxTest::ErrorPrinter >( tmp, argc, argv );
    return status;
}
bool suite_DimensionalSplittingTest_init = false;
#include "/home/vmuser/Tsunami/SWE/src/unit tests/DimensionalSplittingTest.h"

static DimensionalSplittingTest suite_DimensionalSplittingTest;

static CxxTest::List Tests_DimensionalSplittingTest = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_DimensionalSplittingTest( "DimensionalSplittingTest.h", 16, "DimensionalSplittingTest", suite_DimensionalSplittingTest, Tests_DimensionalSplittingTest );

static class TestDescription_suite_DimensionalSplittingTest_testDimensionalSplitting : public CxxTest::RealTestDescription {
public:
 TestDescription_suite_DimensionalSplittingTest_testDimensionalSplitting() : CxxTest::RealTestDescription( Tests_DimensionalSplittingTest, suiteDescription_DimensionalSplittingTest, 33, "testDimensionalSplitting" ) {}
 void runTest() { suite_DimensionalSplittingTest.testDimensionalSplitting(); }
} testDescription_suite_DimensionalSplittingTest_testDimensionalSplitting;

#include <cxxtest/Root.cpp>
const char* CxxTest::RealWorldDescription::_worldName = "cxxtest";
