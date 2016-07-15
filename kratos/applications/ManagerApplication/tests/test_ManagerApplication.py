# import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ManagerApplication import *        #and now our application
from KratosMultiphysics.HenryApplication import *   
   
# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest
   
# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import HenrySystemTests as THenrySystemTests
from SmallTests import HenryUnitaryTests as THenryUnitaryTests
   
## NIGTHLY TESTS
from NightlyTests import HenryHardTests as THenryHardTests
   
def AssambleTestSuites():
     #Populates the test suites to run.
   
     #Populates the test suites to run. At least, it should pupulate the suites:
     #"small", "nighlty" and "all"
   
     #Return
     #------
   
     #suites: A dictionary of suites
     #The set of suites with its test_cases added.
       
     suites = KratosUnittest.KratosSuites
   
     # Create a test suit with the selected tests (Small tests):
     smallSuite = suites['small']
     smallSuite.addTest(THenryUnitaryTests('test_Henry'))
     #smallSuite.addTest(THenrySystemTests('test_Henry'))

   
     # Create a test suit with the selected tests plus all small tests
     nightSuite = suites['nightly']
     nightSuite.addTests(smallSuite)
     #nightSuite.addTest(THenryHardTests('test_Henry'))
   
     # Create a test suit that contains all the tests:
     allSuite = suites['all']
     allSuite.addTests(
     KratosUnittest.TestLoader().loadTestsFromTestCases([THenryUnitaryTests]))
   
     return suites
   
if __name__ == '__main__':
     KratosUnittest.runTests(AssambleTestSuites())
