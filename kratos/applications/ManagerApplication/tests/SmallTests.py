# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

import Kratos_Execute_Henry_Test as Execute_Test
import Kratos_Execute_Henry_UnitaryTest as Execute_Unitarytest

def GetFilePath(fileName):
  return os.path.dirname(__file__) + "/" + fileName

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)
  
class HenrySystemTests(KratosUnittest.TestCase):

    def setUp(self):

      with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
      	self.test = Execute_Test.Kratos_Execute_Test()
      
    def test_Henry(self):
      with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
        self.test.Solve()
    
      
    def tearDown(self):
        pass
    
class HenryUnitaryTests(KratosUnittest.TestCase):

    def setUp(self):
      with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
        self.test = Execute_Unitarytest.Kratos_Execute_Test()
      
    def test_Henry(self):
      with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
        self.test.Solve()
        self.test.AssertVariables(self)
      
    #def tearDown(self):
        #pass

