# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

import Kratos_Execute_Henry_Hard_Test as Execute_Test

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
  
class HenryHardTests(KratosUnittest.TestCase):

    def setUp(self):
      #self.file_name = "henry_Test/henry_test_flow"
      ## Initialize GiD  I/O
      #input_file_name = GetFilePath(self.file_name)
      #parameter_file = open(self.file_name +"_parameters.json",'r')
      #ProjectParameters = Parameters( parameter_file.read())
      ## Creating the model part
      #self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)
      with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
          self.test = Execute_Test.Kratos_Execute_Test()
      
    def test_Henry(self):
      with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
          self.test.Solve()
    
      
    def tearDown(self):
        pass
    
#class DynamicNewmarkTests(KratosUnittest.TestCase):

    #def setUp(self):
      #self.file_name = "dynamic_test/dynamic_newmark_test"
      # Initialize GiD  I/O
      #input_file_name = GetFilePath(self.file_name)
      #parameter_file = open(self.file_name +"_parameters.json",'r')
      #ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      #self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)
      
    #def test_Newmark(self):
      #self.test.Solve()
      
    #def tearDown(self):
        #pass

