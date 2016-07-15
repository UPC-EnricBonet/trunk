#setting the domain size for the problem to be solved
  

#including kratos path
import sys
from KratosMultiphysics import *    #we import the KRATOS  
from KratosMultiphysics.ManagerApplication import *        #and now our application
from KratosMultiphysics.HenryApplication import *        #and now our application


import json
import process_factory
import fractional_iterative_solver
import time as timer


class Kratos_Execute_Test:

    def __init__(self):

        self.domain_size = 2  # 2D problem
        #Create a new empty model part called "ExampleModelPart"
        self.model_part = ModelPart("ExampleModelPart");

        print ("Model part defined: ExampleModelPart")  

        #we import the python file that includes the commands that we need




        #import variables that we will need from solver to our recent created model_part            
        fractional_iterative_solver.AddVariables(self.model_part) 

        # (note that our model part does not have nodes or elements yet) 

        # introducing input & outoput (also postProcess) file name
        output_file_name = "henry/henryUnitaryTest/result_henry"
        input_file_name = "henry/henryUnitaryTest/henry"
        # Mesh built by GID for the postProcess
        gid_mode = GiDPostMode.GiD_PostAscii
        multifile = MultiFileFlag.MultipleFiles
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteElementsOnly
        self.gid_io = GidIO(output_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
        self.model_part_io_fluid = ModelPartIO(input_file_name)
        self.model_part_io_fluid.ReadModelPart(self.model_part)  


        #the buffer size should be set up after the mesh is read for the first time (Transient problem =2,3. Steady problem =1) 
        self.model_part.SetBufferSize(2)



        #################################### STRATEGY (INITIALIZE AND READ) ##############################################################################


       
        # we add the DoFs  
        fractional_iterative_solver.AddDofs(self.model_part)

        print ("Time steps values on each time for unkown=2 (Buffer size)") 

        #creating flow solver (custom)
        self.solverFIS = fractional_iterative_solver.FractionalIterativeSolver(self.model_part,self.domain_size)
        #we import the C++ file that includes the commands that we need
        #import fractional_iterative_strategy

        ## This part is contained is for the FractionalSolver not included at fractional_iterative_solver.py
        self.solverFIS.time_order = 1
        self.solverFIS.echo_level = 3

        print ("Iterative solver create succesfully") 

        self.solverFIS.Initialize()

        #self.Nsteps=1


        print ("Solver inicializate!")    

        #solver.calculateDensityNodes()
        #Arhivo de comparación:
        filePost="henry/henryUnitaryTest/UnitaryTest.txt"
        self.solverFIS.ReadFile(filePost)

        #clean file output for matrix
        open('MatrixSystem.txt', 'w').close()


        ###################################### ELEMENT (INITIALIZE AND READ) ###########################################################   
        self.element = UnitaryTestHenryECU(self.model_part,self.domain_size)
        
  
        ##################################### UTILITIES (INITIALIZE AND READ) ###########################################################
        #self.utilities = UnitaryTestHenryECU(self.model_part,self.domain_size)
       

        ##################################### CONDITIONS (INITIALIZE AND READ) ###########################################################
        #self.condition = UnitaryTestHenryECU(self.model_part,self.domain_size)
        


    def ExecuteFinalizeSolutionStep(self):
        pass


    def Solve(self):        
        #self.solverFIS.Solve()
        self.element.UnitaryTest()#poner UnitaryTestElement()
        #self.utilities.UnitaryTestUtilities()
        #self.condition.UnitaryTestConditions() 
 
    def AssertVariables(self,arg):
        for elem in self.model_part.Elements:
            value1 = elem.GetValue(DENSITY_ELEM)+5.0
            value2 = elem.GetValue(DENSITY_ELEM_TARGET)
            arg.assertAlmostEqual(value1,value2[0],2)



