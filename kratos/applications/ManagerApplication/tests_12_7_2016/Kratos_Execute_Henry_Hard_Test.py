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
        output_file_name = "henry/henryLargeTest/result_henry"
        input_file_name = "henry/henryLargeTest/henry"

        # Mesh built by GID for the postProcess
        gid_mode = GiDPostMode.GiD_PostAscii
        multifile = MultiFileFlag.MultipleFiles
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteElementsOnly
        self.gid_io = GidIO(output_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

        self.model_part_io_fluid = ModelPartIO(input_file_name)
        self.model_part_io_fluid.ReadModelPart(self.model_part)  

        self.model_part.ProcessInfo.SetValue(IS_FLOW_STATIONARY, 1);
        self.model_part.ProcessInfo.SetValue(IS_BUOYANCY, 0);
        self.model_part.ProcessInfo.SetValue(IS_TRANSPORT_STATIONARY, 0);
        self.model_part.ProcessInfo.SetValue(GRAVITY_X, 0.0);
        self.model_part.ProcessInfo.SetValue(GRAVITY_Y, -9.81);

        #Aarchivo = open("MatrixSystem.txt", "r") 
        #Acontenido = archivo.read()
        #Aprint (contenido)


        #the buffer size should be set up after the mesh is read for the first time (Transient problem =2,3. Steady problem =1) 
        self.model_part.SetBufferSize(2)

        # we add the DoFs  
        fractional_iterative_solver.AddDofs(self.model_part)

        print ("Time steps values on each time for unkown=2 (Buffer size)") 

        #creating flow solver (custom)
        self.solver = fractional_iterative_solver.FractionalIterativeSolver(self.model_part,self.domain_size)
        #we import the C++ file that includes the commands that we need
        #import fractional_iterative_strategy

        ## This part is contained is for the FractionalSolver not included at fractional_iterative_solver.py
        self.solver.time_order = 1
        self.solver.echo_level = 3

        print ("Iterative solver create succesfully") 

        self.solver.Initialize()

        self.Nsteps=40
        

        print ("Solver inicializate!")    

        #solver.calculateDensityNodes()
        #Arhivo de comparación:
        filePost="henry/henryLargeTest/SystemTest.txt"
        self.solver.ReadFile(filePost,self.Nsteps)

        mesh_name = 0.0
        self.gid_io.InitializeMesh( mesh_name );
        self.gid_io.WriteMesh((self.model_part).GetMesh());
        self.gid_io.FinalizeMesh()


        self.gid_io.InitializeResults(mesh_name,(self.model_part).GetMesh()) 

        #clean file output for matrix
        open('MatrixSystem_hard.txt', 'w').close()

        print ("Mesh read for the postProcess")


        #if step==0:
	    #for node in model_part.Nodes:
                #node.SetSolutionStepValue(PRESSURE,0,0.0)
                #node.SetSolutionStepValue(CONCENTRATION,0,0.0)
            #if(node.Y < 0.5):
                #node.SetSolutionStepValue(PRESSURE,0,9800.0)
                #node.SetSolutionStepValue(CONCENTRATION,0,0.0)
            #else:
                #node.SetSolutionStepValue(PRESSURE,0,0.0)
                #node.SetSolutionStepValue(CONCENTRATION,0,0.0)




        self.gid_io.WriteNodalResults(PRESSURE,self.model_part.Nodes,0, 0)
        self.gid_io.WriteNodalResults(CONCENTRATION,self.model_part.Nodes,0,0)
        self.gid_io.WriteNodalResults(DENSITY,self.model_part.Nodes,0,1000)
        self.gid_io.WriteNodalResults(DARCY_FLOW_BALANCE,self.model_part.Nodes,0,0)
        self.gid_io.WriteNodalResults(SINKSOURCE_BALANCE,self.model_part.Nodes,0,0)

    def ExecuteFinalizeSolutionStep(self):
        pass
        #for elem in self.model_part.Elements:
            #value = elem.GetSolutionStepValue(DARCY_FLOW_X,0)
            #self.assertAlmostEqual(value,DARCY_FLOW_X_TARGET(self.step, elem))



    def Solve(self):
        Dt=2
        out=0
        out_step=1
        t1 = timer.time()
        for step in range(1,self.Nsteps):
            out=out+1
            time = Dt*step
            print ("new step, time:",time)
    
            self.model_part.CloneTimeStep(time)
            self.solver.Solve()
            self.ExecuteFinalizeSolutionStep()

            if out==out_step:
                out=0
         #      print ("printing step:",step)
                self.gid_io.WriteNodalResults(PRESSURE,self.model_part.Nodes,time,0)
                self.gid_io.WriteNodalResults(CONCENTRATION,self.model_part.Nodes,time,0)
                self.gid_io.WriteNodalResults(DENSITY,self.model_part.Nodes,time,0)
                self.gid_io.WriteNodalResults(DARCY_FLOW_BALANCE,self.model_part.Nodes,time,0)
                self.gid_io.WriteNodalResults(SINKSOURCE_BALANCE,self.model_part.Nodes,time,0)
                self.gid_io.Flush()
    
            t2=timer.time()
            total_time=t2-t1

            print ("total_time", total_time)

        self.gid_io.FinalizeResults()
