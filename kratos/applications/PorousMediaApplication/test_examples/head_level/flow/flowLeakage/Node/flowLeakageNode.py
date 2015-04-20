#setting the domain size for the problem to be solved
domain_size = 2  # 2D problem  

#including kratos path
import sys
from KratosMultiphysics import *    #we import the KRATOS  
from KratosMultiphysics.FlowApplication import *        #and now our application

#Create a new empty model part called "ExampleModelPart"
model_part = ModelPart("ExampleModelPart");

print ("Model part defined: ExampleModelPart")  

#we import the python file that includes the commands that we need
import flow_solver

#import variables that we will need from solver to our recent created model_part            
flow_solver.AddVariables(model_part) 

# (note that our model part does not have nodes or elements yet) 

# introducing input & outoput (also postProcess) file name
output_file_name = "result_flowLeakageNode"
input_file_name = "flowLeakageNode"

# Mesh built by GID for the postProcess
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(output_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)  

#the buffer size should be set up after the mesh is read for the first time (Transient problem =2,3. Steady problem =1) 
model_part.SetBufferSize(2)

 # we add the DoFs  
flow_solver.AddDofs(model_part)

print ("Time steps values on each time for unkown=2 (Buffer size)") 

#creating flow solver (custom)
solver = flow_solver.BasicFlowSolver(model_part,domain_size)
## This part is contained is for the BasicFlowSolver not included at flowSolver.py
solver.time_order = 1
solver.echo_level = 3

print ("Flow solver create succesfully") 

solver.Initialize()

print ("Solver inicializate!")    

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()


gid_io.InitializeResults(mesh_name,(model_part).GetMesh()) 

#clean file output for matrix
open('MatrixSystem.txt', 'w').close()
open('LastMatrixSystem.txt', 'w').close()

print ("Mesh read for the postProcess")

Nsteps=50
Dt=0.1
out=0
out_step=1

gid_io.WriteNodalResults(HEAD_LEVEL,model_part.Nodes,0,0)

import time as timer
t1 = timer.time()
for step in range(1,Nsteps):
    out=out+1
    time = Dt*step
    print "new step, time:",time
    
    model_part.CloneTimeStep(time)

    if step>1:
    	solver.Solve()

    if step==0:
	for node in model_part.Nodes:
                node.SetSolutionStepValue(HEAD_LEVEL,0,0.0) 
    if out==out_step:
   	  out=0
	  print "PRINTING STEP",step 	
   	  
          gid_io.WriteNodalResults(HEAD_LEVEL,model_part.Nodes,time,0)
	  gid_io.Flush()
    
t2=timer.time()
total_time=t2-t1

print "total_time", total_time

gid_io.FinalizeResults()

