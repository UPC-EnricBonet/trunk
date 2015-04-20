#setting the domain size for the problem to be solved
domain_size = 2  # 2D problem  

#including kratos path
import sys
from KratosMultiphysics import *    #we import the KRATOS  
from KratosMultiphysics.PorousMediaApplication import *        #and now our application

#Create a new empty model part called "ExampleModelPart"
model_part = ModelPart("ExampleModelPart");

print ("Model part defined: ExampleModelPart")  

#we import the python file that includes the commands that we need
import fractional_iterative_solver

#import variables that we will need from solver to our recent created model_part            
fractional_iterative_solver.AddVariables(model_part) 

# (note that our model part does not have nodes or elements yet) 

# introducing input & outoput (also postProcess) file name
output_file_name = "result_PressureHorizontal"
input_file_name = "pressureHorizontal"

# Mesh built by GID for the postProcess
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(output_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)  

model_part.ProcessInfo.SetValue(IS_FLOW_STATIONARY, 0);
model_part.ProcessInfo.SetValue(IS_BUOYANCY, 0);
model_part.ProcessInfo.SetValue(IS_TRANSPORT_STATIONARY, 0);
model_part.ProcessInfo.SetValue(GRAVITY_X, 0.0);
model_part.ProcessInfo.SetValue(GRAVITY_Y, 0.0);


#the buffer size should be set up after the mesh is read for the first time (Transient problem =2,3. Steady problem =1) 
model_part.SetBufferSize(2)

 # we add the DoFs  
fractional_iterative_solver.AddDofs(model_part)

print ("Time steps values on each time for unkown=2 (Buffer size)") 

#creating flow solver (custom)
solver = fractional_iterative_solver.FractionalIterativeSolver(model_part,domain_size)
## This part is contained is for the FractionalSolver not included at fractional_iterative_solver.py
solver.time_order = 1
solver.echo_level = 3

print ("Iterative solver create succesfully") 

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

Nsteps=51
Dt=0.2
out=0
out_step=1

gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,0, 0)
gid_io.WriteNodalResults(CONCENTRATION,model_part.Nodes,0,0)
gid_io.WriteNodalResults(DENSITY,model_part.Nodes,0,1000)
#gid_io.WriteNodalResults(DARCY_FLOW_BALANCE,model_part.Nodes,0,0)
#gid_io.WriteNodalResults(SINKSOURCE_BALANCE,model_part.Nodes,0,0)


import time as timer
t1 = timer.time()
for step in range(1,Nsteps):
    out=out+1
    time = Dt*step
    print "new step, time:",time
    
    model_part.CloneTimeStep(time)
    
    solver.Solve()

    if out==out_step:
   	  out=0
	  print "PRINTING STEP",step 	
   	  
          gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
          gid_io.WriteNodalResults(CONCENTRATION,model_part.Nodes,time,0)
          gid_io.WriteNodalResults(DENSITY,model_part.Nodes,time,0)
          #gid_io.WriteNodalResults(DARCY_FLOW_BALANCE,model_part.Nodes,time,0)
          #gid_io.WriteNodalResults(SINKSOURCE_BALANCE,model_part.Nodes,time,0)
	  gid_io.Flush()
    
t2=timer.time()
total_time=t2-t1

print "total_time", total_time

gid_io.FinalizeResults()
