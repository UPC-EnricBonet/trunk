from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#

# setting the domain size for the problem to be solved
domain_size = 2

#
import sys
from KratosMultiphysics import *
from KratosMultiphysics.FlowTransportApplication import *
#from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *


# defining a model part for the fluid
iterative_fluid_model_part = ModelPart("FlowPart")
# defining a model part for the thermal problem
#<****><****><****><****>
iterative_transport_model_part = ModelPart("TransportPart"); 
#<****><****><****><****>

import iterative_fluid_solver
import iterative_transport_solver

#<****><****><****><****> 
#
#
# importing variables
iterative_fluid_solver.AddVariables(iterative_fluid_model_part)
iterative_transport_solver.AddVariables(iterative_transport_model_part)

iterative_fluid_model_part.AddNodalSolutionStepVariable(CONCENTRATION)
iterative_fluid_model_part.AddNodalSolutionStepVariable(POINT_HEAT_SOURCE)
iterative_fluid_model_part.AddNodalSolutionStepVariable(DENSITY)
iterative_fluid_model_part.AddNodalSolutionStepVariable(CONCENTRATION_OLD_IT)
iterative_fluid_model_part.AddNodalSolutionStepVariable(STORAGE_BALANCE)
iterative_fluid_model_part.AddNodalSolutionStepVariable(DARCY_FLOW_BALANCE)
iterative_fluid_model_part.AddNodalSolutionStepVariable(SINKSOURCE_BALANCE)

#solver_module_thermal.AddVariables(thermal_model_part, SolverSettings1)

#<****><****><****><****> 
#now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file  
gid_mode = GiDPostMode.GiD_PostAscii  #we import the python file that includes the commands that we need  
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("FluidTransport",gid_mode,multifile,deformed_mesh_flag,write_conditions)


model_part_io = ModelPartIO("FluidTransport")         # we set the name of the .mdpa file  
model_part_io.ReadModelPart(iterative_fluid_model_part)         # we load the info from the .mdpa 

iterative_fluid_model_part.ProcessInfo.SetValue(IS_FLOW_STATIONARY, 0);
iterative_fluid_model_part.ProcessInfo.SetValue(IS_BUOYANCY, 0);
iterative_fluid_model_part.ProcessInfo.SetValue(IS_TRANSPORT_STATIONARY, 0);
iterative_fluid_model_part.ProcessInfo.SetValue(GRAVITY_X, 0.0);
iterative_fluid_model_part.ProcessInfo.SetValue(GRAVITY_Y, -9.81);


iterative_transport_model_part.ProcessInfo.SetValue(GRAVITY_X, 0.0);
iterative_transport_model_part.ProcessInfo.SetValue(GRAVITY_Y, -9.81);


 # we create a mesh for the postprocess  
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((iterative_fluid_model_part).GetMesh());
gid_io.FinalizeMesh()

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)  
iterative_fluid_model_part.SetBufferSize(3)
iterative_transport_model_part.SetBufferSize(3)

 # we create a mesh for the postprocess  
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((iterative_fluid_model_part).GetMesh());
gid_io.FinalizeMesh()

iterative_fluid_solver.AddDofs(iterative_fluid_model_part)


"""thermal_model_part.Nodes       = fluid_model_part.Nodes         ##make the new model part to use the same nodes as the original model part!
thermal_model_part.ProcessInfo = fluid_model_part.ProcessInfo   ##this makes all of the time informations to coincide between the two model parts 
thermal_model_part.Properties  = fluid_model_part.Properties    ##assign all of the "Material properties"""

   



print ("about to solve!")  
#solver1.Solve()
 


#
#
# Creating the fluid solver

solver_of_fluid = iterative_fluid_solver.IterativeFluidSolver(iterative_fluid_model_part, domain_size)
solver_of_fluid.time_order = 1
solver_of_fluid.echo_level = 3
solver_of_fluid.Initialize()

print ("Solved Fluid!") 
 
model=ConnectivityPreserveModeler()
model.GenerateModelPart(iterative_fluid_model_part,iterative_transport_model_part, "Transport","MassFlow");
solver_of_transport = iterative_transport_solver.IterativeTransportSolver(iterative_transport_model_part,domain_size)
solver_of_transport.time_order = 1
solver_of_transport.echo_level = 3


iterative_transport_solver.AddDofs(iterative_transport_model_part)

solver_of_transport.Initialize()

print ("Solved Transport!") 

#here we set the boundary conditions for the diffusion solver
"""for node in thermal_model_part.Nodes:
   node.SetSolutionStepValue(TEMPERATURE,0,10.0)
   if (node.X<0.1):
      #node.Fix(TEMPERATURE)
      node.SetSolutionStepValue(TEMPERATURE,0,100.0)
   if (node.X>1.99999):
      node.SetSolutionStepValue(POINT_HEAT_SOURCE,0,5.0)
for node in fluid_model_part.Nodes:
   node.SetSolutionStepValue(PRESSURE,0,0.0)
   if (node.X<0.1):
      node.SetSolutionStepValue(PRESSURE,0,100.0)
   if (node.X>1.99999):
      node.SetSolutionStepValue(POINT_FLOW_SOURCE,0,5.0)"""


#### SOLO EN STEADY STATE

"""solver_of_fluid.Solve()
solver_of_thermal.Solve()
#print the values of node and temperature


#for node in fluid_model_part.Nodes:
   #node_pressure = node.GetSolutionStepValue(PRESSURE)
   #print (node.Id,  " ", node.X, " ",node.Y, " ",node.Z," ", PRESSURE," ",node_pressure)

#for node in thermal_model_part.Nodes:
   #node_temperature = node.GetSolutionStepValue(TEMPERATURE)
   #print (node.Id,  " ", node.X, " ",node.Y, " ",node.Z," ", TEMPERATURE," ",node_temperature)



#and we print the results  
gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh()) 
gid_io.WriteNodalResults(TEMPERATURE,fluid_model_part.Nodes,0,0)
gid_io.InitializeResults(mesh_name,(thermal_model_part).GetMesh()) 
gid_io.WriteNodalResults(PRESSURE,thermal_model_part.Nodes,0,0)
gid_io.FinalizeResults()

#since we have already calculated the temp, we can get the mean value
#first the constructor (it could have been called before)
calc_mean=CalculateMeanTemperature(thermal_model_part)
#and we calculate!
calc_mean.Execute()"""

#### SOLO EN STEADY STATE


#### SOLO EN UNSTEADY STATE



						
gid_io.InitializeResults(mesh_name,(iterative_fluid_model_part).GetMesh()) 

Nsteps=51
Dt=0.2
out=0
out_step=1
step = 0


gid_io.WriteNodalResults(PRESSURE,iterative_fluid_model_part.Nodes,0,0)
gid_io.WriteNodalResults(CONCENTRATION,iterative_transport_model_part.Nodes,0,0)
gid_io.WriteNodalResults(DENSITY,iterative_transport_model_part.Nodes,0,1000)
gid_io.WriteNodalResults(DARCY_FLOW_BALANCE,iterative_fluid_model_part.Nodes,0,0)

import time as timer
t1 = timer.time()

for step in range(0, Nsteps):
    time = Dt * step


    iterative_fluid_model_part.CloneTimeStep(time)	

    if step>0:
    	solver_of_fluid.Solve()
        solver_of_transport.Solve()	
        gid_io.WriteNodalResults(CONCENTRATION,iterative_transport_model_part.Nodes,time,0)
	gid_io.WriteNodalResults(PRESSURE,iterative_fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DENSITY,iterative_transport_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DARCY_FLOW_BALANCE,iterative_fluid_model_part.Nodes,time,0)
        gid_io.Flush()

    """if step==0:
	for node in fluid_model_part.Nodes:
                node.SetSolutionStepValue(TEMPERATURE,0,0.0) 
                node.SetSolutionStepValue(PRESSURE,0,0.0)"""
	
 
        

t2=timer.time()
total_time=t2-t1

print ("total_time")
print (total_time)

gid_io.FinalizeResults()
#### SOLO EN UNSTEADY STATE
