
# setting the domain size for the problem to be solved
domain_size = 2

#
from KratosMultiphysics import *
from KratosMultiphysics.FlowTransportApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *



# defining a model part for the fluid
fluid_model_part = ModelPart("FluidPart")
# defining a model part for the thermal problem
#<**><**><**><**>
thermal_model_part = ModelPart("ThermalPart"); 
#<**><**><**><**>

import fluid_solver
import thermal_solver

#<**><**><**><**> 
#
#
# importing variables
fluid_solver.AddVariables(fluid_model_part)
thermal_solver.AddVariables(thermal_model_part)

fluid_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
fluid_model_part.AddNodalSolutionStepVariable(POINT_HEAT_SOURCE)

#solver_module_thermal.AddVariables(thermal_model_part, SolverSettings1)

#<**><**><**><**> 
#now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file  
gid_mode = GiDPostMode.GiD_PostAscii  #we import the python file that includes the commands that we need  
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly


gid_io = GidIO("FlowTransport",gid_mode,multifile,deformed_mesh_flag,write_conditions)


model_part_io = ModelPartIO("FlowTransport")         # we set the name of the .mdpa file 
model_part_io.ReadModelPart(fluid_model_part)         # we load the info from the .mdpa 

 # we create a mesh for the postprocess  
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((fluid_model_part).GetMesh());
gid_io.FinalizeMesh()

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)  
fluid_model_part.SetBufferSize(3)
thermal_model_part.SetBufferSize(3)

 # we create a mesh for the postprocess  
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((fluid_model_part).GetMesh());
gid_io.FinalizeMesh()

fluid_solver.AddDofs(fluid_model_part)


print ("about to solve!")  

# Creating the fluid solver

solver_of_fluid = fluid_solver.StaticPoissonSolver(fluid_model_part, domain_size)
solver_of_fluid.time_order = 1
solver_of_fluid.echo_level = 3
solver_of_fluid.Initialize()

print ("Solved Fluid!") 
 
model=ConnectivityPreserveModeler()
model.GenerateModelPart(fluid_model_part,thermal_model_part, "DiffConv2D","PointHeat");
solver_of_thermal = thermal_solver.StaticPoissonSolver(thermal_model_part,domain_size)
solver_of_thermal.time_order = 1
solver_of_thermal.echo_level = 3


thermal_solver.AddDofs(thermal_model_part)

solver_of_thermal.Initialize()

print ("Solved Transport!") 


#### SOLO EN UNSTEADY STATE

                                                
gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh()) 


time = 1
##output_step = 1
##Dt=0.1
##final_time = 5
gid_io.WriteNodalResults(TEMPERATURE,fluid_model_part.Nodes,0,0)
gid_io.WriteNodalResults(PRESSURE,thermal_model_part.Nodes,0,0)

while(time <= final_time):

    fluid_model_part.CloneTimeStep(time)
    if time>0:
        solver_of_fluid.Solve()
        solver_of_thermal.Solve()
        gid_io.WriteNodalResults(TEMPERATURE,thermal_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
        gid_io.Flush()
    time=time+1

gid_io.FinalizeResults()

f = open("FlowTransport"+".post.lst","w")
f.write("Single\n")
f.write("FlowTransport"+".post.res")
f.close()
#### SOLO EN UNSTEADY STATE
