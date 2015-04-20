##################################################################
##################################################################
#setting the domain size for the problem to be solved


domain_size = 2  # 2D problem  

#including kratos path
from KratosMultiphysics import *    #we import the KRATOS  
from KratosMultiphysics.FlowTransportApplication import *


#defining a model part
thermal_model_part = ModelPart("Example");  #we create a model part  

import thermal_solver
thermal_solver.AddVariables(thermal_model_part)  #from the static_poisson_solver.py we call the function Addvariables so that the model part we have just created has the needed variables 

 # (note that our model part does not have nodes or elements yet) 

 #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file  
gid_mode = GiDPostMode.GiD_PostAscii  #we import the python file that includes the commands that we need  
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("art4",gid_mode,multifile,deformed_mesh_flag,write_conditions)


model_part_io = ModelPartIO("example")             # we set the name of the .mdpa file  
model_part_io.ReadModelPart(thermal_model_part)         # we load the info from the .mdpa 

 # we create a mesh for the postprocess  
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((thermal_model_part).GetMesh());
gid_io.FinalizeMesh()

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)  
thermal_model_part.SetBufferSize(2)

 # we add the DoFs  
thermal_solver.AddDofs(thermal_model_part)
   
#creating a solver object
solver = thermal_solver.StaticPoissonSolver(thermal_model_part,domain_size)
solver.time_order = 1
solver.echo_level = 0
solver.Initialize()


print ("about to solve example flow_transport!")  
solver.Solve()
print ("Solved!")  

#print the values of node and temperature
for node in thermal_model_part.Nodes:
   node_temperature = node.GetSolutionStepValue(TEMPERATURE) 
   print node.Id,  " ", node.X, " ",node.Y, " ",node.Z," ", TEMPERATURE," ",node_temperature

for element in thermal_model_part.Elements:
      print Element," ",element.Id

#and we print the results  
gid_io.InitializeResults(mesh_name,(thermal_model_part).GetMesh()) 
gid_io.WriteNodalResults(TEMPERATURE,thermal_model_part.Nodes,0,0)
gid_io.FinalizeResults()

#since we have already calculated the temp, we can get the mean value
#first the constructor (it could have been called before)
calc_mean=CalculateMeanTemperature(thermal_model_part)
#and we calculate!
calc_mean.Execute()
