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
import flow_solver_advanced

#import variables that we will need from solver to our recent created model_part            
flow_solver_advanced.AddVariables(model_part) 

# (note that our model part does not have nodes or elements yet) 

# introducing input & outoput (also postProcess) file name
input_file_name = "advancedFlow"

# Mesh built by GID for the postProcess
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)
 
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

print ("Mesh read for the postProcess")  

#the buffer size should be set up after the mesh is read for the first time (Transient problem =2,3. Steady problem =1) 
model_part.SetBufferSize(2)

print ("Time steps values on each time for unkown=2 (Buffer size)") 

#creating flow solver (custom)
solver = flow_solver_advanced.AdvancedFlowSolver(model_part,domain_size)
## This part is contained at the class AdvancedFlowSolver at flowSolver.py
#solver.time_order = 1
#solver.echo_level = 3 ##?
#predictor_corrector = False
#ReformDofAtEachIteration = False

print ("Flow solver create succesfully") 

solver.Initialize()

print ("Solver inicializate!")    

# assigning the fluid properties
#permeability = 1.0
#specificStorage = 0.8
headLevel = 0.0
for node in model_part.Nodes:
    #node.SetSolutionStepValue(SPECIFIC_STORAGE, 0, specificStorage); 
    #node.SetSolutionStepValue(PERMEABILITY_WATER, 0, permeability);
    node.SetSolutionStepValue(HEAD_LEVEL, 0, headLevel);

print ("Initials conditions added!") 

#???
#model_part.Properties[0][SPECIFIC_STORAGE] = 0.0   
#model_part.Properties[0][PERMEABILITY_WATER] = 0.0

##no necessary now
# applying a temperature of 100 (Boundary conditions ????)
#for node in model_part.Nodes:
    #if(node.Y > 0.499):
        #node.SetSolutionStepValue(FACE_HEAT_FLUX, 0, 1000.0);

# settings to be changed
Dt = 0.2
Nsteps = 300
time = 0.0
out = 0
step = 0
output_step = 1
Nsteps = 10
for step in range(0, Nsteps):
    time = Dt * step
    model_part.CloneTimeStep(time)

    # solving the fluid problem
    if(step > 1):
        solver.Solve()

    # print the results
    if(out == output_step):
        gid_io.WriteNodalResults(HEAD_LEVEL, model_part.Nodes, time, 0)
        out = 0
    out = out + 1



