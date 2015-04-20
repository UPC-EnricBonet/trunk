##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################

#including kratos path
from KratosMultiphysics import *    #we import the KRATOS  
from KratosMultiphysics.PureDiffusionApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem 


#defining a model part
model_part = ModelPart("PureDiffusionPart");  


#adding of Variables to Model Part
import pure_diffusion_solver
pure_diffusion_solver.AddVariables(model_part)


#now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file  
gid_mode = GiDPostMode.GiD_PostAscii  #we import the python file that includes the commands that we need  
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("PureGid",gid_mode,multifile,deformed_mesh_flag,write_conditions) #output will be the name of the output files

model_part_io = ModelPartIO("PureGid")             # we set the name of the .mdpa file  
model_part_io.ReadModelPart(model_part)                                                 # we load the info from the .mdpa 



# we create a mesh (to be used in the postprocess)
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print model_part

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)  
model_part.SetBufferSize(1)


# we add the DoFs  
pure_diffusion_solver.AddDofs(model_part)
   
#creating a solver object
solver = pure_diffusion_solver.StaticPoissonSolver(model_part,domain_size)
solver.time_order = 1
solver.echo_level = 0
solver.Initialize()


print "about to solve!"  
solver.Solve()
print "Solved!"  


#and we print the results  
gid_io.InitializeResults(mesh_name,(model_part).GetMesh()) 
gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,0,0)
gid_io.FinalizeResults()

f = open("PureGid"+".post.lst","w")
f.write("Single\n")
f.write("PureGid"+"_0.post.res")
f.close()
