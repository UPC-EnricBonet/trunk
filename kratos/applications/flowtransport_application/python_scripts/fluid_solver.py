#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FlowTransportApplication import *

def AddVariables(fluid_model_part):  #this way er only need one command to add all the variables to our problem 
    fluid_model_part.AddNodalSolutionStepVariable(PRESSURE);
    fluid_model_part.AddNodalSolutionStepVariable(POINT_FLOW_SOURCE);

def AddDofs(fluid_model_part):
    for node in fluid_model_part.Nodes:
        node.AddDof(PRESSURE);

    print ("variables for the Poisson solver added correctly")

class StaticPoissonSolver:
    #######################################################################
    def __init__(self,fluid_model_part,domain_size):  #constructor of the class 

        self.fluid_model_part = fluid_model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
        self.poisson_linear_solver =  SkylineLUFactorizationSolver()  #we set the type of solver that we want 
        
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-6,1e-9)  #tolerance for the solver 
        
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = False
        import strategy_python
        self.solver = strategy_python.SolvingStrategyPython(self.fluid_model_part,self.time_scheme,self.poisson_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
      
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)


