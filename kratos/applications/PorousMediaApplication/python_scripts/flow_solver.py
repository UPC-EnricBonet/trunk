#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PorousMediaApplication import *

#Method to add all the variables to our problem: solver to mdpa object
def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(HEAD_LEVEL);

    print ("HEAD_LEVEL variable added correctly")

#Method to add all the unkwons to our problem: solver to mdpa object
def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(HEAD_LEVEL);

    print ("degrees of freedom (HEAD_LEVEL) added correctly")

#Flow solver constructor, arguments: modelpart & domainSize
#Basic Flow solver 
class BasicFlowSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):
        
	self.model_part = model_part
	self.domain_size = domain_size
        
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
	self.basic_linear_solver = SkylineLUFactorizationSolver()

        #definition of the convergence criteria
	#tolerance for the solver
        self.conv_criteria = DisplacementCriteria(1e-5,1e-15)
        
    #######################################################################
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False	
        MoveMeshFlag = False
        
        import strategyBasicFlow
        self.solver = strategyBasicFlow.SolvingStrategyPython(self.model_part,self.time_scheme,self.basic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
	self.solver.SetEchoLevel(3)

    #######################################################################   
    def Solve(self):
	  (self.solver).Solve()
    ####################################################################### 

    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)


