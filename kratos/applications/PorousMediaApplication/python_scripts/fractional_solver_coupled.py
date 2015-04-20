#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PorousMediaApplication import *

#Method to add all the variables to our problem: solver to mdpa object
def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(CONCENTRATION);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    model_part.AddNodalSolutionStepVariable(a);

    print ("PRESSURE &  CONCENTRATION variable added correctly")

#Method to add all the unkwons to our problem: solver to mdpa object
def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(CONCENTRATION);

    print ("degrees of freedom PRESSURE &  CONCENTRATION added correctly")

#Flow solver constructor, arguments: modelpart & domainSize
#Fractional solver 
class FractionalSolverCoupled:
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
        ReformDofSetAtEachStep = True	
        MoveMeshFlag = False
        
        import strategyFractionalCoupled
        self.solver = strategyFractionalCoupled.SolvingStrategyPython(self.model_part,self.time_scheme,self.basic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
	self.solver.SetEchoLevel(3)

    #######################################################################   
    def Solve(self):
        iter=0
        Converged =False
        maxIter = 0
        while(iter <= maxIter and Converged == False):
            self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1)
	    (self.solver).Solve()
            self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 2)
            (self.solver).Solve()
            Converged = (self.solver).CheckGlobalConverge(self.model_part)
            iter=iter+1

    ####################################################################### 

    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)


