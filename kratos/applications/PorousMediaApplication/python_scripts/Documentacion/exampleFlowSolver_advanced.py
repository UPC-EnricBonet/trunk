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


#Advanced Flow solver 
class AdvancedFlowSolver:

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        # assignation of parameters to be used for the strategy
        self.time_order = 1;
        self.prediction_order = 1; #or... predictor_corrector = False
        self.ReformDofAtEachIteration = False;
        self.echo_level = 0

        # definition of the specific solver to be used
        pILUPrecond = ILU0Preconditioner()
        self.linear_solver = BICGSTABSolver(1e-6, 5000, pILUPrecond)
        
        #Falta el timeScheme o se llama desde strategy...?

    def Initialize(self):

        #Esto no vale para nada...?
        self.model_part.ProcessInfo 

        #Call to constructor of our custom strategy in this case (ResidualBasedFlowStrategy)
        self.solver = ResidualBasedFlowStrategy(self.model_part, self.linear_solver, self.ReformDofAtEachIteration, self.time_order, self.prediction_order)
        
        #Esto que significa...?
        (self.solver).SetEchoLevel(self.echo_level)

        print("finished initialization of flow solver & residual based flow strategy")

    def Solve(self):
        #Esto que significa...? (neighbour algoritm)
        ##if(self.ReformDofAtEachIteration): 
        (self.solver).Solve()
