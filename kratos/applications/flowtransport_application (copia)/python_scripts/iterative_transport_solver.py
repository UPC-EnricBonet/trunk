#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FlowTransportApplication import *
CheckForPreviousImport()

def AddVariables(iterative_transport_model_part,config=None):  #this way er only need one command to add all the variables to our problem 
    iterative_transport_model_part.AddNodalSolutionStepVariable(CONCENTRATION);
    iterative_transport_model_part.AddNodalSolutionStepVariable(POINT_HEAT_SOURCE);
    iterative_transport_model_part.AddNodalSolutionStepVariable(DENSITY); 
    iterative_transport_model_part.AddNodalSolutionStepVariable(CONCENTRATION_OLD_IT);
    iterative_transport_model_part.AddNodalSolutionStepVariable(STORAGE_BALANCE); #?????????????????????
    iterative_transport_model_part.AddNodalSolutionStepVariable(DARCY_FLOW_BALANCE);#?????????????????????
    iterative_transport_model_part.AddNodalSolutionStepVariable(SINKSOURCE_BALANCE);#?????????????????????

def AddDofs(iterative_transport_model_part):
    for node in iterative_transport_model_part.Nodes:

        #adding dofs
        node.AddDof(CONCENTRATION);

    print ("variables for the Transport solver added correctly")

class IterativeTransportSolver:
    #######################################################################
    def __init__(self,iterative_transport_model_part,domain_size):  #constructor of the class 

        self.iterative_transport_model_part = iterative_transport_model_part
        self.domain_size = domain_size

        # assignation of parameters to be used
        self.concentrationTol = 0.0000001
        self.maxIterConcentration = 30
        self.timeOrder = 1
        self.predictorCorrector = False

        self.CalculateReactions = False
        self.ReformDofAtEachIteration = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
        self.echoLevel = 1


        self.concentration_linear_solver = SkylineLUFactorizationSolver();



        
    #######################################################################
    def Initialize(self):
       
        self.domain_size = int(self.domain_size)
        self.solverTransportConfiguration = IterativeTransportConfiguration(
            self.iterative_transport_model_part,
            self.concentration_linear_solver,
            self.domain_size)

        self.domain_size = int(self.domain_size)

        self.concentrationTol = float(self.concentrationTol)
        self.maxIterConcentration = int( self.maxIterConcentration)
        self.timeOrder = int(self.timeOrder)
        self.predictorCorrector = bool(self.predictorCorrector)

        self.ReformDofAtEachIteration = bool(self.ReformDofAtEachIteration)
        self.MoveMeshFlag = bool(self.MoveMeshFlag)
        
        self.echoLevel = int(self.echoLevel)

        self.solver = IterativeTransportStrategy(
            self.iterative_transport_model_part,
            self.solverTransportConfiguration,
            self.ReformDofAtEachIteration,
            self.concentrationTol,
            self.maxIterConcentration,
            self.timeOrder,
            self.domain_size,
            self.predictorCorrector,
            self.MoveMeshFlag,
            self.echoLevel)

        self.solver.Check()

        (self.solver).SetEchoLevel(self.echoLevel)
        print("finished initialization of the Fractional Iterative Strategy")


                 
    #######################################################################   
    def Solve(self):
        print("just before solve")
        print(self.iterative_transport_model_part)
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
