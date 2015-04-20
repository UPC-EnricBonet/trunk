#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FlowTransportApplication import *

CheckForPreviousImport()

def AddVariables(iterative_fluid_model_part):  #this way er only need one command to add all the variables to our problem 
    iterative_fluid_model_part.AddNodalSolutionStepVariable(PRESSURE);
    iterative_fluid_model_part.AddNodalSolutionStepVariable(POINT_FLOW_SOURCE);
    iterative_fluid_model_part.AddNodalSolutionStepVariable(DARCY_FLOW_X);
    iterative_fluid_model_part.AddNodalSolutionStepVariable(DARCY_FLOW_Y);
    iterative_fluid_model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
    iterative_fluid_model_part.AddNodalSolutionStepVariable(STORAGE_BALANCE)
    iterative_fluid_model_part.AddNodalSolutionStepVariable(DARCY_FLOW_BALANCE)
    iterative_fluid_model_part.AddNodalSolutionStepVariable(SINKSOURCE_BALANCE)

def AddDofs(iterative_fluid_model_part):
    for node in iterative_fluid_model_part.Nodes:
        node.AddDof(PRESSURE);

    print ("variables for the iterative_fluid_solver added correctly")

class IterativeFluidSolver:
    #######################################################################
    def __init__(self,iterative_fluid_model_part,domain_size):  #constructor of the class 


        self.domain_size = domain_size
        self.iterative_fluid_model_part = iterative_fluid_model_part

        # assignation of parameters to be used
        self.pressureTol = 0.00001
        self.maxIterPressure = 30
        self.timeOrder = 1
        self.predictorCorrector = False

        self.CalculateReactions = False
        self.ReformDofAtEachIteration = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
        self.echoLevel = 1

        # definition of the solver
        self.pressure_linear_solver = SkylineLUFactorizationSolver()


        
    #######################################################################
    def Initialize(self):
        self.domain_size = int(self.domain_size)

        self.solverFluidConfiguration = IterativeFluidConfiguration(
            self.iterative_fluid_model_part,
            self.pressure_linear_solver,
            self.domain_size)

        self.domain_size = int(self.domain_size)

        self.pressureTol = float(self.pressureTol)
        self.maxIterPressure = int( self.maxIterPressure)
        self.timeOrder = int(self.timeOrder)
        self.predictorCorrector = bool(self.predictorCorrector)

        self.ReformDofAtEachIteration = bool(self.ReformDofAtEachIteration)
        self.MoveMeshFlag = bool(self.MoveMeshFlag)
        
        self.echoLevel = int(self.echoLevel)

        self.solver = IterativeFluidStrategy(
            self.iterative_fluid_model_part,
            self.solverFluidConfiguration,
            self.ReformDofAtEachIteration,
            self.pressureTol,
            self.maxIterPressure,
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
        print(self.iterative_fluid_model_part)
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)


