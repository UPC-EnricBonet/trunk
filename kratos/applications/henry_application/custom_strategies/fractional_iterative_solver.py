from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.HenryApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(CONCENTRATION)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
    model_part.AddNodalSolutionStepVariable(CONCENTRATION_OLD_IT)
    model_part.AddNodalSolutionStepVariable(STORAGE_BALANCE)
    model_part.AddNodalSolutionStepVariable(DARCY_FLOW_BALANCE)
    model_part.AddNodalSolutionStepVariable(SINKSOURCE_BALANCE)

    print("variables for the fractional iterative solver added correctly")


def AddDofs(model_part):

    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(PRESSURE)
        node.AddDof(CONCENTRATION)

    print("dofs for the fractional iterative solver added correctly")

class FractionalIterativeSolver:

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        # assignation of parameters to be used
        self.pressureTol = 0.00001
        self.concentrationTol = 0.0000001
        self.maxIterPressure = 30
        self.maxIterConcentration = 30
        self.timeOrder = 1
        self.predictorCorrector = False

        self.CalculateReactions = False
        self.ReformDofAtEachIteration = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
        self.echoLevel = 1

        # definition of the solvers
        self.pressure_linear_solver = SkylineLUFactorizationSolver();
        self.concentration_linear_solver = SkylineLUFactorizationSolver();

        # others
        #pDiagPrecond = DiagonalPreconditioner() X
        #pILUPrecond = ILU0Preconditioner()
        #self.pressure_linear_solver = BICGSTABSolver(1e-3, 5000, pDiagPrecond) X
        #self.concentration_linear_solver = BICGSTABSolver(1e-3, 5000, pDiagPrecond) X
        # self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
        # self.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
    
    def Initialize(self):

        self.domain_size = int(self.domain_size)

        self.solverConfiguration = FractionalIterativeConfiguration(
            self.model_part,
            self.pressure_linear_solver,
            self.concentration_linear_solver,
            self.domain_size)

        self.domain_size = int(self.domain_size)

        self.pressureTol = float(self.pressureTol)
        self.concentrationTol = float(self.concentrationTol)
        self.maxIterPressure = int( self.maxIterPressure)
        self.maxIterConcentration = int(self.maxIterConcentration)
        self.timeOrder = int(self.timeOrder)
        self.predictorCorrector = bool(self.predictorCorrector)

        self.ReformDofAtEachIteration = bool(self.ReformDofAtEachIteration)
        self.MoveMeshFlag = bool(self.MoveMeshFlag)
        
        self.echoLevel = int(self.echoLevel)

        self.solver = FractionalIterativeStrategy(
            self.model_part,
            self.solverConfiguration,
            self.ReformDofAtEachIteration,
            self.pressureTol,
            self.concentrationTol,
            self.maxIterPressure,
            self.maxIterConcentration,
            self.timeOrder,
            self.domain_size,
            self.predictorCorrector,
            self.MoveMeshFlag,
            self.echoLevel)

        self.solver.Check()
        

        (self.solver).SetEchoLevel(self.echoLevel)
        print("finished initialization of the Fractional Iterative Strategy")

    def calculateDensityNodes(self):
        print ("calculateDensity function")
        (self.solver).calculateDensityNodes()

    def Solve(self):
        print("just before solve")
        print(self.model_part)
        (self.solver).Solve()

    def Clear(self):
        (self.solver).Clear()
