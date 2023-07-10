import BoundaryConditions
from BoundaryConditions.CoolingBoundaryCondition import LeftCoolingBoundaryCondition, RightCoolingBoundaryCondition
from BoundaryConditions.FixedValueBoundaryConditions import LeftFixedValueBoundaryCondition, RightFixedValueBoundaryCondition
from BoundaryConditions.HeatFluxBoundaryConditions import LeftHeatFluxBoundaryCondition, RightHeatFluxBoundaryCondition


class HeatPDE:
    # Local Variables
    leftBoundaryConditions = []
    rightBoundaryConditions = []

    # PDE variables
    deltaX = None
    alpha = None
    lateralCoefficientOfCooling = None
    lateralAmbientY = None


    def __init__(self, deltaX):
        self.deltaX = deltaX
        self.alpha = float(input('Enter alpha value: '))  # Future work: Add error handling
        self.lateralCoefficientOfCooling = float(input('Enter lateral coefficient of cooling value: '))  # Future work: Add error handling
        self.lateralAmbientY = float(input('Enter temperature an infinite distance from the sample space: '))  # Future work: Add error handling

        # Initialize boundary conditions
        self.leftBoundaryConditions.append(LeftFixedValueBoundaryCondition(self))
        self.leftBoundaryConditions.append(LeftHeatFluxBoundaryCondition(self))
        self.leftBoundaryConditions.append(LeftCoolingBoundaryCondition(self))

        self.rightBoundaryConditions.append(RightFixedValueBoundaryCondition(self))
        self.rightBoundaryConditions.append(RightHeatFluxBoundaryCondition(self))
        self.rightBoundaryConditions.append(RightCoolingBoundaryCondition(self))


    def GetLeftBoundaryConditions(self):
        return self.leftBoundaryConditions


    def GetRightBoundaryConditions(self, leftBoundaryIndex):
        # Ensures at least one end is a fixed value boundary
        if leftBoundaryIndex == 0:
            return self.rightBoundaryConditions
        else:
            onlyFixed = [self.rightBoundaryConditions[0]]
            return onlyFixed


    # Helper function
    def alphaOverDeltaXSquared(self):
        return self.alpha / (self.deltaX * self.deltaX)


    # Called to create system of ODE's out of the PDE.
    def ODE(self, state, i):
        simplifiedConstants = self.alphaOverDeltaXSquared()
        return simplifiedConstants * (state[i - 1] - 2.0 * state[i] + state[i + 1]) + self.lateralCoefficientOfCooling * (self.lateralAmbientY - state[i])



    # Called to fill jacobian matrix with partial derivatives of ODE with respect or state[i - 1], state[i], and state[i + 1]
    def PartialDerivative(self, jacobianMatrix, state, i):
        simplifiedConstants = self.alphaOverDeltaXSquared()

        jacobianMatrix[i][i - 1] = simplifiedConstants
        jacobianMatrix[i][i] = simplifiedConstants * (- 2.0) - self.lateralCoefficientOfCooling
        jacobianMatrix[i][i + 1] = simplifiedConstants