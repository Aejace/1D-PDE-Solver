import BoundaryConditions
from BoundaryConditions.FixedValueBoundaryConditions import LeftFixedValueBoundaryCondition, RightFixedValueBoundaryCondition


class BatemanBurgersPDE:
    # Local Variables
    leftBoundaryConditions = []
    rightBoundaryConditions = []

    # PDE variables
    deltaX = None
    alpha = None


    def __init__(self, deltaX):
        self.deltaX = deltaX
        self.alpha = float(input('Enter alpha value: '))  # Future work: Add error handling


        # Initialize boundary conditions
        self.leftBoundaryConditions.append(LeftFixedValueBoundaryCondition(self))
        self.rightBoundaryConditions.append(RightFixedValueBoundaryCondition(self))


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
        return (((state[i + 1]**2) - (state[i - 1]**2)) / (4.0 * self.deltaX)) + (simplifiedConstants * (state[i - 1] - 2.0 * state[i] + state[i + 1]))


    # Called to fill jacobian matrix with partial derivatives of ODE with respect or state[i - 1], state[i], and state[i + 1]
    def PartialDerivative(self, jacobianMatrix, state, i):
        simplifiedConstants = self.alphaOverDeltaXSquared()

        jacobianMatrix[i][i - 1] = (state[i - 1] / (2.0 * self.deltaX)) + simplifiedConstants
        jacobianMatrix[i][i] = -2.0 * simplifiedConstants
        jacobianMatrix[i][i + 1] = (state[i + 1]  / (-2.0 * self.deltaX)) + simplifiedConstants