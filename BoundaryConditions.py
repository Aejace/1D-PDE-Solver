class LeftFixedValueBoundaryCondition:
    # Local Variables
    differentialSystem = None

    def __init__(self, differentialSystem):
        self.differentialSystem = differentialSystem


    def ODE(self, state):
        return 0


    def PartialDerivative(self, jacobianMatrix):
        n = len(jacobianMatrix[0])
        for i in range(0, n):
            # left boundary - fixed temperature
            jacobianMatrix[0][i] = 0.0
        return


class RightFixedValueBoundaryCondition:
    # Local Variables
    differentialSystem = None

    def __init__(self, differentialSystem):
        self.differentialSystem = differentialSystem

    def ODE(self, state):
        return 0

    def PartialDerivative(self, jacobianMatrix):
        n = len(jacobianMatrix[0])
        for i in range(0, n):
            # left boundary - fixed temperature
            jacobianMatrix[n - 1][i] = 0.0


class LeftHeatFluxBoundaryCondition:
    # Local Variables
    differentialSystem = None
    leftFlux = None
    simplifiedConstants = None

    def __init__(self, differentialSystem):
        self.differentialSystem = differentialSystem
        self.simplifiedConstants = differentialSystem.alphaOverDeltaXSquared()
        self.leftFlux = float(input('Enter left heat flux value: ')) # Future work: Add error handling


    def ODE(self, state):
        return 2.0 * self.simplifiedConstants * (state[1] - state[0]) - (2.0 / self.differentialSystem.deltaX) * self.leftFlux


    def PartialDerivative(self, jacobianMatrix):
        jacobianMatrix[0][0] = -2.0 * self.simplifiedConstants
        jacobianMatrix[0][1] = 2.0 * self.simplifiedConstants


class RightHeatFluxBoundaryCondition:
    # Local Variables
    differentialSystem = None
    rightFlux = None
    simplifiedConstants = None

    def __init__(self, differentialSystem):
        self.differentialSystem = differentialSystem
        self.simplifiedConstants = differentialSystem.alphaOverDeltaXSquared()
        self.rightFlux = float(input('Enter right heat flux value: ')) # Future work: Add error handling


    def ODE(self, state):
        n = len(state)
        return - 2.0 * self.simplifiedConstants * (state[n - 1] - state[n - 2]) - (2.0 / self.differentialSystem.deltaX) * self.rightFlux


    def PartialDerivative(self, jacobianMatrix):
        n = len(jacobianMatrix[0])
        jacobianMatrix[n - 1][n - 2] = 2.0 * self.simplifiedConstants
        jacobianMatrix[n - 1][n - 1] = -2.0 * self.simplifiedConstants


class LeftCoolingBoundaryCondition:
    # Local Variables
    differentialSystem = None
    simplifiedConstants = None
    leftCoefficientOfCooling = None
    leftAmbientY = None

    def __init__(self, differentialSystem):
        self.differentialSystem = differentialSystem
        self.simplifiedConstants = differentialSystem.alphaOverDeltaXSquared()
        self.leftCoefficientOfCooling = float(input('Enter left coefficient of cooling value: ')) # Future work: Add error handling
        self.leftAmbientY = float(input('Enter temperature an infinite distance from left side of the sample space: '))  # Future work: Add error handling


    def ODE(self, state):
        return 2.0 * self.simplifiedConstants * (state[1] - state[0]) - (2.0 / self.differentialSystem.deltaX) * self.leftCoefficientOfCooling * (state[0] - self.leftAmbientY)


    def PartialDerivative(self, jacobianMatrix):
        jacobianMatrix[0][0] = -2.0 * self.simplifiedConstants - (2.0 / self.differentialSystem.deltaX) * self.leftCoefficientOfCooling
        jacobianMatrix[0][1] = 2.0 * self.simplifiedConstants

class RightCoolingBoundaryCondition:
    # Local Variables
    differentialSystem = None
    simplifiedConstants = None
    rightCoefficientOfCooling = None
    rightAmbientY = None

    def __init__(self, differentialSystem):
        self.differentialSystem = differentialSystem
        self.simplifiedConstants = differentialSystem.alphaOverDeltaXSquared()
        self.rightCoefficientOfCooling = float(input('Enter left coefficient of cooling value: ')) # Future work: Add error handling
        self.rightAmbientY = float(input('Enter temperature an infinite distance from right side of the sample space: '))  # Future work: Add error handling


    def ODE(self, state):
        n = len(state)
        return - 2.0 * self.simplifiedConstants * (state[n - 1] - state[n - 2]) - (2.0 / self.differentialSystem.deltaX) * self.rightCoefficientOfCooling * (state[n - 1] - self.rightAmbientY)


    def PartialDerivative(self, jacobianMatrix):
        n = len(jacobianMatrix[0])
        jacobianMatrix[n - 1][n - 2] = 2.0 * self.simplifiedConstants
        jacobianMatrix[n - 1][n - 1] = -2.0 * self.simplifiedConstants - (2.0 / self.differentialSystem.deltaX) * self.rightCoefficientOfCooling

