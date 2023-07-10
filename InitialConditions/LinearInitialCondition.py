class LinearInitialCondition:
    # Local variables
    leftStartingY = None
    rightStartingY = None

    def __init__(self, differentialSystem):
        self.leftStartingY = float(input("Please enter the left hand side initial value: "))
        self.rightStartingY = float(input("Please enter the right hand side initial value: "))


    def GetValue(self, x):
        return ((self.rightStartingY - self.leftStartingY) * x) + self.leftStartingY