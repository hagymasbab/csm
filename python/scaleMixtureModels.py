
class scaleMixtureModel():

    def __init__(self, fname):
        self.fname = fname

    def generateSyntheticData(N):
        "e"

    def buildStanInput():
        "e"

    def compileStanModel():
        "e"

    def loadStanModel():
        "e"

    def runInference():
        "e"

    def plotResults():
        "e"


class GSM(scaleMixtureModel):
    """docstring for GSM"""

    def __init__(self):
        super(GSM, self).__init__("gsm")


class hierarchicalMixtureModel(scaleMixtureModel):
    """docstring for hierarchicalMixtureModel"""
    def __init__(self, arg):
        super(hierarchicalMixtureModel, self).__init__()
        self.arg = arg

    def buildPriorCovariance():
        "e"

    def buildPriorMean():
        "e"


class CSM(hierarchicalMixtureModel):
    """docstring for CSM"""

    def __init__(self):
        super(CSM, self).__init__("csm")


class CSM2(hierarchicalMixtureModel):
    """docstring for CSM2"""
    def __init__(self):
        super(CSM2, self).__init__("csm2")


class MSM(hierarchicalMixtureModel):
    """docstring for ClassName"""
    def __init__(self):
        super(MSM, self).__init__("msm")
