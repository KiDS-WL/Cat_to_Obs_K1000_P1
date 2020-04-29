# 28/11/2019, B. M. Giblin
# Code to read in parameter file for shear ratio test
# and extract important info (MICE vs K1000, true/estimated P(z), mag on/off)

class Get_Input:
    def __init__(self, paramfile):
        self.paramfile = paramfile
        self.paraminput = open(self.paramfile).read()

    def Source_Type(self):
        # Check if sources are 'K1000' or 'MICE2_KV450'
        return self.paraminput.split('Source_Type = ')[-1].split(' ')[0].split('\n')[0].split('\t')[0]

    def Lens_Type(self):
        # Check if lenses are 'BOSS_data' or something else
        return self.paraminput.split('Lens_Type = ')[-1].split(' ')[0].split('\n')[0].split('\t')[0]

    def Random_Type(self):
        # Check if randoms are 'BOSS_random' or something else
        return self.paraminput.split('Random_Type = ')[-1].split(' ')[0].split('\n')[0].split('\t')[0]

    def Mag_OnOff(self):
        # Check if mag is 'on' or 'off' for MICE mocks
        return self.paraminput.split('Mag_OnOff = ')[-1].split(' ')[0].split('\n')[0].split('\t')[0]

    def Pz_TrueEstimated(self):
        # Check if P(z) is 'True' or 'Estimated' for MICE mocks
        return str(self.paraminput.split('Pz_TrueEstimated = ')[-1].split(' ')[0].split('\n')[0].split('\t')[0])

    def SN(self):
        # Check if there is shape noise included in the mock shears
        return str(self.paraminput.split('SN = ')[-1].split(' ')[0].split('\n')[0].split('\t')[0])
