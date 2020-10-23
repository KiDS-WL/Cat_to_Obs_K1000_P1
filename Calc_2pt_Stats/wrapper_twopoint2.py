

    ##############################
    ##  wrapper_twopoint2.py    ##
    ##  Chieh-An Lin            ##
    ##  Version 2020.02.19      ##
    ##############################


import numpy as np
import scipy.interpolate as itp
import astropy.io.fits as fits
import twopoint
import wrapper_twopoint as wtp


###############################################################################
## Functions related to mean & covariance

class TwoPointBuilder:
    
    def __init__(self,
                  nbTomoN=2, nbTomoG=5, 
                  N_theta=9, theta_min=0.5, theta_max=300,
                  N_ell=8, ell_min=100, ell_max=1500,
                  nbModes=5,
                  prefix_Flinc=None,
                  prefix_CosmoSIS=None,
                  verbose=True):
    
        ## File prefixes
        self.prefix_Flinc    = 'data/mockFootprint/' if prefix_Flinc is None else prefix_Flinc
        self.prefix_CosmoSIS = 'data/mockFootprint/milvus/cosmosis/' if prefix_CosmoSIS is None else prefix_CosmoSIS
        
        ## Customize the above for your own inference; but don't touch the below
        ########################################################################
        
        ## Tomographic bins
        self.nbTomoN = nbTomoN
        self.nbTomoG = nbTomoG
        
        ## Define angular bin parameters
        self.N_theta   = N_theta
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.N_ell     = N_ell
        self.ell_min   = ell_min
        self.ell_max   = ell_max
        self.nbModes   = nbModes
        
        self.verbose    = verbose
        self.mean       = None
        self.cov        = None
        self.kernelList = None
        
        self.setAngBins()
        self.setNbPairs()
        return
    
    ## Initialize
    def setAngBins(self):
        bAndC_theta = np.logspace(np.log10(self.theta_min), np.log10(self.theta_max), 2*self.N_theta+1)
        self.theta     = bAndC_theta[1::2] ## [arcmin]
        self.bin_theta = bAndC_theta[0::2]
        bAndC_ell = np.logspace(np.log10(self.ell_min), np.log10(self.ell_max), 2*self.N_ell+1)
        self.ell       = bAndC_ell[1::2]
        self.bin_ell   = bAndC_ell[0::2]
        self.nArr      = np.arange(self.nbModes) + 1.0
        return
    
    def setNbPairs(self):
        self.__nbPairsNN  = self.nbTomoN * (self.nbTomoN+1) // 2
        self.__nbPairsNG  = self.nbTomoN * self.nbTomoG
        self.__nbPairsGG  = self.nbTomoG * (self.nbTomoG+1) // 2
        self.__pairListNN = [(i, j) for i in range(self.nbTomoN) for j in range(i, self.nbTomoN)]
        self.__pairListNG = [(i, j) for i in range(self.nbTomoN) for j in range(self.nbTomoG)]
        self.__pairListGG = [(i, j) for i in range(self.nbTomoG) for j in range(i, self.nbTomoG)]
        return
    
    ## Mean
    def makeTomoAngDict(self):
        labConv = wtp.LabelConvention()
        tomoAngDict = { ## Don't touch the order of this list
            labConv.w:       [self.__pairListNN, self.N_theta, self.theta],
            labConv.gamma_t: [self.__pairListNG, self.N_theta, self.theta], 
            labConv.gamma_x: [self.__pairListNG, self.N_theta, self.theta], 
            labConv.xi_p:    [self.__pairListGG, self.N_theta, self.theta], 
            labConv.xi_m:    [self.__pairListGG, self.N_theta, self.theta], 
            labConv.P_nn:    [self.__pairListNN, self.N_ell,   self.ell], 
            labConv.P_ne_E:  [self.__pairListNG, self.N_ell,   self.ell], 
            labConv.P_ne_B:  [self.__pairListNG, self.N_ell,   self.ell], 
            labConv.P_ee_E:  [self.__pairListGG, self.N_ell,   self.ell], 
            labConv.P_ee_B:  [self.__pairListGG, self.N_ell,   self.ell], 
            labConv.E_n:     [self.__pairListGG, self.nbModes, self.nArr], 
            labConv.B_n:     [self.__pairListGG, self.nbModes, self.nArr]
        }
        
        assert len(tomoAngDict) == len(labConv.kernelTypeDict)
        return tomoAngDict
    
    def _makeMean_none(self, statsTag):
        labConv = wtp.LabelConvention()
        tomoAngDict = self.makeTomoAngDict()
        statsList = statsTag.split('+')
        statsList_complete = tomoAngDict.keys()
        
        for stats in statsList:
            if stats not in statsList_complete:
                raise ValueError('\"%s\" not allowed' % statsTag)
        
        stock = []
        
        for stats, line in tomoAngDict.items():
            pairList = line[0]
            N_ang    = line[1]
          
            if stats in statsList:
                stock += [0.0] * len(pairList) * N_ang
        return stock
    
    def _loadNpyDataMat_Flinc(self, aves, statsTag, randTag, verbose=True):
        labConv = wtp.LabelConvention()
        statsList = statsTag.split('+')
        stock = []
        
        for stats in statsList:
            stats_c = labConv.defaultToCustomStatsTag(stats)
            name = '%s%s/MFP_combDataMat/dataMat_%s_%s_full.npy' % (self.prefix_Flinc, aves, stats_c, randTag)
            data = np.load(name)
            stock.append(data)
            if verbose == True:
                print('Loaded \"%s\"' % name)
        
        ## Cut
        NList = [data.shape[0] for data in stock]
        N_min = min(NList)
        stock = [data[:N_min] for data in stock]
        stock = np.hstack(stock)
        
        if verbose == True:
            print('N = %s' % stock.shape[0])
        return stock
    
    def _makeMean_Flinc(self, aves, statsTag, verbose=True):
        data = self._loadNpyDataMat_Flinc(aves, statsTag, 'signal', verbose=verbose)
        mean = np.mean(data, axis=0)
        return mean
    
    def _interpolateCF(self, x, y):
        theta = self.theta
        inter = itp.interp1d(np.log10(x), x*y, bounds_error=False, fill_value='extrapolate')
        CF    = inter(np.log10(theta)) / theta
        return theta, CF

    def _loadAsciiMean_CosmoSIS(self, stats, i, j, xArr, verbose=True):
        labConv = wtp.LabelConvention()
        
        if stats == labConv.w:
            name = '%sgalaxy_xi/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            yArr = loadAscii(name, verbose=verbose)
            theta, w = self._interpolateCF(xArr, yArr)
            return w
        
        if stats == labConv.gamma_t:
            name = '%sgalaxy_shear_xi/bin_%d_%d.txt' % (self.prefix_CosmoSIS, i+1, j+1)
            yArr = loadAscii(name, verbose=verbose)
            theta, gamma_t = self._interpolateCF(xArr, yArr)
            return gamma_t
          
        if stats == labConv.gamma_x:
            gamma_x = [0] * self.N_theta
            return gamma_x
        
        if stats == labConv.xi_p:
            name = '%sxi_binned_plus/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            xi_p = loadAscii(name, verbose=verbose)
            return xi_p
          
        if stats == labConv.xi_m:
            name = '%sxi_binned_minus/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            xi_m = loadAscii(name, verbose=verbose)
            return xi_m
          
        if stats == labConv.P_nn:
            name = '%sbandpower_galaxy/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            P_nn = loadAscii(name, verbose=verbose)
            return P_nn
        
        if stats == labConv.P_ne_E:
            name = '%sbandpower_galaxy_shear/bin_%d_%d.txt' % (self.prefix_CosmoSIS, i+1, j+1)
            P_ne_E = loadAscii(name, verbose=verbose)
            return P_ne_E
        
        if stats == labConv.P_ne_B:
            P_ne_B = [0] * self.N_ell
            return P_ne_B
        
        if stats == labConv.P_ee_E:
            name   = '%sbandpower_shear_e/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            P_ee_E = loadAscii(name, verbose=verbose)
            return P_ee_E
          
        if stats == labConv.P_ee_B:
            P_ee_B = [0] * self.N_ell
            return P_ee_B
          
        if stats == labConv.E_n:
            name = '%scosebis/bin_%d_%d.txt' % (self.prefix_CosmoSIS, j+1, i+1)
            E_n  = loadAscii(name, verbose=verbose)
            return E_n
        
        if stats == labConv.B_n:
            B_n = [0] * self.nbModes
            return B_n
        
        return None
    
    def _makeMean_CosmoSIS(self, statsTag, verbose=True):
        labConv = wtp.LabelConvention()
        tomoAngDict = self.makeTomoAngDict()
        statsList = statsTag.split('+')
        statsList_complete = tomoAngDict.keys()
        
        for stats in statsList:
            if stats not in statsList_complete:
                raise ValueError('\"%s\" not allowed' % statsTag)
        
        name  = '%sshear_xi_plus/theta.txt' % self.prefix_CosmoSIS
        xArr  = loadAscii(name, verbose=verbose) * (60.0 * 180.0 / np.pi) ## [arcmin]
        stock = []
        
        for stats, line in tomoAngDict.items():
            pairList = line[0]
          
            if stats in statsList:
                for i, j in pairList:
                    value = self._loadAsciiMean_CosmoSIS(stats, i, j, xArr, verbose=verbose)
                    stock.append(value)
        
        stock = np.concatenate(stock)
        return stock
    
    def setMean(self, tag, name=None, statsTag=None, verbose=True):
        if tag is None or tag == 'none':
            if statsTag == None:
                print('No mean')
                return None
        
            self.mean = self._makeMean_none(statsTag)
            print('Set mean to dummy values')
            return
        elif tag == 'Flinc':
            self.mean = self._makeMean_Flinc(name, statsTag, verbose=verbose) ## Consider `name` as `aves`
        elif tag == 'CosmoSIS':
            self.mean = self._makeMean_CosmoSIS(statsTag, verbose=verbose)
        elif tag == 'variable':
            self.mean = name ## Consider `name` as the variable which contains the mean
        
        ## Otherwise, consider `name` as the path to the file
        ## which contains the mean data vector
        else:
            try:
                if name[-4:] == '.npy':
                    self.mean = np.load(name).flatten()
                elif name[-4:] == '.fit' or name[-5:] == '.fits':
                    self.mean = fits.getdata(name, 1).field(0)
                else:
                    self.mean = np.loadtxt(name).flatten()
                if verbose:
                    print('Loaded \"%s\"' % name)
            except:
                raise OSError('\"%s\" not found' % name)
        return
    
    ## Covariance
    def _makeCov_Flinc(self, aves, statsTag, verbose=True):
        data = self._loadNpyDataMat_Flinc(aves, statsTag, 'obs', verbose=verbose)
        cov  = np.cov(data, rowvar=0, ddof=1)
        return cov
    
    def _makePairIndex_2PCF(self, tomo1, tomo2, sign, order=-1):
        if order > 0:
            beginGG = 0
            beginNG = self.__nbPairsGG * 2
            beginNN = self.__nbPairsGG * 2 + self.__nbPairsNG
            split   = self.nbTomoG
            tomo1G  = tomo1 - 0
            tomo1N  = tomo1 - split
            tomo2G  = tomo2 - 0
            tomo2N  = tomo2 - split
            
            isGG    = asInt(tomo2 < split)
            isNG    = asInt((tomo1 < split) * (tomo2 >= split))
            isNN    = asInt(tomo1 >= split)
            ind     = isGG * (beginGG + pairIndex(self.nbTomoG, tomo1G, tomo2G) + self.__nbPairsGG * sign)
            ind    += isNG * (beginNG + tomo2N + self.nbTomoN * tomo1G) ## Vary lens tomo bins first
            ind    += isNN * (beginNN + pairIndex(self.nbTomoN, tomo1N, tomo2N))
          
        else:
            beginNN = 0
            beginNG = self.__nbPairsNN
            beginGG = self.__nbPairsNN + self.__nbPairsNG
            split   = self.nbTomoN
            tomo1N  = tomo1 - 0
            tomo1G  = tomo1 - split
            tomo2N  = tomo2 - 0
            tomo2G  = tomo2 - split
            
            isNN    = asInt(tomo2 < split)
            isNG    = asInt((tomo1 < split) * (tomo2 >= split))
            isGG    = asInt(tomo1 >= split)
          
        ind  = isNN * (beginNN + pairIndex(self.nbTomoN, tomo1N, tomo2N))
        ind += isNG * (beginNG + tomo2G + self.nbTomoG * tomo1N) ## Vary source tomo bins first
        ind += isGG * (beginGG + pairIndex(self.nbTomoG, tomo1G, tomo2G) + self.__nbPairsGG * sign)
        return ind
    
    def _makePairIndex_BP(self, tomo1, tomo2, order=-1):
        if order > 0:
            beginGG = 0
            beginNG = self.__nbPairsGG
            beginNN = self.__nbPairsGG + self.__nbPairsNG
            split   = self.nbTomoG
            tomo1G  = tomo1 - 0
            tomo1N  = tomo1 - split
            tomo2G  = tomo2 - 0
            tomo2N  = tomo2 - split
            
            isGG    = asInt(tomo2 < split)
            isNG    = asInt((tomo1 < split) * (tomo2 >= split))
            isNN    = asInt(tomo1 >= split)
            ind     = isGG * (beginGG + pairIndex(self.nbTomoG, tomo1G, tomo2G))
            ind    += isNG * (beginNG + tomo2N + self.nbTomoN * tomo1G) ## Vary lens tomo bins first
            ind    += isNN * (beginNN + pairIndex(self.nbTomoN, tomo1N, tomo2N))
          
        else:
            beginNN = 0
            beginNG = self.__nbPairsNN
            beginGG = self.__nbPairsNN + self.__nbPairsNG
            split   = self.nbTomoN
            tomo1N  = tomo1 - 0
            tomo1G  = tomo1 - split
            tomo2N  = tomo2 - 0
            tomo2G  = tomo2 - split
            
            isNN    = asInt(tomo2 < split)
            isNG    = asInt((tomo1 < split) * (tomo2 >= split))
            isGG    = asInt(tomo1 >= split)
            
        ind  = isNN * (beginNN + pairIndex(self.nbTomoN, tomo1N, tomo2N))
        ind += isNG * (beginNG + tomo2G + self.nbTomoG * tomo1N) ## Vary source tomo bins first
        ind += isGG * (beginGG + pairIndex(self.nbTomoG, tomo1G, tomo2G))
        return ind
    
    def _covListToMatrix_2PCF(self, data, cleanNaN=True, CTag='tot'):
        if cleanNaN:
            ind = np.isnan(data)
            data[ind] = 0.0
        
        tomoA1 = asInt(data[0]) - 1
        tomoA2 = asInt(data[1]) - 1
        tomoB1 = asInt(data[6]) - 1
        tomoB2 = asInt(data[7]) - 1
        signA  = asInt(data[4])
        signB  = asInt(data[10])
        binA   = asInt(data[5])
        binB   = asInt(data[11])
        if CTag == 'tot':
            value = data[12:].sum(axis=0)
        elif CTag == 'CNG2':
            value = data[12] + data[13] + data[14]
        elif CTag == 'CNG':
            value = data[12] + data[13]
        else:
            value = data[12]
        
        ## Make index
        indA  = self._makePairIndex_2PCF(tomoA1, tomoA2, signA, order=-1)
        indB  = self._makePairIndex_2PCF(tomoB1, tomoB2, signB, order=-1)
        indA  = binA + self.N_theta * indA
        indB  = binB + self.N_theta * indB
        d_tot = indA.max() + 1
        
        ## Fill the other triangle
        cov = np.zeros((d_tot, d_tot), dtype=float)
        cov[indA, indB] = value
        ind = np.arange(d_tot, dtype=int)
        cov[ind, ind] *= 0.5
        cov += cov.T
        return cov
    
    def _covListToMatrix_BP(self, data, cleanNaN=True, CTag='tot'):
        if cleanNaN:
            ind = np.isnan(data)
            data[ind] = 0.0
        
        tomoA1 = asInt(data[0]) - 1
        tomoA2 = asInt(data[1]) - 1
        tomoB1 = asInt(data[5]) - 1
        tomoB2 = asInt(data[6]) - 1
        binA   = asInt(data[4])
        binB   = asInt(data[9])
        if CTag == 'tot':
            value = data[10:].sum(axis=0)
        else:
            value = data[10]
        
        ## Make index
        indA  = self._makePairIndex_BP(tomoA1, tomoA2, order=-1)
        indB  = self._makePairIndex_BP(tomoB1, tomoB2, order=-1)
        indA  = binA + self.N_ell * indA
        indB  = binB + self.N_ell * indB
        d_tot = indA.max() + 1
        
        ## Fill the other triangle
        cov = np.zeros((d_tot, d_tot), dtype=float)
        cov[indA, indB] = value
        ind = np.arange(d_tot, dtype=int)
        cov[ind, ind] *= 0.5
        cov += cov.T
        return cov
    
    def _get2PCFIndexForCut(self, statsTag):
        labConv = wtp.LabelConvention()
        wTh = True if labConv.w in statsTag else False
        gT  = True if labConv.gamma_t  in statsTag else False
        xiP = True if labConv.xi_p in statsTag else False
        xiM = True if labConv.xi_m in statsTag else False
        ind = [wTh]*self.N_theta*self.__nbPairsNN + [gT]*self.N_theta*self.__nbPairsNG + [xiP]*self.N_theta*self.__nbPairsGG + [xiM]*self.N_theta*self.__nbPairsGG
        return ind
    
    def _getBPIndexForCut(self, statsTag):
        labConv = wtp.LabelConvention()
        statsList = statsTag.split('+')
        Pnn  = True if labConv.P_nn in statsTag else False
        PneE = True if labConv.P_ne_E in statsTag else False
        PneB = True if labConv.P_ne_B in statsTag else False
        PeeE = True if labConv.P_ee_E in statsTag else False
        PeeB = True if labConv.P_ee_B in statsTag else False
        
        if PneE or PeeE == True:
            ind = [Pnn]*self.N_ell*self.__nbPairsNN + [PneE]*self.N_ell*self.__nbPairsNG + [PeeE]*self.N_ell*self.__nbPairsGG
        else:
            ind = [Pnn]*self.N_ell*self.__nbPairsNN + [PneB]*self.N_ell*self.__nbPairsNG + [PeeB]*self.N_ell*self.__nbPairsGG
        return ind
    
    def _getCOSEBIIndexForCut(self, statsTag):
        labConv = wtp.LabelConvention()
        En  = True if labConv.E_n in statsTag else False
        Bn  = True if labConv.B_n in statsTag else False
        ind = [En]*self.nbModes*self.__nbPairsGG + [Bn]*self.nbModes*self.__nbPairsGG
        return ind
    
    @classmethod
    def getCategory(cls, statsTag):
        labConv = wtp.LabelConvention()
        statsList = statsTag.split('+')
        is2PCF    = False
        isBPE     = False
        isBPB     = False
        isCOSEBI  = False
        
        for stats in statsList:
            if stats in [labConv.w, labConv.gamma_t, labConv.gamma_x, labConv.xi_p, labConv.xi_m]:
                is2PCF   = is2PCF or True
            elif stats in [labConv.P_ne_E, labConv.P_ee_E]:
                isBPE    = isBPE or True
            elif stats in [labConv.P_ne_B, labConv.P_ee_B]:
                isBPB    = isBPB or True
            elif stats in [labConv.E_n, labConv.B_n]:
                isCOSEBI = isCOSEBI or True
        
        category = 1*int(is2PCF) + 2*int(isBPE) + 4*int(isBPB) + 8*int(isCOSEBI)
        return category

    def _makeCov_list(self, name, statsTag, cleanNaN=True, CTag='tot', verbose=True):
        covList  = loadAscii(name, verbose=verbose)
        category = TwoPointBuilder.getCategory(statsTag)
        
        if category == 1:
            cov = self._covListToMatrix_2PCF(covList, cleanNaN=cleanNaN, CTag=CTag)
            ind = self._get2PCFIndexForCut(statsTag)
        elif category in [0, 2, 4]:
            cov = self._covListToMatrix_BP(covList, cleanNaN=cleanNaN, CTag=CTag)
            ind = self._getBPIndexForCut(statsTag)
        elif category == 8:
            raise NotImplementedError('Reading COSEBI cov from list format not implemented')
        else:
            raise ValueError('statsTag = \"%s\" not allowed' % statsTag)
        
        cov = cov[ind].T[ind].T 
        return cov
    
    def setCov(self, tag, name=None, statsTag=None, cleanNaN=True, CTag='tot', verbose=True):
        if tag is None or tag == 'none':
            print('No covariance')
            return
        elif tag == 'Flinc':
            self.cov = self._makeCov_Flinc(name, statsTag, verbose=verbose) ## Consider `name` as `aves`
        elif tag == 'list':
            self.cov = self._makeCov_list(name, statsTag, cleanNaN=cleanNaN, CTag=CTag, verbose=verbose)
        elif tag == 'variable':
            self.cov = name ## Consider `name` as the variable which contains the covariance
      
        ## Otherwise, consider `name` as the path to the file
        ## which contains the covariance matrix
        else:
            try:
                if name[-4:] == '.npy':
                    self.cov = np.load(name)
                elif name[-4:] == '.fit' or name[-5:] == '.fits':
                    self.cov = fits.getdata(name, 1)
                else:
                    self.cov = np.loadtxt(name)
                if verbose:
                    print('Loaded \"%s\"' % name)
            except:
                raise OSError('\"%s\" not found' % name)
        return
    
    ## n(z)
    def _makeKernel(self, name, nOfZNameList, nGalList, sigmaEpsList):
        if len(nOfZNameList) == 0:
            return None
        
        nOfZName = nOfZNameList[0]
        if nOfZName[-4:] == '.npy':
            nOfZList = [np.load(nOfZName) for nOfZName in nOfZNameList]
        elif nOfZName[-4:] == '.fit' or nOfZName[-5:] == '.fits':
            nOfZList = [fits.getdata(nOfZName, 1) for nOfZName in nOfZNameList]
            nOfZList = [[data.field(0), data.field(1)] for data in nOfZList]
        else:
            nOfZList = [loadAscii(nOfZName, verbose=self.verbose) for nOfZName in nOfZNameList]
        
        zArr     = nOfZList[0][0]
        nOfZList = [nOfZ[1][:-1] for nOfZ in nOfZList]
        z_lower  = zArr[:-1]
        z_upper  = zArr[1:]
        z_middle = 0.5 * (z_lower + z_upper)
        kernel   = twopoint.NumberDensity(name, z_lower, z_middle, z_upper, nOfZList, ngal=nGalList, sigma_e=sigmaEpsList)
        return kernel
      
    def setNOfZ(self, nOfZNameList, nGalList=None, sigmaEpsList=None):
        if nOfZNameList is None or nOfZNameList == 'none':
            print('No n(z)')
            return
        
        nbTomo = self.nbTomoN + self.nbTomoG
        
        ## Assert
        if nbTomo == len(nOfZNameList):
            nOfZNameListN = nOfZNameList[:self.nbTomoN]
            nOfZNameListG = nOfZNameList[self.nbTomoN:]
        else:
            raise AssertionError('Bad length of nOfZNameList')
        
        ## Assert
        if nGalList is None:
            nGalListN = None
            nGalListG = None
        elif nbTomo == len(nGalList):
            nGalListN = nGalList[:self.nbTomoN]
            nGalListG = nGalList[self.nbTomoN:]
        else:
            raise AssertionError('Bad length of nGalList')
        
        ## Assert
        if sigmaEpsList is None:
            pass
        elif len(sigmaEpsList) == self.nbTomoG:
            pass
        elif len(sigmaEpsList) == nbTomo:
            sigmaEpsList = sigmaEpsList[self.nbTomoN:]
        else:
            raise AssertionError('Bad length of sigmaEpsList')
        
        ## Make
        labConv = wtp.LabelConvention()
        kernelN = self._makeKernel(labConv.lens, nOfZNameListN, nGalListN, None)
        kernelG = self._makeKernel(labConv.source, nOfZNameListG, nGalListG, sigmaEpsList)
        
        if kernelN is None:
            self.kernelList = [kernelG]
        elif kernelG is None:
            self.kernelList = [kernelN]
        else:
            self.kernelList = [kernelN, kernelG]
        return
    
    ## Build up & save
    def _makeTwoPoint_withCov(self, labConv, statsTag_c):
        tomoAngDict = self.makeTomoAngDict()
        statsList_c  = statsTag_c.split('+')
        statsList_c_complete = labConv.kernelTypeDict.keys()
      
        for stats_c in statsList_c:
            if stats_c not in statsList_c_complete:
                raise ValueError('\"%s\" not allowed' % statsTag_c)
    
        statsNameDict = {}
        scBuilder = twopoint.SpectrumCovarianceBuilder()
        binInd    = 0
      
        for stats_c, line1, line2 in zip(labConv.kernelTypeDict.keys(), labConv.kernelTypeDict.values(), tomoAngDict.values()):
            ker1, ker2, type1, type2, unit = line1
            pairList, N_ang, angle         = line2
            
            if stats_c in statsList_c:
                statsNameDict[(ker1, ker2, type1, type2)] = stats_c
                for i, j in pairList:
                    for angInd in range(N_ang):
                        x = angle[angInd]
                        y = self.mean[binInd]
                        binInd += 1
                        scBuilder.add_data_point(ker1, ker2, type1, type2, i+1, j+1, x, angInd+1, y)
    
        ## Make TP
        scBuilder.set_names(statsNameDict)
        spectra, cov_info = scBuilder.generate(self.cov, 'arcmin')
        
        TP = wtp.TwoPointWrapper.from_spectra(spectra, kernels=self.kernelList, covmat_info=cov_info) ## kernels & covmat_info can be None
        return TP
    
    def _makeTwoPoint_withoutCov(self, labConv, statsTag_c):
        tomoAngDict = self.makeTomoAngDict()
        statsList_c  = statsTag_c.split('+')
        statsList_c_complete = labConv.kernelTypeDict.keys()
        
        for stats_c in statsList_c:
            if stats_c not in statsList_c_complete:
                raise ValueError('\"%s\" not allowed' % statsTag_c)
        
        spectra = []
        binInd  = 0
        
        for stats_c, line1, line2 in zip(labConv.kernelTypeDict.keys(), labConv.kernelTypeDict.values(), tomoAngDict.values()):
            ker1, ker2, type1, type2, unit = line1
            pairList, N_ang, angle         = line2
            
            if stats_c in statsList_c:
                sBuilder = wtp.SpectrumBuilder()
              
                for i, j in pairList:
                    value   = self.mean[binInd:binInd+N_ang]
                    binInd += N_ang
                    sBuilder.addTomo(i, j, angle, value)
          
                spec = sBuilder.makeSpectrum(stats_c, (type1, type2), unit, kernels=(ker1, ker2)) ## kernels can be None
                spectra.append(spec)
        
        ## Make
        TP = wtp.TwoPointWrapper.from_spectra(spectra, kernels=self.kernelList, covmat_info=None) ## kernels & covmat_info can be None
        return TP
    
    def _makeTwoPoint_onlyNOfZ(self):
        TP = wtp.TwoPointWrapper.from_spectra([], kernels=self.kernelList, covmat_info=None) ## spectra = [], kernels & covmat_info can be None
        return TP
    
    def makeTwoPoint(self, labConv, statsTag_c):
        if self.mean is None and self.cov is None:
          TP = self._makeTwoPoint_onlyNOfZ()
        elif self.cov is None:
          TP = self._makeTwoPoint_withoutCov(labConv, statsTag_c)
        else:
          TP = self._makeTwoPoint_withCov(labConv, statsTag_c)
        return TP
    
    @classmethod
    def addCovToTwoPoint(cls, TP, statsList, dimList, cov):
        cov_info = twopoint.CovarianceMatrixInfo('COVMAT', statsList, dimList, cov)
        TP2 = wtp.TwoPointWrapper.from_spectra(TP.spectra, kernels=TP.kernels, covmat_info=cov_info)
        return TP2
        
###############################################################################
## Auxiliary functions

def loadAscii(name, sep=None, cmt='#', verbose=True):
    data = np.loadtxt(name, comments=cmt, delimiter=sep)
    if verbose is True:
        print('Loaded \"%s\"' % name)
    return data.T

def asInt(a, copy=True):
    if np.isscalar(a):
        return int(a)
    if type(a) is list:
        return np.array(a, dtype=int)
    return a.astype(int, copy=copy)
  
def pairIndex(N, i, j):
    ind  = (N + (N+1-i)) * i // 2
    ind += j - i
    return ind

def getPrefix():
  return 'data/KiDS/kcap/Input_mean_cov_nOfZ/'

###############################################################################
## Main function snippet

def saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict={},
        meanTag=None, meanName=None,
        covTag=None, covName=None,
        nOfZNameList=None, nGalList=None, sigmaEpsList=None,
        saveName=None
    ):
    
    labConv = wtp.LabelConvention()
    
    ## Custom
    TPBuilder = TwoPointBuilder(
        nbTomoN=nbTomoN, nbTomoG=nbTomoG,
        N_theta=N_theta, theta_min=theta_min, theta_max=theta_max,
        N_ell=N_ell, ell_min=ell_min, ell_max=ell_max,
        nbModes=nbModes,
        prefix_Flinc=prefix_Flinc,
        prefix_CosmoSIS=prefix_CosmoSIS
    )
    
    ## Labels
    statsList, scArgs = labConv.makeScaleCutsArgs(scDict)
    if statsList is not None:
        statsTag_c = labConv.defaultToCustomStatsTag('+'.join(statsList))
        statsTag   = labConv.customToDefaultStatsTag(statsTag_c)
    else:
        statsTag_c = None
        statsTag   = None
        
    ## Make
    TPBuilder.setMean(meanTag, name=meanName, statsTag=statsTag)
    TPBuilder.setCov(covTag, name=covName, statsTag=statsTag)
    TPBuilder.setNOfZ(nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList)
    TP = TPBuilder.makeTwoPoint(labConv, statsTag_c)
    
    ## Cut
    TP.cutScales(cutCross=scArgs[0], statsTag_tomoInd_tomoInd_list=scArgs[1], statsTag_binIndList_dict=scArgs[2])
    TP.keepScales(statsTag_tomoInd1_tomoInd2__angMin_angMax_dict=scArgs[3], statsTag__angMin_angMax_dict=scArgs[4])
    
    ## Save
    TP.to_fits(saveName, overwrite=True, clobber=True)
    print('Saved \"%s\"' % saveName)
    return

###############################################################################
## Save 1 - Flinc & others

def makeMilvusNOfZ(sourceOnly=False):
    prefix = 'data/mockFootprint/milvus/MFP_selection/'
    nOfZNameList = [
        '%snOfZ_hist_BOSSA_tomo0.dat' % prefix, 
        '%snOfZ_hist_BOSSA_tomo1.dat' % prefix, 
        '%snOfZ_hist_KiDSVDA_tomo0.dat' % prefix, 
        '%snOfZ_hist_KiDSVDA_tomo1.dat' % prefix, 
        '%snOfZ_hist_KiDSVDA_tomo2.dat' % prefix, 
        '%snOfZ_hist_KiDSVDA_tomo3.dat' % prefix, 
        '%snOfZ_hist_KiDSVDA_tomo4.dat' % prefix
    ]
    nGalList     = [0.016702, 0.016725,   0.85651964, 1.5625674, 2.2357926, 1.5223054, 1.384634]
    sigmaEpsList = [0.284605] * 5
    
    if sourceOnly:
        nOfZNameList = nOfZNameList[2:]
        nGalList     = nGalList[2:]
    
    return nOfZNameList, nGalList, sigmaEpsList

def makeBucerosNOfZ(sourceOnly=False):
    prefix = 'data/mockFootprint/buceros/MFP_selection/'
    nOfZNameList = [
        '%snOfZ_hist_BOSS_2dFLenS_combined_tomo0.dat' % prefix, 
        '%snOfZ_hist_BOSS_2dFLenS_combined_tomo1.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo0.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo1.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo2.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo3.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_tomo4.dat' % prefix
    ]
    nGalList     = [0.010839, 0.011899,   0.856520, 1.562567, 2.235793, 1.522305, 1.384634]
    sigmaEpsList = [0.284965, 0.278972, 0.288525, 0.279802, 0.289495]
    
    if sourceOnly:
        nOfZNameList = nOfZNameList[2:]
        nGalList     = nGalList[2:]
    
    return nOfZNameList, nGalList, sigmaEpsList

def makeBucerosBroadNOfZ(sourceOnly=False):
    prefix = 'data/mockFootprint/buceros/MFP_selection/'
    nOfZNameList = [
        '%snOfZ_hist_BOSS_2dFLenS_combined_broad_tomo0.dat' % prefix, 
        '%snOfZ_hist_BOSS_2dFLenS_combined_broad_tomo1.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_broad_tomo0.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_broad_tomo1.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_broad_tomo2.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_broad_tomo3.dat' % prefix, 
        '%snOfZ_hist_KiDS_varDepth_A_broad_tomo4.dat' % prefix
    ]
    nGalList     = [0.010839, 0.011899,   0.856520, 1.562567, 2.235793, 1.522305, 1.384634]
    sigmaEpsList = [0.284965, 0.278972, 0.288525, 0.279802, 0.289495]
    
    if sourceOnly:
        nOfZNameList = nOfZNameList[2:]
        nGalList     = nGalList[2:]
    
    return nOfZNameList, nGalList, sigmaEpsList

def makeKV450NOfZ():
    prefix = 'data/KiDS/KV-450_shear_CF/'
    nOfZNameList = [
        '%snofz/Nz_DIR_z0.1t0.3.asc' % prefix, 
        '%snofz/Nz_DIR_z0.3t0.5.asc' % prefix, 
        '%snofz/Nz_DIR_z0.5t0.7.asc' % prefix, 
        '%snofz/Nz_DIR_z0.7t0.9.asc' % prefix, 
        '%snofz/Nz_DIR_z0.9t1.2.asc' % prefix
    ]
    nGalList = [1]*5
    sigmaEpsList = [0.28]*5
    return nOfZNameList, nGalList, sigmaEpsList

def makeCOSEBICov1():
    name = 'data/KiDS/KV-450_COSEBIs/Covariance_nMaximum_20_0.50_300.00_nBins5_NoiseJustForNoise.ascii'
    data1 = loadAscii(name)
    name = 'data/KiDS/KV-450_COSEBIs/Cov_SSC_thKV450_nell_500_20_0.50-300.00_nBins_5.ascii'
    data2 = loadAscii(name)
    name = 'data/KiDS/KV-450_COSEBIs/Covariance_nMaximum_20_0.50_300.00_nBins5_sigma_m_0.0200.ascii'
    data3 = loadAscii(name)
    
    cov = data1 + data2 + data3
    return cov

def saveFitsTwoPoint_NOfZ_milvus():
    nOfZNameList, nGalList, sigmaEpsList = makeMilvusNOfZ()
    saveName = '%stwoPoint_None_mean_None_cov_None_nOfZ_milvus.fits' % getPrefix()
    
    saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict={},
        meanTag=None, meanName=None,
        covTag=None, covName=None,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def saveFitsTwoPoint_NOfZ_buceros():
    nOfZNameList, nGalList, sigmaEpsList = makeBucerosNOfZ()
    saveName = 'twoPoint_None_mean_None_cov_None_nOfZ_buceros.fits'
    #saveName = '%stwoPoint_None_mean_None_cov_None_nOfZ_buceros.fits' % getPrefix()
    
    saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict={},
        meanTag=None, meanName=None,
        covTag=None, covName=None,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def saveFitsTwoPoint_COSEBI_KV450():
    scDict = {'use_stats': 'En'.lower()}
    meanName = 'data/KiDS/KV-450_COSEBIs/En_nBins_5_0.50-300.00.ascii'
    cov = makeCOSEBICov1()
    nOfZNameList, nGalList, sigmaEpsList = makeKV450NOfZ()
    saveName = '%stwoPoint_En_mean_KV450_cov_KV450_nOfZ_KV450.fits' % getPrefix()
    
    saveFitsTwoPoint(
        nbTomoN=0, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=20,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict=scDict,
        meanTag='file', meanName=meanName,
        covTag='variable', covName=cov,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def saveFitsTwoPoint_En_theory():
    scDict = {'use_stats': 'En'.lower()}
    covName  = 'data/mockFootprint/otis/MFP_for_others/cov_th_En_obs.dat'
    nOfZNameList, nGalList, sigmaEpsList = makeMilvusNOfZ(sourceOnly=True)
    saveName = '%stwoPoint_En_mean_theoryNoiseFree_cov_theoryOtis_nOfZ_milvus.fits' % getPrefix()
    
    saveFitsTwoPoint(
        nbTomoN=0, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS='data/mockFootprint/milvus/cosmosis/',
        scDict=scDict,
        meanTag='CosmoSIS', meanName=None,
        covTag='file', covName=covName,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

###############################################################################
## Save 2 - BPCovTest

def saveFitsTwoPoint_BPCovTest():
    loadPrefix = 'data/mockFootprint/'
    savePrefix = getPrefix()
    scDict = {'use_stats': 'PeeE'.lower()}
    covNameList  = [
        '%svultur/MFP_for_others/cov_th_PeeE_obs.dat' % loadPrefix,
        '%szosterops/MFP_for_others/cov_th_PeeE_obs.dat' % loadPrefix,
        '%svultur/MFP_for_others/cov_sim_PeeE_obs.dat' % loadPrefix,
        '%szosterops/MFP_for_others/cov_sim_PeeE_obs.dat' % loadPrefix
    ]
    nOfZNameList, nGalList, sigmaEpsList = makeMilvusNOfZ(sourceOnly=True)
    saveNameList = [
        '%stwoPoint_PeeE_mean_theoryNoiseFree_cov_theoryVultur_nOfZ_milvus.fits' % savePrefix,
        '%stwoPoint_PeeE_mean_theoryNoiseFree_cov_theoryZosterops_nOfZ_milvus.fits' % savePrefix,
        '%stwoPoint_PeeE_mean_theoryNoiseFree_cov_simVultur_nOfZ_milvus.fits' % savePrefix,
        '%stwoPoint_PeeE_mean_theoryNoiseFree_cov_simZosterops_nOfZ_milvus.fits' % savePrefix
    ]
    
    for covName, saveName in zip(covNameList, saveNameList):
        saveFitsTwoPoint(
            nbTomoN=0, nbTomoG=5,
            N_theta=9, theta_min=0.5, theta_max=300,
            N_ell=8, ell_min=100, ell_max=1500,
            nbModes=5,
            prefix_Flinc=None,
            prefix_CosmoSIS='data/mockFootprint/milvus/cosmosis/',
            scDict=scDict,
            meanTag='CosmoSIS', meanName=None,
            covTag='file', covName=covName,
            nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
            saveName=saveName
        )
    return

def saveFitsTwoPoint_copySameMean():
    scDict = {'use_stats': 'PeeE'.lower()}
    loadPrefix = 'data/mockFootprint/'
    savePrefix = getPrefix()
    mean = wtp.TwoPointWrapper.from_fits('%stwoPoint_PeeE_mean_theoryNoise1_cov_theoryVultur_nOfZ_milvus.fits' % savePrefix, covmat_name=None).makeMeanVector()
    covNameList  = [
        '%szosterops/MFP_for_others/cov_th_PeeE_obs.dat' % loadPrefix,
        '%svultur/MFP_for_others/cov_sim_PeeE_obs.dat' % loadPrefix,
        '%szosterops/MFP_for_others/cov_sim_PeeE_obs.dat' % loadPrefix
    ]
    nOfZNameList, nGalList, sigmaEpsList = makeMilvusNOfZ(sourceOnly=True)
    saveNameList = [
        '%stwoPoint_PeeE_mean_theoryNoise1_cov_theoryZosterops_nOfZ_milvus.fits' % savePrefix,
        '%stwoPoint_PeeE_mean_theoryNoise1_cov_simVultur_nOfZ_milvus.fits' % savePrefix,
        '%stwoPoint_PeeE_mean_theoryNoise1_cov_simZosterops_nOfZ_milvus.fits' % savePrefix
    ]
    
    for covName, saveName in zip(covNameList, saveNameList):
        saveFitsTwoPoint(
            nbTomoN=0, nbTomoG=5,
            N_theta=9, theta_min=0.5, theta_max=300,
            N_ell=8, ell_min=100, ell_max=1500,
            nbModes=5,
            prefix_Flinc=None,
            prefix_CosmoSIS='data/mockFootprint/milvus/cosmosis/',
            scDict=scDict,
            meanTag='variable', meanName=mean,
            covTag='file', covName=covName,
            nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
            saveName=saveName
        )
    return

###############################################################################
## Save 3 - GGLTest

def saveAvgLensNOfZ():
    prefix = 'data/KiDS/KiDS-1000_GGL_LensCats_N_of_Z/N_of_Z/'
    name1  = '%sBOSS_n_of_z1.txt' % prefix
    name2  = '%sBOSS_n_of_z2.txt' % prefix
    data1  = loadAscii(name1)
    data2  = loadAscii(name2)
    zArr   = np.concatenate([data1[0][:-1], data2[0]])
    nOfZ_BK_1 = np.concatenate([data1[1][:-1], np.zeros_like(data2[1], dtype=float)]) 
    nOfZ_BK_2 = np.concatenate([np.zeros_like(data1[1][:-1], dtype=float), data2[1]])
    
    for i, nOfZ in enumerate([nOfZ_BK_1, nOfZ_BK_2]):
        name3 = '%snOfZ_BOSS_overlap_z%d.dat' % (prefix, i+1)
        f = open(name3, 'w')
        f.write('## n(z) for BOSS & KiDS overlap\n')
        for z, n in zip(zArr, nOfZ):
            f.write(' %.6f  %.6f\n' % (z, n))
        f.close()
        print('Saved \"%s\"' % name3)
    
    name1  = '%s2dFLenS_n_of_z1.txt' % prefix
    name2  = '%s2dFLenS_n_of_z2.txt' % prefix
    data1  = loadAscii(name1)
    data2  = loadAscii(name2)
    zArr   = np.concatenate([data1[0][:-1], data2[0]])
    nOfZ_2K_1 = np.concatenate([data1[1][:-1], np.zeros_like(data2[1], dtype=float)]) 
    nOfZ_2K_2 = np.concatenate([np.zeros_like(data1[1][:-1], dtype=float), data2[1]])
    
    for i, nOfZ in enumerate([nOfZ_2K_1, nOfZ_2K_2]):
        name3 = '%snOfZ_2dFLenS_overlap_z%d.dat' % (prefix, i+1)
        f = open(name3, 'w')
        f.write('## n(z) for 2dFLenS & KiDS overlap\n')
        for z, n in zip(zArr, nOfZ):
            f.write(' %.6f  %.6f\n' % (z, n))
        f.close()
        print('Saved \"%s\"' % name3)
    
    wgtList = [[4.631, 1.987], [5.302, 2.074]]
    
    for i, wgt, nOfZ_BK, nOfZ_2K in zip(range(2), wgtList, [nOfZ_BK_1, nOfZ_BK_2], [nOfZ_2K_1, nOfZ_2K_2]):
        nOfZ  = (wgt[0]*nOfZ_BK + wgt[1]*nOfZ_2K) / (wgt[0] + wgt[1])
        name3 = '%snOfZ_avgLens_z%d.dat' % (prefix, i+1)
        f = open(name3, 'w')
        f.write('## Averaged lens n(z) between BOSS & 2dFLenS\n')
        f.write('## Hand-calculating weights\n')
        for z, n in zip(zArr, nOfZ):
            f.write(' %.6f  %.6f\n' % (z, n))
        f.close()
        print('Saved \"%s\"' % name3)
    return

def showQuick(i):
    import matplotlib.pyplot as plt
    fig = plt.gcf()
    fig.clf()
    
    prefix = 'data/KiDS/KiDS-1000_GGL_LensCats_N_of_Z/N_of_Z/'
    wgtList = [[1.987, 4.631], [2.074, 5.302]]
    
    name1 = '%snOfZ_BOSS_overlap_z%d.dat' % (prefix, i+1)
    name2 = '%snOfZ_2dFLenS_overlap_z%d.dat' % (prefix, i+1)
    name3 = '%snOfZ_avgLens_z%d.dat' % (prefix, i+1)
    
    data1 = loadAscii(name1)
    data2 = loadAscii(name2)
    data3 = loadAscii(name3)
    
    plt.plot(data1[0], data1[1])
    plt.plot(data2[0], data2[1])
    plt.plot(data3[0], data3[1])
    return

def makeBOSSAnd2dFLenSNOfZ(ind=0):
    prefix1 = 'data/KiDS/KiDS-1000_GGL_LensCats_N_of_Z/N_of_Z/'
    prefix2 = 'data/mockFootprint/milvus/MFP_selection/'
    nOfZNameList = [
        '%snOfZ_BOSS_overlap_z1.dat' % prefix1, 
        '%snOfZ_BOSS_overlap_z2.dat' % prefix1, 
        '%snOfZ_2dFLenS_overlap_z1.dat' % prefix1, 
        '%snOfZ_2dFLenS_overlap_z2.dat' % prefix1, 
        '%snOfZ_avgLens_z1.dat' % prefix1, 
        '%snOfZ_avgLens_z2.dat' % prefix1
    ]
    nOfZNameList = nOfZNameList[2*ind:2*ind+2] +  [
        '%snOfZ_hist_KiDSVDA_tomo0.dat' % prefix2, 
        '%snOfZ_hist_KiDSVDA_tomo1.dat' % prefix2, 
        '%snOfZ_hist_KiDSVDA_tomo2.dat' % prefix2, 
        '%snOfZ_hist_KiDSVDA_tomo3.dat' % prefix2, 
        '%snOfZ_hist_KiDSVDA_tomo4.dat' % prefix2
    ]
    nGalList     = [0, 0] + [0.85651964, 1.5625674, 2.2357926, 1.5223054, 1.384634]
    sigmaEpsList = [0.284605] * 5
    return nOfZNameList, nGalList, sigmaEpsList

def saveFitsTwoPoint_NOfZ_GGLTest():
    saveNameList = [
        '%stwoPoint_None_mean_None_cov_None_nOfZ_BOSSOverlapKiDSVDA.fits' % getPrefix(), 
        '%stwoPoint_None_mean_None_cov_None_nOfZ_2dFLenSOverlapKiDSVDA.fits' % getPrefix(), 
        '%stwoPoint_None_mean_None_cov_None_nOfZ_avgLensKiDSVDA.fits' % getPrefix()
    ]
    
    for i, saveName in enumerate(saveNameList):
        nOfZNameList, nGalList, sigmaEpsList = makeBOSSAnd2dFLenSNOfZ(i)
        
        saveFitsTwoPoint(
            nbTomoN=2, nbTomoG=5,
            N_theta=9, theta_min=0.5, theta_max=300,
            N_ell=8, ell_min=100, ell_max=1500,
            nbModes=5,
            prefix_Flinc=None,
            prefix_CosmoSIS=None,
            scDict={},
            meanTag=None, meanName=None,
            covTag=None, covName=None,
            nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
            saveName=saveName
        )
    return

###############################################################################
## Save 4 - Methodology

def saveFitsTwoPoint_methodology_theory():
    loadPrefix  = 'data/mockFootprint/'
    savePrefix  = getPrefix()
    scDict      = {'use_stats': 'PneE PeeE'.lower()}
    covNameList = [
        '%sbuceros/cosmokids/thps_cov_kids1000_buceros_bp_apo_simple_paircount_obs_bandpower_E_apod_list.dat' % loadPrefix,
        '%segretta/cosmokids/thps_cov_kids1000_egretta_bp_apo_obs_bandpower_E_apod_list.dat' % loadPrefix,
        '%segretta/cosmokids/thps_cov_kids1000_egretta_bp_apo_obs_bandpower_E_apod_list_mbias.dat' % loadPrefix,
    ]
    nOfZNameList, nGalList, sigmaEpsList = makeBucerosBroadNOfZ()
    saveNameList = [
        '%stwoPoint_PneE+PeeE_mean_None_cov_theoryBuceros_nOfZ_bucerosBroad.fits' % savePrefix,
        '%stwoPoint_PneE+PeeE_mean_None_cov_theoryEgretta_nOfZ_bucerosBroad.fits' % savePrefix,
        '%stwoPoint_PneE+PeeE_mean_None_cov_theoryEgrettaMCorr_nOfZ_bucerosBroad.fits' % savePrefix,
    ]
    
    for covName, saveName in zip(covNameList, saveNameList):
        saveFitsTwoPoint(
            nbTomoN=2, nbTomoG=5,
            N_theta=9, theta_min=0.5, theta_max=300,
            N_ell=8, ell_min=100, ell_max=1500,
            nbModes=5,
            prefix_Flinc=None,
            prefix_CosmoSIS=None,
            scDict=scDict,
            meanTag=None, meanName=None,
            covTag='list', covName=covName,
            nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
            saveName=saveName
        )
    return

def saveFitsTwoPoint_methodology_sim():
    savePrefix = getPrefix()
    scDict     = {'use_stats': 'PneE PeeE'.lower()}
    avesList   = ['buceros', 'egretta']
    nOfZNameList, nGalList, sigmaEpsList = makeBucerosBroadNOfZ()
    saveNameList = [
        '%stwoPoint_PneE+PeeE_mean_None_cov_simBuceros_nOfZ_bucerosBroad.fits' % savePrefix,
        '%stwoPoint_PneE+PeeE_mean_None_cov_simEgretta_nOfZ_bucerosBroad.fits' % savePrefix
    ]
    
    for aves, saveName in zip(avesList, saveNameList):
        saveFitsTwoPoint(
            nbTomoN=2, nbTomoG=5,
            N_theta=9, theta_min=0.5, theta_max=300,
            N_ell=8, ell_min=100, ell_max=1500,
            nbModes=5,
            prefix_Flinc=None,
            prefix_CosmoSIS=None,
            scDict=scDict,
            meanTag=None, meanName=None,
            covTag='Flinc', covName=aves,
            nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
            saveName=saveName
        )
    return

###############################################################################
## Sandbox

def saveFitsTwoPoint_test():
    scDict    = {'use_stats': 'PneE PeeE'.lower()}
    #scDict   = {'use_stats': 'PneE'.lower(), 'cut_pair_pnee': '1+1 2+1 2+2 2+3', 'keep_ang_peee': '100 1500', 'keep_ang_pnee': '100 1070'}
    covName  = 'data/mockFootprint/buceros/cosmokids/thps_cov_kids1000_buceros_bp_apo_simple_paircount_obs_bandpower_E_apod_list.dat'
    nOfZNameList, nGalList, sigmaEpsList = makeBucerosBroadNOfZ()
    saveName = 'twoPoint_test.fits'
    
    saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict=scDict,
        meanTag='CosmoSIS', meanName=None,
        covTag='list', covName=covName,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def makeKV450MeanAndCov():
    prefix = 'data/KiDS/KV-450_shear_CF/'
    
    name1  = '%sdata_vector/KV450_reweight_3x4x4_v2_good_xipm_mcor_5bin.dat' % prefix
    data1  = loadAscii(name1)
    data_p = data1[1:].T[:9].T.flatten()
    data_m = data1[1:].T[9:].T.flatten()
    mean   = np.concatenate([data_p, data_m])
    
    #WARNING Not sure about the convention
    name2 = '%scovariance/cov_analytic_montepython_mcorr.txt' % prefix
    cov   = loadAscii(name2)
    ind   = np.arange(2*9*15, dtype=int).reshape(2, 9, 15)
    np.swapaxes(ind, 1, 2)
    ind   = ind.flatten()
    cov   = cov[ind].T[ind]
    return mean, cov

def saveFitsTwoPoint_KV450():
    scDict = {
        'use_stats': 'xiP xiM'.lower(),
        'keep_ang_xiP'.lower(): '0.5 75',
        'keep_ang_xiM'.lower(): '4 300'
    }
    mean, cov = makeKV450MeanAndCov()
    nOfZNameList, nGalList, sigmaEpsList = makeKV450NOfZ()
    saveName = '%stwoPoint_xiPM_mean_KV450Data_cov_ForKV450_nOfZ_KV450Data.fits' % getPrefix()
    
    saveFitsTwoPoint(
        nbTomoN=0, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict=scDict,
        meanTag='variable', meanName=mean,
        covTag='variable', covName=cov,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def saveFitsTwoPoint_list():
    scDict   = {'use_stats': 'PneE PeeE'.lower()}
    covName  = 'data/mockFootprint/buceros/cosmokids/thps_cov_kids1000_buceros_bp_apo_simple_paircount_obs_bandpower_E_noap_list.dat'
    nOfZNameList, nGalList, sigmaEpsList = makeBucerosNOfZ()
    saveName = 'twoPoint_BP_mean_theoryNoiseFree_cov_theoryBuceros_nOfZ_buceros.fits' #% getPrefix()
    
    saveFitsTwoPoint(
        nbTomoN=2, nbTomoG=5,
        N_theta=9, theta_min=0.5, theta_max=300,
        N_ell=8, ell_min=100, ell_max=1500,
        nbModes=5,
        prefix_Flinc=None,
        prefix_CosmoSIS=None,
        scDict=scDict,
        meanTag='CosmoSIS', meanName=None,
        covTag='list', covName=covName,
        nOfZNameList=nOfZNameList, nGalList=nGalList, sigmaEpsList=sigmaEpsList,
        saveName=saveName
    )
    return

def saveAsciiFakeCovariance():
    name = 'data/mockFootprint/milvus/cosmokids/thps_cov_kids1000_milvus_obs_complex_list.dat'
    data = loadAscii(name)
    
    ind = data[0] == 1
    data[0][~ind] += 1
    ind = data[1] == 1
    data[1][~ind] += 1
    ind = data[6] == 1
    data[6][~ind] += 1
    ind = data[7] == 1
    data[7][~ind] += 1
    
    name = 'fakeCovariance.dat'
    f = open(name, 'w')
    for line in data.T:
        f.write('%d  %d  %d  %d  %d  %d    %d  %d  %d  %d  %d  %d    % .10e  % .10e  % .10e\n' %\
            (line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14]))
    
    f.close()
    print('Saved \"%s\"' % name)
    return

def saveFitsTwoPoint_fakeCTerm():
    prefix = '../KCAP_scale_cuts/kv450_data/'
    name3  = '%ssystematics/KV450_ALL_c12_treecorr.out' % prefix
    data3  = loadAscii(name3)
    cTerm  = data3[3]
    
    tpType4 = twopoint.Types.galaxy_position_fourier    ## GPF
    tpType5 = twopoint.Types.galaxy_shear_emode_fourier ## GEF
    tpType6 = twopoint.Types.galaxy_shear_bmode_fourier ## GBF
    TPBuilder = TwoPointBuilder()
    
    sBuilder = wtp.SpectrumBuilder()
    for i in range(5):
        for j in range(i, 5):
            sBuilder.addTomo(i, j, TPBuilder.theta, cTerm)
    spec1 = sBuilder.makeSpectrum('2d_cterm', (tpType4, tpType4), None)
    
    sBuilder = wtp.SpectrumBuilder()
    for i in range(5):
        for j in range(i, 5):
            sBuilder.addTomo(i, j, TPBuilder.theta, cTerm)
    spec2 = sBuilder.makeSpectrum('cterm_cos', (tpType5, tpType5), None)
    
    sBuilder = wtp.SpectrumBuilder()
    for i in range(5):
        for j in range(i, 5):
            sBuilder.addTomo(i, j, TPBuilder.theta, cTerm)
    spec3 = sBuilder.makeSpectrum('cterm_sin', (tpType6, tpType6), None)
    
    TP = wtp.TwoPointWrapper.from_spectra([spec1, spec2, spec3], kernels=None, covmat_info=None)
    
    name = '%stwoPoint_fakeCTerm.fits' % getPrefix()
    TP.to_fits(name, overwrite=True, clobber=True)
    print('Saved \"%s\"' % name)
    return

###############################################################################
## Check

def initialize():
    prefix = getPrefix()
    #name = '%stwoPoint_None_mean_None_cov_None_nOfZ_milvus.fits' % prefix
    #name = '%stwoPoint_En_mean_KV450_cov_KV450_nOfZ_KV450.fits' % prefix
    #name = '%stwoPoint_PeeE_mean_theoryNoiseFree_cov_theoryVultur_nOfZ_milvus.fits' % prefix
    #name = '%stwoPoint_PneE+PeeE_mean_None_cov_theoryBuceros_nOfZ_bucerosBroad.fits' % prefix
    name = 'twoPoint_test.fits'
    return name

def printTwoPointHDU(name=None, ind=1):
    name = initialize() if name is None else name
    hdr  = fits.getheader(name, ind)
    data = fits.getdata(name, ind)
    
    print()
    print(hdr.tostring(sep='\n'))
    print(data)
    return 

def _checkExtensions_fromFile(hdrList):
    str_ = [hdr['EXTNAME'] for hdr in hdrList]
    print('Name of extensions:')
    print('  %s' % '  '.join(str_))
    return
  
def _checkCovariance_fromFile(hdrList):
    try:
        extNameList = [hdr['EXTNAME'] for hdr in hdrList]
    except:
        print('  One of the extensions doesn\'t have \"EXTNAME\" key.')
        print('  That is not normal. Exit.')
        return
    
    print('Check covariance:')
    
    try:
        ind = extNameList.index('COVMAT')
        print('  Passed.')
        hdr   = hdrList[ind]
        count = 0
        stock = []
        print()
        print('Dimension of the covariance:')
        
        while count >= 0:
            try:
                stock.append(hdr['STRT_%d' % count])
                count += 1
            except:
                count = -1
        count = len(stock)
        stock.append(hdr['NAXIS1'])
        
        for i in range(count):
            print('  %8s = %3d' % (hdr['NAME_%d' % i], stock[i+1]-stock[i]))
        print('     Total = %3d' % stock[-1])

    except:
        print('  Failed. Skipped.')
    return

def _checkSpectra_fromFile(hdrList):
    print('Dimension of the data vector:')
    print('     stats   dim  tomo_1  tomo_2  N_ang')
    
    for hdr in hdrList:
        try:
            print('  %8s   %3d      %2d      %2d     %2d' % (hdr['EXTNAME'], hdr['NAXIS2'], hdr['N_ZBIN_1'], hdr['N_ZBIN_2'], hdr['N_ANG']))
        except:
            continue
    return

def _checkKernels_fromFile(hdrList, dataList):
    print('Redshift ranges:')
    print('     kernel     z_min     z_max')
    
    for hdr, data in zip(hdrList, dataList):
        try:
            print('  %9s  %8f  %8f' % (hdr['EXTNAME'], data['Z_LOW'][0], data['Z_HIGH'][-1]))
        except:
            continue
    return 

def _checkNGal_fromFile(hdrList):
    print('Galaxy number densities:')
    
    for hdr in hdrList:
        try:
            dummy = hdr['NGAL_1']
            count = 1
            print('  %s' % hdr['EXTNAME'])
            while count > 0:
                try:
                    print('    NGAL%d = %8f' % (count, hdr['NGAL_%d' % count]))
                    count += 1
                except:
                    count = -1
        except:
            continue
    return

def printTwoPoint_fromFile(name=None):
    name    = initialize() if name is None else name
    HDUList = fits.open(name)
    
    ## Check default HDU
    print()
    print('Check default HDU:')
    
    if 'SIMPLE' in HDUList[0].header:
        HDUList = HDUList[1:]
        print('  Passed.')
    else:
        print('  No default HDU.')
        print('  Means that this file was not generated in the standard way.')
        print('  Will continue.')
    
    hdrList  = [HDU.header for HDU in HDUList]
    dataList = [HDU.data for HDU in HDUList]
    
    print()
    _checkExtensions_fromFile(hdrList)
    print()
    _checkCovariance_fromFile(hdrList)
    print()
    _checkSpectra_fromFile(hdrList)
    print()
    _checkKernels_fromFile(hdrList, dataList)
    print()
    _checkNGal_fromFile(hdrList)
    return 

def _printSpectrum(spectrum):
    print('name = %s' % spectrum.name)
    print('tomoInd1 = %s' % spectrum.bin1)
    print('tomoInd2 = %s' % spectrum.bin2)
    print('pair = %s' % spectrum.bin_pairs)
    ##print(spectrum.type1)
    ##print(spectrum.type2)
    print('nOfZ1 = %s' % spectrum.kernel1)
    print('nOfZ2 = %s' % spectrum.kernel2)
    print('angInd = %s' % spectrum.angular_bin)
    #print('ang = %s' % spectrum.angle)
    ##print(spectrum.angle_min)
    ##print(spectrum.angle_max)
    print('mean = %s' % spectrum.value)
    ##print(spectrum.npairs)
    ##print(spectrum.varxi)
    ##print(spectrum.windows)
    #print('std = %s' % spectrum.error)
    ##print(spectrum.metadata)
    #print('unit = %s' % spectrum.angle_unit)
    #print('binInd = %s' % spectrum.dv_index)
    ##print(spectrum.extra_cols)
    return

def _printKernel(kernel):
    print('name = %s' % kernel.name)
    print('zlow = %s' % kernel.zlow)
    print('z = %s' % kernel.z)
    print('zhigh = %s' % kernel.zhigh)
    print('nbNOfZ = %s' % kernel.nbin)
    print('nbZBins = %s' % kernel.nsample)
    print('nArr = %s' % kernel.nzs)
    print('n_gal = %s' % kernel.ngal)
    print('sigma_e = %s' % kernel.sigma_e)
    return

def _printCovMatInfo(covMatInfo):
    pass
    #print('name = %s' % covMatInfo.name)
    print('names = %s' % covMatInfo.names)
    print('lengths = %s' % covMatInfo.lengths)
    print('starts = %s' % covMatInfo.starts)
    print('covmat = %s' % covMatInfo.covmat)
    #print('diagonal = %s' % covMatInfo.diagonal)
    return

def printTwoPoint(TP, mean=True, cov=True, nOfZ=True):
    if mean:
        print()
        print('Spectra:')
        for spectrum in TP.spectra:
            print()
            _printSpectrum(spectrum)
    
    if cov:
        print()
        print('Covariance:')
        if hasattr(TP, 'covmat_info') and TP.covmat_info is not None:
            print()
            _printCovMatInfo(TP.covmat_info)
            print('Direct cov.shape = %s' % str(TP.covmat.shape))
        else:
            print()
            print('Did not find `covmat_info` attribute')
    
    if nOfZ:
        print()
        print('Kernels:')
        for kernel in TP.kernels:
            print()
            _printKernel(kernel)
    
    ##print(TP.windows)
    ##print(TP._spectrum_index)
    return

def printTwoPoint_fromObj(name=None, mean=True, cov=True, nOfZ=True):
    name = initialize() if name is None else name
    try:
        TP = wtp.TwoPointWrapper.from_fits(name, covmat_name='COVMAT')
    except:
        TP = wtp.TwoPointWrapper.from_fits(name, covmat_name=None)
    printTwoPoint(TP, mean=mean, cov=cov, nOfZ=nOfZ)
    return

def _checkExtensions_unitary(HDUList1, HDUList2):
    str1 = [HDU.header['EXTNAME'] for HDU in HDUList1]
    str2 = [HDU.header['EXTNAME'] for HDU in HDUList2]
    
    if len(str1) != len(str2):
        print()
        print('No, because extensions are different.')
        print()
        print('Extensions from File 1 are')
        print('  %s' % '  '.join(str1))
        print('Extensions from File 2 are')
        print('  %s' % '  '.join(str2))
        return 1
    
    for n1, n2 in zip(str1, str2):
        if n1 != n2:
            print()
            print('No, because extensions are different.')
            print()
            print('Extensions from File 1 are')
            print('  %s' % '  '.join(str1))
            print('Extensions from File 2 are')
            print('  %s' % '  '.join(str2))
            return 1
    return 0

def _checkCovariance_unitary(HDU1, HDU2):
    cov1 = HDU1.data
    cov2 = HDU2.data
    
    if cov1.shape != cov2.shape:
        print()
        print('No, because covariance sizes are different.')
        print()
        print('The size of the covariance from File 1 is')
        print('  %s' % str(cov1.shape))
        print('The size of the covariance from File 2 is')
        print('  %s' % str(cov2.shape))
        return 1
    
    nonZero1 = cov1 != 0
    nonZero2 = cov2 != 0
    
    if np.any(nonZero1 != nonZero2):
        print()
        print('No, because covariances are different.')
        print()
        print('I checked the positions of zeros and some are different.')
        print('You have to find where they are yourself though.')
        return 1
    
    value1 = cov1[nonZero1]
    value2 = cov2[nonZero2]
    threshold = 1e-6
    ind = np.fabs(value1 / value2 - 1) > threshold
    
    if ind.sum() > 0:
        print()
        print('No, because covariances are different.')
        print()
        print('I checked the non-zero elements and some differ by more than %g%%.' % (threshold*100))
        print('You have to find where they are yourself though.')
        return 1
        
    return 0

def _checkSpectrum_unitary(HDU1, HDU2):
    extName = HDU1.header['EXTNAME']
    
    if HDU1.header['NAXIS2'] != HDU2.header['NAXIS2']:
        print()
        print('No, because vector lengths are different in \"%s\".' % extName)
        print()
        print('The length of \"%s\" from File 1 is' % extName)
        print('  %d' % HDU1.header['NAXIS2'])
        print('The length of \"%s\" from File 2 is' % extName)
        print('  %d' % HDU2.header['NAXIS2'])
        return 1
    
    if HDU1.header['TFIELDS'] != HDU2.header['TFIELDS']:
        print()
        print('No, because the numbers of columns are different in \"%s\".' % extName)
        print()
        print('The number of columns in \"%s\" from File 1 is' % extName)
        print('  %d' % HDU1.header['TFIELDS'])
        print('The number of columns in \"%s\" from File 2 is' % extName)
        print('  %d' % HDU2.header['TFIELDS'])
        return 1
    
    binA1   = HDU1.data.field(0)
    binA2   = HDU2.data.field(0)
    binB1   = HDU1.data.field(1)
    binB2   = HDU2.data.field(1)
    
    if np.any(binA1 != binA2) or np.any(binA1 != binA2):
        print()
        print('No, because tomo bin indices are different in \"%s\".' % extName)
        print()
        print('You have to find where they are yourself though.')
        print('It can very likely be an issue of mis-ordering.')
        return 1
        
    angBin1 = HDU1.data.field(2)
    angBin2 = HDU2.data.field(2)
    
    if np.any(angBin1 != angBin2):
        print()
        print('No, because angular bin indices are different in \"%s\".' % extName)
        print()
        print('You have to find where they are yourself though.')
        print('It can very likely be an issue of mis-ordering.')
        return 1
        
    value1 = HDU1.data.field(3)
    value2 = HDU2.data.field(3)
    threshold = 1e-6
    ind = np.fabs(value1 / value2 - 1) > threshold
    
    if ind.sum() > 0:
        print()
        print('No, because mean values are different in \"%s\".' % extName)
        print()
        print('Some differ by more than %g%%.' % (threshold*100))
        print('You have to find where they are yourself though.')
        print('It can very likely be an issue of mis-ordering.')
        return 1
        
    ang1 = HDU1.data.field(4)
    ang2 = HDU2.data.field(4)
    threshold = 1e-6
    ind = np.fabs(ang1 / ang2 - 1) > threshold
    
    if ind.sum() > 0:
        print()
        print('No, because angle/scale/mode values are different in \"%s\".' % extName)
        print()
        print('Some differ by more than %g%%.' % (threshold*100))
        print('You have to find where they are yourself though.')
        print('It can very likely be an issue of mis-ordering.')
        return 1
        
    return 0

def _checkKernel_unitary(HDU1, HDU2):
    extName = HDU1.header['EXTNAME']
    
    if HDU1.header['NAXIS2'] != HDU2.header['NAXIS2']:
        print()
        print('No, because the numbers of bins are different in \"%s\".' % extName)
        print()
        print('The number of bins in \"%s\" from File 1 is' % extName)
        print('  %d' % HDU1.header['NAXIS2'])
        print('The number of bins in \"%s\" from File 2 is' % extName)
        print('  %d' % HDU2.header['NAXIS2'])
        return 1
    
    if HDU1.header['TFIELDS'] != HDU2.header['TFIELDS']:
        print()
        print('No, because the numbers of columns are different in \"%s\".' % extName)
        print()
        print('The number of columns in \"%s\" from File 1 is' % extName)
        print('  %d' % HDU1.header['TFIELDS'])
        print('The number of columns in \"%s\" from File 2 is' % extName)
        print('  %d' % HDU2.header['TFIELDS'])
        return 1
    
    for i in range(HDU1.header['TFIELDS']):
        value1 = HDU1.data.field(i)
        value2 = HDU2.data.field(i)
        nonZero1 = value1 != 0
        nonZero2 = value2 != 0
        
        if np.any(nonZero1 != nonZero2):
            print()
            print('No, because bin or pdf values of Column %d (starting with 0) are different in \"%s\".' % (i, extName))
            print()
            print('I checked the positions of zeros and some are different.')
            print('You have to find where they are yourself though.')
            return 1
        
        value1 = value1[nonZero1]
        value2 = value2[nonZero2]
        threshold = 1e-6
        ind = np.fabs(value1 / value2 - 1) > threshold
        
        if ind.sum() > 0:
            print()
            print('No, because bin or pdf values of Column %d (starting with 0) are different in \"%s\".' % (i, extName))
            print()
            print('I checked the non-zero elements and some differ by more than %g%%.' % (threshold*100))
            print('You have to find where they are yourself though.')
            return 1
        
    return 0

def unitaryTest(name1, name2):
    print('File 1 = \"%s\"' % name1)
    print('File 2 = \"%s\"' % name2)
    print('Are they identical?')
    
    HDUList1 = fits.open(name1)
    HDUList2 = fits.open(name2)
    
    ## Check default HDU
    if 'SIMPLE' in HDUList1[0].header:
        HDUList1 = HDUList1[1:]
    if 'SIMPLE' in HDUList2[0].header:
        HDUList2 = HDUList2[1:]
    
    if _checkExtensions_unitary(HDUList1, HDUList2):
        return
    
    for HDU1, HDU2 in zip(HDUList1, HDUList2):
        if 'COVMAT' in HDU1.header['EXTNAME']:
            if _checkCovariance_unitary(HDU1, HDU2):
                return
        elif 'nz' in HDU1.header['EXTNAME'].lower():
            if _checkKernel_unitary(HDU1, HDU2):
                return
        else:
            if _checkSpectrum_unitary(HDU1, HDU2):
                return
    
    print()
    print('Yes, they passed all tests!')
    return

def printUnitaryTest():
    name1 = initialize()
    #name2 = initialize()
    name2 = '%stwoPoint_xiPM_mean_KV450Data_cov_KV450_nofz_KV450Data.fits' % getPrefix()
    unitaryTest(name1, name2)
    return
    
###############################################################################

