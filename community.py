import numpy as np


class community:
    """Class for a community of auxotrophs, based on model described in van Vliet et al, PLoS Comp Bio 2020"""

    def __init__(self, c_up_T=1, c_up_P=1, c_prod_T=1, c_prod_P=1, mu_WT=None, mu_dT=None, mu_dP=None, **kwargs):
        """Initialize the community with a set of parameters
        
            c_up_T: factor that scales uptake rate of trp in dTrp
            c_up_P: factor that scales uptake rate of pro in dPro
            c_prod_T: factor that scales production rate of trp in dPro
            c_prod_P: factor that scales production rate of pro in dTrp
            mu_WT: growth rate of WT in h^-1
            mu_dT: growth rate of dTrp mutant in h^-1
            mu_dP: growth rate of dPro mutant in h^-1          
        
        """
        
        self.par_names = ['upPro', 'upTrp', 'DiffP', 'DiffT', 'lPro', 'lTrp', 'ic', 'mu0', 'rho', 'rho2d', 'cell_l', 'cell_w', 'beta']
          
        # parameters for the model, see S1 Table in the paper
        mu0 = 1.29/3600 
        self.upPro = 2.04 
        self.upTrp = 24.05 
        self.DiffP = 879 
        self.DiffT = 659 
        self.lPro = 1.59e-5
        self.lTrp = 6.04e-7 
        self.ic = 20 
        self.mu0 = mu0 if mu_WT is None else mu_WT/3600
        self.rho = 0.65 
        self.rho2d = 0.22 
        self.cell_l = 5.2 
        self.cell_w = 0.68 
        self.beta = 0.88
        
        #set scaling factors
        self.c_up_T = c_up_T
        self.c_up_P = c_up_P
        self.c_prod_T = c_prod_T
        self.c_prod_P = c_prod_P
        self.mu_dT = mu0 if mu_dT is None else mu_dT/3600
        self.mu_dP = mu0 if mu_dP is None else mu_dP/3600
        
        #set other parameters
        self.set_parameters(kwargs)
        
        #calc community properties
        self.calc_interaction_parameters()
        self.calc_eq_prop() 
        
        return None
            
    ## model initialization functions
    def set_parameters(self, dict):
        '''set model parameters from dictionary
        input:
            dict: dictionary with model parameters
        '''
        #assign custom variables
        for key, value in dict.items():
            if key in self.par_names:
                setattr(self, key, value)    
            else:
                Warning(("%s is an invalid parameter name, ignoring this" % key))  
        return None  
    
    
    def calc_interaction_parameters(self):
        #calculate interaction parameters of the community
        #DP depends on Pro, producer of Pro is DT, non-producer is DP
        self.muP = self.calc_mu_max(up_auxo = self.upPro * self.c_up_P, 
                                    leak_auxo = self.lPro, 
                                    up_prod = self.upPro, 
                                    leak_prod = self.lPro * self.c_prod_P,
                                    mu_auxo = self.mu_dP
                                    )
        
        self.rangeP = self.calc_int_range(up_auxo = self.upPro * self.c_up_P, 
                                          leak_auxo = self.lPro, 
                                          up_prod = self.upPro, 
                                          leak_prod = self.lPro * self.c_prod_P,
                                          D = self.DiffP,
                                          mu_auxo = self.mu_dP
                                          )
        
        self.nbP = self.calc_int_nb(self.rangeP)
        
        
        # self.sector_widthP = self.sector_width(up_auxo = self.upPro * self.c_up_P, 
        #                                   leak_auxo = self.lPro, 
        #                                   up_prod = self.upPro, 
        #                                   leak_prod = self.lPro * self.c_prod_P,
        #                                   D = self.DiffP,
        #                                   mu_auxo = self.mu_dP,
        #                                   rho = self.rho)
                                                  
        
        #DT depends on Trp, producer of Trp is DP, non-producer is DT
        self.muT = self.calc_mu_max(up_auxo = self.upTrp * self.c_up_T, 
                                    leak_auxo = self.lTrp, 
                                    up_prod = self.upTrp, 
                                    leak_prod = self.lTrp * self.c_prod_T,
                                    mu_auxo = self.mu_dT)
        
        self.rangeT = self.calc_int_range(up_auxo = self.upTrp * self.c_up_T, 
                                          leak_auxo = self.lTrp, 
                                          up_prod = self.upTrp, 
                                          leak_prod = self.lTrp * self.c_prod_T,
                                          D = self.DiffT,
                                          mu_auxo = self.mu_dT)        
        
        self.nbT = self.calc_int_nb(self.rangeT)
        
        # self.sector_widthT = self.sector_width(up_auxo = self.upTrp * self.c_up_T,
        #                                     leak_auxo = self.lTrp, 
        #                                     up_prod = self.upTrp, 
        #                                     leak_prod = self.lTrp * self.c_prod_T,
        #                                     D = self.DiffT,
        #                                     mu_auxo = self.mu_dT,
        #                                     rho = self.rho)
        
        return None
    
    
    def calc_mu_max(self, up_auxo=1, up_prod=1, leak_auxo=1, leak_prod=1, mu_auxo=None):
        """Calculate the maximum growth rate of the cell

        Following equation 10 & 11 in S1 Text
        
        Growth rate is returned in units of h^-1 
        
        """
        
        #as we have over-producers we can no longer assume that leakage is negligible compared to growth we thus need exact equations
        mu_auxo = self.mu0 if mu_auxo is None else mu_auxo 
        e_max = leak_prod * self.ic / (leak_prod + up_prod) #eq 8 in S1 Text
        i_max = self.int_of_ext(e_max, up=up_auxo, leak=leak_auxo, mu0=mu_auxo) #eq 7 in S1 Text
        mu_max = self.mu_of_int(i_max, mu0=mu_auxo) #eq 9 in S1 Text

        # theta = ((leak_prod * self.ic)/(2 * mu_auxo)) * ((up_auxo + leak_auxo)/(up_prod + leak_prod))
        # mu_max = 3600 * mu_auxo * theta * (np.sqrt(1 + 2 / theta) - 1)        
        return mu_max
    
    
    def calc_int_range(self, up_auxo=1, up_prod=1, leak_auxo=1, leak_prod=1, D=1, mu_auxo=None):
        '''calculate the interaction range of the cell
        
        Following equation 17 in S1 Text. Note there was an error made when simplifying eq 17 from Mathematica to Latex, that has been corrected here.
        
        Range is returned in units of um
        
        '''
        
        #as we have over-producers we can no longer assume that leakage is negligible compared to growth we thus need exact equations nr 17 in S1 Text
        r0_auxo = self.calc_r0(up_auxo, leak_auxo, D, self.rho)
        r0_prod = self.calc_r0(up_prod, leak_prod, D, self.rho)
        delta = leak_prod*self.ic * 2*r0_auxo * (up_auxo + leak_auxo) / (2*mu_auxo * (r0_prod + r0_auxo) * (up_prod + leak_prod))
        rel_leak = leak_auxo / (mu_auxo * delta)
        c1 = delta * (rel_leak + 1)

        num = delta * mu_auxo * (2*rel_leak - 1) * np.sqrt(4*delta + c1**2) - delta**2 * mu_auxo * (4/delta + 1 + 3 * rel_leak - 2 * rel_leak**2)  
        denom = mu_auxo * delta * (2 * delta * rel_leak * (rel_leak - 1) - 1)
            
        range = self.beta * r0_auxo * np.log(num / denom)   
            
        return range
    
    
    def calc_int_nb(self, range):
        '''convert interaction range to number of neighbors
        
        Following equation 21 in S1 Text
        
        '''
        range = (2 * range * (self.cell_l - self.cell_w) + np.pi * (range + self.cell_w/2)**2 - np.pi * (self.cell_w/2)**2) * self.rho2d
        return np.round(range)
        
    def calc_eq_prop(self):
        '''calculate equilibrium properties of the community
        
        Equilibrium frequency (eq) is calculated following equation 16 in S2 Text
        Equilibrium frequency in well mixed system (eq_wm) is calculated following equation 20 in S2 Text
        Clustering of DT (clusteringT) and DP (clusteringP) is calculated following equation 17 and 18 in S2 Text
        Community growth rate relative to growth rate in well mixed system (rel_growth) is calculated following equation 29 & 30 in S2 Text
        
        '''
        self.eq = (self.muT * (self.nbT - 2)/self.nbT + (self.muT/self.nbT - self.muP/self.nbP)) / (self.muT * (self.nbT - 2)/self.nbT + self.muP * (self.nbP - 2)/self.nbP)
        self.eq_wm = self.muT / (self.muT + self.muP)      
        
        self.clusteringT = (self.nbT - 2) / (self.nbT - 1)
        self.clusteringP = (self.nbP - 2) / (self.nbP - 1)

        pT = self.eq
        pP = 1-self.eq
        pPT = self.clusteringT * pP 
        pTP = self.clusteringP * pT
        growth = pT * pPT * self.muT + pP * pTP * self.muP
        growth_wm = self.eq_wm * (1-self.eq_wm) * self.muT + (1-self.eq_wm) * self.eq_wm * self.muP
                
        self.rel_growth_wt = growth / (3600 * self.mu0)        
        self.rel_growth_wm = growth / growth_wm
        
        return None
    
    def calc_growth_profiles(self, x_vec=np.linspace(0,10,100)):
        ''' Calculates growth profiles of the community, assuming the two types occupy different sides of the community.
        
        DT on left, x<0; DP on right, x>0
        '''
        
        #Calculate growth profile of DT
        lTrp_prod = self.lTrp * self.c_prod_T # production rate of Trp in DP
        upTrp_auxo = self.upTrp * self.c_up_T # uptake rate of Trp in DT
        ext_Trp = self.ext_profile_auxo(x_vec, 
                                        up_prod=self.upTrp, 
                                        leak_prod=lTrp_prod, 
                                        r0_prod=self.calc_r0(self.upTrp, lTrp_prod, self.DiffT, self.rho), 
                                        r0_auxo=self.calc_r0(self.upTrp, self.lTrp, self.DiffT, self.rho))
        
        int_Trp = self.int_of_ext(ext_Trp, up=upTrp_auxo, leak=self.lTrp, mu0=self.mu_dT)
        self.mu_profile_T = self.mu_of_int(int_Trp, mu0=self.mu_dT) 
        self.x_T = -x_vec
        
        #Calculate growth profile of DP
        lPro_prod = self.lPro * self.c_prod_P # production rate of Pro in DT
        upPro_auxo = self.upPro * self.c_up_P # uptake rate of Pro in DP
        ext_Pro = self.ext_profile_auxo(x_vec, 
                                        up_prod=self.upPro, 
                                        leak_prod=lPro_prod, 
                                        r0_prod=self.calc_r0(self.upPro, self.lPro, self.DiffP, self.rho), 
                                        r0_auxo=self.calc_r0(upPro_auxo, self.lPro, self.DiffP, self.rho))
        
        int_Pro = self.int_of_ext(ext_Pro, up=upPro_auxo, leak=self.lPro, mu0=self.mu_dP)
        
        self.mu_profile_P = self.mu_of_int(int_Pro, mu0=self.mu_dP)
        self.x_P = x_vec

        return None
    
    def report_properties(self):
        '''report properties of the community'''
        print(f"Frequency dT = {self.eq:.2f}")
        print(f"Clustering dT = {self.clusteringT:.2f}, Clustering dP = {self.clusteringP:.2f}")
        print(f"Growth relative to WT = {self.rel_growth_wt:.2f}")        
        print(f"Growth defect spatial community = {self.rel_growth_wm:.2f}")        
        print(f"mu_max dP = {self.muP:.2f}, mu_max dT = {self.muT:.2f}")
        print(f"range dP = {self.rangeP:.2f}um, range dT = {self.rangeT:.2f}um, range dP/dT = {self.rangeP/self.rangeT:.2f}")
        #print(f"sector width dP = {self.sector_widthP:.2f}um, sector width dT = {self.sector_widthT:.2f}um, sector width dP/dT = {self.sector_widthP/self.sector_widthT:.2f}")
        
        return None
    
    def ext_profile_prod(self, x_vec, leak_prod=1, up_prod=1, r0_prod=1, r0_auxo=1):
        ''' Calculate external concentration profile in area of producer
        
        Following equation 14 in S1 Text
        
        '''
        ext = self.ic * leak_prod * (r0_auxo - r0_prod * np.exp(x_vec/r0_prod) + r0_prod) / \
                ((leak_prod + up_prod) * (r0_auxo + r0_prod))
        return ext
        
    def ext_profile_auxo(self, x_vec, leak_prod=1, up_prod=1, r0_prod=1, r0_auxo=1):
        ''' Calculate external concentration profile in area of auxotroph
        
        Following equation 14 in S1 Text
        
        '''        
        ext = self.ic * leak_prod * r0_auxo * np.exp(-x_vec/r0_auxo) / \
                ((leak_prod + up_prod) * (r0_auxo + r0_prod))
        return ext     
    
   
        
    def calc_r0(self, up, leak, D, rho):
        ''' Calculate length scale of diffusion, r0
        
        Following equation 12 in S1 Text
        
        '''
        
        effective_D = D * (1 - rho)**2 / (rho * (1 + rho / 2))
        return np.sqrt( (effective_D) / (up + leak))
    
    def int_of_ext(self, ext, up=1, leak=1, mu0=1):
        ''' Convert external concentration to internal concentration
        
        Following equation 7 in S1 Text
        
        '''
        intC = ((up + leak) * ext - leak + np.sqrt(((up + leak)*ext + leak)**2 + 4 * (up + leak) * mu0 * ext) ) / (2 * (mu0 + leak))
        return intC
    
    def mu_of_int(self, intC, mu0=1):
        ''' Calculate growth rate from internal concentration
        
        Following Monod equation
        
        Output is in units of h^-1
        
        '''
        return 3600 * (mu0 * intC) / (1 + intC)    
    
    def ext_of_x_sector(self, x, L, up_auxo=1, up_prod=1, leak_auxo=1, leak_prod=1, D=1, rho=1):
        r0_prod = self.calc_r0(up_prod, leak_prod, D, rho)
        r0_auxo = self.calc_r0(up_auxo, leak_auxo, D, rho)
        
        num = (np.cosh(x/r0_auxo)*self.ic*leak_prod*r0_auxo) 
        denom = ((np.cosh(L/r0_auxo)*r0_auxo + np.sinh(L/r0_auxo)*r0_prod)*(leak_prod + up_prod))

        return num/denom

    def int_of_x_sector(self, x, L, up_auxo=1, up_prod=1, leak_auxo=1, leak_prod=1, D=1, rho=1, mu_auxo=1):
        ext = self.ext_of_x_sector(x, L, up_auxo=up_auxo, up_prod=up_prod, leak_auxo=leak_auxo, leak_prod=leak_prod, D=D, rho=rho)
        int = self.int_of_ext(ext, up=up_prod, leak=leak_prod, mu0=mu_auxo)
        return int

    def sector_width(self, up_auxo=1, up_prod=1, leak_auxo=1, leak_prod=1, D=1, rho=1, mu_auxo=1):
        
        L = np.linspace(0,1000,100001)
        
        int0 = self.int_of_x_sector(0, L, up_auxo=up_auxo, up_prod=up_prod, leak_auxo=leak_auxo, leak_prod=leak_prod, D=D, rho=rho, mu_auxo=mu_auxo)
        intL = self.int_of_x_sector(L, L, up_auxo=up_auxo, up_prod=up_prod, leak_auxo=leak_auxo, leak_prod=leak_prod, D=D, rho=rho, mu_auxo=mu_auxo)

        mu_at0 = self.mu_of_int(int0, mu0=mu_auxo)
        mu_atL =self.mu_of_int(intL, mu0=mu_auxo)

        rel_mu_at0 = mu_at0/mu_atL
        
        sector_width = L[np.argmin(np.abs(rel_mu_at0 - 0.2))]

        return sector_width
    
    def calc_growth_profiles_sectors(self):
        ''' Calculates growth profiles of the community, assuming sector growth
        
        DT on left, x<0; DP on right, x>0
        '''
        
        #Calculate growth profile of DT
        upTrp_auxo = self.upTrp * self.c_up_T # uptake rate of Trp in DT
        x_vec = np.linspace(0,self.sector_widthT,100)
        ext_Trp = self.ext_of_x_sector(x_vec, self.sector_widthT,
                                        leak_auxo = self.lTrp, 
                                            up_prod = self.upTrp, 
                                            leak_prod = self.lTrp * self.c_prod_T,
                                            D = self.DiffT,
                                            rho = self.rho)
        
        ext_Trp = np.flip(ext_Trp)
         
        int_Trp = self.int_of_ext(ext_Trp, up=upTrp_auxo, leak=self.lTrp, mu0=self.mu_dT)
        self.mu_profile_sector_T = self.mu_of_int(int_Trp, mu0=self.mu_dT) 
        self.x_sector_T = -x_vec
        
        #Calculate growth profile of DP
        upPro_auxo = self.upPro * self.c_up_P # uptake rate of Pro in DP
        x_vec = np.linspace(0,self.sector_widthP,100)
        ext_Pro = self.ext_of_x_sector(x_vec, self.sector_widthP,
                                          up_auxo = self.upPro * self.c_up_P, 
                                          leak_auxo = self.lPro, 
                                          up_prod = self.upPro, 
                                          leak_prod = self.lPro * self.c_prod_P,
                                          D = self.DiffP,
                                          rho = self.rho)
        
        ext_Pro = np.flip(ext_Pro)
        
        int_Pro = self.int_of_ext(ext_Pro, up=upPro_auxo, leak=self.lPro, mu0=self.mu_dP)
        
        self.mu_profile_sector_P = self.mu_of_int(int_Pro, mu0=self.mu_dP)
        self.x_sector_P = x_vec