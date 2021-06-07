import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import os

def pot_func(phi,V):
    pot=V*(1.+np.cos(phi))
    return pot

def pot_funcS(phi,V):
    pot=V*(1.+np.sin(phi))
    return pot

def Vij(p1,p2,g):
    #make abstraction of phi value (if possible)
    V12=g*(np.cos(p1-p2)-3.*np.cos(p1)*np.cos(p2))
    return V12

def Vijphi_rel(dphi,g):
    Vij=(-g/2.)*(np.cos(dphi))
    return Vij

def VijPhi_CM(dphi,g):
    Vij=(-g*3./2.)*(np.cos(dphi))
    return Vij

def pot_matrix(size):
    V_mat=np.zeros((size,size),float)
    for m in range(size):
        for mp in range(size):
            if m==(mp+1):
                V_mat[m,mp]=.5
            if m==(mp-1):
                V_mat[m,mp]=.5
    return V_mat

def gen_prob_dist(Ng,rho_phi):
    p = np.zeros((Ng, Ng, Ng),float)
    # Normalize:
    P_norm = np.zeros((Ng,Ng),float)
    for i0 in range(Ng):
        for i1 in range(Ng):
            di01=i0 - i1
            if di01 < 0:
                di01+=Ng
            for i2 in range(Ng):
                di12= i1- i2
                if di12 < 0:
                    di12 +=Ng
                p[i0,i1,i2]=rho_phi[di01]*rho_phi[di12]
                P_norm[i0,i2] += p[i0,i1,i2]
    for i0 in range(Ng):
        for i1 in range(Ng):
            for i2 in range(Ng):
                p[i0,i1,i2]=p[i0,i1,i2]/P_norm[i0,i2]
    return p

def gen_prob_dist_end(Ng,rho_phi):
    p = np.zeros((Ng, Ng),float)
    # Normalize:
    P_norm = np.zeros(Ng,float)
    for i0 in range(Ng):
        for i1 in range(Ng):
            di01=i0 - i1
            if di01 < 0:
                di01+=Ng
            p[i0,i1]=rho_phi[di01]
            P_norm[i0] += p[i0,i1]
    for i0 in range(Ng):
        for i1 in range(Ng):
            p[i0,i1]=p[i0,i1]/P_norm[i0]
    return p

def calculateOrientationalCorrelations(p1, p2):
    return np.cos(p1-p2)

class ChoiceMC(object):
    def __init__(self, m_max, P, g, MC_steps, N, Nskip=100, Nequilibrate=0, PIGS=False, T=1, B=1, V0=0., potentialField='transverse'):
        """
        Creates a ChoiceMC object. This object can be used to generate the density
        matrices and performed the PIMC method based on the inputs. The results
        are stored as attributes of the object.
        
        Example:
        PIMC = ChoiceMC(m_max=50, P=9, g=1, MC_steps=1000, N=3, PIGS=True, Nskip=100, Nequilibrate=100)
        This sets up the system using 50 grid points, 9 beads, an interaction strength of 1, 3 rotors,
        path integral ground state enabled, 100 skip steps and 100 equilibrium steps.

        Parameters
        ----------
        m_max : int
            The maximum size of the free rotor eigenstate basis used to construct
            the rotational density matrix.
        P : int
            The number of beads to use in the path-integral.
        g : float
            Interaction strength between the rotors.
        MC_steps : int
            The number of steps to use in the Monte Carlo method.
        N : int
            The number of rotors to be simulated.
        Nskip : int, optional
            The number of steps to skip when saving the trajectory. The default is 100.
        Nequilibrate : int, optional
            The number of steps to skip before the average properties are accumulated 
            to allow for system equilibration. The default is 0.
        PIGS : bool, optional
            Enables path-integral ground state calculations. The default is False.
        T : float, optional
            The system temperature. The default is 1.
        B : float, optional
            The rotational constant for the rotors in 1/cm. The default is 1.
        V0 : float, optional
            The external potential field for the system. The default is 0.
        potentialField : string, optional
            The type of external potential field for the system. The default is transverse.

        Returns
        -------
        All inputs mentioned above are stored as attributes in the system.
        self.beta: float
            The beta value based on the system temperature.
        self.tau: float
            The tau value for the path integral method based on the beta value
            and the number of beads.
        self.Ngrid: int
            The number of steps to discretize the angle phi into.
        self.delta_phi: float
            The change in phi between the discretized phi values.
        self.potFunc: function handle
            The function handle that matches the desired potential function.
        """
        
        # Extracting information from kwargs for extra arguments
        self.Nskip = Nskip
        self.Nequilibrate = Nequilibrate
        self.PIGS = PIGS
        self.T = T
        self.B = B
        self.V0 = V0
        # Setting which potential function will be used
        if potentialField == "transverse":
            self.potFunc = pot_funcS
        elif potentialField == 'parallel':
            self.potFunc = pot_func
        else:
            raise Exception("Unrecognized potential model, allowed models are transverse or parallel")
        
        self.beta = 1./self.T
        self.P = P
        self.tau = self.beta/float(self.P)
        
        self.m_max = m_max
        self.Ngrid = 2 * m_max + 1
        self.delta_phi = 2. * np.pi / float(self.Ngrid)
        self.g = g
        self.MC_steps = MC_steps
        self.N = N
        
        # Creating a folder to save the output data
        self.path = os.path.join(os.getcwd(), "ChoiceMC_P" + str(P) + "_N" + str(N) + "_g" + str(round(g,3)) + "_MCSteps" + str(MC_steps)+"_V" + str(V0) + "_mMax" + str(m_max))
        try:
            os.mkdir(self.path)
        except FileExistsError:
            pass
        
        # Throwing a warning if the center bead will not be an integer value when pigs is enabled
        if self.P % 2 == 0 and self.PIGS:
            raise Warning("PIGS is enabled and an even number of beads was input, it is recommended to use an odd number")
    
    def createFreeRhoMarx(self):
        """
        Creates the free density matrix using the Marx method.
        This function will overwrite the self.rho_phi attribute of the object
        used during the Monte Carlo method. The Marx method is most accurate, and
        is used as the default in runMC.
        
        Returns
        -------
        self.free_rho_marx: numpy array
            Nx2 numpy array with the phi value in the 1st column and free density 
            matrix values in the 2nd column
        self.rho_phi: numpy array
            Nx1 numpy array with density matrix values, used in the runMC function
        """
        # Creating the free rotor density matrix using the MARX method, most accurate
        rho_phi=np.zeros(self.Ngrid,float)
        rho_marx_out=open(os.path.join(self.path,'rhofree_marx'),'w')
        self.free_rho_marx = np.zeros((self.Ngrid, 2),float)
        for i in range(self.Ngrid):
            dphi=float(i) * self.delta_phi
            integral = 0.
            for m in range(self.m_max):
                integral += np.exp(-1./(4.*self.tau*self.B)*(dphi+2.*np.pi*float(m))**2)
            for m in range(1,self.m_max):
                integral+=np.exp(-1./(4.*self.tau*self.B)*(dphi+2.*np.pi*float(-m))**2)
            integral*=np.sqrt(1./(4.*np.pi*self.B*self.tau))
            rho_phi[i]=integral
            rho_marx_out.write(str(dphi)+' '+str(integral)+'\n')
            self.free_rho_marx[i,:] = [dphi, integral]
        rho_marx_out.close()
        
        # Overwrites the current rho_phi to match the marx method
        self.rho_phi = rho_phi
        return self.free_rho_marx.copy()        
    
    def createFreeRhoSOS(self):
        """
        Creates the free density matrix using the SOS method.
        This function will overwrite the self.rho_phi attribute of the object
        used during the Monte Carlo method.
        
        Returns
        -------
        self.free_rho_sos: numpy array
            Nx2 numpy array with the phi value in the 1st column and free density 
            matrix values in the 2nd column
        self.rho_phi: numpy array
            Nx1 numpy array with density matrix values, used in the runMC function
        """
        # Creating the free rotor density matrix using the SOS method
        rho_phi=np.zeros(self.Ngrid,float)
        rhofree_sos_out = open(os.path.join(self.path,'rhofree_sos'),'w')
        self.free_rho_sos = np.zeros((self.Ngrid, 2),float)
        for i in range(self.Ngrid):
            dphi = float(i) * self.delta_phi
            integral = 0.
            for m in range(1,self.m_max):
                integral += (2. * np.cos(float(m) * dphi)) * np.exp(-self.tau * self.B * m**2)
            integral = integral / (2.*np.pi)
            integral = integral + 1./(2.*np.pi)
            rho_phi[i]=np.fabs(integral)
            rhofree_sos_out.write(str(dphi)+' '+str(rho_phi[i])+'\n')
            self.free_rho_sos[i,:] = [dphi, rho_phi[i]]
        rhofree_sos_out.close()
        
        # Overwrites the current rho_phi to match the sos method
        self.rho_phi = rho_phi
        return self.free_rho_sos.copy()
    
    def createFreeRhoPQC(self):
        """
        Creates the free density matrix using the PQC method.
        This function will overwrite the self.rho_phi attribute of the object
        used during the Monte Carlo method.
        
        Returns
        -------
        self.free_rho_pqc: numpy array
            Nx2 numpy array with the phi value in the 1st column and free density 
            matrix values in the 2nd column
        self.rho_phi: numpy array
            Nx1 numpy array with density matrix values, used in the runMC function
        """
        # Creating the free rotor density matrix using the PQC method
        rho_phi_pqc=np.zeros(self.Ngrid,float)
        rho_pqc_out=open(os.path.join(self.path,'rhofree_pqc'),'w')
        
        self.free_rho_pqc = np.zeros((self.Ngrid, 2),float)
        for i in range(self.Ngrid):
            dphi=float(i) * self.delta_phi
            rho_phi_pqc[i]=np.sqrt(1./(4.*np.pi*self.B*self.tau))*np.exp(-1./(2.*self.tau*self.B)*(1.-np.cos(dphi)))
            rho_pqc_out.write(str(dphi)+' '+str(rho_phi_pqc[i])+'\n')
            self.free_rho_pqc[i,:] = [dphi, rho_phi_pqc[i]]
        rho_pqc_out.close()
        
        # Overwrites the current rho_phi to match the pqc method
        self.rho_phi = rho_phi_pqc
        return self.free_rho_pqc.copy()
    
    def createRhoSOS(self):
        """
        Performs the calculation of the density matrix using the SOS method.
        Stores the values outlined below in the object and also returns the
        final density matrix. This function will overwrite the self.rho_phi 
        attribute of the object used during the Monte Carlo method.
        
        Returns
        -------
        self.rho_sos: numpy array
            Nx2 numpy array with the phi value in the 1st column and density 
            matrix values in the 2nd column
        self.rho_phi: numpy array
            Nx1 numpy array with density matrix values, used in the runMC function
        self.Z_sos: float
            Partition function calculated by the SOS method
        self.A_sos: float
            Helmholtz energy calculated by the SOS method
        self.E0_sos: float
            Ground state energy calculated by the SOS method
        self.E0_PIGS_sos: float
            Ground state energy calculated using PIGS and the SOS method
        """
        # 1 body Hamiltonian
        V = self.V0 * pot_matrix(2*self.m_max+1)
        H = V.copy()
        
        for m in range(self.Ngrid):
            m_value = -self.m_max+m
            H[m,m] = self.B * float(m_value**2) + self.V0 # constant potential term on diagonal
        evals, evecs = np.linalg.eigh(H)
        
        rho_mmp=np.zeros((self.Ngrid,self.Ngrid), float)
        
        Z_exact = 0.  #sum over state method
        for m in range(self.Ngrid):
            Z_exact += np.exp(-self.beta * evals[m])
            for mp in range(self.Ngrid):
                for n in range(self.Ngrid):
                    rho_mmp[m,mp]+=np.exp(-self.beta * evals[n]) * evecs[m,n] * evecs[mp,n]
        
        self.Z_sos = Z_exact
        self.A_sos = -(1./self.beta)*np.log(Z_exact)
        self.E0_sos = evals[0]
        if self.PIGS == True:
            Z_exact_pigs=rho_mmp[self.m_max,self.m_max]
            rho_dot_V_mmp=np.dot(rho_mmp,H)
            E0_pigs_sos=rho_dot_V_mmp[self.m_max,self.m_max]
            self.E0_PIGS_sos = E0_pigs_sos/Z_exact_pigs
            
        print('Z (sos) = ', self.Z_sos)
        print('A (sos) = ', self.A_sos)
        print('E0 (sos) =', self.E0_sos)
        if self.PIGS == True:
            print('E0 (pigs sos) =', self.E0_PIGS_sos)
        print(' ')
        
        # <phi|m><m|n> exp(-beta E n) <n|m'><m'|phi>
        rho_sos_out=open(os.path.join(self.path,'rho_sos'),'w')
        
        #built basis
        psi_m_phi=np.zeros((self.Ngrid,self.Ngrid),float)
        for i in range(self.Ngrid):
            for m in range(self.Ngrid):
                m_value =- self.m_max+m
                psi_m_phi[i,m] = np.cos(i * self.delta_phi*m_value) / np.sqrt(2.*np.pi)
                
        psi_phi=np.zeros((self.Ngrid,self.Ngrid),float)
        for i in range(self.Ngrid):
            for n in range(self.Ngrid):
                for m in range(self.Ngrid):
                    psi_phi[i,n] += evecs[m,n] * psi_m_phi[i,m]
        
        self.rho_sos = np.zeros((self.Ngrid, 3),float)
        for i in range(self.Ngrid):
            rho_exact=0.   
            for n in range(self.Ngrid):
                rho_exact+=np.exp(-self.beta*evals[n])*(psi_phi[i,n]**2)
            rho_exact/=(Z_exact)
            rho_sos_out.write(str(i * self.delta_phi)+ ' '+str(rho_exact)+' '+str(psi_phi[i,0]**2)+'\n')
            self.rho_sos[i,:] = [i * self.delta_phi, rho_exact, psi_phi[i,0]**2]
        rho_sos_out.close()
        self.rho_phi = self.rho_sos[:,1]
        return self.rho_sos.copy()
    
    def createRhoNMM(self):
        """
        Creates the density and free density matrices by the NMM method. This 
        function will overwrite the self.rho_phi attribute of the object used 
        during the Monte Carlo method.
        
        Returns
        -------
        self.potential: numpy array
            Nx2 array containing the potential for the rotors
        self.rho_nmm: numpy array
            Nx2 numpy array with the phi value in the 1st column and density 
            matrix values in the 2nd column
        self.free_rho_nmm: numpy array
            Nx2 numpy array with the phi value in the 1st column and free density 
            matrix values in the 2nd column
        self.rho_phi: numpy array
            Nx1 numpy array with density matrix values, used in the runMC function
        self.Z_nmm:float
            Partition function calculated by the NMM method
        self.E0_nmm: float
            Ground state energy calculated by the NMM method
        """
        rho_free=np.zeros((self.Ngrid,self.Ngrid),float)
        rho_potential=np.zeros(self.Ngrid,float)
        potential=np.zeros(self.Ngrid,float)
        for i in range(self.Ngrid):
            potential[i] = self.potFunc(float(i)*self.delta_phi,self.V0)
            rho_potential[i]=np.exp(-(self.tau/2.)*potential[i])
            for ip in range(self.Ngrid):
                integral=0.
                dphi=float(i-ip)*self.delta_phi
                for m in range(self.m_max):
                    integral+=np.exp(-1./(4.*self.tau*self.B)*(dphi+2.*np.pi*float(m))**2)
                for m in range(1,self.m_max):
                    integral+=np.exp(-1./(4.*self.tau*self.B)*(dphi+2.*np.pi*float(-m))**2)
                integral*=np.sqrt(1./(4.*np.pi*self.B*self.tau))
                rho_free[i,ip]=integral
        
        self.potential = np.zeros((self.Ngrid, 2),float)
        #output potential to a file
        potential_out=open(os.path.join(self.path,'V'),'w')
        for i in range(self.Ngrid):
                potential_out.write(str(float(i)*self.delta_phi)+' '+str(potential[i])+'\n')
                self.potential[i,:] = [float(i)*self.delta_phi, potential[i]]
        potential_out.close()
        
        # construct the high temperature density matrix
        rho_tau=np.zeros((self.Ngrid,self.Ngrid),float)
        for i1 in range(self.Ngrid):
                for i2 in range(self.Ngrid):
                        rho_tau[i1,i2]=rho_potential[i1]*rho_free[i1,i2]*rho_potential[i2]
        
        # form the density matrix via matrix multiplication        
        rho_beta=rho_tau.copy()
        for k in range(self.P-1):
                rho_beta=self.delta_phi*np.dot(rho_beta,rho_tau)
        
        E0_nmm=0.
        rho_dot_V=np.dot(rho_beta,potential)
        Z0=0. # pigs pseudo Z
        rho_nmm_out=open(os.path.join(self.path,'rho_nmm'),'w')
        Z_nmm=rho_beta.trace()*self.delta_phi # thermal Z
        Z_free_nmm = rho_free.trace()*self.delta_phi
        
        self.rho_nmm = np.zeros((self.Ngrid, 2), float)
        self.free_rho_nmm = np.zeros((self.Ngrid, 2), float)
        for i in range(self.Ngrid):
            E0_nmm += rho_dot_V[i]
            rho_nmm_out.write(str(i*self.delta_phi)+ ' '+str(rho_beta[i,i]/Z_nmm)+'\n')
            self.rho_nmm[i,:] = [i*self.delta_phi, rho_beta[i,i]/Z_nmm]
            self.free_rho_nmm[i,:] = [i*self.delta_phi, rho_free[i,i]/Z_free_nmm]
            for ip in range(self.Ngrid):
                Z0 += rho_beta[i,ip]
        self.rho_phi = self.rho_nmm[:,1]
        rho_nmm_out.close()
        E0_nmm/=Z0
        
        print('Z (tau) = ',Z_nmm)
        print('E0 (tau) = ',E0_nmm)
        self.Z_nmm = Z_nmm
        self.E0_nmm = E0_nmm
        E0_vs_tau_out=open('Evst','a')
        E0_vs_tau_out.write(str(self.tau)+' '+str(E0_nmm)+'\n')
        E0_vs_tau_out.close()
        print(' ')
    
    def createRhoVij(self):
        """
        Creates the rhoVij matrix for the nearest neighbour interactions of the system

        Returns
        -------
        self.rhoVij: numpy array
            NxN array containing the probability density based on the interaction
            potential between nearest neighbours
        self.rhoV: numpy array
            Nx1 array containing the probability density based on the on-site interactions

        """
        # potential rho
        self.rhoV = np.zeros((self.Ngrid),float)
        for i_new in range(self.Ngrid):
            self.rhoV[i_new] = np.exp(-self.tau * (self.potFunc(float(i_new)*self.delta_phi, self.V0)))
        # rho pair
        self.rhoVij = np.zeros((self.Ngrid,self.Ngrid),float)
        for i in range(self.Ngrid):
            for j in range(self.Ngrid):
                self.rhoVij[i,j] = np.exp(-self.tau * (Vij(i*self.delta_phi, j*self.delta_phi, self.g)))
    
    def runMC(self, averagePotential = True, averageEnergy = True, orientationalCorrelations = True):
        """
        Performs the monte carlo integration to simulate the system.
        
        Parameters
        ----------
        averagePotential : bool, optional
            Enables the tracking and calculation of the average potential. 
            The default is True.
        averageEnergy : bool, optional
            Enables the tracking and calculation of the average energy. 
            The default is True.
        orientationalCorrelations : bool, optional
            Enables the tracking and calculation of the orientational correllations. 
            The default is True.

        Returns
        -------
        self.V_MC: float
            The resultant average system potential.
        self.E_MC: float
            The resultant average system energy.
        self.eiej_MC: float
            The resultant average system orientational correlation.
        self.histo: numpy array
            Nx6 array containing:
                1st column: The angle phi
                2nd column: PIMC Histogram
                3rd column: Middle bead histogram
                4th column: Left bead histogram
                5th column: Right bead histogram
                6th column: Initial overall histogram
        
        Outputs
        -------
        histo_A_P_N: Nx2 txt file
            A saved version of the self.histo outlined above.
        traj_A: Nx3 dat file
            The phi values of the left, right and middle beads in columns 1, 2 
            and 3 respectively, with each row corresponding to a specific rotor
            during a specific MC step.

        """
        histo_L=np.zeros(self.Ngrid,float)
        histo_R=np.zeros(self.Ngrid,float)
        histo_middle=np.zeros(self.Ngrid,float)
        histo_pimc=np.zeros(self.Ngrid,float)
        histo_initial=np.zeros(self.Ngrid,float)
        
        if not hasattr(self, 'rho_phi'):
            self.createFreeRhoMarx()
            
        if not hasattr(self, 'rhoVij'):
            self.createRhoVij()
        
        p_dist=gen_prob_dist(self.Ngrid, self.rho_phi)
        p_dist_end = gen_prob_dist_end(self.Ngrid, self.rho_phi) if self.PIGS == True else None
        
        path_phi=np.zeros((self.N,self.P),int) ## i  N => number of beads
        for i in range(self.N):
            for p in range(self.P): # set path at potential minimum
                # Initial conditions
                # All rotors have angle of pi
                #path_phi[i,p]=int(self.Ngrid/2)
                # All rotors have angle of 0
                #path_phi[i,p]=0
                # Rotors have random angles
                path_phi[i,p]=np.random.randint(self.Ngrid)
                histo_initial[path_phi[i,p]]+=1.
            
        traj_out=open(os.path.join(self.path,'traj_A.dat'),'w')
        
        # recommanded numpy random number initialization
        rng = default_rng()
        
        print('start MC')
        
        # Initializing the observables to be tracked
        V_average = 0. if averagePotential == True else None
        E_average = 0. if averageEnergy == True else None
        eiej_corr = 0. if orientationalCorrelations == True else None
            
        P_middle = int((self.P-1)/2)
        prob_full=np.zeros(self.Ngrid,float)
        
        for n in range(self.MC_steps):
            # Total and potential energy estimators
            V_total=0.
            E_total=0.
            eiej=0.
            for i in range(self.N):
                for p in range(self.P): 
                    p_minus=p-1
                    p_plus=p+1
                    if (p_minus<0):
                        p_minus+=self.P
                    if (p_plus>=self.P):
                        p_plus-=self.P  
                    
                    # kinetic action, links between beads
                    if self.PIGS==True:
                        # This uses a split open path of beads for PIGS
                        if p==0:
                            # Special case for leftmost bead
                            for ip in range(self.Ngrid):
                                prob_full[ip]=p_dist_end[ip,path_phi[i,p_plus]]
                        elif p==(self.P-1):
                            # Special case for rightmost bead
                            for ip in range(self.Ngrid):
                                prob_full[ip]=p_dist_end[path_phi[i,p_minus],ip]
                        elif (p!=0 and p!= (self.P-1)):
                            # All other beads
                            for ip in range(self.Ngrid):
                                prob_full[ip]=p_dist[path_phi[i,p_minus],ip,path_phi[i,p_plus]]
                    else:
                        # Regular kinetic interactions between beads, periodic conditions
                        for ip in range(self.Ngrid):
                            prob_full[ip]=p_dist[path_phi[i,p_minus],ip,path_phi[i,p_plus]]                         
                    
                    # Local on site interaction with the potential field
                    if self.V0 != 0.:
                        for ip in range(self.Ngrid):
                           prob_full[ip]*=self.rhoV[ip]
        
                    # NN interactions and PBC(periodic boundary conditions)
                    if (i<(self.N-1)):
                        # Interaction with the rotor to the right
                        for ir in range(len(prob_full)):
                            prob_full[ir]*=self.rhoVij[ir,path_phi[i+1,p]]
                    if (i>0):
                        # Interaction with the rotor to the left
                        for ir in range(len(prob_full)):
                            prob_full[ir]*=self.rhoVij[ir,path_phi[i-1,p]]
                    if (i==0):
                        # Periodic BC for the leftmost rotor
                        for ir in range(len(prob_full)):
                            prob_full[ir]*=self.rhoVij[ir,path_phi[self.N-1,p]]
                    if (i==(self.N-1)):
                        # Periodic BC for the rightmost rotor
                        for ir in range(len(prob_full)):
                            prob_full[ir]*=self.rhoVij[ir,path_phi[0,p]]
        
                    # Normalize
                    norm_pro=0.
                    for ir in range(len(prob_full)):
                        norm_pro+=prob_full[ir]
                    for ir in range(len(prob_full)):
                        prob_full[ir]/=norm_pro
                    index=rng.choice(self.Ngrid,1, p=prob_full)
                    # Rejection free sampling
                    path_phi[i,p] = index
        
                    histo_pimc[path_phi[i,p]]+=1.
                
                if (n%self.Nskip==0):
                    traj_out.write(str(path_phi[i,0]*self.delta_phi)+' ')
                    traj_out.write(str(path_phi[i,self.P-1]*self.delta_phi)+' ')
                    traj_out.write(str(path_phi[i,P_middle]*self.delta_phi)+' ') #middle bead
                    traj_out.write('\n')
                    
                histo_L[path_phi[i,0]]+=1.
                histo_R[path_phi[i,self.P-1]]+=1.
                histo_middle[path_phi[i,P_middle]]+=1.
                
                if n >= self.Nequilibrate and averagePotential == True:
                    # External field
                    V_total += self.potFunc(float(path_phi[i,P_middle])*self.delta_phi,self.V0)
                    # Nearest neighbour interactions
                    if (i<(self.N-1)):
                        # Double check and think more about this, we only look at right neighbour?
                        V_total += Vij(path_phi[i,P_middle]*self.delta_phi, path_phi[i+1,P_middle]*self.delta_phi, self.g)
                    elif (i==(self.N-1)):
                        V_total += Vij(path_phi[i,P_middle]*self.delta_phi, path_phi[0,P_middle]*self.delta_phi, self.g)
                    
                if n >= self.Nequilibrate and averageEnergy == True:
                    # External field
                    E_total += self.potFunc(float(path_phi[i,0])*self.delta_phi,self.V0)
                    E_total += self.potFunc(float(path_phi[i,self.P-1])*self.delta_phi,self.V0)
                    # Nearest neighbour interactions
                    if (i<(self.N-1)):
                        # Double check and think more about this, we only look at right neighbour?
                        E_total += Vij(path_phi[i,0]*self.delta_phi, path_phi[i+1,0]*self.delta_phi, self.g)
                        E_total += Vij(path_phi[i,self.P-1]*self.delta_phi, path_phi[i+1,self.P-1]*self.delta_phi, self.g)
                    elif (i==(self.N-1)):
                        E_total += Vij(path_phi[i,0]*self.delta_phi, path_phi[0,0]*self.delta_phi, self.g)
                        E_total += Vij(path_phi[i,self.P-1]*self.delta_phi, path_phi[0,self.P-1]*self.delta_phi, self.g)
                
                if n >= self.Nequilibrate and orientationalCorrelations == True:
                    if (i<(self.N-1)):
                        eiej += calculateOrientationalCorrelations(path_phi[i,P_middle]*self.delta_phi, path_phi[i+1,P_middle]*self.delta_phi)
                    if (i==(self.N-1)):
                        eiej += calculateOrientationalCorrelations(path_phi[i,P_middle]*self.delta_phi, path_phi[0,P_middle]*self.delta_phi)
                        
            if n >= self.Nequilibrate and averagePotential == True:
                V_average += V_total
            if n >= self.Nequilibrate and averageEnergy == True:
                E_average += E_total
            if n >= self.Nequilibrate and orientationalCorrelations == True:
                eiej_corr += eiej
                
        traj_out.close()
        
        if averagePotential == True:
            print('<V> = ',V_average/(self.MC_steps-self.Nequilibrate)/self.N)
            self.V_MC = V_average/(self.MC_steps-self.Nequilibrate)/self.N
        if averageEnergy == True:
            print('<E> = ',E_average/(self.MC_steps-self.Nequilibrate)/2./self.N)
            self.E_MC = E_average/(self.MC_steps-self.Nequilibrate)/2./self.N
        if orientationalCorrelations == True:
            print('<ei.ej> = ',eiej_corr/(self.MC_steps-self.Nequilibrate)/self.N)
            self.eiej_MC = eiej_corr/(self.MC_steps-self.Nequilibrate)/self.N

        self.histo = np.zeros((self.Ngrid,6))
        histo_out=open(os.path.join(self.path,'histo_A_P'+str(self.P)+'_N'+str(self.N)),'w')
        for i in range(self.Ngrid):
          histo_out.write(str(i*self.delta_phi) + ' ' +
                          str(histo_pimc[i]/(self.MC_steps*self.N*self.P)/self.delta_phi) + ' ' +
                          str(histo_middle[i]/(self.MC_steps*self.N)/self.delta_phi) + ' ' +
                          str(histo_L[i]/(self.MC_steps*self.N)/self.delta_phi) + ' ' +
                          str(histo_R[i]/(self.MC_steps*self.N)/self.delta_phi) + ' ' +
                          str(histo_initial[i]/(self.N*self.P)/self.delta_phi)+'\n')
          self.histo[i,:] = [i*self.delta_phi, 
                            histo_pimc[i]/(self.MC_steps*self.N*self.P)/self.delta_phi,
                            histo_middle[i]/(self.MC_steps*self.N)/self.delta_phi,
                            histo_L[i]/(self.MC_steps*self.N)/self.delta_phi,
                            histo_R[i]/(self.MC_steps*self.N)/self.delta_phi,
                            histo_initial[i]/(self.N*self.P)/self.delta_phi]
          
        histo_out.close()                
    
    def runMCReplica(self):
        """
        Performs the monte carlo integration to simulate the system with entanglement considered. This employs
        the replica trick and the extended ensemble to do so.
        
        ############################# TO ADD: Ratio trick

        Returns
        -------
        self.S2: float
            The resultant second Renyi entropy
        self.histo: numpy array
            Nx6 array containing:
                1st column: The angle phi
                2nd column: PIMC Histogram
                3rd column: Middle bead histogram
                4th column: Left bead histogram
                5th column: Right bead histogram
                6th column: Initial overall histogram
        
        Outputs
        -------
        histo_A_P_N: Nx2 txt file
            A saved version of the self.histo outlined above.
        traj_A: Nx3 dat file
            The phi values of the left, right and middle beads in columns 1, 2 
            and 3 respectively, with each row corresponding to a specific rotor
            during a specific MC step.

        """
        histo_L=np.zeros(self.Ngrid,float)
        histo_R=np.zeros(self.Ngrid,float)
        histo_middle=np.zeros(self.Ngrid,float)
        histo_pimc=np.zeros(self.Ngrid,float)
        histo_initial=np.zeros(self.Ngrid,float)
        
        #############################################################
        # Add check for even number of rotors, if the ratio trick isn't being employed
        
        if not hasattr(self, 'rho_phi'):
            self.createFreeRhoMarx()
            
        if not hasattr(self, 'rhoVij'):
            self.createRhoVij()
        
        if not self.PIGS:
            raise Exception("PIGS must be enabled to run runMCReplica, please create a choiceMC object with this enabled")
        
        p_dist=gen_prob_dist(self.Ngrid, self.rho_phi)
        p_dist_end=gen_prob_dist_end(self.Ngrid, self.rho_phi)
        
        path_phi=np.zeros((self.N,self.P),int) ## i  N => number of beads
        for i in range(self.N):
            for p in range(self.P): # set path at potential minimum
                # Initial conditions
                # All rotors have angle of pi
                #path_phi[i,p]=int(self.Ngrid/2)
                # All rotors have angle of 0
                #path_phi[i,p]=0
                # Rotors have random angles
                path_phi[i,p]=np.random.randint(self.Ngrid)
                histo_initial[path_phi[i,p]]+=1.
            
        traj_out=open(os.path.join(self.path,'traj_A.dat'),'w')
        
        # recommanded numpy random number initialization
        rng = default_rng()
        
        print('start MC')
        
        # Counter for the number of MC steps in the swapped and unswapped configuration
        N_swapped = 0
        N_unswapped = 0
        swapped = False
        
        # Index of the partition in the chain
        N_partition = self.N // 2
        
        # Creating a replica of the original path
        path_phi_replica = path_phi.copy()
        
        # The indeces of the beads in the middle and to the left of the middle of the chain
        P_middle = int((self.P-1)/2)
        P_midLeft = P_middle - 1
        
        prob_full=np.zeros(self.Ngrid,float)
        prob_full_replica=np.zeros(self.Ngrid,float)
        
        for n in range(self.MC_steps):
            # Entanglement estimators
            N_swapped += swapped
            N_unswapped += (not swapped)
            
            # As the rotors are looped through, the only ones that can partake in the swapped/unswapped
            # configuration are the rotors in the "A" partition
            for i in range(self.N):
                for p in range(self.P): 
                    p_minus=p-1
                    p_plus=p+1
                    if (p_minus<0):
                        p_minus+=self.P
                    if (p_plus>=self.P):
                        p_plus-=self.P  
                    
                    # kinetic action
                    if p==0:
                        # Special case on the left end of the chain
                        for ip in range(self.Ngrid):
                            prob_full[ip]=p_dist_end[ip,path_phi[i,p_plus]]
                            prob_full_replica[ip]=p_dist_end[ip,path_phi_replica[i,p_plus]]
                    elif p==(self.P-1):
                        # Special case on the right end of the chain
                        for ip in range(self.Ngrid):
                            prob_full[ip]=p_dist_end[path_phi[i,p_minus],ip]
                            prob_full_replica[ip]=p_dist_end[path_phi_replica[i,p_minus],ip]
                    elif (p==P_midLeft) and swapped and i < N_partition:
                        # Special case for extended ensemble for the bead to the left of the middle
                        # This only applies to the "A" partition
                        for ip in range(self.Ngrid):
                            prob_full[ip]=p_dist[path_phi[i,p_minus],ip,path_phi_replica[i,p_plus]]
                            prob_full_replica[ip]=p_dist[path_phi_replica[i,p_minus],ip,path_phi[i,p_plus]]
                    elif (p==P_middle) and swapped and i < N_partition:
                        # Special case for extended ensemble on the middle bead
                        # This only applies to the "A" partition
                        for ip in range(self.Ngrid):
                            prob_full[ip]=p_dist[path_phi_replica[i,p_minus],ip,path_phi[i,p_plus]]
                            prob_full_replica[ip]=p_dist[path_phi[i,p_minus],ip,path_phi_replica[i,p_plus]]
                    elif (p!=0 and p!=(self.P-1)):
                        # Interactions for non-swapped and non-end beads
                        for ip in range(self.Ngrid):
                            prob_full[ip]=p_dist[path_phi[i,p_minus],ip,path_phi[i,p_plus]]
                            prob_full_replica[ip]=p_dist[path_phi_replica[i,p_minus],ip,path_phi_replica[i,p_plus]]                         
                    
                    # Local on site interaction with the external potential field
                    if self.V0 != 0.:
                        for ip in range(self.Ngrid):
                           prob_full[ip]*=self.rhoV[ip]
                           prob_full_replica[ip]*=self.rhoV[ip]
                       
                    # NN interactions and PBC(periodic boundary conditions)
                    
                    if (i<(self.N-1)):
                        # Interaction with right neighbour
                        if (p==P_middle) and swapped and i == (N_partition-1):
                            # Swaps the right interaction for the middle bead of the rotor at the partition on the A side
                            for ir in range(len(prob_full)):
                                prob_full[ir]*=self.rhoVij[ir,path_phi_replica[i+1,p]]
                                prob_full_replica[ir]*=self.rhoVij[ir,path_phi[i+1,p]]
                        else:
                            for ir in range(len(prob_full)):
                                prob_full[ir]*=self.rhoVij[ir,path_phi[i+1,p]]
                                prob_full_replica[ir]*=self.rhoVij[ir,path_phi_replica[i+1,p]]
                    if (i>0):
                        # Interaction with left neighbour
                        if (p==P_middle) and swapped and i == N_partition:
                            # Swaps the left interaction for the middle bead of the rotor at the partition on the B side
                            for ir in range(len(prob_full)):
                                prob_full[ir]*=self.rhoVij[ir,path_phi_replica[i-1,p]]
                                prob_full_replica[ir]*=self.rhoVij[ir,path_phi[i-1,p]]
                        else:
                            for ir in range(len(prob_full)):
                                prob_full[ir]*=self.rhoVij[ir,path_phi[i-1,p]]
                                prob_full_replica[ir]*=self.rhoVij[ir,path_phi_replica[i-1,p]] 
                    if (i==0):
                        # Periodic BC for the leftmost rotor
                        if (p==P_middle) and swapped:
                            # Swaps the left interaction for the middle bead of the leftmost rotor
                            for ir in range(len(prob_full)):
                                prob_full[ir]*=self.rhoVij[ir,path_phi_replica[self.N-1,p]]
                                prob_full_replica[ir]*=self.rhoVij[ir,path_phi[self.N-1,p]]
                        else:
                            for ir in range(len(prob_full)):
                                prob_full[ir]*=self.rhoVij[ir,path_phi[self.N-1,p]]
                                prob_full_replica[ir]*=self.rhoVij[ir,path_phi_replica[self.N-1,p]]
                    if (i==(self.N-1)):
                        # Periodic BC for the rightmost rotor
                        if (p==P_middle) and swapped:
                            # Swaps the left interaction for the middle bead of the rightmost rotor
                            for ir in range(len(prob_full)):
                                prob_full[ir]*=self.rhoVij[ir,path_phi_replica[0,p]]
                                prob_full_replica[ir]*=self.rhoVij[ir,path_phi[0,p]]
                        else:
                            for ir in range(len(prob_full)):
                                prob_full[ir]*=self.rhoVij[ir,path_phi[0,p]]
                                prob_full_replica[ir]*=self.rhoVij[ir,path_phi_replica[0,p]]
        
                    #normalize
                    norm_pro=0.
                    norm_pro_replica=0.
                    for ir in range(len(prob_full)):
                        norm_pro+=prob_full[ir]
                        norm_pro_replica+=prob_full_replica[ir]
                    for ir in range(len(prob_full)):
                        prob_full[ir]/=norm_pro
                        prob_full_replica[ir]/=norm_pro_replica
                    index=rng.choice(self.Ngrid,1, p=prob_full)
                    index_replica=rng.choice(self.Ngrid,1, p=prob_full_replica)
                    
                    # Rejection free sampling
                    path_phi[i,p] = index
                    path_phi_replica[i,p] = index_replica
                    
                    histo_pimc[path_phi[i,p]]+=1.
                    
                if (n%self.Nskip==0):
                    traj_out.write(str(path_phi[i,0]*self.delta_phi)+' ')
                    traj_out.write(str(path_phi[i,self.P-1]*self.delta_phi)+' ')
                    traj_out.write(str(path_phi[i,P_middle]*self.delta_phi)+' ') #middle bead
                    traj_out.write('\n')
                    
                histo_L[path_phi[i,0]]+=1.
                histo_R[path_phi[i,self.P-1]]+=1.
                histo_middle[path_phi[i,P_middle]]+=1.
            
            # Metropolis critereon
            
            ######################################################
            # This needs to be double checked for double counting
            
            # The interaction with the external potential field is ignored, as this will
            # be the same for both the swapped and unswapped ensembles
            # Any unchanged interactions (kinetic or potential) between the swapped and
            # unswapped configurations are not calculated
            
            rhoUnswapped = 1.
            rhoSwapped = 1.
            
            # Kinetic contribution from the swapped interactions in the A partition
            ##########################################################
            # This needs to be double checked, should this be multiplying the middle and midleft bead?
            for i in range(N_partition):
                # Middle bead
                rhoUnswapped *= p_dist[path_phi[i,P_midLeft],path_phi[i,P_middle],path_phi[i,P_middle+1]]
                rhoUnswapped *= p_dist[path_phi_replica[i,P_midLeft],path_phi_replica[i,P_middle],path_phi_replica[i,P_middle+1]]
                rhoSwapped *= p_dist[path_phi_replica[i,P_midLeft],path_phi[i,P_middle],path_phi[i,P_middle+1]]
                rhoSwapped *= p_dist[path_phi[i,P_midLeft],path_phi_replica[i,P_middle],path_phi_replica[i,P_middle+1]]
                # Bead to the left of the middle
                rhoUnswapped *= p_dist[path_phi[i,P_midLeft-1],path_phi[i,P_midLeft],path_phi[i,P_middle]]
                rhoUnswapped *= p_dist[path_phi_replica[i,P_midLeft-1],path_phi_replica[i,P_midLeft],path_phi_replica[i,P_middle]]
                rhoSwapped *= p_dist[path_phi[i,P_midLeft-1],path_phi[i,P_midLeft],path_phi_replica[i,P_middle]]
                rhoSwapped *= p_dist[path_phi_replica[i,P_midLeft-1],path_phi_replica[i,P_midLeft],path_phi[i,P_middle]]
                
            # Potential contribution, this only impacts the middle bead
            # Swapped interactions for the periodic BCs
            rhoUnswapped *= self.rhoVij[path_phi[0,P_middle],path_phi[self.N-1,P_middle]]
            rhoUnswapped *= self.rhoVij[path_phi_replica[0,P_middle],path_phi_replica[self.N-1,P_middle]]
            rhoSwapped *= self.rhoVij[path_phi_replica[0,P_middle],path_phi[self.N-1,P_middle]]
            rhoSwapped *= self.rhoVij[path_phi[0,P_middle],path_phi_replica[self.N-1,P_middle]]
            # Swapped interactions at the partition between A and B
            rhoUnswapped *= self.rhoVij[path_phi[N_partition-1,P_middle],path_phi[N_partition,P_middle]]
            rhoUnswapped *= self.rhoVij[path_phi_replica[N_partition-1,P_middle],path_phi_replica[N_partition,P_middle]]
            rhoSwapped *= self.rhoVij[path_phi_replica[N_partition-1,P_middle],path_phi[N_partition,P_middle]]
            rhoSwapped *= self.rhoVij[path_phi[N_partition-1,P_middle],path_phi_replica[N_partition,P_middle]]
                    
            # Determing if we should be sampling the swapped or unswapped distribution
            if swapped:
                ratio = rhoUnswapped/rhoSwapped
                if ratio > 1:
                    swapped = False
                elif ratio > np.random.uniform():
                    swapped = False
            elif not swapped:
                ratio = rhoSwapped/rhoUnswapped
                if ratio > 1:
                    swapped = True
                elif ratio > np.random.uniform():
                    swapped = True  
            
        traj_out.close()
        
        
        self.S2 = -np.log(N_swapped/N_unswapped)
        print('S2 = ', str(self.S2))

        self.histo = np.zeros((self.Ngrid,6))
        histo_out=open(os.path.join(self.path,'histo_A_P'+str(self.P)+'_N'+str(self.N)),'w')
        for i in range(self.Ngrid):
            histo_out.write(str(i*self.delta_phi) + ' ' +
                            str(histo_pimc[i]/(self.MC_steps*self.N*self.P)/self.delta_phi) + ' ' +
                            str(histo_middle[i]/(self.MC_steps*self.N)/self.delta_phi) + ' ' +
                            str(histo_L[i]/(self.MC_steps*self.N)/self.delta_phi) + ' ' +
                            str(histo_R[i]/(self.MC_steps*self.N)/self.delta_phi) + ' ' +
                            str(histo_initial[i]/(self.N*self.P)/self.delta_phi)+'\n')
            self.histo[i,:] = [i*self.delta_phi, 
                              histo_pimc[i]/(self.MC_steps*self.N*self.P)/self.delta_phi,
                              histo_middle[i]/(self.MC_steps*self.N)/self.delta_phi,
                              histo_L[i]/(self.MC_steps*self.N)/self.delta_phi,
                              histo_R[i]/(self.MC_steps*self.N)/self.delta_phi,
                              histo_initial[i]/(self.N*self.P)/self.delta_phi]
          
        histo_out.close()
        
    def plotRho(self, rhoList):
        """
        Plots and saves the specified density matrices.
        
        Parameters
        ----------
        rhoList : List
            A list containing the names of the density matrices to be plotted.
            Allowed choices: rho_sos, free_rho_sos, free_rho_pqc, free_rho_marx, rho_nmm, free_rho_nmm.
        
        Outputs
        -------
        DensityMatrices.png: png file
            Image of the plot of the density matrices.
        """
        
        fig, (ax, ax_table) = plt.subplots(1, 2, gridspec_kw={"width_ratios": [5, 1]}, figsize=(8,5))
        for rho in rhoList:
            if hasattr(self, rho):
                exec("ax.plot(self." + rho + "[:,0],self." + rho + "[:,1], label='" + rho + "')")
            else:
                raise Warning(rho + " has not been generated, make sure that the name is correct and that it has been generated")
            
        ax.set_xlabel('Phi')
        ax.set_ylabel('Density Matrix')
        ax.set_title('Density Matrices')
        ax.legend()
        text = np.array([["Temperature", str(self.T) + " K"],
                        ["Number of\nGrid Points", str(self.Ngrid)],
                        ["Number of\nBeads", str(self.P)],
                        ["Number of\nMC Steps", str(self.MC_steps)],
                        ["Number of\nRotors", str(self.N)],
                        ["Interaction\nStrength", str(round(self.g,3))],
                        ["Rotational\nConstant", str(self.B)],
                        ["Skip Steps", str(self.Nskip)],
                        ["Equilibration\nSteps", str(self.Nequilibrate)],
                        ["PIGS", str(self.PIGS)],
                        ["Potential (V0)", str(self.V0)]])
        ax_table.axis("off")
        table = ax_table.table(cellText=text, 
                               loc='center',
                               colWidths=[1.5, 1],
                               cellLoc='center')
        table.set_fontsize(9)
        table.scale(1,2.0)
        fig.tight_layout()
        fig.savefig(os.path.join(self.path, "DensityMatrices.png"))
    
    def plotHisto(self, *args):
        """
        Plots and saves the specified histograms.
        
        Parameters
        ----------
        *args : String
            Keywords to determine which histograms to plot, for the overall 
            PIMC, left, middle and right beads and the initial conditions.
            Allowed choices: PIMC, left, middle, right, initial.
        
        Outputs
        -------
        labels_Histograms.png: png file
            Image of the plots of the histograms specified.
        """
        if not hasattr(self, "histo"):
            raise Warning("The histograms do not exist, please execute ChoiceMC.runMC()")
            return
        
        histo_dict = {"PIMC" : 1, 
                      "middle" : 2, 
                      "left" : 3, 
                      "right" : 4, 
                      "initial" : 5}
        fig_label = ""
        fig, (ax, ax_table) = plt.subplots(1, 2, gridspec_kw={"width_ratios": [5, 1]}, figsize=(8,5))
        for hist in args:
            exec("ax.plot(self.histo[:,0],self.histo[:," + str(histo_dict[hist]) + "], label='" + hist + "')")
            fig_label += hist + "_"   
        ax.set_xlabel('Phi')
        ax.set_ylabel('Probability Density')
        ax.set_title('Histograms')
        ax.legend()
        text = np.array([["Temperature", str(self.T) + " K"],
                        ["Number of\nGrid Points", str(self.Ngrid)],
                        ["Number of\nBeads", str(self.P)],
                        ["Number of\nMC Steps", str(self.MC_steps)],
                        ["Number of\nRotors", str(self.N)],
                        ["Interaction\nStrength", str(round(self.g,3))],
                        ["Rotational\nConstant", str(self.B)],
                        ["Skip Steps", str(self.Nskip)],
                        ["Equilibration\nSteps", str(self.Nequilibrate)],
                        ["PIGS", str(self.PIGS)],
                        ["Potential (V0)", str(self.V0)]])
        ax_table.axis("off")
        table = ax_table.table(cellText=text, 
                               loc='center',
                               colWidths=[1.5, 1],
                               cellLoc='center')
        table.set_fontsize(9)
        table.scale(1,2.0)
        fig.tight_layout()
        fig.savefig(os.path.join(self.path, fig_label + "Histograms.png"))