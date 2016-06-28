import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
from pylab import *
from decimal import Decimal

# Uncomment the following line if you wish to see the entire matrix C in the output
np.set_printoptions(threshold=np.nan)


# Start to gather the single particle data
spdata = np.genfromtxt ('spdata.dat')

ind = spdata[:,2]
n = spdata[:,3]
l = spdata[:,4]
j = spdata[:,5]
mj = spdata[:,6]
tz = spdata[:,7]


# Decide whether you want to analyze protons or neutrons
IsospinChoice = raw_input("Do you wish to consider protons, neutrons, or both? \n")

if IsospinChoice in ['n', 'N', 'neutron', 'Neutron', 'NEUTRON', 'neutrons', 'Neutrons', 'NEUTRONS']:
    twoiso = 1
    ind_n = ind[np.where(tz == twoiso)]
    numParticles = 8
    print 'Neutrons selected'
elif IsospinChoice in ['p', 'P', 'proton', 'Proton', 'PROTON', 'protons', 'Protons', 'PROTONS']:
    twoiso = -1
    ind_n = ind[np.where(tz == twoiso)]
    numParticles = 8
    print 'Protons selected'
else:
    ind_n = ind
    numParticles = 16
    print 'Both selected'

numStates = size(ind_n)


# Decide whether you want to include the two-body interaction or not
includetwobody = raw_input("Do you wish to include the two-body interaction? Y/N \n")

if includetwobody in ['y', 'Y', 'yes', 'Yes', 'YES']:
    useV = 1
else:
    useV = 0

# If you do decide to include two-body, this will build a 4-dimensional array
# in which to store your matrix elements
if useV == 1 and numParticles == 8:
    twobody = np.genfromtxt ('twobody.dat')
    V0 = np.zeros((numStates,numStates,numStates,numStates))
    for alpha in ind_n: #np.linspace(68, 68, num=1):#
        alphas = twobody[np.where(twobody[:,0] == alpha)]
        alpha_ind = np.where(ind_n==alpha)[0]
        for beta in ind_n: #np.linspace(49, 49, num=1):#
            betas = alphas[np.where(alphas[:,2] == beta)]#and
            beta_ind = np.where(ind_n==beta)[0]
            for gamma in ind_n:
                gammas = betas[np.where(betas[:,1] == gamma)]
                gamma_ind = np.where(ind_n==gamma)[0]
                for delta in ind_n:
                    V_agbd = gammas[np.where(gammas[:,3] == delta)]
                    delta_ind = np.where(ind_n==delta)[0]
                    if size(V_agbd) == 5:
                        V0[alpha_ind[0],gamma_ind[0],beta_ind[0],delta_ind[0]] = Decimal(V_agbd.item(0,4))
                    elif size(V_agbd) != 0:
                        print "Something is wrong with your filter"
elif useV == 1 and numParticles == 16:
    V0 = np.zeros((numStates,numStates,numStates,numStates))
    twobody = np.genfromtxt ('twobody.dat')
    with open("twobody.dat", "r") as infile:
        for line in infile:
            number = line.split()
            a = int(number[0]) - 1
            b = int(number[1]) - 1
            c = int(number[2]) - 1
            d = int(number[3]) - 1
            #print a, b, c, d, float(l[4])
            V0[a][b][c][d] = Decimal(number[4])
else:
    twobody = np.zeros((5,5))
    V0 = np.zeros((numStates,numStates,numStates,numStates))

h = np.zeros((numStates,numStates))

E = np.zeros((numStates, 1))
EPrev = np.zeros((numStates, 1))

C0 = np.identity(numStates)
C = C0
Cstar = np.transpose(np.conjugate(C))
rho = np.zeros((numStates,numStates))
for gamma in ind_n:
    gamma_int = np.where(ind_n==gamma)[0]
    gamma_ind = gamma_int[0]
    for delta in ind_n:
        delta_int = np.where(ind_n==delta)[0]
        delta_ind = delta_int[0]
        for i in np.linspace(0, numParticles-1, num=numParticles):
            rho[gamma_ind,delta_ind] = Decimal(rho[gamma_ind,delta_ind]+ Cstar[i,gamma_ind]*C[delta_ind,i])

if useV==1:
    h0 = np.zeros((numStates,numStates))
    for alpha in ind_n: #np.linspace(68, 68, num=1):#
        alphas = twobody[np.where(twobody[:,0] == alpha)]
        #print alphas
        for beta in ind_n:
            betas = alphas[np.where(alphas[:,2] == beta)]#and
            #print betas
            h0[np.where(ind_n==alpha),np.where(ind_n==alpha)] = 10*(
                2*n[np.where(ind==alpha)]+l[np.where(ind==alpha)]+1.5)
else:
    h0 = 10*(2*n[np.where(tz == twoiso)]+l[np.where(tz == twoiso)]+1.5)*np.identity(numStates)
        
EigValues, EigVectors = LA.eigh(h0)
EPrev = EigValues # Initialize the energy to the HO energies
print 'Energy eigenvalues initialized'
#print rho

nMax = 20

for nIter in np.linspace(1, nMax, num=nMax):
    V = np.zeros((numStates,numStates))
    for alpha in ind_n: #np.linspace(68, 68, num=1):#
        alpha_int = np.where(ind_n==alpha)[0]
        alpha_ind = alpha_int[0]
        for beta in ind_n: #np.linspace(49, 49, num=1):#
            beta_int = np.where(ind_n==beta)[0]
            beta_ind = beta_int[0]
            for gamma in ind_n:
                gamma_int = np.where(ind_n==gamma)[0]
                gamma_ind = gamma_int[0]
                for delta in ind_n:
                    delta_int = np.where(ind_n==delta)[0]
                    delta_ind = delta_int[0]
                    V0_agbd = V0[alpha_ind,gamma_ind,beta_ind,delta_ind]
                    rho_gd = rho[gamma_ind,delta_ind]
                    V[alpha_ind,beta_ind] = Decimal( V[alpha_ind,beta_ind] + rho_gd*V0_agbd )
                                       

    h = h0 + V
    #print h
    EigValues, EigVectors = LA.eigh(h)
    
    # sort eigenvectors and eigenvalues
    permute = EigValues.argsort()
    E = EigValues#[permute]
    C = EigVectors#[:,permute]
    Cstar = np.transpose(np.conjugate(C))
    rho = np.zeros((numStates,numStates))
    for gamma in ind_n:
        gamma_int = np.where(ind_n==gamma)[0]
        gamma_ind = gamma_int[0]
        for delta in ind_n:
            delta_int = np.where(ind_n==delta)[0]
            delta_ind = delta_int[0]
            for i in np.linspace(0, numParticles-1, num=numParticles):
                rho[gamma_ind,delta_ind] = Decimal(rho[gamma_ind,delta_ind]+ Cstar[i,gamma_ind]*C[delta_ind,i])
    diff = E-EPrev
    print 'After iteration', nIter, ', E='
    print E
    print 'diff=', np.linalg.norm(diff)
    if np.linalg.norm(diff) < 10**(-5):
        print "The calculation converged after ", nIter, "iterations."
        print "The converged energy eigenvalues are"
        print "\n", EigValues[permute], "\n"
#        print "\n and the corresponding eigenvectors are"
#        print "\n", EigVectors[:,permute], "\n"
        break
    EPrev = E
else:
    print "The calculation did not converge after", nMax, "iterations."
    print "The non-converged energy eigenvalues are"
    print "\n", EigValues[permute]
    print "\n and the corresponding eigenvectors are"
    print "\n", EigVectors[:,permute], "\n"

