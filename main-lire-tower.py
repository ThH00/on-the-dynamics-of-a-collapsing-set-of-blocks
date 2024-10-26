import numpy as np
import math
import time
import os
import argparse

start_time = time.time()

n = 5
mu_val = .3
output_path = "G:\\My Drive\\Research\\08-Stability of Stacked Objects\\Code-2024-02\\stacked-blocks-2d\\outputs\\max-overhang-new"

# nondimensionalization parameters
l_nd = 1                    # m, length nondimensionalization paramter
m_nd = 1                    # kg, mass nondimensionalization parameter
a_nd = 9.81                 # m/(s**2), acceleration nondimensionalization parameter
t_nd = np.sqrt(l_nd/a_nd)   # s, time nondimensionalization parameter

# quantites following this line are nondimensional

ndof = 3*(n-1)              # total number of degress of freedom
gr = 9.81/a_nd              # gravitational acceleration

# block dimensions
m = np.ones(n)/m_nd         # block mass
w = 1*np.ones(n)/l_nd     # block width
h = 0.2*np.ones(n) /l_nd    # block height
mu = mu_val*np.ones(n)         # friction coefficient

dtime = 0.01
ntime = 1000
tf = dtime*ntime
t = np.linspace(0,tf,ntime)          # time array

xbb = np.zeros(ntime)
xbbdot = np.zeros(ntime)
xbbddot = np.zeros(ntime)

# constraint count
ng = 0
ngamma = 0
nN = 2*(n-1)
nF = n-1
nX = 3*ndof+3*ng+3*ngamma+3*nN+2*nF

# fixed basis vectors
Ex = np.array([1,0])
Ey = np.array([0,1])

# generalized alpha parameters
MAXITERn = 20           # maximum number of newton iterations
r = 0.3
rho_inf = 0.5
rho_infinity_initial = rho_inf
# eq. 72
alpha_m = (2*rho_inf-1)/(rho_inf+1)
alpha_f = rho_inf/(rho_inf+1)
gama = 0.5+alpha_f-alpha_m
beta = 0.25*(0.5+gama)**2

# coefficients of restitution
eN = 0                  # normal coefficient of restitution
eF = 0                  # friction coefficient of retitution

# mass matrix (constant)
M_diagonal = np.zeros(ndof)
for i in range(1,n):
    M_diagonal[(i-1)*3:i*3] = np.array([m[i], m[i], m[i]/12*(w[i]**2+h[i]**2)])
M = np.diag(M_diagonal.flatten())

# applied forces (weight)
force = np.zeros(ndof)
for i in range(1,n):
    force[(i-1)*3:i*3] = np.array([0,-m[i]*gr,0])

def remove_last_line(file_path):
    # Read all lines from the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Rewrite the file without the last line
    with open(file_path, 'w') as file:
        file.writelines(lines[:-1])

def get_gN(q,u,a):
    """Calculate the normal contact constraint"""

    q_copy = q
    u_copy = u
    a_copy = a

    qbb = np.array([xbb[iter],h[0]/2,0])
    ubb = np.array([xbbdot[iter],0,0])
    abb = np.array([xbbddot[iter],0,0])
    
    q = np.concatenate((qbb,q),axis=None)
    u = np.concatenate((ubb,u),axis=None)
    a = np.concatenate((abb,a),axis=None)

    x, y, theta = get_xyt(q)
    xdot, ydot, thetadot = get_xyt(u)
    xddot, yddot, thetaddot = get_xyt(a)

    ex = np.zeros((n,2))
    ey = np.zeros((n,2))
    
    exdot = np.zeros((n,2))
    eydot = np.zeros((n,2))
    
    exddot = np.zeros((n,2))
    eyddot = np.zeros((n,2))
    
    r = np.zeros((n,2))
    v = np.zeros((n,2))
    a = np.zeros((n,2))

    ra = np.zeros((n,2))
    rb = np.zeros((n,2))
    rc = np.zeros((n,2))
    rd = np.zeros((n,2))
    
    va = np.zeros((n,2))
    vb = np.zeros((n,2))
    vc = np.zeros((n,2))
    vd = np.zeros((n,2))
    
    aa = np.zeros((n,2))
    ab = np.zeros((n,2))
    ac = np.zeros((n,2))
    ad = np.zeros((n,2))
    
    for i in range(n):  # range of blocks
        ex[i,:] = np.array([np.cos(theta[i]), np.sin(theta[i])])
        ey[i,:] = np.array([-np.sin(theta[i]), np.cos(theta[i])])

        exdot[i,:] = thetadot[i]*ey[i,:]
        eydot[i,:] = -thetadot[i]*ex[i,:]

        exddot[i,:] = thetaddot[i]*ey[i,:]-thetadot[i]**2*ex[i,:]
        eyddot[i,:] = -thetaddot[i]*ex[i,:]-thetadot[i]**2*ey[i,:]

        r[i,:] = np.array([x[i],y[i]])
        v[i,:] = np.array([xdot[i],ydot[i]])
        a[i,:] = np.array([xddot[i],yddot[i]])

        ra[i,:] = r[i,:]+w[i]/2*ex[i,:]+h[i]/2*ey[i,:]
        rb[i,:] = r[i,:]-w[i]/2*ex[i,:]+h[i]/2*ey[i,:]
        rc[i,:] = r[i,:]-w[i]/2*ex[i,:]-h[i]/2*ey[i,:]
        rd[i,:] = r[i,:]+w[i]/2*ex[i,:]-h[i]/2*ey[i,:]

        va[i,:] = v[i,:]+w[i]/2*exdot[i,:]+h[i]/2*eydot[i,:]
        vb[i,:] = v[i,:]-w[i]/2*exdot[i,:]+h[i]/2*eydot[i,:]
        vc[i,:] = v[i,:]-w[i]/2*exdot[i,:]-h[i]/2*eydot[i,:]
        vd[i,:] = v[i,:]+w[i]/2*exdot[i,:]-h[i]/2*eydot[i,:]

        aa[i,:] = a[i,:]+w[i]/2*exddot[i,:]+h[i]/2*eyddot[i,:]
        ab[i,:] = a[i,:]-w[i]/2*exddot[i,:]+h[i]/2*eyddot[i,:]
        ac[i,:] = a[i,:]-w[i]/2*exddot[i,:]-h[i]/2*eyddot[i,:]
        ad[i,:] = a[i,:]+w[i]/2*exddot[i,:]-h[i]/2*eyddot[i,:]

    gNa = np.zeros(n)
    gNb = np.zeros(n)
    gNc = np.zeros(n)
    gNd = np.zeros(n)

    dgNa_dq = np.zeros((n,3*n))
    dgNb_dq = np.zeros((n,3*n))
    dgNc_dq = np.zeros((n,3*n))
    dgNd_dq = np.zeros((n,3*n))

    gNa_dot = np.zeros(n)
    gNb_dot = np.zeros(n)
    gNc_dot = np.zeros(n)
    gNd_dot = np.zeros(n)

    gNa_ddot = np.zeros(n)
    gNb_ddot = np.zeros(n)
    gNc_ddot = np.zeros(n)
    gNd_ddot = np.zeros(n)

    gammaF_allcases = np.zeros((n,4))
    dgammaF_dq_allcases = np.zeros((n,4,3*n))
    gammaFdot_allcases = np.zeros((n,4))
    
    gN = np.zeros(nN)
    gNdot = np.zeros(nN)
    gNddot = np.zeros(nN)
    WN = np.zeros((nN,3*n))

    corners = np.zeros(nN,dtype=np.int8)

    for i in range(1,n):
        gNa[i] = np.dot(rc[i,:]-ra[i-1,:],ey[i,:])
        gNb[i] = np.dot(rc[i,:]-rb[i-1,:],ey[i,:])
        gNc[i] = np.dot(rc[i,:]-rb[i-1,:],ey[i-1,:])
        gNd[i] = np.dot(rd[i,:]-ra[i-1,:],ey[i-1,:])

        dgNa_dq[i,3*(i-1):3*(i+1)] = np.array([np.sin(theta[i]), -np.cos(theta[i]),-h[i-1]*np.sin(theta[i] - theta[i-1])/2 - w[i-1]*np.cos(theta[i] - theta[i-1])/2,-np.sin(theta[i]), np.cos(theta[i]),h[i-1]*np.sin(theta[i] - theta[i-1])/2 + w[i-1]*np.cos(theta[i] - theta[i-1])/2 - x[i]*np.cos(theta[i]) + x[i-1]*np.cos(theta[i]) - y[i]*np.sin(theta[i]) + y[i-1]*np.sin(theta[i])])
        dgNb_dq[i,3*(i-1):3*(i+1)] = np.array([np.sin(theta[i]), -np.cos(theta[i]),-h[i-1]*np.sin(theta[i] - theta[i-1])/2 + w[i-1]*np.cos(theta[i] - theta[i-1])/2,-np.sin(theta[i]), np.cos(theta[i]),h[i-1]*np.sin(theta[i] - theta[i-1])/2 - w[i-1]*np.cos(theta[i] - theta[i-1])/2 - x[i]*np.cos(theta[i]) + x[i-1]*np.cos(theta[i]) - y[i]*np.sin(theta[i]) + y[i-1]*np.sin(theta[i])])
        dgNc_dq[i,3*(i-1):3*(i+1)] = np.array([np.sin(theta[i-1]), -np.cos(theta[i-1]),-h[i]*np.sin(theta[i] - theta[i-1])/2 + w[i]*np.cos(theta[i] - theta[i-1])/2 - x[i]*np.cos(theta[i-1]) + x[i-1]*np.cos(theta[i-1]) - y[i]*np.sin(theta[i-1]) + y[i-1]*np.sin(theta[i-1]),-np.sin(theta[i-1]), np.cos(theta[i-1]),h[i]*np.sin(theta[i] - theta[i-1])/2 - w[i]*np.cos(theta[i] - theta[i-1])/2])
        dgNd_dq[i,3*(i-1):3*(i+1)] = np.array([np.sin(theta[i-1]), -np.cos(theta[i-1]),-h[i]*np.sin(theta[i] - theta[i-1])/2 - w[i]*np.cos(theta[i] - theta[i-1])/2 - x[i]*np.cos(theta[i-1]) + x[i-1]*np.cos(theta[i-1]) - y[i]*np.sin(theta[i-1]) + y[i-1]*np.sin(theta[i-1]),-np.sin(theta[i-1]), np.cos(theta[i-1]),h[i]*np.sin(theta[i] - theta[i-1])/2 + w[i]*np.cos(theta[i] - theta[i-1])/2])

        gNa_dot[i] = np.dot(vc[i,:]-va[i-1,:],ey[i,:])+np.dot(rc[i,:]-ra[i-1,:],eydot[i,:])
        gNb_dot[i] = np.dot(vc[i,:]-vb[i-1,:],ey[i,:])+np.dot(rc[i,:]-rb[i-1,:],eydot[i,:])
        gNc_dot[i] = np.dot(vc[i,:]-vb[i-1,:],ey[i-1,:])+np.dot(rc[i,:]-rb[i-1,:],eydot[i-1,:])
        gNd_dot[i] = np.dot(vd[i,:]-va[i-1,:],ey[i-1,:])+np.dot(rd[i,:]-ra[i-1,:],eydot[i-1,:])

        gNa_ddot[i] = np.dot(ac[i,:]-aa[i-1,:],ey[i,:])+2*np.dot(vc[i,:]-va[i-1,:],eydot[i,:])+np.dot(rc[i,:]-ra[i-1,:],eyddot[i,:])
        gNb_ddot[i] = np.dot(ac[i,:]-ab[i-1,:],ey[i,:])+2*np.dot(vc[i,:]-vb[i-1,:],eydot[i,:])+np.dot(rc[i,:]-rb[i-1,:],eyddot[i,:])
        gNc_ddot[i] = np.dot(ac[i,:]-ab[i-1,:],ey[i-1,:])+2*np.dot(vc[i,:]-vb[i-1,:],eydot[i-1,:])+np.dot(rc[i,:]-rb[i-1,:],eyddot[i-1,:])
        gNd_ddot[i] = np.dot(ad[i,:]-aa[i-1,:],ey[i-1,:])+2*np.dot(vd[i,:]-va[i-1,:],eydot[i-1,:])+np.dot(rd[i,:]-ra[i-1,:],eyddot[i-1,:])

        gammaF_allcases[i,0] = np.dot(va[i-1,:]-vc[i,:],ex[i,:])    # a
        gammaF_allcases[i,1] = np.dot(vb[i-1,:]-vc[i,:],ex[i,:])    # b
        gammaF_allcases[i,2] = np.dot(va[i-1,:]-vc[i,:],ex[i-1,:])  # c
        gammaF_allcases[i,3] = np.dot(va[i-1,:]-vd[i,:],ex[i-1,:])  # d

        dgammaF_dq_allcases[i,0,3*(i-1):3*(i+1)] = np.array([np.cos(theta[i]), np.sin(theta[i]),-h[i-1]*np.cos(theta[i] - theta[i-1])/2 + w[i-1]*np.sin(theta[i] - theta[i-1])/2,-np.cos(theta[i]), -np.sin(theta[i]), -h[i]/2])
        dgammaF_dq_allcases[i,1,3*(i-1):3*(i+1)] = np.array([np.cos(theta[i]), np.sin(theta[i]),-h[i-1]*np.cos(theta[i] - theta[i-1])/2 - w[i-1]*np.sin(theta[i] - theta[i-1])/2,-np.cos(theta[i]), -np.sin(theta[i]), -h[i]/2])
        dgammaF_dq_allcases[i,2,3*(i-1):3*(i+1)] = np.array([np.cos(theta[i-1]), np.sin(theta[i-1]), -h[i-1]/2, -np.cos(theta[i-1]),-np.sin(theta[i-1]),-h[i]*np.cos(theta[i] - theta[i-1])/2 - w[i]*np.sin(theta[i] - theta[i-1])/2])
        dgammaF_dq_allcases[i,3,3*(i-1):3*(i+1)] = np.array([np.cos(theta[i-1]), np.sin(theta[i-1]), -h[i-1]/2, -np.cos(theta[i-1]),-np.sin(theta[i-1]),-h[i]*np.cos(theta[i] - theta[i-1])/2 + w[i]*np.sin(theta[i] - theta[i-1])/2])

        gammaFdot_allcases[i,0] = np.dot(va[i-1,:]-vc[i,:],exdot[i,:])+np.dot(aa[i-1,:]-ac[i,:],ex[i,:])        # a
        gammaFdot_allcases[i,1] = np.dot(vb[i-1,:]-vc[i,:],exdot[i,:])+np.dot(ab[i-1,:]-ac[i,:],ex[i,:])        # b
        gammaFdot_allcases[i,2] = np.dot(va[i-1,:]-vc[i,:],exdot[i-1,:])+np.dot(aa[i-1,:]-ac[i,:],ex[i-1,:])    # c
        gammaFdot_allcases[i,3] = np.dot(va[i-1,:]-vd[i,:],exdot[i-1,:])+np.dot(aa[i-1,:]-ad[i,:],ex[i-1,:])    # d

        if np.dot(rb[i-1,:]-rd[i,:],ex[i,:])>0 or np.dot(ra[i-1,:]-rc[i,:],ex[i,:])<0:
            # complete horizontal detachement
            corners[2*(i-1)] = 4
            corners[2*(i-1)+1] = 4
        else:
            if np.dot(rb[i-1,:]-rc[i,:],ex[i,:])>0:
                # b = 1
                gN[2*(i-1)] = gNb[i]
                gNdot[2*(i-1)] = gNb_dot[i]
                gNddot[2*(i-1)] = gNb_ddot[i]
                WN[2*(i-1),:] = dgNb_dq[i,:]
                corners[2*(i-1)] = 1
            else:
                # c = 2
                gN[2*(i-1)] = gNc[i]
                gNdot[2*(i-1)] = gNc_dot[i]
                gNddot[2*(i-1)] = gNc_ddot[i]
                WN[2*(i-1),:] = dgNc_dq[i,:]
                corners[2*(i-1)] = 2
            if np.dot(ra[i-1,:]-rd[i,:],ex[i,:])<0:
                # a = 0
                gN[2*(i-1)+1] = gNa[i]
                gNdot[2*(i-1)+1] = gNa_dot[i]
                gNddot[2*(i-1)+1] = gNa_ddot[i]
                WN[2*(i-1)+1,:] = dgNa_dq[i,:]
                corners[2*(i-1)+1] = 0
            else:
                # d = 3
                gN[2*(i-1)+1] = gNd[i]
                gNdot[2*(i-1)+1] = gNd_dot[i]
                gNddot[2*(i-1)+1] = gNd_ddot[i]
                WN[2*(i-1)+1,:] = dgNd_dq[i,:]
                corners[2*(i-1)+1] = 3
          

    WN = np.transpose(WN[:,3:3*n])    # Remove derivatives wrt coordinates of first block

    q = q_copy
    u = u_copy
    a = a_copy

    return gN, gNdot, gNddot, WN, corners, gammaF_allcases, dgammaF_dq_allcases, gammaFdot_allcases

def get_R(X,prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF,*index_sets):
    """Calculates the residual"""

    [prev_a,_,_,_,_,_,_,_,_,_,prev_lambdaN,_,prev_lambdaF] = get_X_components(prev_X)
    [a,U,Q,Kappa_g,Lambda_g,lambda_g,Lambda_gamma,lambda_gamma,
     KappaN,LambdaN,lambdaN,LambdaF,lambdaF] = get_X_components(X)
    
    # AV - Auxiliary Variables [abar, lambdaNbar, lambdaFbar]
    prev_abar = prev_AV[0:ndof]
    prev_lambdaNbar = prev_AV[ndof:ndof+nN]
    prev_lambdaFbar = prev_AV[ndof+nN:ndof+nN+nF]

    # auxiliary variables update
    # eq. 49
    abar = (alpha_f*prev_a+(1-alpha_f)*a-alpha_m*prev_abar)/(1-alpha_m)
    # eq. 96
    lambdaNbar = (alpha_f*prev_lambdaN+(1-alpha_f)*lambdaN-alpha_m*prev_lambdaNbar)/(1-alpha_m)
    # eq. 114
    lambdaFbar = (alpha_f*prev_lambdaF+(1-alpha_f)*lambdaF-alpha_m*prev_lambdaFbar)/(1-alpha_m)

    AV = np.concatenate((abar,lambdaNbar,lambdaFbar),axis=None)

    # velocity update (73)
    u = prev_u+dtime*((1-gama)*prev_abar+gama*abar)+U
    # position update (73)
    q = prev_q+dtime*prev_u+dtime**2/2*((1-2*beta)*prev_abar+2*beta*abar)+Q

    # bilateral constraints at position level
    g = np.zeros((ng))
    gdot = np.zeros((ng))
    gddot = np.zeros((ng))
    Wg = np.zeros((ndof,ng))

    # bilateral constraints at velocity level
    gamma = np.zeros((ngamma))
    gammadot = np.zeros((ngamma))
    Wgamma = np.zeros((ndof,ngamma))

    # normal gap distance constraints and some frictional quantities
    gN, gNdot, gNddot, WN, corners, gammaF_allcases, dgammaF_dq_allcases, gammaFdot_allcases = get_gN(q,u,a)

    # eq. 44
    ksiN = gNdot+eN*prev_gNdot
    # discrete normal percussion eq. 95
    PN = LambdaN+dtime*((1-gama)*prev_lambdaNbar+gama*lambdaNbar)
    # eq. 102
    Kappa_hat_N = KappaN+dtime**2/2*((1-2*beta)*prev_lambdaNbar+2*beta*lambdaNbar)

    # slip speed frictional constraints
    gammaF = np.zeros((nF))
    gammaFdot = np.zeros((nF))
    WF = np.zeros((ndof,nF))
    
    if index_sets == ():
        A = np.zeros((nN))
        B = np.zeros((nN))
        C = np.zeros((nN))

        for i in range(nN):
            if corners[i] != 4 and r*gN[i] - Kappa_hat_N[i] <=0:
                A[i] = 1
                if r*ksiN[i]-PN[i] <= 0:
                    B[i] = 1
                    if r*gNddot[i]-lambdaN[i] <= 0:
                        C[i] = 1

        # A = np.ones((nN))
        # B = np.ones((nN))
        # C = np.ones((nN))
        
    else:
        A = index_sets[0]
        B = index_sets[1]
        C = index_sets[2]
        D = index_sets[3]
        E = index_sets[4]

        flag_slip_check = False

        for i in range(nN):
            if corners[i] == 4 and (A[i]!=0 or B[i]!=0 or C[i]!=0):
                flag_slip_check = True
                A[i] = 0
                B[i] = 0
                C[i] = 0
        
        if flag_slip_check == True:
            for i in range(nF):
                if A[2*i] == 0 and A[2*i+1] == 0:
                    D[i] = 0
                    E[i] = 0


    # Assigning the friction force based on the contact region
    for i in range(nF):
        if A[2*i] == 1:
            # if first contact is closed
            gammaF[i] = gammaF_allcases[i+1,corners[2*i]]
            gammaFdot[i] = gammaFdot_allcases[i+1,corners[2*i]]
            WF[:,i] = np.transpose(dgammaF_dq_allcases[i+1,corners[2*i],3:3*n])
        else:
            # the second contact is closed (A[2*i+1] == 1) OR both contacts are open
            gammaF[i] = gammaF_allcases[i+1,corners[2*i+1]%4]
            gammaFdot[i] = gammaFdot_allcases[i+1,corners[2*i+1]%4]
            WF[:,i] = np.transpose(dgammaF_dq_allcases[i+1,corners[2*i+1]%4,3:3*n])

    # eq. 48
    ksiF = gammaF+eN*prev_gammaF
    # eq. 113
    PF = LambdaF+dtime*((1-gama)*prev_lambdaFbar+gama*lambdaFbar)

    if index_sets == ():
        D = np.zeros(nF)
        E = np.zeros(nF)

        for i in range(nF):
            if A[2*i] == 1 or A[2*i+1] == 1:
                if np.abs(r*ksiF[i]-PF[i])<=mu[i+1]*(PN[2*i]+PN[2*i+1]):
                    # D-stick
                    D[i] = 1
                    # if B[2*i] == 1 or B[2*i+1] == 1:
                    if np.abs(r*gammaFdot[i]-lambdaF[i])<=mu[i+1]*(lambdaN[2*i]+lambdaN[2*i+1]):
                        # E-stick
                        E[i] = 1

        # D = np.ones(nF)
        # E = np.ones(nF)

    R_LambdaN = np.zeros(nN)
    R_lambdaN = np.zeros(nN)
    R_KappaN = np.zeros(nN)
    for i in range(nN):
        if A[i] == 1:
            R_KappaN[i] = gN[i]
        else:            
            R_KappaN[i] = Kappa_hat_N[i]
        if B[i] == 1:
            R_LambdaN[i] = ksiN[i]
        else:
            R_LambdaN[i] = PN[i]
        if C[i] == 1:
            R_lambdaN[i] = gNddot[i]
        else:
            R_lambdaN[i] = lambdaN[i]

    R_LambdaF = np.zeros(nF)
    R_lambdaF = np.zeros(nF)
    for i in range(nF):
        if A[2*i] == 1 or A[2*i+1] == 1:
            if D[i] == 1:
                # D-stick            
                R_LambdaF[i] = ksiF[i]
                if E[i] == 1:
                    # E_stick
                    R_lambdaF[i] = gammaFdot[i]
                else:
                    # E_slip
                    R_lambdaF[i] = lambdaF[i]+mu[i+1]*(lambdaN[2*i]+lambdaN[2*i+1])*gammaFdot[i]/np.abs(gammaFdot[i])
            else:
                # D_slip 
                R_LambdaF[i] = PF[i]+mu[i+1]*(PN[2*i]+PN[2*i+1])*ksiF[i]/np.abs(ksiF[i])
                R_lambdaF[i] = lambdaF[i]+mu[i+1]*(lambdaN[2*i]+lambdaN[2*i+1])*gammaF[i]/np.abs(gammaF[i])
        else: # no touch
            R_LambdaF[i] = PF[i]
            R_lambdaF[i] = lambdaF[i]
        
        
    Rs = np.concatenate(([M@a-force-Wg@lambda_g-Wgamma@lambda_gamma-WN@lambdaN-WF@lambdaF],
               [M@U-Wg@Lambda_g-Wgamma@Lambda_gamma-WN@LambdaN-WF@LambdaF],
               [M@Q-Wg@Kappa_g-WN@KappaN-dtime/2*(Wgamma@Lambda_gamma+WF@LambdaF)],
               g,
               gdot,
               gddot,
               gamma,
               gammadot),axis=None)
    
    Rc = np.concatenate((R_KappaN, R_LambdaN, R_lambdaN, R_LambdaF, R_lambdaF),axis=None)
    
    R = np.concatenate([Rs, Rc],axis=None)

    if index_sets == ():
        print(f"A={A}")
        print(f"B={B}")
        print(f"C={C}")
        print(f"D={D}")
        print(f"E={E}")
        return R, AV, q, u, gNdot, gammaF, A, B, C, D, E
    else:
        return R, AV, q, u, gNdot, gammaF

def get_R_J(X,prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF,*fixed_contact):
    global A_save, B_save, C_save, D_save, E_save
    global iter

    epsilon = 1e-6
    fixed_contact_regions = False

    if fixed_contact != ():
        fixed_contact = fixed_contact[0]
        fixed_contact_regions = True
        A = fixed_contact[0:nN]
        B = fixed_contact[nN:2*nN]
        C = fixed_contact[2*nN:3*nN]
        D = fixed_contact[3*nN:3*nN+nF]
        E = fixed_contact[3*nN+nF:3*nN+2*nF]
        R, AV, q, u, gNdot, gammaF =  get_R(X,prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF, A, B, C, D, E)
    else:
        R, AV, q, u, gNdot, gammaF, A, B, C, D, E = get_R(X,prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF)
        contacts_nu = np.concatenate((A,B,C,D,E),axis=None)

    # Initializing the Jacobian
    J = np.zeros((nX,nX))
    I = np.identity(nX)

    A_save[:,iter] = A
    B_save[:,iter] = B
    C_save[:,iter] = C
    D_save[:,iter] = D
    E_save[:,iter] = E

    # Constructing the Jacobian column by column
    for i in range(nX):
        # print(i)
        R_plus_epsilon,_,_,_,_,_ = get_R(X+epsilon*I[:,i],prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF, A, B, C, D, E)
        J[:,i] = (R_plus_epsilon-R)/epsilon

    if fixed_contact_regions:
        return R, AV, q, u, gNdot, gammaF, J
    else:
        return R, AV, q, u, gNdot, gammaF, J, contacts_nu

def update(prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF,*fixed_contact):
    """Takes components at time t and return values at time t+dt"""
    
    nu = 0
    X = prev_X
    
    if fixed_contact != ():
        fixed_contact = fixed_contact[0]
        fixed_contact_regions = True
    else:
        fixed_contact_regions = False

    try:
        if fixed_contact_regions == True:
            R, AV, q, u, gNdot, gammaF, J = get_R_J(X,prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF,fixed_contact)
        else:
            R, AV, q, u, gNdot, gammaF, J, contacts_nu = get_R_J(X,prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF)
            contacts = np.zeros((MAXITERn+1,3*nN+2*nF))
            contacts[nu,:] = contacts_nu
        norm_R = np.linalg.norm(R,np.inf)
        print(f"nu = {nu}")
        print(f"norm(R) = {norm_R}")

        while np.abs(np.linalg.norm(R,np.inf))>10**(-10) and nu<MAXITERn:
            # Newton Update
            X = X-np.linalg.solve(J,R)
            # Calculate new EOM and residual
            nu = nu+1
            if fixed_contact_regions:
                R, AV, q, u, gNdot, gammaF, J = get_R_J(X,prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF,fixed_contact)
            else:
                R, AV, q, u, gNdot, gammaF, J, contacts_nu = get_R_J(X,prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF)
                contacts[nu,:] = contacts_nu
            norm_R = np.linalg.norm(R,np.inf)
            print(f"nu = {nu}")
            print(f"norm(R) = {norm_R}")
        if nu == MAXITERn:
            print(f"No Convergence for nu = {nu} at rho_inf = {rho_inf}")
            raise TimeoutError

    except TimeoutError as e:
        if fixed_contact_regions is False:
            # if unique contact regions were already determined, don't recalculate them
            unique_contacts = np.unique(contacts, axis=0)
            do_not_unpack = True    
            # because if the number of contact regions is 6 which is the original number
            # of outputs of update, each row of unique contacts will be assinged as an output variable
            return unique_contacts, do_not_unpack
        return 
    except np.linalg.LinAlgError as e:
        print(e)
        # increment rho_inf        
        update_rho_inf()
        # calling function recursively
        update(prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF,fixed_contact)
    except Exception as e:
        raise
    
    return X,AV,q,u,gNdot,gammaF

def update_rho_inf():
    global rho_inf, alpha_m, alpha_f, gama, beta 
    rho_inf = rho_inf+0.05  #0.01
    print(rho_inf)
    if np.abs(rho_inf - rho_infinity_initial) < 0.001:
        print("possibility of infinite loop")
        raise Exception
    if rho_inf > 1.001:
        rho_inf = 0
    # eq. 72
    alpha_m = (2*rho_inf-1)/(rho_inf+1)
    alpha_f = rho_inf/(rho_inf+1)
    gama = 0.5+alpha_f-alpha_m
    beta = 0.25*(0.5+gama)**2

def get_X_components(X):
    a = X[0:ndof]
    U = X[ndof:2*ndof]
    Q = X[2*ndof:3*ndof]
    Kappa_g = X[3*ndof:3*ndof+ng]
    Lambda_g = X[3*ndof+ng:3*ndof+2*ng]
    lambda_g = X[3*ndof+2*ng:3*ndof+3*ng]
    Lambda_gamma = X[3*ndof+3*ng:3*ndof+3*ng+ngamma]
    lambda_gamma = X[3*ndof+3*ng+ngamma:3*ndof+3*ng+2*ngamma]
    Kappa_N = X[3*ndof+3*ng+2*ngamma:3*ndof+3*ng+2*ngamma+nN]
    Lambda_N = X[3*ndof+3*ng+2*ngamma+nN:3*ndof+3*ng+2*ngamma+2*nN]
    lambda_N = X[3*ndof+3*ng+2*ngamma+2*nN:3*ndof+3*ng+2*ngamma+3*nN]
    Lambda_F = X[3*ndof+3*ng+2*ngamma+3*nN:3*ndof+3*ng+2*ngamma+3*nN+nF]
    lambda_F = X[3*ndof+3*ng+2*ngamma+3*nN+nF:3*ndof+3*ng+2*ngamma+3*nN+2*nF]
    return a,U,Q,Kappa_g,Lambda_g,lambda_g,Lambda_gamma,lambda_gamma,\
        Kappa_N,Lambda_N,lambda_N,Lambda_F,lambda_F

def get_xyt(q):
    x = np.zeros(n)
    y = np.zeros(n)
    theta = np.zeros(n)
    for i in range(n):
        x[i] = q[3*i]
        y[i] = q[3*i+1]
        theta[i] = q[3*i+2]
    return x, y, theta

def increment_saved_arrays():
    global q_save, u_save, X_save, gNdot_save, gammaF_save, AV_save
    q_save_addition = np.tile(q_save[leaves_counter,:,:],(1,1,1))
    q_save = np.vstack((q_save,q_save_addition))
    u_save_addition = np.tile(u_save[leaves_counter,:,:],(1,1,1))
    u_save = np.vstack((u_save,u_save_addition))
    X_save_addition = np.tile(X_save[leaves_counter,:,:],(1,1,1))
    X_save = np.vstack((X_save,X_save_addition))
    gNdot_save_addition = np.tile(gNdot_save[leaves_counter,:,:],(1,1,1))
    gNdot_save = np.vstack((gNdot_save,gNdot_save_addition))
    gammaF_save_addition = np.tile(gammaF_save[leaves_counter,:,:],(1,1,1))
    gammaF_save = np.vstack((gammaF_save,gammaF_save_addition))
    AV_save_addition = np.tile(AV_save[leaves_counter,:,:],(1,1,1))
    AV_save = np.vstack((AV_save,AV_save_addition))

def solve(iter_start):
    global q_save, u_save, X_save, gNdot_save, gammaF_save, AV_save
    global leaves_counter
    global iter

    fixed_contact_regions = False
    increment_leaves = True

    prev_X = X_save[leaves_counter,:,iter_start-1]
    prev_AV = AV_save[leaves_counter,:,iter_start-1]
    prev_q = q_save[leaves_counter,:,iter_start-1]
    prev_u = u_save[leaves_counter,:,iter_start-1]
    prev_gNdot = gNdot_save[leaves_counter,:,iter_start-1]
    prev_gammaF = gammaF_save[leaves_counter,:,iter_start-1]

    for iter in range(iter_start,ntime):
        print(f"iteration {iter}")

        current_time = time.time()
        if current_time-start_time>1800:
            raise Exception

        try:
            X,AV,q,u,gNdot,gammaF = update(prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF)

            prev_X = X
            prev_AV = AV
            prev_q = q
            prev_u = u
            prev_gNdot = gNdot
            prev_gammaF = gammaF

            q_save[leaves_counter,:,iter] = prev_q
            u_save[leaves_counter,:,iter] = prev_u
            X_save[leaves_counter,:,iter] = prev_X
            gNdot_save[leaves_counter,:,iter] = prev_gNdot
            gammaF_save[leaves_counter,:,iter] = prev_gammaF
            AV_save[leaves_counter,:,iter] = prev_AV
        except ValueError as e:
            unique_contacts,_ = update(prev_X,prev_AV,prev_q,prev_u,prev_gNdot,prev_gammaF)
            solve_bifurcation(iter,unique_contacts)
            increment_leaves = False
            break   # this break is important 
        except RuntimeWarning:
            if iters_per_oscillation<160:
                iters_per_oscillation = iters_per_oscillation+10
            else:
                raise Exception
        except Exception as e:
            raise

    if increment_leaves == True:

        increment_saved_arrays()

        q_save_shape = np.shape(q_save)
        leaves_counter = leaves_counter + 1
        print(f'leaves_counter = {leaves_counter}')

        if leaves_counter>50:
            raise Exception

    return

global bif_counter  # used in bifurcation log
bif_counter = 0

def solve_bifurcation(iter_bif,*fixed_contact_region_params):
    global q_save, u_save, X_save, gNdot_save, gammaF_save, AV_save
    global leaves_counter
    global iter
    global bif_counter
    bif_counter +=1

    global increment_leaves
    increment_leaves = True # I think this is unneccessary

    # f.write(f'Running solve_bifurcations at iter_bif {iter_bif} and leaf {leaves_counter}\n') 
    
    # fixed_contact_regions = True
    unique_contacts = fixed_contact_region_params[0]
    n_unique_contacts = np.shape(unique_contacts)[0]

    # f.write(f'The number of unique contacts is {n_unique_contacts}\n')

    nonconvergence_counter = 0

    for k in range(n_unique_contacts):
        iter = iter_bif

        # print(f"k = {k}")
        # g.write("     |"*bif_counter)
        # g.write(f'__ {k+1} of {n_unique_contacts}  ')
        
        try:
            fixed_contact  = unique_contacts[k,:]
            X,AV,q,u,gNdot,gammaF = update(X_save[leaves_counter,:,iter_bif-1],AV_save[leaves_counter,:,iter_bif-1],
                                       q_save[leaves_counter,:,iter_bif-1],u_save[leaves_counter,:,iter_bif-1],
                                       gNdot_save[leaves_counter,:,iter_bif-1],gammaF_save[leaves_counter,:,iter_bif-1],
                                       fixed_contact)

            prev_X = X
            prev_AV = AV
            prev_q = q
            prev_u = u
            prev_gNdot = gNdot
            prev_gammaF = gammaF

            q_save[leaves_counter,:,iter_bif] = prev_q
            u_save[leaves_counter,:,iter_bif] = prev_u
            X_save[leaves_counter,:,iter_bif] = prev_X
            gNdot_save[leaves_counter,:,iter_bif] = prev_gNdot
            gammaF_save[leaves_counter,:,iter_bif] = prev_gammaF
            AV_save[leaves_counter,:,iter_bif] = AV

            # f.write(f'{k}-th unique contact convergence successfull\n')            

            solve(iter_bif+1)   

        except TypeError as e:       
            # make a provision for if we always passed and never converged.
            # f.write(f'{k}-th unique contact convergence unsuccessfull\n') 
            # g.write('unsuccessful\n')   
            nonconvergence_counter = nonconvergence_counter+1
            # f.write(f'nonconvergence_counter = {nonconvergence_counter}\n')
            # num_bif_contacts[leaves_counter,3] = num_bif_contacts[leaves_counter,3]-1
            if nonconvergence_counter == n_unique_contacts:
                # exception raised when None of the fixed contact regions converged
                # raise Exception
                global MAXITERn
                nonconvergence_counter = 0
                if MAXITERn < 10:  # INCOMPLETE, TO BE FIXED
                    # try to increase number of iterations
                    MAXITERn = 200
                    solve_bifurcation(iter_bif,unique_contacts)
                else:
                    try:
                        update_rho_inf()
                        solve_bifurcation(iter_bif,unique_contacts)
                    except:
                        # we cannot update rho_inf anymore
                        # we need to abandon this leaf
                        # g.write(f'bifurcation convergence failed\n')
                        pass
                # solve_bifurcation(iter_bif,unique_contacts) # maybe wrong, remove
            else:
                pass

    bif_counter = bif_counter-1
    # increment_leaves = False

    return 

# initial normal force
lambdaN0 = np.zeros(nN)
lambdaN0[nN-2] = 0
lambdaN0[nN-1] = m[n-1]*gr
for i in range(n-1,1,-1):
    lambdaN0[2*(i-1)-2] = 0
    lambdaN0[2*(i-1)-1] = m[i-1]*gr+lambdaN0[2*(i-1)+1]

X0 = np.zeros(nX)
X0[3*ndof+3*ng+2*ngamma+2*nN:3*ndof+3*ng+2*ngamma+3*nN] = lambdaN0

# initial auxiliary variable
nAV = ndof+nN+nF
AV0 = np.zeros(nAV)
AV0[ndof:ndof+nN] = lambdaN0

# initial position
x0 = np.zeros(n)
y0 = np.zeros(n)
# the blocks are initially stacked
y0[0] = h[0]/2
for i in range(1,n):
    y0[i] = y0[i-1]+h[i-1]/2+h[i]/2
    x0[i] = x0[i-1]+w[i]/(2*(n-i))
x0[1] = x0[1]+w[1]/100
theta0 = np.zeros(n)
# assembling position coordinates in q0
q0 = np.zeros(3*n)
for i in range(1,n):
    q0[3*i:3*i+3] = [x0[i],y0[i],theta0[i]]

# initial velocity
u0 = np.zeros(3*n)  # starting from rest

# initial normal gap speeds
gNdot0 = np.zeros(nN)   # starting from rest

# inital slip speeds
gammaF0 = np.zeros(nF)  # starting from rest

prev_X = X0
prev_AV = AV0
prev_q = q0[3:3*n]
prev_u = u0[3:3*n]
prev_gNdot = gNdot0
prev_gammaF = gammaF0

q_save = np.zeros((1,ndof,ntime))
u_save = np.zeros((1,ndof,ntime))
X_save = np.zeros((1,nX,ntime))
gNdot_save = np.zeros((1,nN,ntime))
gammaF_save = np.zeros((1,nF,ntime))
AV_save = np.zeros((1,ndof+nN+nF,ntime))

q_save[0,:,0] = prev_q
u_save[0,:,0] = prev_u
X_save[0,:,0] = prev_X
gNdot_save[0,:,0] = prev_gNdot
gammaF_save[0,:,0] = prev_gammaF
AV_save[0,:,0] = prev_AV

A_save = np.zeros((nN,ntime))
B_save = np.zeros((nN,ntime))
C_save = np.zeros((nN,ntime))
D_save = np.zeros((nF,ntime))
E_save = np.zeros((nF,ntime))

# parameters for bifurcations
fixed_contact_regions = False   # COMMENT THIS AND SEE IF NECESSARY

leaves_counter = 0
# running the code!

try:
    solve(1)
finally:
    # removing last added leaf_counter
    leaves_counter = leaves_counter-1
    # removing the last added block to q_save
    q_save = q_save[0:np.shape(q_save)[0]-1,:,:]        
    end_time = time.time()
    execution_time = end_time-start_time
    print("Execution time:", execution_time, "seconds")

    block0 = np.stack((xbb,h[0]/2*np.ones((ntime)),np.zeros((ntime))))
    block0_tiled = np.tile(block0,(np.shape(q_save)[0],1,1))
    q_save_total = np.concatenate((block0_tiled,q_save),axis=1)

    import scipy.io
    file_name = str(f'{output_path}/q.mat')
    scipy.io.savemat(file_name,dict(q=q_save_total))

    file_name_test = str(f'{output_path}/test.mat')
    scipy.io.savemat(file_name_test,dict(test=u_save))

print("Saved.")