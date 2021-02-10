import SAM
from aircraft import CP1

from math import sin, cos, asin, sqrt

def power_required(plane, V, air):

    q = .5*air.rho*V**2*plane.S
    CL = plane.W/q
    CD = plane.CD0 + plane.K*CL**2
    D = q*CD
    PR = D*V

    # # alternative calculation
    # num = sqrt(2*plane.W/plane.S/air.rho)
    # den =  CL**1.5/CD
    # PR = plane.W*(num/den)

    return PR

def min_power_required(plane, air):
    V_PRmin = sqrt(2*plane.W/plane.S/air.rho/
                    sqrt(3*plane.CD0/plane.K))
    return V_PRmin

def power_available(plane, air):
    if plane.engine == 'prop':
        PA = plane.PAs*(air.rho/SAM.air_sea.rho)
    return PA

def climb_rate(plane, V, air):
    PR = power_required(plane, V, air)
    PA = power_available(plane, air)
    RC = (PA - PR)/plane.W
    gmm = asin(RC/V)
    return RC, gmm

def climb_rate_max(plane, hs):

    # convert hs into list if its scaler
    try:
        n = len(hs)
    except:
        n = 1
        h = [hs]

    RCmaxs = []
    for h in hs:
        air = SAM.tropo(h)
        WS = plane.W/plane.S
        if plane.engine = 'jet':
            TW = plane.F_s*(air.rho/SAM.air_sea.rho)/plane.W
            WS = plane.W/plane.S
            temp1 = WS*TW/(3*air.rho*plane.C_D_0)
            temp2 = 1+sqrt(1+12*plane.K*plane.C_D_0/TW**2)
            V_rcmax = sqrt(temp1*temp2)
            q_rcmax = 0.5*air.rho*V_rcmax**2
            RCmax = V_rcmax*(TW - q_rcmax*plane.C_D_0/WS - plane.K*WS/q_rcmax)
        if plane.engine = 'prop':
            PA = plane.Ps*sqrt(air.rho/SAM.air_sea.rho)
            PRmin = plane.W*sqrt(plane.K*plane.C_D_0)*4/sqrt(3)*
                    sqrt(2*WS/sqrt(3*plane.C_D_0/plane.K)/air.rho)
            RCmax = (PA - PRmin)/plane.W
        RCmaxs.append(RCmax)

    # return scaler if input hs is scalar (n=1)
    if n == 1:
        return RCmaxs[0]
    else:
        return RCmaxs

def hodograph(V, gmm):
    Vx = [v*cos(g) for v, g in zip(V, gmm)]
    Vy = [v*sin(g) for v, g in zip(V, gmm)]
    return Vx, Vy

def power_vs_v(plane, air, Vs):

    PR, PA, RC, Gmm = [], [], [], []

    for v in Vs:
        pr = power_required(plane, v, air)
        pa = power_available(plane, air)
        rc, g = climb_rate(plane, v, air)

        PR.append(pr)
        PA.append(pa)
        RC.append(rc)
        Gmm.append(g)

    return PR, PA, RC, Gmm

def CL_CD_vs_h(plane, V, h):
    # standard atmosphere model
    air = SAM.tropo(h)
    q = .5*air.rho*V**2*plane.S
    CL = plane.W/q
    CD = plane.CD0 + plane.K*CL**2
    return CL/CD

def CLsqrt_CD_vs_h(plane, V, h):
    # standard atmosphere model
    air = SAM.tropo(h)
    q = .5*air.rho*V**2*plane.S
    CL = plane.W/q
    CD = plane.CD0 + plane.K*CL**2
    return sqrt(CL)/CD

def endurance(plane, air):
    CL_CD_max = 1/(2*sqrt(plane.K*plane.C_D_0))
    W0 = plane.W
    W1 = W0 - plane.Wf
    if plane.engine == 'jet':
        E = 1/plane.c_t*CL_CD_max*log(W0/W1)/3600
    # TODO
    if plane.engine == 'prop' :
        air = SAM.tropo(h)
        E = 0
    return E

def range(plane, air):
    CLsqrt_CD_max = 3/4/sqrt(sqrt(3*plane.K*plane.C_D_0**3))
    W0 = plane.W
    W1 = W0 - plane.Wf
    air = SAM.tropo(h)
    if plane.engine == 'jet':
        R = 2*sqrt(2/(air.rho*plane.S))/plane.c_t*
            CLsqrt_CD_max*(sqrt(W0) - sqrt(W1))
    # TODO
    if plane.engine == 'prop':
        R = 0

    return

# TODO 1:
# # maximum angle of climb
# poly_s = [1/2 *rho_s^2* C_D_0/WS 0 0 1/2*rho_s*PAs/W -2*K*WS]
# poly_h = [1/2 *rho_h^2* C_D_0/WS 0 0 1/2*rho_h*(PAs * sqrt(rho_h / rho_s))/W -2*K*WS]
# Vgamma_s = roots(poly_s)
# Vgamma_h = roots(poly_h)
# Vgammamax_s_est = 2*K*WS / (1/2*rho_s*PAs/W)
# Vgammamax_h_est = 2*K*WS / (1/2*rho_h*(PAs * sqrt(rho_h / rho_s))/W)

# q_gamma_s = 1/2*rho_s*Vgammamax_s_est^2
# CL_gamma_s = W / (q_gamma_s *S)
# CD_gamma_s = C_D_0 + K*CL_gamma_s^2
# D_gamma_s = q_gamma_s*S*CD_gamma_s
# RC_gamma_s = (PAs - D_gamma_s*Vgammamax_s_est)/W
# gammamax_s = asin(RC_gamma_s / Vgammamax_s_est) * 180/pi

# q_gamma_h = 1/2*rho_s*Vgammamax_h_est^2
# CL_gamma_h = W / (q_gamma_h *S)
# CD_gamma_h = C_D_0 + K*CL_gamma_h^2
# D_gamma_h = q_gamma_h*S*CD_gamma_h
# RC_gamma_h = (PAs * sqrt(rho_h / rho_s) - D_gamma_h*Vgammamax_h_est)/W
# gammamax_h = asin(RC_gamma_h / Vgammamax_h_est) * 180/pi


# TODO 2:
# # get a sense of corresponding speed with respect to P_A
# P = PAs * sqrt(rho_h / rho_s)
# v = 0.1:1:100
# p1 = [.5*rho_h*S*C_D_0 0 0 -P K*S*WS^2/(.5*rho_h)]
# # p2 = [.5*rho_h*S*C_D_0 0 -P 0 K*S*WS^2/(.5*rho_h)]
# f1 = polyval(p1,v)
# # f2 = polyval(p2,v)

# figure(2)
# plot(v,f1)
# title(['speed profile under P_A =', num2str(P)])
# # plot(v,f1,v,f2)
# htype = findobj(gcf,'type','line')
# set(htype,'linewidth',2)
  