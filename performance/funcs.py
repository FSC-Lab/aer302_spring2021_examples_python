import SAM
from aircraft import CP1

from math import sin, cos, asin

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

# # minimum power from graphics
# [PRmin_s,iMin_s] = min(PR_s)
# VPRmin_s = 10+iMin_s*1
# [PRmin_h,iMin_h] = min(PR_h)
# VPRmin_h = (10+iMin_h*1)* sqrt(rho_s / rho_h)

# # minimum power from calculation (analytical)
# VPRmin_s_calc = sqrt( 2*WS / rho_s / sqrt(3*C_D_0/K) )
# VPRmin_h_calc = sqrt( 2*WS / rho_h / sqrt(3*C_D_0/K) )

# # maximum rate of climb
# VRC_max_s_calc = VPRmin_s_calc
# VRC_max_h_calc = VPRmin_h_calc
# C_D_rcmax = 4*C_D_0 # C_D_0 = 1/3 K CL^2
# RCmax_s = (PAs - PRmin_s)/W
# RCmax_h = (PAs * sqrt(rho_h / rho_s) - PRmin_h)/W


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




# figure(1)
# plot(Vs,PR_s,'-',Vs,PA_s,'-',Vh,PR_h,'--',Vh,PA_h,'--',VPRmin_s,PRmin_s,'o',VPRmin_h,PRmin_h,'*')
# axis([0 250 0 2e06])
# # plot(V,Pr_s,'-',V,Pa_s,'-',Vh,Pr_h,'--',Vh,Pa_h,'--')
# #plot(V,Pr_s,V,Pa_s,Vh,Pr_h,Vh,Pa_h)
# title('Power Required and Available for CP1')
# xlabel(' velocity (m/s)')
# ylabel(' power (N m/s)')
# legend('PR at sea level','PA at sea level',['PR at altitude ' num2str(h) ' m'],['PA at altitude ' num2str(h) ' m'])
# htype = findobj(gcf,'type','line')
# set(htype,'linewidth',2)
# # plotdlg

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

# figure(3)
# plot(Vh_s,Vv_s,'-',Vh_h,Vv_h,'--')
# title(' Hodograph for [CP1]')
# legend('sea-level', '6000m')
# axis([0 100 0 10]) 
# xlabel(' Vh (m/s)')
# ylabel('Vv (m/s)')
# htype = findobj(gcf,'type','line')
# set(htype,'linewidth',2)
# # plotdlg  

  