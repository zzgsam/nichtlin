import sympy as sp
import symbtools as st
import numpy as np
from matplotlib import pyplot as plt
from sympy.utilities.lambdify import lambdify
def sys_schweberkorper(x,u):
    f_temp = sp.Matrix([x[1], k_L / m * ((k_V * (x[2] + eta_0) - A_B * x[1]) / A_SP) ** 2 - g, -1 / T_M * x[2]])
    g_temp = sp.Matrix([0, 0, k_M / T_M])
    sys_temp = f_temp + g_temp * u
    return sys_temp

# def input_schweberkorper(x):
#     input_temp=sp.Matrix([-350 * x[0] - 22.6832957040952 * x[1] ** 3 + 1.44408050744045 * x[1] ** 2 * x[2] - 28.3609827456799 * x[1] ** 2 - 0.0306447020041123 * x[1] * x[2] ** 2 + 0.952932309414394 * x[1] * x[2] - 20.0782834036154 * x[1] + 0.000216769484360984 * x[2] ** 3 - 0.00745037397804703 * x[2] ** 2 - 2.69465238109664 * x[2] + 424.843845417481) / (-0.0777354129346506 * x[1] + 0.00164961617602003 * x[2] + 0.0522378455739677)])
#     return input_temp
if __name__ == '__main__':
    x1, x2, x3= xx = sp.symbols('x1:4')
    u = sp.Symbol("u_pwm")
    print(xx,u)
    T=0.001
    T_matrix=sp.Matrix([T])
    #A_B,A_SP,m,g,T_M,k_M,k_V,k_L,eta_0,h_ub, h_lb,eta_ub,eta_lb,T,N_pred,state_r,input_r,h,h_p,eta,u_pwm = sp.symbols('A_B,A_SP,m,g,T_M,k_M,k_V,k_L,eta_0,h_ub, h_lb,eta_ub,eta_lb,T,N_pred,state_r,input_r,h,h_p,eta,u_pwm')
    A_B= 2.8274e-3  # [m**2]
    A_SP = 0.4299e-3  # [m**2]
    m = 2.8e-3  # [kg]
    g = 9.81  # [m/(s**2)]
    T_M = 0.57  # [s]
    k_M = 0.31  # [s**-1]
    k_V = 6e-5  # [m**3]
    k_L = 2.18e-4  # [kg/m]
    eta_0 = 1900 / 60  # / 60 * 2 * pi
    h_r,h_pr,eta_r=sp.symbols('h_r,h_pr,eta_r')
    h_r=0.8
    h_r0=0.6
    T_tra_value=4
    f_sys = sp.Matrix([x2, k_L / m * ((k_V * (x3 + eta_0) - A_B * x2) / A_SP)**2 - g, -1 / T_M * x3 ])
    print(f_sys)
    g_sys = sp.Matrix([0, 0,  k_M / T_M ])
    print(g)
    sys=f_sys+g_sys*u
    y=h=x1
    Lfh = st.lie_deriv(h, f_sys, xx)
    Lf2h = st.lie_deriv(h, f_sys, xx, 2)
    Lf3h = st.lie_deriv(h, f_sys, xx, 3)
    print(Lf3h)

    Lgh = st.lie_deriv(h, g_sys, xx)
    LgLfh = st.lie_deriv(Lfh, g_sys, xx)
    LgLf2h=st.lie_deriv(Lf2h, g_sys, xx)
    print("LgLf2h",LgLf2h)

    a0,a1,a2,a3=alpha=sp.symbols('alpha0:4')
    u_input=sp.Symbol('u_{input}')
    w=sp.Symbol('w') #Führungsgröße
    print(alpha,u_input)

##y_ref_trajectory
    rg = 3  # relative_grade
    y_end = h_r
    t = sp.Symbol('t')
    T_tra = sp.Symbol('T_tra')
    y_poly = sp.Symbol('y_poly')
    y_poly = 0
    for i in range(rg + 1):
        y_poly = y_poly + (h_r-h_r0) * sp.factorial(2 * rg + 1) / (sp.factorial(rg) * T_tra ** (2 * rg + 1)) * (
                    (-1) ** (rg - i) / (sp.factorial(i) * sp.factorial(rg - i) * (2 * rg - i + 1)) * T_tra ** i * t ** (
                        2 * rg - i + 1))
    print(y_poly)
    #first variable is period, the second one is current sampling time
    y_poly_func=lambdify([T_tra,t],y_poly)
    print(sp.diff(y_poly, t, 1).subs([(t, 0)]))##check if the trajectory satify all requirements
    y_r=y_poly+h_r0
    y_r_p=sp.diff(y_r,t)
    y_r_pp=sp.diff(y_r,t,2)
    y_r_ppp=sp.diff(y_r,t,3)
    print(y_r_p)
    print(y_r_pp)
    print(y_r_ppp)

    u_input=1/(LgLf2h)*(a0*(y_r-h)+a1*(y_r_p-Lfh)+a2*(y_r_pp - Lf2h)+a3*(y_r_ppp- Lf3h))

    s=sp.Symbol('s')
    nenner=(s+1)*(s+1)*(s+1)
    print(sp.expand(nenner))
    #print(nenner.subs([(a0,10),(a1,2),(a2,3)]))
    #sp.expand(sp.solve(nenner.subs([(a0,10),(a1,2),(a2,3)]),s)[0])

    u_input=sp.simplify(sp.expand(u_input.subs([(T_tra,T_tra_value),(a3,1),(a2,7),(a1,11),(a0,5)])))
    print(u_input)
    u_input_func = lambdify([t, x1,x2,x3],u_input)

    u_input_T=1/(LgLf2h)*(a0*(h_r-h)+a1*(-Lfh)+a2*( - Lf2h)+a3*(- Lf3h)).subs([(a3,1),(a2,7),(a1,11),(a0,5)])
    print(u_input_T)
    u_input_T_func = lambdify([x1,x2,x3],u_input_T)
    xx_act = sp.Matrix([h_r0, 0, 48])
    xx_print =xx_act[0]
    h_r_print=h_r0
    T_print = 0
    u_act=0
    T_act=0
    print(sys_schweberkorper(xx_act,u_act))
    print( xx_act + sys_schweberkorper(xx_act, u_act)*T_matrix )
    for i in range(5000):
        xx_act = xx_act + sys_schweberkorper(xx_act, u_act)*T_matrix
        # if xx_act[0] <= 0:
        #     xx_act[0] = 0
        # if xx_act[0] <= 0 and xx_act[1] <= 0:
        #     xx_act[1] = 0
        xx_print=np.append(xx_print,xx_act[0])
        T_act += T

        #print(T_act)
        T_print=np.append(T_print,[T_act])
        if T_act > T_tra_value:
            u_act = u_input_T_func(xx_act[0], xx_act[1], xx_act[2])
            h_r_print=np.append(h_r_print,h_r0+y_poly_func(T_tra_value,T_tra_value))
        else:
            u_act = u_input_func(T_act,xx_act[0],xx_act[1],xx_act[2])
            h_r_print = np.append(h_r_print, h_r0+y_poly_func(T_tra_value, T_act))
        # if xx_act[0] <= 0:
        #     xx_act[0] = 0
        #     u_act=255
        # if xx_act[0] <= 0 and xx_act[1] <= 0:
        #     xx_act[1] = 0
        #     u_act = 255

        print(u_act)

    #print(xx_print)
    plt.plot(T_print,xx_print)
    plt.plot(T_print, h_r_print)
    plt.show()

