import sympy as sp
import symbtools as st
import numpy as np
from matplotlib import pyplot as plt
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
    T=0.1
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

    Trans=sp.Matrix([h,Lfh,Lf2h])
    print(Trans)

    zz = st.symb_vector("z1:4")

    res = sp.solve(Trans-zz, xx)[0]
    x_z_beziehung = st.lzip(xx, res)
    print("res",res)
    print("x_z_beziehung",x_z_beziehung)


    z_dot_tmp = Trans.jacobian(xx)*(f_sys + g_sys*u)

    z_dot = z_dot_tmp.subs(x_z_beziehung)
    z_dot.simplify()
    z_dot
    print(z_dot)

    a0,a1,a2,a3=alpha=sp.symbols('alpha0:4')
    u_input=sp.Symbol('u_{input}')
    w=sp.Symbol('w') #Führungsgröße
    print(alpha,u_input)

    u_input=1/(LgLf2h)*(a0*w-(a0*h+a1*Lfh+a2*Lf2h+a3*Lf3h))
    u_input.subs([(w,0.8)])

    s=sp.Symbol('s')
    nenner=(s+1)*(s+1)*(s+1)
    print(sp.expand(nenner))
    #print(nenner.subs([(a0,10),(a1,2),(a2,3)]))
    #sp.expand(sp.solve(nenner.subs([(a0,10),(a1,2),(a2,3)]),s)[0])

    u_input=sp.simplify(sp.expand(u_input.subs([(w,h_r),(a3,1),(a2,7),(a1,11),(a0,5)])))
    print(u_input)

    xx_act = sp.Matrix([0, 0, 0])
    xx_print =xx_act[0]
    T_print = 0
    u_act=0
    T_act=0
    print(sys_schweberkorper(xx_act,u_act))
    print( xx_act + sys_schweberkorper(xx_act, u_act)*T_matrix )
    for i in range(100):
        xx_act = xx_act + sys_schweberkorper(xx_act, u_act)*T_matrix
        if xx_act[0] <= 0:
            xx_act[0] = 0
        if xx_act[0] <= 0 and xx_act[1] <= 0:
            xx_act[1] = 0
        xx_print=np.append(xx_print,xx_act[0])
        T_act += T
        #print(T_act)
        T_print=np.append(T_print,[T_act])
        u_act = u_input.subs([(x1, xx_act[0]), (x2, xx_act[1]), (x3, xx_act[2])])
        if xx_act[0] <= 0:
            xx_act[0] = 0
            u_act=255
        if xx_act[0] <= 0 and xx_act[1] <= 0:
            xx_act[1] = 0
            u_act = 255

        print(u_act)

    #print(xx_print)
    plt.plot(T_print,xx_print)
    plt.show()