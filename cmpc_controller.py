import numpy as np

import cvxpy as cp


def calc_Jacobian(x, u, param):


    L_f = param["L_f"]

    L_r = param["L_r"]

    dt   = param["h"]


    psi = x[2]

    v   = x[3]

    delta = u[1]

    a   = u[0]


    # Jacobian of the system dynamics

    A = np.zeros((4, 4))

    B = np.zeros((4, 2))


    #############################################################################

    #                    TODO: Implement your code here                         #

    #############################################################################


    # A = ...

    # B = ...


    #given the formula in TrajectoryTracking.ipynb for A

    A[0, 0] = 1

    A[0, 1] = 0

    A[0, 2] = -dt * v * np.sin(psi + np.arctan((L_r * np.arctan(delta)) / (L_f + L_r)))

    A[0, 3] = dt * np.cos(psi + np.arctan((L_r * np.arctan(delta)) / (L_f + L_r)))

    A[1, 0] = 0

    A[1, 1] = 1

    A[1, 2] = dt * v * np.cos(psi + np.arctan((L_r * np.arctan(delta)) / (L_f + L_r)))

    A[1, 3] = dt * np.sin(psi + np.arctan((L_r * np.arctan(delta)) / (L_f + L_r)))

    A[2, 0] = 0

    A[2, 1] = 0

    A[2, 2] = 1

    A[2, 3] = (dt * np.arctan(delta)) / (((L_r**2 * np.arctan(delta)**2) / (L_f + L_r)**2 + 1)**(1/2) * (L_f + L_r))

    A[3, 0] = 0

    A[3, 1] = 0

    A[3, 2] = 0

    A[3, 3] = 1



    #similar job for B

    B[0, 0] = 0

    B[0, 1] = -(dt * L_r * v * np.sin(psi + np.arctan((L_r * np.arctan(delta)) / (L_f + L_r)))) / ((delta**2 + 1) * ((L_r**2 * np.arctan(delta) ** 2) / (L_f + L_r) ** 2 + 1) * (L_f + L_r))

    B[1, 0] = 0

    B[1, 1] = (dt * L_r * v * np.cos(psi + np.arctan((L_r * np.arctan(delta)) / (L_f + L_r)))) / ((delta ** 2 + 1) * ((L_r**2 * np.arctan(delta) ** 2) / (L_f + L_r) ** 2 + 1) * (L_f + L_r))

    B[2, 0] = 0

    B[2, 1] = (dt * v) / ((delta ** 2 + 1) * ((L_r ** 2 * np.arctan(delta) ** 2) / (L_f + L_r)**2 + 1) ** (3/2) * (L_f + L_r))

    B[3, 0] = dt

    B[3, 1] = 0


    #############################################################################

    #                            END OF YOUR CODE                               #

    #############################################################################


    return [A, B]


def LQR_Controller(x_bar, u_bar, x0, param):

    len_state = x_bar.shape[0]  #21

    len_ctrl  = u_bar.shape[0]  #20

    dim_state = x_bar.shape[1]  #4

    dim_ctrl  = u_bar.shape[1]  #2


    n_u = len_ctrl * dim_ctrl  #40

    n_x = len_state * dim_state #84

    n_var = n_u + n_x #124



    n_eq  = dim_state * len_ctrl # dynamics

    n_ieq = dim_ctrl * len_ctrl  # input constraints


    

    #############################################################################

    #                    TODO: Implement your code here                         #

    #############################################################################


    # define the parameters

    # Q = np.eye(4)  * ...

    # R = np.eye(2)  * ...

    # Pt = np.eye(4) * ...

    Q = np.eye(4) * 5

    R = np.eye(2) * 1

    Pt = np.eye(4) * 5


    # define the cost function

    # P = ...

    # q = ...

    '''

    know that num_x = num_u + 1, so sN-Pt is correspond to last value of x; sk-Q is first num_u of x;

    uk-R is num_u

    create a matrix P that has 4xnum_x + 2xnum_u in dimension.

    In diagonal, the first num_ux4 is Q, the num_u + 1 x4 is Pt, the last num_ux2 is R

    q is still zero because no existence of single variable

    '''

    P = np.eye(n_var) #124 = 20*4 + 4 + 20*2


    #for first Q, num_u * 4, dim_state = 4

    for i in range(len_ctrl):

        P[i * dim_state:(i+1) * dim_state, i * dim_state:(i+1) * dim_state] = Q

    

    #for Pt in num_u + 1

    P[len_ctrl * dim_state:(len_ctrl+1) * dim_state, len_ctrl * dim_state:(len_ctrl+1) * dim_state] = Pt


    #for rest of num_u in R

    for i in range(len_ctrl):

        P[n_x + i*dim_ctrl:n_x + (i+1)*dim_ctrl, n_x + i*dim_ctrl:n_x + (i+1)*dim_ctrl] = R

    

    

    # define the constraints

    # A = ...

    # b = ...


    #for each Ak*delta_sk + Bk*delta_uk = delta_sk1, there will be len_ctrl (20) equations

    #for each equation, the Ab matrix need to incorporate 4x4A and 4x2B, so A.shape[0] should be 4x20 = 80

    #because A@x = b, for x in shape of n_var, A.shape[1] = n_var, shape of b is n_eq


    #when adding second constraint

    new_eq = dim_state * (len_ctrl + 1)

    A = np.zeros((n_eq, n_var))

    b = np.zeros(n_eq)

    


    for i in range(len_ctrl):

        #get Jacobian

        Ak, Bk = calc_Jacobian(x_bar[i, :], u_bar[i, :], param)


        #in Ax=b, so Ak*delta_sk + Bk*delta_uk - delta_sk1 = 0

        #Assuming A is a two part matrix, the first 80x84 is for delta_xk in diagonal, the rest 80x40 is for delta_uk in diagonal

        #so that each xk and uk will be multiplied by their corresponding Ak and Bk


        A[i * dim_state:(i+1) * dim_state, i * dim_state:(i+1) * dim_state] = Ak

        A[i * dim_state:(i+1) * dim_state, (i+1) * dim_state:(i+2) * dim_state] = -1 * np.eye(dim_state)

        A[i * dim_state:(i+1) * dim_state, n_x + i * dim_ctrl:n_x + (i+1) * dim_ctrl] = Bk


    #second constraint at the end of Ab matrix

    #delta_s_init = x0 - x_bar[0, :]

    #for the first s0, the A = 1 so 1*delta_s0 = delta_s_init

    #A[84, 124] [80]



    #adding second constraint separately to constraint

    #A[n_eq:, :dim_state] = np.eye(dim_state)

    #b[n_eq:] = delta_s_init



    #[Ak -1 0 0 0 0 0 0 0 0 ... Bk 0 00 0 0 ....]

    #[0 Ak -1...] (20*4, 124)

    #[np.eye 0 00 0 0 0 0 0 0 0 00], b = [delta_s_init  0 0 0 0 0 ...]


    # Define and solve the CVXPY problem.

    # x = cp.Variable(n_var)

    # prob = cp.Problem(...)

    # prob.solve(verbose=False, max_iter = 10000)

    x = cp.Variable(n_var)


    #x - x_bar 21*4, u_bar 20*2

    #x [84] + [40]


    #cost function is 1/2xTPx

    objective = cp.Minimize(1/2 * cp.quad_form(x, P))


    #constraint is Ax=b

    constraints = [A @ x == b]

    constraints.append(x[0:dim_state] == (x0 - x_bar[0, :]))


    prob = cp.Problem(objective, constraints)

    prob.solve(verbose = False, max_iter = 10000)



    #############################################################################

    #                            END OF YOUR CODE                               #

    #############################################################################


    u_act = x.value[n_x:n_x + dim_ctrl] + u_bar[0, :]

    return u_act

    


def CMPC_Controller(x_bar, u_bar, x0, param):

    len_state = x_bar.shape[0]

    len_ctrl  = u_bar.shape[0]

    dim_state = x_bar.shape[1]

    dim_ctrl  = u_bar.shape[1]

    

    n_u = len_ctrl * dim_ctrl

    n_x = len_state * dim_state

    n_var = n_u + n_x


    n_eq  = dim_state * len_ctrl # dynamics

    n_ieq = dim_ctrl * len_ctrl # input constraints


    a_limit = param["a_lim"]

    delta_limit = param["delta_lim"]

    

    #############################################################################

    #                    TODO: Implement your code here                         #

    #############################################################################

    

    # define the parameters

    # Q = np.eye(4)  * ...

    # R = np.eye(2)  * ...

    # Pt = np.eye(4) * ...g

    Q = np.eye(4) * 20

    R = np.eye(2) * 1

    Pt = np.eye(4) * 20

    

    # define the cost function

    # P = ...

    # q = ...

    

    #the same as LQR

    P = np.eye(n_var)


    #for first Q, num_u * 4, dim_state = 4

    for i in range(len_ctrl):

        P[i * dim_state:(i+1) * dim_state, i * dim_state:(i+1) * dim_state] = Q

    

    #for Pt in num_u + 1

    P[len_ctrl * dim_state:(len_ctrl+1) * dim_state, len_ctrl * dim_state:(len_ctrl+1) * dim_state] = Pt


    #for rest of num_u in R

    for i in range(len_ctrl):

        P[n_x + i*dim_ctrl:n_x + (i+1)*dim_ctrl, n_x + i*dim_ctrl:n_x + (i+1)*dim_ctrl] = R

    

    

    # define the constraints

    # A = ...

    # b = ...

    # G = ...

    # ub = ...

    # lb = ...

    


    new_eq = dim_state * (len_ctrl + 1)

    A = np.zeros((n_eq, n_var))

    b = np.zeros(n_eq)


    #define umin and umax

    u_min = np.zeros(2)

    u_max = np.zeros(2)

    u_min[0] = a_limit[0]

    u_min[1] = delta_limit[0]

    u_max[0] = a_limit[1]

    u_max[1] = delta_limit[1]

    

    A_ie_lb = np.zeros((len_ctrl * dim_ctrl, n_var))

    b_ie_lb = np.zeros(len_ctrl * dim_ctrl)

    A_ie_ub = np.zeros((len_ctrl * dim_ctrl, n_var))

    b_ie_ub = np.zeros(len_ctrl * dim_ctrl)

    #the same as LQR

    for i in range(len_ctrl):

        Ak, Bk = calc_Jacobian(x_bar[i, :], u_bar[i, :], param)



        A[i * dim_state:(i+1) * dim_state, i * dim_state:(i+1) * dim_state] = Ak

        A[i * dim_state:(i+1) * dim_state, (i+1) * dim_state:(i+2) * dim_state] = -1 * np.eye(dim_state)

        A[i * dim_state:(i+1) * dim_state, n_x + i * dim_ctrl:n_x + (i+1) * dim_ctrl] = Bk


        #define the ub and lb

        lb = u_min - u_bar[i, :]

        ub = u_max - u_bar[i, :]


        #define Ab matrix for ub and lb

        #lb

        A_ie_lb[i * dim_ctrl:(i+1) * dim_ctrl, n_x + i * dim_ctrl:n_x + (i+1) * dim_ctrl] = np.eye(dim_ctrl)

        b_ie_lb[i * dim_ctrl:(i+1) * dim_ctrl] = lb


        #ub

        A_ie_ub[i * dim_ctrl:(i+1) * dim_ctrl, n_x + i * dim_ctrl:n_x + (i+1) * dim_ctrl] = np.eye(dim_ctrl)

        b_ie_ub[i * dim_ctrl:(i+1) * dim_ctrl] = ub


        


    #delta_s_init = x0 - x_bar[0, :]

    #A[n_eq:, :dim_state] = np.eye(dim_state)

    #b[n_eq:] = delta_s_init


    # Define and solve the CVXPY problem.

    # x = cp.Variable(n_var)

    # prob = cp.Problem(...)

    # prob.solve(verbose=False, max_iter = 10000)

    x = cp.Variable(n_var)


    #cost function is 1/2xTPx

    objective = cp.Minimize(1/2 * cp.quad_form(x, P))


    #constraint is Ax=b

    constraints = [A @ x == b]

    constraints.append(x[0:dim_state] == (x0 - x_bar[0, :]))



    #inequality constraint

    constraints.append(A_ie_lb @ x >= b_ie_lb)

    constraints.append(A_ie_ub @ x <= b_ie_ub) 



    prob = cp.Problem(objective, constraints)

    prob.solve(verbose = False, max_iter = 10000)


    #############################################################################

    #                            END OF YOUR CODE                               #

    #############################################################################

    

    u_act = x.value[n_x:n_x + dim_ctrl] + u_bar[0, :]

    return u_act
