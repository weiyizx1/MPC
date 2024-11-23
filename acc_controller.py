import numpy as np


def ACC_Controller(t, x, param):

    vd = param["vd"]

    v0 = param["v0"]

    m = param["m"]

    Cag = param["Cag"]

    Cdg = param["Cdg"]


    # cost function and constraints for the QP

    P = np.zeros((2, 2))

    q = np.zeros([2, 1])

    A = np.zeros([5, 2])

    b = np.zeros([5])

    

    #############################################################################

    #                    TODO: Implement your code here                         #

    #############################################################################


    # set the parameters

    # lam = ...

    # alpha = ...

    # w = ...

    lam = 10

    alpha = 0.15

    w = 1000000




    # construct the cost function

    # P = ...

    # q = ...

    #x = [Fw; delta], xTx = [Fw, delta][Fw; delta] = Fw^2 + delta^2; required is Fw^2 + wdelta^2

    #no existence of single [Fw, delta], only Fw^2 and delta^2 exist, so q = [0]

    #[Fw, delta][1 0; 0 w] = [Fw, wdelta]; [Fw, wdelta][Fw; delta] = Fw^2 + wdelta^2

    P[0, 0] = 2

    P[0, 1] = 0

    P[1, 0] = 0

    P[1, 1] = 2 * w

    

    # construct the constraints

    # A = ...

    # b = ...



    #define v and D

    D = x[0]

    v = x[1]



    #define h: h = (v-vd)^2/2

    h = (v - vd) ** 2 / 2

    #define B

    B = D - 1/2 * (v0 - v)**2 / Cdg - 1.8 * v

    #for first constraint

    #if i = (v - vd)/m, j = lamdba*h, then iFw - delta <= -j; when in Ax <= b

    #A = [i -1]; b = -j

    A[0, 0] = (v - vd) / m

    A[0, 1] = -1

    b[0] = -lam * h



    #for second constraint

    '''

    rewrite as -1/m*(1.8+(v-v0)/cdg)*Fw >= -alphaB - (v0-v)

    set i = 1/m*(1.8+(v-v0)/cdg)*Fw, iFw <= alpha*B + (v0-v)

    '''

    A[1, 0] = 1/m * (1.8 + (v - v0) / Cdg)

    A[1, 1] = 0

    b[1] = alpha * B + (v0 - v)


    #for third constraint

    '''

    rewrite it as -Fw/m <= cdg

    '''

    A[2, 0] = -1/m

    A[2, 1] = 0

    b[2] = Cdg


    #for fourth constraint

    '''

    Fw/m <= cag

    '''

    A[3, 0] = 1/m

    A[3, 1] = 0

    b[3] = Cag

    


    #for fifth constraint

    '''

    rewrite as -delta <= 0

    '''

    A[4, 0] = 0

    A[4, 1] = -1

    b[4] = 0



    #############################################################################

    #                            END OF YOUR CODE                               #

    #############################################################################

    

    return A, b, P, q
