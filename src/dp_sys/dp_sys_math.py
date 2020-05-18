import numpy as np
import math

from src.dp_sys.utils import map2pi


def A(th1, th2, al1, al2, m1, m2, l1, l2):
    a11 = -m2 * l1 * l2 * al2 * math.sin(th1 - th2)
    a12 = m2 * l1 * l2 * al2 * math.sin(th1 - th2)
    a21 = -m2 * l1 * l2 * al1 * math.sin(th1 - th2)
    a22 = m2 * l1 * l2 * al1 * math.sin(th1 - th2)
    return np.array([
        [a11, a12],
        [a21, a22]
    ])


def B_inv(th1, th2, al1, al2, m1, m2, l1, l2):
    b11 = (m1 + m2) * l1 * l1
    b12 = m2 * l1 * l2 * math.cos(th1 - th2)
    b21 = m2 * l1 * l2 * math.cos(th1 - th2)
    b22 = m2 * l2 * l2

    v = (1 / (b11 * b22 - b12 * b21)) * np.array([
        [b22, -b12],
        [-b21, b11]
    ])

    return v

def partial_th1(th1, th2, al1, al2, m1, m2, l1, l2, g):
    return -m2 * l1 * l2 * al1 * al2 * math.sin(th1 - th2) - g * (m1 + m2) * l1 * math.sin(th1)


def partial_th2(th1, th2, al1, al2, m1, m2, l1, l2, g):
    return m2 * l1 * l2 * al1 * al2 * math.sin(th1 - th2) - g * m2 * l2 * math.sin(th2)


def F(th1, th2, al1, al2, m1, m2, l1, l2, g):
    I = np.eye(2, 2)
    Z = np.zeros((2, 2))
    mat_B_inv = B_inv(th1, th2, al1, al2, m1, m2, l1, l2)
    mat_A = A(th1, th2, al1, al2, m1, m2, l1, l2)
    par_th1 = partial_th1(th1, th2, al1, al2, m1, m2, l1, l2, g)
    par_th2 = partial_th2(th1, th2, al1, al2, m1, m2, l1, l2, g)
    tb = np.matmul(-mat_B_inv, mat_A)
    scale_mat = np.matmul(np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 0, tb[0][0], tb[0][1]],
        [0, 0, tb[1][0], tb[1][1]],
    ]), np.array([th1, th2, al1, al2]).transpose())

    b = np.matmul(mat_B_inv, np.array([par_th1, par_th2]).transpose())
    offset_mat = np.array([0, 0, b[0], b[1]])

    return scale_mat + offset_mat


def runge_kutt_4(th1, th2, al1, al2, m1, m2, l1, l2, g, dt):
    state = [th1, th2, al1, al2]
    # stage 1
    state1 = state
    k1 = F(state1[0], state1[1], state1[2], state1[3], m1, m2, l1, l2, g) * dt
    # stage 2
    state2 = k1 / 2 + state
    k2 = F(state2[0], state2[1], state2[2], state2[3], m1, m2, l1, l2, g) * dt

    # stage 3
    state3 = k2 / 2 + state
    k3 = F(state3[0], state3[1], state3[2], state3[3], m1, m2, l1, l2, g) * dt

    # stage 4
    state4 = k3 + state
    k4 = F(state4[0], state4[1], state4[2], state4[3], m1, m2, l1, l2, g) * dt

    # calc y
    state_dt = state + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6

    # 限制在 -pi 到 pi 之间
    state_dt[0] = map2pi(state_dt[0])
    state_dt[1] = map2pi(state_dt[1])

    return state_dt
