import numpy as np
import scipy.sparse
import scipy.sparse.linalg

def solve_qp(R, Q, q, A, B):
    """
    minimize (1/2) u^T R u + (1/2)x^T Q x - q^T x
    subject to Ax + Bu = 0
    """
    BRB = B @ scipy.sparse.linalg.spsolve(R, B.T)
    A_sys  = scipy.sparse.bmat([[Q, A.T], [A, -BRB]], format='csc')
    b_sys  = np.concatenate([q, np.zeros(A.shape[0])], axis=0)
    x_sys  = scipy.sparse.linalg.spsolve(A_sys, b_sys)
    X_star = x_sys[0:A.shape[1]]
    return X_star

def make_newton_solver(Q, A, BRB, G, w):
    N_x = A.shape[1]
    N_z = w.shape[0]
    w2inv = w**(-2)
    W2inv = scipy.sparse.spdiags(w2inv, [0], N_z, N_z)
    S     = Q + (G.T @ W2inv @ G)
    sys_A = scipy.sparse.bmat([[S, A.T], [A, -BRB]], format='csc')
    sys_B = scipy.sparse.linalg.splu(sys_A)
    def solver(rx, ry, rz, rs):
        rz2 = rz - (w * rs)
        rx2 = rx + G.T @ (w2inv * rz2)
        sys_x = sys_B.solve(np.concatenate([rx2, ry]))
        dx = sys_x[0:N_x]
        dy = sys_x[N_x:]
        dz = w2inv * (G @ dx - rz2)
        ds = w * (rs - w * dz)
        return dx, dy, dz, ds
    return solver

def calc_alpha(z, s, dz, ds):
    x  = np.concatenate([z, s], axis=0)
    dx = np.concatenate([dz, ds], axis=0)
    cond = (dx < 0)
    if cond.shape[0] == 0:
        return 1.0
    else:
        return min(1.0, 0.95 * np.min(- x[cond] / dx[cond]))

def solve_qp_with_inequality(R, Q, q, A, B, G, h):
    """
    minimize (1/2) u^T R u + (1/2)x^T Q x - q^T x
    subject to Ax + Bu = 0
               Gx + s  = h, s >= 0
    """
    N_x = A.shape[1]
    N_y = A.shape[0]
    N_z = G.shape[0]
    N_s = N_z

    x = np.zeros((N_x, ))
    y = np.zeros((N_y, ))
    z = np.ones((N_z, ))
    s = np.ones((N_z, ))
    
    X_ZERO = np.zeros_like(x)
    Y_ZERO = np.zeros_like(y)
    Z_ZERO = np.zeros_like(z)
    BRB = B @ scipy.sparse.linalg.spsolve(R, B.T)

    success  = False
    LOOP_MAX = 20
    for loop in range(LOOP_MAX):
        rx  = - ( (Q @ x) - q + (A.T @ y) + (G.T @ z) )
        ry  = - ( (A @ x) - (BRB @ y) )
        rz  = - ( (G @ x) + s - h )
        gap = np.dot(z, s)
        terminal_cond = ( (np.sqrt(np.mean(rx**2)) < 1e-10) and
                          (np.sqrt(np.mean(ry**2)) < 1e-10) and
                          (np.sqrt(np.mean(rz**2)) < 1e-10) and
                          (gap / N_z < 1e-10) )
        if terminal_cond:
            success = True
            break

        sqrt_s = np.sqrt(s)
        sqrt_z = np.sqrt(z)
        w = sqrt_s / sqrt_z
        lambda_ = sqrt_s * sqrt_z
        newton_solver = make_newton_solver(Q, A, BRB, G, w)

        dx1, dy1, dz1, ds1 = newton_solver(rx, ry, rz, rs=-lambda_)
        alpha = calc_alpha(z, s, dz1, ds1)
        mu_prev   = gap / N_z
        mu_hat    = np.dot(z + alpha*dz1, s + alpha*ds1) / N_z
        mu_target = mu_prev * (mu_hat / mu_prev)**3

        rs = (mu_target - (dz1 * ds1)) / lambda_
        dx2, dy2, dz2, ds2 = newton_solver(X_ZERO, Y_ZERO, Z_ZERO, rs)

        dx = dx1 + dx2
        dy = dy1 + dy2
        dz = dz1 + dz2
        ds = ds1 + ds2
        alpha = calc_alpha(z, s, dz, ds)
        x += alpha * dx
        y += alpha * dy
        z += alpha * dz
        s += alpha * ds

    if not success:
        raise RuntimeError('Not solved')
    X_star = x
    # print(np.max(G @ X_star))
    return X_star