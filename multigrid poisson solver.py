import numpy as np
import matplotlib.pyplot as plt

def jacobi_iteration(u, f, h):
    """
    Performs one iteration of the Jacobi method for solving the Poisson equation
    u_xx + u_yy = f on a square domain [0,1] x [0,1] with Dirichlet boundary conditions
    u = 0 on the boundary.

    Parameters:
    u (ndarray): current solution
    f (ndarray): right-hand side
    h (float): grid spacing

    Returns:
    ndarray: updated solution
    """
    m, n = u.shape
    u_new = np.zeros((m, n))

    for i in range(1, m-1):
        for j in range(1, n-1):
            u_new[i, j] = 0.25 * (u[i-1, j] + u[i+1, j] + u[i, j-1] + u[i, j+1] - h**2 * f[i, j])

    return u_new

def residual(u, f, h):
    """
    Computes the residual r = f - Au for the Poisson equation u_xx + u_yy = f
    on a square domain [0,1] x [0,1] with Dirichlet boundary conditions u = 0 on the boundary.

    Parameters:
    u (ndarray): current solution
    f (ndarray): right-hand side
    h (float): grid spacing

    Returns:
    ndarray: residual vector
    """
    m, n = u.shape
    r = np.zeros((m, n))

    for i in range(1, m-1):
        for j in range(1, n-1):
            r[i, j] = f[i, j] - (u[i-1, j] + u[i+1, j] + u[i, j-1] + u[i, j+1] - 4*u[i, j]) / h**2

    return r

def restrict(f):
    # restriction operator from fine to coarse grid
    # assumes indices 0,2,4,... are the coarse points on the fine grid

    n = f.shape[0] // 2
    f_coarse = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            f_coarse[i, j] = (
                f[2*i, 2*j] +
                f[2*i+1, 2*j] +
                f[2*i, 2*j+1] +
                f[2*i+1, 2*j+1]
            ) / 4

    return f_coarse



def prolong(u):
    """
    Prolongs the function u from a coarser grid to a finer grid using linear interpolation.

    Parameters:
    u (ndarray): function on coarser grid

    Returns:
    ndarray: function on finer grid
    """
    m, n = u.shape
    mf = (m-1)*2 + 1
    nf = (n-1)*2 + 1
    uf = np.zeros((mf, nf))

    for i in range(1, m):
        for j in range(1, n):
            uf[2*i-1, 2*j-1] = u[i, j]
            uf[2*i-1, 2*j] = 0.5 * (u[i, j] + u[i, j-1])
            uf[2*i-1, 2*j+1] = u[i, j-1]
            uf[2*i, 2*j-1] = 0.5 * (u[i, j] + u[i-1, j])
            uf[2*i, 2*j] = 0.25 * (u[i, j] + u[i, j-1] + u[i-1, j] + u[i-1, j-1])
            uf[2*i, 2*j+1] = 0.5 * (u[i-1, j-1] + u[i-1, j])

    # boundary conditions
    uf[0, :] = uf[-1, :] = uf[:, 0] = uf[:, -1] = 0

    return uf

def v_cycle(u, f, h, nu1, nu2, max_level):
    """
    Performs a V-cycle to solve the Poisson equation u_xx + u_yy = f on a square domain
    [0,1] x [0,1] with Dirichlet boundary conditions u = 0 on the boundary.

    Parameters:
    u (ndarray): current solution
    f (ndarray): right-hand side
    h (float): grid spacing on current level
    nu1 (int): number of pre-smoothing iterations
    nu2 (int): number of post-smoothing iterations
    max_level (int): maximum number of levels for the multigrid method

    Returns:
    ndarray: updated solution
    """
    m, n = u.shape

    # check if we've reached the coarsest grid
    if m == 3:
        # solve exactly on the coarsest grid
        u = np.linalg.solve(np.array([[4, -1, 0], [-1, 4, -1], [0, -1, 4]]) / h**2, f[1, 1]*np.ones((3, 3)))
        return u

    # pre-smoothing
    for i in range(nu1):
        u = jacobi_iteration(u, f, h)

    # compute the residual
    r = residual(u, f, h)

    # restrict the residual to a coarser grid
    rc = restrict(r)

    # solve the coarse problem exactly
    uc = np.zeros((rc.shape[0], rc.shape[1]))
    for i in range(max_level-1):
        uc = v_cycle(uc, rc, h*2**(i+1), nu1, nu2, max_level-1)

    # interpolate the coarse solution back to the current grid
    uf = np.zeros((m, n))
    mf, nf = uc.shape
    for i in range(1, mf-1):
        for j in range(1, nf-1):
            uf[2*i-1, 2*j-1] = uc[i, j]
            uf[2*i-1, 2*j] = 0.5 * (uc[i, j] + uc[i, j-1])
            uf[2*i-1, 2*j+1] = uc[i, j-1]
            uf[2*i, 2*j-1] = 0.5 * (uc[i, j] + uc[i-1, j])
            uf[2*i, 2*j] = 0.25 * (uc[i, j] + uc[i, j-1] + uc[i-1, j] + uc[i-1, j-1])
            uf[2*i, 2*j+1] = 0.5 * (uc[i-1, j-1] + uc[i-1, j])

    # boundary conditions
    uf[0, :] = uf[-1, :] = uf[:, 0] = uf[:, -1] = 0

    # correct
    for i in range(nu2):
        u = jacobi_iteration(u, f, h)

    return u

def w_cycle(u, f, h, nu1, nu2, max_level, num_cycles):
    # perform num_cycles V-cycles at each level of the hierarchy
    for k in range(max_level-1):
        for i in range(num_cycles):
            u = v_cycle(u, f, h, nu1, nu2, k+1)
        # interpolate the solution to the next level
        u_fine = u.copy()
        u = restrict(u_fine, k+1)
        f = restrict(f, k+1)

    # solve at the coarsest level using a direct solver
    u = direct_solve(f, h)

    # perform num_cycles V-cycles at each level of the hierarchy in reverse order
    for k in range(max_level-2, -1, -1):
        # interpolate the solution from the coarser level
        u_coarse = u.copy()
        u = prolong(u_coarse, k+1)
        # correct the solution using the interpolated error
        u += v_cycle(np.zeros_like(u), f - apply_operator(u, h, k+1), h, nu1, nu2, k+1)
        for i in range(num_cycles-1):
            u += v_cycle(np.zeros_like(u), f - apply_operator(u, h, k+1), h, nu1, nu2, k+1)
    return u



# set up the problem
n = 129  # number of grid points in each direction
h = 1 / (n-1)  # grid spacing
x = np.linspace(0, 1, n)
X, Y = np.meshgrid(x, x)
f = -2 * (np.pi**2) * np.sin(np.pi*X) * np.sin(np.pi*Y)  # right-hand side
u = np.zeros((n, n))  # initial guess

# solve using multigrid with W-cycle
nu1 = 2  # number of pre-smoothing iterations
nu2 = 2  # number of post-smoothing iterations
max_level = int(np.log2(n-1)) + 1  # maximum number of levels
num_cycles = 2  # number of V-cycles to perform at each level
u = v_cycle(u, f, h, nu1, nu2, max_level)

# plot the solution
plt.imshow(u, cmap='jet', origin='lower', extent=[0, 1, 0, 1])
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.show()