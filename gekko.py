# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 15:18:15 2023

@author: ISA
"""

from gekko import GEKKO
import numpy as np

m = GEKKO()

m.time = np.linspace(0,500,500+1)

beta = m.MV(value=1, ub=20, lb=0)
beta.status = 1


N = m.Var(value=1)
I0 = m.Var(value=0.01)
R0 = m.Var(value=0)
U0 = m.Var(value=N - I0 - R0)
J0 = m.Var(value=I0)
Lf0 = m.Var(value=0)
Ls0 = m.Var(value=0)
gamma = m.Var(value=365/75)
mu = m.Var(value=1/80)
muTB= m.Var(value=1/6)
sigma= m.Var(value=1/6)
rho = m.Var(value=0.03)
u= m.Var(value=0.88)
v= m.Var(value=0.083)
w = m.Var(value=0.0006)
lamda = m.Var(beta * I0)
clamda = m.Var(0.2 * lamda)
b =  m.Var((mu * (U0 + Lf0 + Ls0 + R0)) + (muTB * I0))


m.Equation(U0.dt()== b - ((lamda + mu) * U0))
m.Equation(Lf0.dt()== (lamda*U0) + ((clamda)*(Ls0 + R0)) - ((u + v + mu) * Lf0))
m.Equation(Ls0.dt()== (u * Lf0) - ((w + clamda + mu) * Ls0))
m.Equation(I0.dt()== w*Ls0 + v*Lf0 - ((gamma + muTB + sigma) * I0) + (rho * R0))
m.Equation(R0.dt()== ((gamma + sigma) * I0) - ((rho + clamda + mu) * R0))
m.Equation(J0.dt()== w*Ls0 + v*Lf0 + (rho * R0))

m.Minimize(J0[-1] - 300/100000)

m.options.IMODE = 6
m.solve()


