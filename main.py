""""
Group: Chaya Mizrachi ID: 214102584,
        Yael Siman-Tov ID:325181295,
        Linoy Nisim Pur ID: 324029685
Source Git: lihiSabag https://github.com/lihiSabag/Numerical-Analysis-2023.git
GitHub of this project: https://github.com/YaelSimanTov/Analiza-nomerit-work3.git
 """




import numpy as np
import sympy as sp
import math
from sympy.utilities.lambdify import lambdify
from math import e

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[90m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    # Background colors:
    GREYBG = '\033[100m'
    REDBG = '\033[101m'
    GREENBG = '\033[102m'
    YELLOWBG = '\033[103m'
    BLUEBG = '\033[104m'
    PINKBG = '\033[105m'
    CYANBG = '\033[106m'


x = sp.symbols('x')


def max_steps(a, b, err):

    startPoint=a
    endPoint=b
    exp = pow(10, -10)
    if endPoint - startPoint == 0:
        return 100
    s=int(((-1) * math.log(exp / (endPoint - startPoint), e)) / math.log(2, e))
    return s


def newton_rephson(P, test_range, epsilon):
    """
    This function calculates polynom root (P) by using Newton's formula: xr+1= xr-(f(x) / f'(x))
    :param P: Polynom
    :type P: function
    :param test_range: A range of values in which we will check if a root exists
    :type test_range: list with 2 values
    :param epsilon: Another stopping condition - when the difference of x between the 2 iterations is less than or equal
     to epsilon
    :type epsilon:float
    :return: None
    """

    p_derivative = sp.diff(P, x)  # p_derivative = p'(x)
    p_derivative = lambdify(x, p_derivative)  # convert to an equivalent NumPy function to perform calculations
    P = lambdify(x, P) # convert to an equivalent NumPy function to perform calculations

    if len(test_range) == 0:
        x0 = 0.5
        max_step = 50

    elif len(test_range) == 1:
        x0 = test_range[0]


    else:
        x0 = (test_range[-1] + test_range[0]) / 2  # x0 is the first value that are guessed and placed in a Newton's formula
        max_step = max_steps(test_range[0], test_range[-1], epsilon)

    if P(x0) == 0:
        #print('The root is ', x0, '. The number of iterations is ', count_iteration)
        return x0 , 0

    for i in range(0, max_step + 1):

        """
        if abs(p_derivative(x0)) == 0 or abs(p_derivative(x0)) < 1:

            print('The convergence condition is not met')
            return
        """

        if abs(p_derivative(x0)) < epsilon:
            print('The convergence condition (|f\'(x)| < 1) is not met at x =', x0, ', therefore in the range:',
                  test_range, 'there may not be roots or the method may converge slowly.\n\t')
            break


        x_root = x0 - ((P(x0)) / (p_derivative(x0)))  # x value of the next iteration, Newton's formula
        if abs(P(x_root)) < epsilon or abs(x_root-x0) <= epsilon:

           #print('The root is ', round(x_root, 6), '\n f( ', x_root, ')= ',0, '\nThe number of iterations is ', i+1)
           return x_root, i

        x0=x_root

    return None, 0

def secant_method(f, test_range, TOL=1e-6, N=50):
    f = lambdify(x, f) # convert to an equivalent NumPy function to perform calculations

    x0=test_range[0]
    x1=test_range[1]
    #print("{:<10} {:<15} {:<15} {:<15}".format("Iteration", "xo", "x1", "p"))
    for i in range(1,N+1):
        if f(x1) - f(x0) == 0:
            print( " 'Secant' method cannot continue!\n\t")
            break

        p = x0 - ((x1 - x0) / (f(x1) - f(x0))) * f(x0)

        if abs(p - x1) < TOL:
            return p , i # Procedure completed successfully
        #print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f}".format(i, x0, x1,p))
        x0 = x1
        x1 = p
    return None, 0




if __name__ == '__main__':
    # if you want to try other functions please change here f

    f = x ** 3 - x - 1
    #f = -2 * sp.log(x) + 2
    #f= sp.sin(x ** 2) - x ** 3 - 1
    #f = x**2-4*x
    #f = x**3-sp.cos(x)
    #f = x ** 4 + x ** 3 - 3 * x ** 2
    #f= sp.sin(x ** 2) - x ** 3 - 1

    epsilon = 0.0001
    domain = [-5, 6]
    root = None

    print("Finding roots of the equation f(X) :\n\t")
    choice = int(input(
        "Which method Which methods do you want? \n\t1.Newton Raphson\n\t2.Secant Method\n"))
    if choice == 1:
        root, num_iter= newton_rephson(f, domain, epsilon)
        if root is None:
           print("Try found a root using 'Secant' method !\n\t")
           root, num_iter= secant_method(f, domain, epsilon)

    elif choice == 2:
        #root, num_iter= secant_rephson(f, domain, epsilon)
        root, num_iter= secant_method(f, domain, epsilon )
    else:
        print(bcolors.FAIL, "Invalid input", bcolors.ENDC)

    if root is not None:
        print('The root is ', round(root, 6), '\n f( ', root, ')= ', 0, '\nThe number of iterations is ', num_iter)

    elif choice == 1 or choice == 2:
        print("The equation f(x) has not an approximate root \n\t")