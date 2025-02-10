#Purpose: To give an accurate solution to the Kepler Equation
#By Shayaan Hooda - https://github.com/crazystuffxyz
#Please keep this here in the code, or give contribution. Thanks!
#Also, please accept my poor coding skills, I just recently started learning Python

import mpmath as mp
final_dps1=10100 #Default accuracy
#Main function that does the same as the desmos link (https://www.desmos.com/calculator/xcjxasjdfr)
def solve_for_x_approximation(y, m, final_dps=final_dps1, max_iterations=10000):
    previous_accuracy = mp.inf
    #Initial guess based on an "accurate" approximation
    initial_guess = mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)+mp.mpf(m)*mp.sin(mp.mpf(y)))))))))))))))))))))))))))))))))))))))))
    rough_dps = 50
    with mp.workdps(rough_dps):
        x = mp.mpf(initial_guess)
        y_low = mp.mpf(y)
        m_low = mp.mpf(m)
        tol_low = mp.mpf('1e-40')
        for i in range(max_iterations):
            fx = x - m_low * mp.sin(x) - y_low
            dfx = 1 - m_low * mp.cos(x)
            if mp.fabs(dfx) < tol_low:
                break
            x_new = x - fx/dfx
            if mp.fabs(x_new - x) < tol_low:
                break
            x = x_new
    with mp.workdps(final_dps):
        x = mp.mpf(x)
        y_high = mp.mpf(y)
        m_high = mp.mpf(m)
        tol_high = mp.mpf(f'1e-{final_dps - 100}')
        for i in range(max_iterations): # Iterate by using Newton's Method
            fx = x - m_high * mp.sin(x) - y_high #f(x)
            dfx = 1 - m_high * mp.cos(x) # f'(x)
            if mp.fabs(dfx) < tol_high:
                print("Derivative too small; Newton's method may fail to converge. Returning x value.")
                return x
            x_new = x - fx/dfx
            accuracy = mp.ceil(-mp.log10(mp.fabs((x_new-m_high*mp.sin(x_new))-mp.mpf(y)))) #Accuracy based on base 10
            if accuracy == accuracy:
                print(f"On iteration {i+1} recieved accuracy of {accuracy}")
            if accuracy > final_dps:
                print(f"Converged in {i+1} high-precision iterations.")
                final_dps1 = final_dps
                return x_new
            if accuracy == previous_accuracy:
                #Sometimes, the code caps at an accuracy. This will see when it caps and if it does, it will return the maximum accuracy. Don't worry, the accuracy is usually above 10000 or near it.
                print(f"Reached maximum accuracy for values ({accuracy}) in {i+1} high-precision iterations.")
                final_dps1 = final_dps
                return x_new
            if accuracy != accuracy:
                # If accuracy is NaN
                print(f"Error finding value. You probably put invalid numbers (such as infinity, strings, complex numbers, etc.) for the values, or the program doesn't work for the numbers. The numbers you put for m is {m} and you put {y} for y.")
                break
            previous_accuracy = accuracy
            x = x_new
    if accuracy == accuracy:
        print("Max iterations reached; the solution may not be accurate.")
        return x

# Main Code:
if __name__ == "__main__":
    y_value = -0.4
    m_value = 1
    solution_x = solve_for_x_approximation(y_value, m_value)
    if solution_x:
        print(f"The solution x is approximately {mp.nstr(solution_x, final_dps1)}")
