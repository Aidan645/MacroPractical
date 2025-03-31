import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

def linear_function(x,a,b):
    return a*x +b

def r_squared(x,y,popt):
    residuals = y - linear_function(x,*popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y-np.mean(y))**2)
    r_squared = 1-(ss_res/ss_tot)
    return r_squared


labels = ["Methyl 3-phenylpropionate", "4-phenyl-2-butanone", "2'-phenylethyl Acetate"]

def show_plots():
    for i in range(len(labels)):
        fig,ax = plt.subplots(1,1, figsize = (5,5))
        data = pd.read_csv(f"{labels[i]}.csv")
        x = data[f"c(M)"]
        y = data["Area"]
        popt, _ = curve_fit(linear_function,x,y)
        r_sq = r_squared(x,y,popt)
        
        ax.scatter(x,y, label = f"Experimental data")
        ax.plot(x,y, linestyle = "dashed")
        ax.plot(x,linear_function(x,*popt),label = f"Fitted Calibration Curve, {popt[0]:.05e}x + {popt[1]:.05e}, R^2 = {r_sq:.02f}")
        ax.set(title = f"Calibration Curve Data for {labels[i]}", xlabel = "c(M)", ylabel = "Area")
        ax.legend()
        plt.tight_layout()
        plt.show()

def print_coeffs():
    for i in range(len(labels)):
        data = pd.read_csv(f"{labels[i]}.csv")
        x = data[f"c(M)"]
        y = data["Area"]
        popt, _ = curve_fit(linear_function,x,y)
        r_sq = r_squared(x,y,popt)
        
        print(f"{labels[i]}: \t{'' if i!=1 else '        '} f(x) = {popt[0]}x \t+ {popt[1]}, \t r^2 = {r_sq}")

print_coeffs()