import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm
import pandas as pd
import matplotlib as mpl

amounts_data = pd.read_excel("data3/Agustina_German.xlsx")
dd_student_name = amounts_data["Student"].values[0]
dd_amounts_data = amounts_data[amounts_data["Student"] == dd_student_name]
dd_polymer = dd_amounts_data["Polymer"]
print(dd_polymer)



class ramanspectrum():
    x = []
    y = []
    
    x_unit = "wavenumbers"
    y_unit = "counts"
    
    def __init__(self,x,y) -> None:
        
        self.x = np.array(x, dtype=np.float64)
        self.y = np.array(y, dtype= np.float64)
        
        pass
    
    
    
class macroramanspectrum():
    
    

    
    def __init__(self,name:str, kinetic_path = "", thermal_path = "") -> None:
        self.name = name
        if kinetic_path != "" :
            self.data = pd.read_csv(kinetic_path,delimiter = ";",dtype=str,decimal=",")
            self.kinetic_data_process(self.data)
        if thermal_path != "" :
            self.data3 = pd.read_csv(thermal_path,delimiter = ";",dtype=float,decimal=",")
            self.thermal_data_process(self.data3)
        pass    
        
    def thermal_data_process(self,dataframe):
        self.thermal_spectra = {}
        self.labels = ["Glass Spectrum","Before heating", "After heating to 80 C", "After heating to 100 C", "After heating to 120 C"]
        for i in range(len(self.labels)+1):
            if i == 0: x = dataframe.iloc[:,i].values
            else: self.thermal_spectra[self.labels[i-1]] = ramanspectrum(x,dataframe.iloc[:,i].values)
        
        pass 
    
    def kinetic_data_process(self, df: pd.DataFrame):
        self.kinetic_spectra = []
        df.iloc[:, 0] = df.iloc[:, 0].str.replace(',', '.', regex=False)

        # Convert first column to float with 4 decimal places
        df.iloc[:, 0] = pd.to_numeric(df.iloc[:, 0], errors='coerce').round(4)

        # Convert the entire DataFrame to float
        df = df.apply(lambda x: pd.to_numeric(x, errors='coerce'))

        for i in range(len(df.keys())):
            if i == 0: x = df.iloc[:,i].values
            else: self.kinetic_spectra.append(ramanspectrum(x,df.iloc[:,i].values))
            
        
        
        # baseline
        
        
        for i in range(len(self.kinetic_spectra)):
            id = self.find_nearest_id(self.kinetic_spectra[i].x,1080)
            id2 = self.find_nearest_id(self.kinetic_spectra[i].x,1770)
            a = (-self.kinetic_spectra[i].y[id]+self.kinetic_spectra[i].y[id2])/(-self.kinetic_spectra[i].x[id]+self.kinetic_spectra[i].x[id2])
            b = self.kinetic_spectra[i].y[id]
            self.kinetic_spectra[i].y -= -a*(self.kinetic_spectra[0].x[id]-self.kinetic_spectra[i].x) +b
        #normalise
        id_norm = self.find_nearest_id(self.kinetic_spectra[0].x, 1155 )
        for i in range(len(self.kinetic_spectra)):
            self.kinetic_spectra[i].y *= self.kinetic_spectra[0].y[id_norm]
            self.kinetic_spectra[i].y /= float(self.kinetic_spectra[i].y[id_norm])
        n=1
        for i in range(len(self.kinetic_spectra)):
            n = max(np.max(self.kinetic_spectra[i].y),n)
        for i in range(len(self.kinetic_spectra)):
            self.kinetic_spectra[i].y /= n
            
            
    def display_kinetic_plot(self, ax):
        colours = plt.colormaps["viridis"](np.linspace(1,0,len(self.kinetic_spectra)))
        for i in range(len(self.kinetic_spectra)):
            cur_spec = self.kinetic_spectra[i]
            ax.plot(cur_spec.x, cur_spec.y,color = colours[i])
        ax.set(xlabel = "Wavenumbers [1/cm]", ylabel = "Normalised Intensity", title = f"Kinetic Raman Spectra of Ratio {self.name}")
    
    @staticmethod
    def find_nearest_id(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx
    
    def reaction_profile_plot(self, ax, min, max, c = False):
        self.integrals = []
        
        
        
        for i in range(len(self.kinetic_spectra)):
            cur_spec = self.kinetic_spectra[i]
            self.min_id = self.find_nearest_id(cur_spec.x, min)
            self.max_id = self.find_nearest_id(cur_spec.x, max)
            self.integrals.append(np.sum(cur_spec.y[self.min_id:self.max_id]-cur_spec.y[self.max_id]))
        
        self.integrals /= self.integrals[0]
        
        conversion = 1-(self.integrals/self.integrals[0])
        ax.set(ylabel = f"{'Peak integral' if not c else 'Conversion Factor'}", xlabel = "Time [minute]", title = f"{'conversion' if c else 'reaction'} profile (using integral of peaks between 1620 and 1660 wavenumbers 1/cm)")    
        start = 4.4 if self.name == "B" else 0 
        if c:
            ax.plot(np.linspace(start,start+90,len(self.kinetic_spectra)),conversion, label = f"Conversion Profile of Ratio {self.name}")
        else:
            ax.plot(np.linspace(start,start+90,len(self.kinetic_spectra)),self.integrals, label = f"Reaction Profile of Ratio {self.name}")
        pass
    
    
    
    def thermal_plot(self, ax):
        colours = plt.colormaps["PuRd"](np.linspace(0.2,0.8,len(self.thermal_spectra.keys())))
        i = 0
        labels = []
        c80 = self.thermal_spectra[self.labels[2]]
        c100 = self.thermal_spectra[self.labels[3]]
        c120 = self.thermal_spectra[self.labels[4]]
        c0 = self.thermal_spectra[self.labels[1]]
        cglass = self.thermal_spectra[self.labels[0]]
        
        specs = [c0,c80,c100,c120]
        
        for spec in specs:
            spec.y -= cglass.y
        
        id = self.find_nearest_id(c0.x,1525)
        for i in range(len(specs)):
            cur_spec = specs[i]
            cur_spec.y -=  c0.y[id]
            
            ax.plot(cur_spec.x, cur_spec.y,color = colours[i], label = self.labels[i+1])
            ax.legend()
            i+=1
        ax.set(xlabel = "Wavenumbers [1/cm]", ylabel = "Intensity [-]", title = f"Thermal Raman Spectra of Ratio {self.name}")
        pass


fig,ax = plt.subplots(1,1, figsize = (8,8), dpi=70)

letters = ["A","B", "C", "D", "E", "F", "G", "H", "I"]

dataes = [f"data3/S{letter}_kinetic.csv" for letter in letters]



def kinetic(letter):
    macroB = macroramanspectrum(letter, kinetic_path=f"data3/S{letter}_kinetic.csv")
    macroB.display_kinetic_plot(ax)


def thermal():
    macrob = macroramanspectrum("E", kinetic_path="", thermal_path="data3/SE_Thermal.csv")
    macrob.thermal_plot(ax)
    ax.set(xlim=(1530,1680))
    fig.tight_layout()
    
def profiles(c=False):
    macros = []
    for i in range(len(dataes)):
        macros.append(macroramanspectrum(letters[i],kinetic_path=dataes[i]))
        macros[i].reaction_profile_plot(ax,1620,1660,c=c)
    ax.legend()
    
def allthermal(ax):
    df = pd.read_csv("data3/All_120_after.csv",delimiter = ";",dtype=float,decimal=",")
    spectra: list[ramanspectrum] = []
    for i in range(len(df.columns)):
        if i == 0: x = df.iloc[:,i]
        else: spectra.append(ramanspectrum(x,df.iloc[:,i]))
    spectra = spectra[:-1]
    
    for i in range(len(spectra)):
       spectrum = spectra[i]
       ax.plot(spectrum.x,spectrum.y, label = f"Ratio {letters[i]}")
    ax.set(xlabel = f"Wavenumbers [1/cm]", ylabel = "Intensity", title = "Raman Spectra after 120 C")
    ax.legend()
    
def single_profile():
    macrob = macroramanspectrum("E", kinetic_path="data3/SE_kinetic.csv")
    macrob.reaction_profile_plot(ax,1620,1660,c=True)
    ax.legend()

single_profile()

# norm = mpl.colors.Normalize(vmin=90, vmax=0) 
# # creating ScalarMappable 
# sm = plt.cm.ScalarMappable(cmap=plt.colormaps["viridis_r"], norm=norm) 
# sm.set_array([]) 
# plt.colorbar(sm,ax=ax,ticks=np.linspace(90, 0, 5),label = "Time [M]")

plt.tight_layout()
plt.show()
    
    