"""
Name  : The Method of Evaluation of the Diffuseness of the Wavefield (EDWav Method)
Author: Bo YANG (seism.yang@foxmail.com)
Latest: 2024-03-16   

Copyright (C) 2024 Bo YANG
""" 

__Version__ = '1.1.0'
__Testing_Environment__ = 'Python 3.11'

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import warnings
from numba import jit


class EDWav():
    """ 
    Procedure for evaluating the diffuseness proxy of seismic waveform

    A simple example:
        st=obspy.read()
        wv=st[0].data
        FS=st[0].stats.sampling_rate

        obj=EDWav.evaluate(wv,FS) # evaluation
        obj.plot()                # show result

    Tip: The program uses 'numba@jit' speed up.
         There may be many warnings when calling the program for the first time and it runs slowly. 
         Users can ignore it and it will not appear in the next times.
    """

    def __init__(self):
        """
        Properties of 'EDWav' class
        """
        self.A=0        # Condition A
        self.B=0        # Condition B
        self.C=0        # Condition C
        self.proxy=[]   # Diffuseness proxies of condition A, B and C
        self.freqs=[]   # Frequency of the three conditions
        self.stats=[]   # Some information related to the evaluation results

    def evaluate(wave:np.ndarray=[],FS=1,FB=[],FR=100,NT=1,SF=0.05):
        """
        Function introduction: evaluate the diffuseness proxy of input seismic waveform segment.
        
        Input parameters introduction:
        wave    single channel seismic waveform data, for example: st=obspy.read('test.sac'); wave=st[0].data;
        FS      waveform data sampling rate;
        FB      frequency range to be evaluated, for example: [1, 100], unit: Hz;
        FR      the frequency resolution of the result
        NT      number of multitapers;
        SF      a parameter of sRMS.

        Output parameters introduction: EDWav object instance.
        """
        obj=EDWav()
        NT =int(NT)
        FR =int(FR)

        # Number of waveform sampling points
        npts=len(wave)

        # Frequency range to be analyzed
        if len(FB) != 2:
            FB=[0, FS/2]
        
        # Calculate the minimum length of the window according to the 'FR' and 'FB' entered
        len_of_win=math.ceil(FR*FS/(FB[1]-FB[0])) 
        if npts < len_of_win:
            raise Exception("The waveform length is too small, end the evaluation!")
        
        # The number of windows
        num_of_win=math.floor(npts/len_of_win)    
        if num_of_win < 30:
            warnings.warn("Depending on the settings entered, the final result may not be accurate. Recommended input waveform length npts npts>" + str(math.ceil(30*len_of_win)))

        # Length of window in final
        len_of_win=math.floor(len_of_win+(npts-num_of_win*len_of_win)/num_of_win)

        # Frequencies that this window can analyze
        freqs=np.arange(0, FS/2+FS/len_of_win, FS/len_of_win)
        if freqs[-1]>FS/2:
            freqs=np.delete(freqs,-1)

        # Frequencies to be analyzed
        fin=np.zeros(2,dtype=int)
        fin[0]=len(np.where(freqs<FB[0])[0])
        fin[1]=len(np.where(freqs<=FB[1])[0])
        obj.freqs=freqs[fin[0]:fin[1]]

        # taper
        tapers, weight = sinusoidal_tapers(len_of_win,NT)

        # Split waveform
        wave_seg=np.zeros((len_of_win,num_of_win))
        for i in range(num_of_win):
            wave_seg[:,i]=wave[i*len_of_win:(i+1)*len_of_win]

        np.seterr(invalid='ignore')

        # Start to calculate the three conditions A, B and C of the wave-field.
        for i in range(NT):
            temp = np.fft.rfft(wave_seg*np.expand_dims(tapers[:,i], axis=1), axis=0) / len_of_win
            temp[1:-1,:] = 2*temp[1:-1,:]
            wave_fft = temp[fin[0]:fin[1],:]

            E_power = np.mean(np.abs(wave_fft)**2, axis=1)
            EEpower = E_power * E_power.reshape(-1, 1)

            A = np.abs(np.mean(wave_fft, axis=1)) **2 / E_power
            tempB = 0; tempC = 0
            for j in range(num_of_win):
                tempB = tempB + wave_fft[:,j] * wave_fft[:,j].reshape(-1, 1)
                tempC = tempC + wave_fft[:,j] * np.conj(wave_fft[:,j].reshape(-1, 1))
            B = np.abs(tempB/num_of_win) **2 / EEpower
            C = np.abs(tempC/num_of_win) **2 / EEpower

            obj.A = obj.A + weight[i]*A
            obj.B = obj.B + weight[i]*B
            obj.C = obj.C + weight[i]*C

        # Calculate the diffuseness proxies of the three conditions
        obj = obj.sRMS(SF)

        np.seterr(invalid='warn')

        obj.stats = {'len_of_win'  : len_of_win, \
                     'num_of_win'  : num_of_win, \
                     'num_of_taper': NT,         \
                     'SF_of_sRMS'  : SF}

        return obj

    @jit
    def sRMS(self,SF):
        """ 
        Calculate the sRMS function values of conditions A, B and C at SF.
        """ 
        A=self.A
        B=self.B
        C=self.C-np.eye(len(A))

        N = math.floor(len(A))
        N2= math.ceil(N/2)
        s = math.ceil(N*SF)

        meanA=np.mean(A)
        meanB=np.mean(B)
        meanC=np.mean(C)

        wA=np.zeros(N)
        wB=np.zeros((N,N))
        wC=np.zeros((N,N))
        proxy=np.zeros(3)

        if s==N:
            wA=1
            wB=1
            wC=1
        else:
            for i in range(N):
                indexi=np.arange(i-s,i+s+1)
                fin=[0,0]
                fin[0]=len(np.where(indexi<0)[0])
                fin[1]=len(np.where(indexi<N)[0])
                indexi=indexi[fin[0]:fin[1]]

                wA[i]=np.mean(A[indexi])

                for j in range(i,N):
                    indexj=np.arange(j-s,j+s+1)
                    fin=[0,0]
                    fin[0]=len(np.where(indexj<0)[0])
                    fin[1]=len(np.where(indexj<N)[0])
                    indexj=indexj[fin[0]:fin[1]]

                    wB[i,j]=np.mean(B[indexi[0]:indexi[-1]+1,indexj[0]:indexj[-1]+1])
                    wC[i,j]=np.mean(C[indexi[0]:indexi[-1]+1,indexj[0]:indexj[-1]+1])
                    wB[j,i]=wB[i,j]
                    wC[j,i]=wC[i,j]

            wA=wA/meanA
            wB=wB/meanB
            wC=wC/meanC

        proxy[0]=np.sqrt(np.mean((wA * A) **2))
        proxy[1]=np.sqrt(np.mean((wB * B) **2))
        proxy[2]=np.sqrt(np.mean((wC * C) **2))
        
        self.proxy = proxy
        return self

    def result(self):
        """ 
        Plot the calculated diffuseness evaluation result
        """
        x  , y    = np.meshgrid(self.freqs, self.freqs)
        fig, axes = plt.subplots(1, 3, figsize=(11, 3))
        
        h0=axes[0].plot(self.freqs, self.A, linestyle='-', color='k')
        axes[0].set_title('Cond. A   $P_A$ = ' + str(np.around(self.proxy[0], 4)) )
        axes[0].set_xlabel('Freq. (Hz)')
        axes[0].set_xlim((np.min(self.freqs),np.max(self.freqs))) 
        axes[0].set_ylim((0,1)) 
        axes[0].set_aspect(1./axes[0].get_data_ratio(), adjustable='box')

        h1=axes[1].pcolor(x, y, self.B, vmin=0, vmax=1)
        axes[1].invert_yaxis()
        axes[1].set_title('Cond. B   $P_B$ = ' + str(np.around(self.proxy[1], 4)))
        axes[1].set_xlabel('Freq. (Hz)')
        axes[1].set_ylabel('Freq. (Hz)')
        axes[1].set_aspect('equal', 'box')

        h2=axes[2].pcolor(x, y, self.C, vmin=0, vmax=1)
        axes[2].invert_yaxis()
        axes[2].set_title('Cond. C   $P_C$ = ' + str(np.around(self.proxy[2], 4)))
        axes[2].set_xlabel('Freq. (Hz)')
        axes[2].set_ylabel('Freq. (Hz)')
        axes[2].set_aspect('equal', 'box')
        
        cax = add_right_cax(axes[2], pad=0.02, width=0.015)
        fig.colorbar(h2, cax=cax)

        plt.show()



def sinusoidal_tapers(len_of_win,num_of_taper):
    """ 
    Sinusoidal tapers
    """

    points=np.arange(1,len_of_win+1)
    tapers=np.zeros((len_of_win,num_of_taper))
    weight=np.zeros(num_of_taper)
    for i in range(num_of_taper):
        tapers[:,i]=np.sqrt(2/(len_of_win+1))*np.sin((np.pi*(i+1)*points)/(len_of_win+1))
        weight[i]=1/num_of_taper

    return tapers, weight

def add_right_cax(ax, pad, width):
    """ 
    Add a 'cax' of the same height to the right of an 'ax'.
    'pad' is the distance between 'cax' and 'ax', and 'width' is the width of 'cax'.
    """ 
    axpos = ax.get_position()
    caxpos = mtransforms.Bbox.from_extents(
        axpos.x1 + pad,
        axpos.y0,
        axpos.x1 + pad + width,
        axpos.y1
    )
    cax = ax.figure.add_axes(caxpos)

    return cax
