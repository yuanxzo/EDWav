"""
Name  : The Method of Evaluation of the Diffuseness of the Wavefield (EDWav Method)
Author: Bo YANG (seism.yang@foxmail.com)
Reference: Yang, B., Meng, H., Gu, N., Liu, X., Chen, X., & Ben‐Zion, Y. (2024). 
           A frequency domain methodology for quantitative evaluation of diffuse wavefield with applications to seismic imaging. 
           Journal of Geophysical Research: Solid Earth, 129, e2024JB028895. https://doi.org/10.1029/2024JB028895
Latest: 2024-05-14   

Copyright (C) 2024 Bo YANG
""" 

__Version__ = '2.0.0'
__Testing_Environment__ = 'Python==3.11.7  numba==0.58.1'    # If the run fails, consider using these specified versions

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
        obj.result()              # show result

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
        Function introduction: evaluate the diffuseness of input seismic waveform segment.
        
        Input parameters introduction:
        - wave :    single channel seismic waveform data, for example: st=obspy.read('test.sac'); wave=st[0].data;
        - FS   :    waveform data sampling rate;
        - FB   :    frequency range to be evaluated, for example: [1, 100], unit: Hz;
        - FR   :    the frequency resolution of the result;
        - NT   :    number of multitapers;
        - SF   :    a parameter of sRMS, 0<=SF<=1. If SF==-1, then the sRMS is not calculated

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
            raise Exception("The waveform length is too short, end the evaluation!")
        
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
        nfin=len(obj.freqs)
               
        # taper
        tapers, weight = sinusoidal_tapers(len_of_win,NT)

        # Split waveform
        wave_seg=np.zeros((len_of_win,num_of_win))
        for i in range(num_of_win):
            wave_seg[:,i]=wave[i*len_of_win:(i+1)*len_of_win]

        np.seterr(invalid='ignore')

        # Start to calculate the three conditions A, B and C of the waveform.
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

        # Suppressing side lobe effect of sinusoidal tapers
        # (Note that this is only an empirical operation. You can choose to comment on the following lines of code to reject this operation.)
        Cw=np.ones((nfin*nfin,1))
        CC=np.logspace(1,0,NT+1)
        for k in range(1,NT+1):
            Cw[k:-1:nfin+1]=CC[k-1]
        Cw=Cw.reshape(nfin,nfin)*Cw.reshape(nfin,nfin).T
        obj.C=obj.C/Cw
               
        # Calculate the diffuseness proxies of the three conditions
        if SF>=0 and SF<=1:
            obj = obj.sRMS(SF)
        else:
            obj.proxy=[np.nan, np.nan, np.nan]
        
        np.seterr(invalid='warn')

        obj.stats = {'len_of_win'  : len_of_win, \
                     'num_of_win'  : num_of_win, \
                     'num_of_taper': NT,         \
                     'SF_of_sRMS'  : SF}

        return obj

    def evaluate_sliding(wave:np.ndarray=[], FS=1, FB=[], FR=100, NT=1, SF=0.05, NW=30, DW=None, NM=1):
        """
        Function introduction:
        - Evaluate sliding window analysis on an input seismic waveform.

        Parameters:
        - wave (np.ndarray): Single channel seismic waveform data. Example: st = obspy.read('test.sac'); wave = st[0].data.
        - FS (int):          Sampling rate of the waveform data.
        - FB (list):         Frequency range to be evaluated, e.g., [1, 100] in Hz.
        - FR (int):          Frequency resolution of the result.
        - NT (int):          Number of multitapers.
        - SF (float):        A parameter of sRMS, 0 <= SF <= 1. If SF == -1, the sRMS is not calculated.
        - NW (int):          Number of windows. Default is 30.
        - DW (int):          Window step size. Default is equal to NW.
        - NM (int):          Whether to eliminate the result dimension? Default is yes (1).

        Returns:
        - proxy (np.ndarray): Proxy array.
        - stats (dict): Dictionary containing various statistics:
            - 'freqs': Frequencies analyzed.
            - 'len_of_win': Length of each window.
            - 'num_of_win': Number of windows.
            - 'num_of_sld': Number of sliding windows.
            - 'num_of_taper': Number of tapers.
            - 'SF_of_sRMS': Scaling factor for sRMS.
            - 'NW': Number of windows.
            - 'DW': Window step size.
            - 'FB': Frequency band analyzed.

        Raises:
        - Exception: If the length of the sliding window is greater than the length of the waveform.
        """
        obj=EDWav()
        NT = int(NT)
        FR = int(FR)
        NW = int(NW)
        if DW is None:
            DW = NW
        DW = int(DW)
        NM = int(NM)

        # Number of waveform sampling points
        npts=len(wave)
        wave = np.array(wave).reshape(-1, 1)

        # Frequency range to be analyzed
        if len(FB) != 2:
            FB=[0, FS/2]

        # Calculate the minimum length of the window according to the 'FR' and 'FB' entered
        len_of_win=math.ceil(FR*FS/(FB[1]-FB[0])) 
        if npts < len_of_win * NW:
            raise ValueError("The length of the sliding window is greater than the length of the waveform, end the evaluation!")

       # Frequencies that this window can analyze
        freqs=np.arange(0, FS/2+FS/len_of_win, FS/len_of_win)
        if freqs[-1]>FS/2:
            freqs = freqs[:-1]

        # Frequencies to be analyzed
        fin = (freqs >= FB[0]) & (freqs <= FB[1])
        obj.freqs=freqs[fin]
        nfin = len(obj.freqs)

        # Number of windows
        num_of_win = npts // len_of_win

        # Number of sliding windows
        sindex = np.arange(0, num_of_win - NW + 1, DW)
        num_of_sld = len(sindex)

        # Length of window in final
        len_of_win = int(np.floor(len_of_win + (npts - num_of_win * len_of_win) / num_of_win))
        
        # taper
        tapers, weight = sinusoidal_tapers(len_of_win, NT)

        # split waveform
        wave_seg = np.zeros((len_of_win, num_of_win))
        for i in range(num_of_win):
            wave_seg[:, i] = wave[i * len_of_win:(i + 1) * len_of_win,0]

        # obtain the spectrum of all windows under the action of different tapers
        wave_fft = np.zeros((nfin, num_of_win, NT), dtype=complex)
        for i in range(NT):
            wave_tmp = np.fft.rfft(wave_seg*np.expand_dims(tapers[:,i], axis=1), axis=0) / len_of_win
            wave_tmp[1:-1,:] =  wave_tmp[1:-1,:]*2
            wave_fft[:, :, i] = wave_tmp[fin,:]

        np.seterr(invalid='ignore')
        
        # the denominator of the three conditions
        E_power = np.ones((nfin, NT))
        EEpower = np.ones((nfin, nfin, NT))
        if NM == 1:
            for i in range(NT):
                E_power[:,i] = np.median(np.abs(wave_fft[:, :, i]) ** 2, axis=1)
                EEpower[:,:,i] = E_power[:,i] * E_power[:,i].reshape(-1, 1)

        # begin sliding evaluation
        proxy = np.zeros(npts)
        stats = {'stack': np.zeros(npts)}
        for k in range(num_of_sld):

            index = [sindex[k], sindex[k] + NW - 1]
            indexs = np.arange(index[0] * len_of_win, (index[1] + 1)* len_of_win)
            
            obj.A = 0
            obj.B = 0
            obj.C = 0

            for i in range(NT):
                wave_use = wave_fft[:, index[0]:index[1] + 1, i]

                A = np.abs(np.mean(wave_use, axis=1)) ** 2 / E_power[:,i]
                tempB = 0
                tempC = 0
                for j in range(NW):
                    tempB += wave_use[:, j] * wave_use[:, j].reshape(-1, 1)
                    tempC += wave_use[:, j] * np.conj(wave_use[:, j].reshape(-1, 1))
                B = np.abs(tempB / NW) ** 2 / EEpower[:,:,i]
                C = np.abs(tempC / NW) ** 2 / EEpower[:,:,i]

                obj.A += weight[i] * A
                obj.B += weight[i] * B
                obj.C += weight[i] * C

            # adjust the result
            if np.any(np.isnan(obj.A)) or np.any(np.isnan(obj.B)) or np.any(np.isnan(obj.C)):
                proxy[indexs] += 0
                stats['stack'][indexs] += 0
            else:
                # Suppressing side lobe effect of sinusoidal tapers
                # (Note that this is only an empirical operation. You can choose to comment on the following lines of code to reject this operation.)
                Cw = np.ones((nfin, nfin))
                CC = np.logspace(1, 0, NT + 1)
                for kk in range(2, NT + 2):
                    Cw[kk:, kk:] = CC[kk - 1]
                Cw = Cw * Cw.T
                obj.C = obj.C / Cw

                # claculate the diffuseness proxies of the three conditions
                obj = obj.sRMS(SF)

                # save the result
                proxy[indexs] += np.mean(obj.proxy)  # 忽略obj.C
                stats['stack'][indexs] += 1

        # return the result
        proxy = proxy / stats['stack']

        stats['freqs'] = obj.freqs
        stats['len_of_win'] = len_of_win
        stats['num_of_win'] = num_of_win
        stats['num_of_sld'] = num_of_sld
        stats['num_of_taper'] = NT
        stats['SF_of_sRMS'] = SF
        stats['NW'] = NW
        stats['DW'] = DW
        stats['FB'] = FB

        np.seterr(invalid='warn')
        
        return proxy, stats

           
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
        # Considering the similarity between conditions B and C, it is possible to make the proxy of condition C 
        # equal to that of condition B to avoid the influence of diagonal sidelobes on the proxy of condition C
        # (Note that this is only an empirical operation. You can choose to comment on the following line of code to reject this operation.)
        proxy[2]=proxy[1]
        
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
    weight=np.arange(num_of_taper,0,-1)
    weight=weight/np.sum(weight)
    for i in range(num_of_taper):
        tapers[:,i]=np.sqrt(2/(len_of_win+1))*np.sin((np.pi*(i+1)*points)/(len_of_win+1))

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
