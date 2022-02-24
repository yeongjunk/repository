#!/usr/bin/python
import glob
import math
import os
import numpy as np
import time
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.mlab import griddata
start_time = time.time()

corr = 1 # 1: uncorrelated, random, 2: anti-correlated, random, 3 correlated: AAH 

path = "./"

x_text = 0.9
y_text = 0.07
color_abc = 'white'

### hamiltonian parameters #############################

lwdth = 0.2 # line width in plots

f_size_1 = 28 #fontsize
lab_size_1 = 22 #labelsize
legendsize = 20

#color_list=['deepskyblue', 'red', 'black', 'indigo', 'lightgreen', 'forestgreen','saddlebrown', 'darkgreen', 'green', 'blue', 'red']
color_list=['black', 'darkviolet', 'navy', 'blue', 'dodgerblue', 'cyan', 'lightgreen', 'forestgreen','saddlebrown', 'darkgreen', 'green']

cmap = mpl.cm.get_cmap('gist_heat')

# https://matplotlib.org/api/markers_api.html
marker_list = ['s', 'o', 'X', '*', 'v', '.', 'D', '1', '2', '3', '4', 'x', '8', 'P', 'H']

for zzz in range(15):
    
    ######################## Plotting  parameters ####################################
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True

    #plt.rcParams.update({'lines.markeredgewidth': 2})

    #plt.rcParams['lines.markeredgewidth'] = 2

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    nrows = 1
    ncols = 1
    fig, axs = plt.subplots(nrows,ncols,figsize=(13.0, 13.0*np.sqrt(2.0)/2.0))#was before np.sqrt(2.0)/2
    
    jMAX = 2

    for j in range(jMAX):
        
        t0 = 1.0

        i_max = 100000

        start = 0

        list_of_E = [np.absolute(t0)] #[-2.0, 2.0, 2.0001, 2.01]

        E_FB = np.absolute(t0)
        
        ThetaI = 1.0*np.pi*(random.random()-0.5)
        phi_plus = 2.0*np.pi*(random.random()-0.5)
        
        Theta2 = np.pi*(random.random()-0.5) #0.78615 #2*np.pi*(random.random()-0.5)
        phi1 = 2.0*np.pi*(random.random()-0.5)
        phi2 = 2.0*np.pi*(random.random()-0.5)

        z = np.cos(Theta2)*(np.cos(phi1)+np.sin(phi1)*1.0j)
        w = np.sin(Theta2)*(np.cos(phi2)+np.sin(phi2)*1.0j)

        tI = - np.sin(2.0*ThetaI)*(np.cos(phi_plus)+np.sin(phi_plus)*1.0j)

        print ThetaI, Theta2, phi1,phi2
        print z, w, t0, np.conj(t0)

        for i_E, E in enumerate(list_of_E):
            xi_W = []
            xi_W2 = []
            xi_std = []
            
            color = cmap((j+1.0)/(jMAX+1.0))
            
            #list_of_Ws = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.14, 0.2, 0.28, 0.4, 0.56, 0.8, 1.0, 1.2, 1.4, 1.7, 2.0, 2.4, 2.8, 3.4, 4.0, 4.8, 5.6, 6.8, 8.0, 9.0, 10.0, 14.0, 20.0, 28.0, 40.0, 56.0, 100.0]
            
            list_of_Ws = list(np.logspace(-5, -1, num = 15, endpoint=False)) + list(np.logspace(-1, 1, num = 20, endpoint=False)) + list(np.logspace(1, 3, num = 5)) + list(np.linspace(0.5*E_FB*5., 1.5*E_FB*5., num = 31,endpoint=True))
            #list_of_Ws = list(np.logspace(-5, -1, num = 5, endpoint=False)) + list(np.logspace(-1, 1, num = 5, endpoint=False)) + list(np.logspace(1, 3, num = 5)) + list(np.linspace(0.5*E_FB*5., 1.5*E_FB*5., num = 6,endpoint=True))
            
            list_of_Ws.sort()
            
            #print list_of_Ws
            
            for j_W, W in enumerate(list_of_Ws):
                print "W is ", W
                
                R_0 = 1.0
                
                xi = []
                
                xi_p = 0.0
                
                for i in range(1, i_max+1):
                    mu_a = W*(random.random()-0.5)
                    if corr == 1:
                        mu_b = W*(random.random()-0.5) # uncorrelated
                    elif corr == 2:
                        mu_b = -mu_a
                    else:
                        mu_b = 0.0
                    
                    e_p =  t0*np.cos(2.0*ThetaI) + (mu_a*z*np.conj(z) + mu_b*w*np.conj(w))
                    e_f = -t0*np.cos(2.0*ThetaI) + (mu_b*z*np.conj(z) + mu_a*w*np.conj(w))
                    
                    t_j = (mu_a-mu_b)*w*np.conj(z)
                    
                    R_1 = (E - e_p - t_j/R_0)/tI
                    R_0 = (E - e_f - np.conj(tI)/R_1)/np.conj(t_j)
                    
                    xi_p += np.log(np.absolute(R_1)) + np.log(np.absolute(R_0))
                                
                    xi.append(xi_p/(2.0*i+0.0))
                    
                    if i%100000 < 1:
                        print i, xi_p/(2*i+0.0)
                        print R_1, R_0

                print "xi to -1 and xi are ", np.mean(xi), 1.0/np.mean(xi)
                
                        #print xi
                print np.mean(xi), np.mean(xi[start:]), np.std(xi), np.std(xi[start:])
                
                xi_W.append(1.0/np.mean(xi[start:]))
                xi_W2.append(1.0/np.absolute(np.mean(xi[start:])))
                xi_std.append(np.std(xi[start:])/np.absolute(np.mean(xi[start:]))*np.absolute(np.mean(xi[start:])))
            
            path2 = "./data/"
            if corr == 1:
                filename = "xi_Carlo-Creutz_exp_E_%.6f_tI_%.5f_Theta2_%.5f_phi1_%.5f_phi2_%.5f_W_%.5f.dat" %(E, np.absolute(tI), Theta2, phi1, phi2, W) 
            elif corr == 2:
                filename = "xi_Carlo-Creutz_anti_exp_E_%.6f_tI_%.5f_Theta2_%.5f_phi1_%.5f_phi2_%.5f_W_%.5f.dat" %(E, np.absolute(tI), Theta2, phi1, phi2, W) 
            else:
                filename = "xi_Carlo-Creutz_other_exp_E_%.6f_tI_%.5f_Theta2_%.5f_phi1_%.5f_phi2_%.5f_W_%.5f.dat" %(E, np.absolute(tI), Theta2, phi1, phi2, W) 
            f_out = open(path2 + filename, 'w')

            # https://stackoverflow.com/questions/37377366/writing-lists-of-data-to-a-text-file-column-by-column
            for val in zip(list_of_Ws, xi_W, xi_W2, xi_std):
                f_out.write('{}\t{}\t{}\n'.format(val[0], val[1], val[2], val[3]))                
            f_out.close()

            # https://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html
            # https://stackoverflow.com/questions/7601334/how-to-set-the-line-width-of-error-bar-caps-in-matplotlib
            #y_err_low  = np.array(L_M) - np.array(std_L_M)
            #y_err_high = np.array(L_M) + np.array(std_L_M)
            
            #(_, caps, _) = axs.errorbar(list_of_Ws, L_M, yerr = [y_err_low, y_err_high], color = color_list[i_M], capsize=4, elinewidth=2, fmt='o', marker = marker_list[i_M], markersize = marker_size, linestyle = '-', label = r'$M = %lu$' %(M), zorder = 1)
            
            y_err_high = np.array(xi_std)
            
            (_, caps, _) = axs.errorbar(list_of_Ws, xi_W2, yerr = y_err_high, linestyle = 'none', marker = marker_list[j], markersize = 10, fillstyle = 'none', color = color, label = r'$E = %.1f, \Theta_I = %.3f, \Theta_2 = %.3f, \varphi_1 = %.3f, \varphi_2 = %.3f$'%(E, ThetaI, Theta2, phi1, phi2))
            
            for cap in caps:
                #cap.set_color('red')
                cap.set_markeredgewidth(2)
                
            #axs.plot(list_of_Ws, xi_W, linestyle = 'none', marker = marker_list[2*i_E], markersize = 7, fillstyle = 'none', color = color)
            #axs.plot(list_of_Ws, xi_W2, linestyle = 'none', marker = marker_list[i_E%(len(list_of_E))], markersize = 10, fillstyle = 'none', color = color, label = r'$E = %.6f$'%(E))
            
            lam_mean = np.mean(xi_W2[:14])
            
            i_max = np.argmax(xi_W2)
            
            axs.text(0.3, 0.45-(j)*0.2, r'$\xi_{0} = %.4f \pm %.4f, \, \xi_{max} = %.4f$' %(lam_mean, np.std(xi_W2[:14]), max(xi_W2)), transform=axs.transAxes, fontsize=f_size_1-6, color = color, horizontalalignment='center')
            axs.text(0.85, 0.45-(j)*0.2, r'$W_{max} = %.4f$' %(list_of_Ws[i_max]), transform=axs.transAxes, fontsize=f_size_1-6, color = color, horizontalalignment='center')
            axs.text(0.2, 0.35-(j)*0.2, r'$\xi_{max} = %.4f \xi_{0}$' %(max(xi_W2)/lam_mean), transform=axs.transAxes, fontsize=f_size_1-6, color = color, horizontalalignment='center')
            axs.text(0.65, 0.35-(j)*0.2, r'$W_{max} = %.4f E_{FB}$' %(list_of_Ws[i_max]/E_FB), transform=axs.transAxes, fontsize=f_size_1-6, color = color, horizontalalignment='center')

    xlim_min = 8*np.power(10.0,-6)
    xlim_max = 1.125*np.power(10.0,3)

    x_ticks = list(np.logspace(-5, 3, num = 9, endpoint=True))
    x_minor = list(np.logspace(-5, 3, num = 9*10, endpoint=True))

    axs.set_xlim([xlim_min, xlim_max])
    axs.set_xticks(x_ticks)
    axs.set_xticks(x_minor, minor=True)

    ylim_min = np.power(10.0, -1)
    ylim_max = np.power(10.0, 1)

    y_ticks = list(np.logspace(-1, 1, num = 3, endpoint=True))
    y_minor = list(np.logspace(-1, 1, num = 3*10, endpoint=True))

    axs.set_ylim([ylim_min, ylim_max])
    axs.set_yticks(y_ticks)
    axs.set_yticks(y_minor, minor=True)

    axs.legend(fontsize=legendsize, loc='best', numpoints = 1, ncol=1)

    axs.set_xlabel(r'$W$', fontsize = f_size_1)
    axs.set_ylabel(r'$\xi$', rotation = 0, fontsize = f_size_1, labelpad = 15)

    axs.set_xscale('log')
    axs.set_yscale('log')

    axs.plot([4.*E,4.*E], [ylim_min, ylim_max], linestyle = '-.', linewidth=1, color = 'gray')
    axs.plot([5.*E,5.*E], [ylim_min, ylim_max], linestyle = '-', linewidth=1.5, color = 'black')
    axs.plot([6.*E,6.*E], [ylim_min, ylim_max], linestyle = '-.', linewidth=1, color = 'gray')

    print E

    print("--- %s seconds ---" % (time.time() - start_time))  

    ############### Outputting figure ###########################
    output_path='./'
    if corr == 1:
        #pdf_file_name = "Diamond-xi_uncorrelated_Re_t0_%.5f_Theta2_%.5f_phi1_%.5f_phi2_%.5f_W_%.5f.dat" %(np.real(t0), Theta2, phi1, phi2, W) 
        pdf_file_name = "Carlo-Creutz-xi_uncorrelated_r1_%lu_t0_01.dat" %(zzz)
    elif corr == 2:
        pdf_file_name = 'Carlo-Creutz-xi_anti_correlated_t1_%.5f_t2_%.5f_t1t_%.5f_t2t_%.5f_BOTH.pdf' %(t1, t2, t1t, t2t)
    else:
        pdf_file_name = 'Carlo-Creutz-xi_other_t1_%.5f_t2_%.5f_t1t_%.5f_t2t_%.5f_BOTH.pdf' %(t1, t2, t1t, t2t)
    pdf_file_name = os.path.join(output_path, pdf_file_name)
    plt.savefig(pdf_file_name, dpi=600, facecolor='w', edgecolor='w', orientation='portrait', format="pdf", bbox_inches='tight')

    #plt.show()
