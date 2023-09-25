

# plotting elements
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
mpl.rcParams['font.weight'] = 'medium';mpl.rcParams['xtick.labelsize'] = 10;mpl.rcParams['ytick.labelsize'] = 10


class parameter_figure():
    def __init__(self,pltkeys):

        self.fig =  plt.figure(figsize=(12.5,3.5),facecolor='white')

        xmin = 0.04
        ymin = 0.12
        dx = 0.18
        dy = 0.42
        xbuf = 0.06

        # upper row
        self.ax1 = self.fig.add_axes([xmin+0*(dx+xbuf),ymin+1*dy,dx,dy]) # f
        self.ax2 = self.fig.add_axes([xmin+1*(dx+xbuf),ymin+1*dy,dx,dy]) # Lx
        self.ax3 = self.fig.add_axes([xmin+2*(dx+xbuf),ymin+1*dy,dx,dy]) # Ly
        self.ax4 = self.fig.add_axes([xmin+3*(dx+xbuf),ymin+1*dy,dx,dy]) # Lz

        self.ax5 = self.fig.add_axes([xmin+0*(dx+xbuf),ymin+0*dy,dx,dy]) # alpha
        self.ax6 = self.fig.add_axes([xmin+1*(dx+xbuf),ymin+0*dy,dx,dy]) # sigmax
        self.ax7 = self.fig.add_axes([xmin+2*(dx+xbuf),ymin+0*dy,dx,dy]) # sigmay
        self.ax8 = self.fig.add_axes([xmin+3*(dx+xbuf),ymin+0*dy,dx,dy]) # sigmaz

        self.axlist = [self.ax1,self.ax2,self.ax3,self.ax4,self.ax5,self.ax6,self.ax7,self.ax8]

        self.cwheel = ['blue','black','red']

        self.pltkeys = pltkeys

    def _add_data(self,compname,ikey,radval,median,lo,hi):

        if compname=='bar':
            self.axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='black')
            self.axlist[ikey].scatter(radval,median,marker='X',edgecolor='none',facecolor='black')
        elif compname=='disc':
            self.axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='blue')
            self.axlist[ikey].scatter(radval,median,marker='X',edgecolor='none',facecolor='blue')
        elif compname=='knot':
            self.axlist[ikey].plot([radval,radval],[median+lo,median+hi],color='red')
            self.axlist[ikey].scatter(radval,median,marker='X',edgecolor='none',facecolor='red')

    def _frame_figure(self):

        for ax in self.axlist:
            ax.tick_params(axis="both",direction="in",which="both")

        self.ax1.axis([0.,5.,0.,1.1])
        self.ax2.axis([0.,5.,-50,50.])
        self.ax3.axis([0.,5.,-50,50.])
        self.ax4.axis([0.,5.,-800,100.])
        self.ax5.axis([0.,5.,0.,90.])
        self.ax6.axis([0.,5.,0.,220.])
        self.ax7.axis([0.,5.,0.,220.])
        self.ax8.axis([0.,5.,0.,220.])

        for ax in [self.ax1,self.ax2,self.ax3,self.ax4]:
            ax.set_xticklabels(())

        self.ax6.set_xlabel('radius (kpc)',x=1.3)

        for ikey,key in enumerate(self.pltkeys):
            self.axlist[ikey].set_ylabel(key)

    def _print_figure(self,modeltag,appendix):

        plt.savefig('figures/fitvalues_{0}{1}.png'.format(modeltag,appendix),dpi=300)
