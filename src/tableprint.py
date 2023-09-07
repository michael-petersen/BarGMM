

import numpy as np


class table_print():
    def __init__(self,format='terminal',outputfile=None):
        """
        format:
            'csv'
            'terminal'
            'markdown'
            'tex'

        #

        """
        self.format = format
        self.outputfile = outputfile

        if ((format == 'csv') | (format == 'markdown') | (format == 'tex')) & (outputfile == None):
            # throw error
            print('we need a file output!')

        if ((format == 'csv') | (format == 'markdown') | (format == 'tex')):
            self.f = outputfile
        else:
            self.f = None


    def print_header(self):

        if self.format=='terminal':
            #self.f = None
            print('{0:4s}{1:3s}{2:9s}{3:9s} {4:9s} {5:9s} {6:9s} {7:9s} {8:9s} {9:9s}'.format(' ',' ','f','Lx','Ly','Lz','alpha','sx','sy','sz'))
        if self.format=='csv':
            print('comp,binmin,binmax,starmed,starmed-,starmed+,f,f-,f+,Lx,Lx-,Lx+,Ly,Ly-,Ly+,Lz,Lz-,Lz+,alpha,alpha-,alpha+,sigmax,sigmax-,sigmax+,sigmay,sigmay-,sigmay+,sigmaz,sigmaz-,sigmaz+,',file=self.f)
        if self.format=='markdown':
            #self.f = open(self.outputfile,'w')
            print('|comp|radii| f | L<sub>x</sub> | L<sub>y</sub> | L<sub>z</sub> | angle | w<sub>x</sub> | w<sub>y</sub> | w<sub>z</sub> |',file=self.f)
            print('|---|---|---| ---| --- | ---| --- | --- | --- | --- |',file=self.f)
        if self.format=='tex':
            #self.f = open(self.outputfile,'w')
            print('comp &radii & f & $L_x$& $L_y$ & $L_z$ & $\alpha$ & $\sigma_x$ & $\sigma_y$ & $\sigma_z$ \\\\',file=self.f)

    def print_column1(self,comptag,binprefacs,minrad,maxrad,irad,cnum,Stars,criteria):

        if self.format=='terminal':
            print('{0:2.1f}-{1:2.1f}'.format(np.round(binprefacs[irad]*minrad,1),np.round(binprefacs[irad]*maxrad,1)))

        if self.format=='markdown':
            # print to markdown
            print('|'+comptag[minrad][cnum],end='|',file=self.f)
            print('{0:2.1f}-{1:2.1f}'.format(np.round(binprefacs[irad]*minrad,1),np.round(binprefacs[irad]*maxrad,1)),end='|',file=self.f)

        if self.format=='csv':
            # print to csv
            print(comptag[minrad][cnum],end=',',file=self.f)
            print('{0:2.1f},{1:2.1f}'.format(np.round(binprefacs[irad]*minrad,1),np.round(binprefacs[irad]*maxrad,1)),end=',',file=self.f)
            print('{0:3.2f},{1:3.2f},{2:3.2f}'.format(np.nanmedian(Stars['R'][criteria]),np.nanpercentile(Stars['R'][criteria],14.),np.nanpercentile(Stars['R'][criteria],86.)),end=',',file=self.f)

        if self.format=='tex':
            # print to latex
            print(comptag[minrad][cnum],end='&',file=self.f)
            print('{0:2.1f}-{1:2.1f}'.format(np.round(binprefacs[irad]*minrad,1),np.round(binprefacs[irad]*maxrad,1)),end='&',file=self.f)

    def print_key_column(self,median,lo,hi,rounding=1):

        if self.format=='terminal':
            print(' {0:9.4f}'.format(median),end='')

        if self.format=='markdown':
            print('{0}<sup>+{1}</sup><sub>-{2}</sub>'.format(np.round(median,rounding),np.round(np.abs(lo),rounding),np.round(np.abs(hi),rounding)),end='|',file=self.f)

        if self.format=='tex':
            print('{0}^{{+{1}}}_{{-{2}}}'.format(np.round(median,rounding),np.round(np.abs(lo),rounding),np.round(np.abs(hi),rounding)),end='&',file=self.f)

        if self.format=='csv':
            print('{0},-{2},+{1}'.format(np.round(median,rounding),np.round(np.abs(lo),rounding),np.round(np.abs(hi),rounding)),end=',',file=self.f)

    def print_columnX(self,comptag,minrad,cnum):

        if self.format=='terminal':
            print(' '+comptag[minrad][cnum])

        if self.format=='markdown':
            print('',file=self.f)

        if self.format=='tex':
            print('\\\\',file=self.f)

        if self.format=='csv':
            print('',file=self.f)
