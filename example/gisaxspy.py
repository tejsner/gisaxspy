# new gisaxspy header
import numpy as np
import matplotlib.colors as mpcl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# conversions (base is nanometer)
micrometer = 1e3
nanometer = 1
milimeter = 1e6
angstrom = 0.1
degree = 2*np.pi/360

class gisaxsExperiment:
    wavelength = 1.155*angstrom
    SD = 1845.0*milimeter
    alpha_i = 0.14*degree
    db_x = 534
    db_y = 53
    pixelsize = 46.9*micrometer
    clim = (10.0,3e3)
    qylim = [-0.6, 0.6]
    qzlim = [0, 0.7]
    dpi = 80

    time_format = 'minutes'  # options are: 'minutes', 'seconds'

    phi_plot = True

    def __init__(self, wli_file, gisaxs_file, tifpath, offset=0):
        self.wliData = np.genfromtxt(wli_file, delimiter=',')
        self.gisaxsData = np.loadtxt(gisaxs_file, delimiter=',')
        self.tifpath = tifpath
        self.Ddry = min(self.wliData[:,1])
        self.wli_offset = offset

        # set time to start at t=0 and add the WLI offset
        self.gisaxsData[:,1] = self.gisaxsData[:,1] - self.gisaxsData[0,1]
        self.wliData[:,0] = self.wliData[:,0] - self.wliData[0,0] + self.wli_offset

        # create the lists containing the finalized data
        self.gisaxs_num = []
        self.time = []
        self.thickness = []
        self.samx = []

        # sync data
        self.syncData()

    def getLeadingZeros(self, imnum):
        num_len = self.tifpath.count('#')
        num_str = str(imnum)
        return (num_len - len(num_str))

    def showWLI(self):
        x = np.array(self.time)
        y = np.array(self.thickness)
        plt.xlabel('ERROR: time_format set incorrectly. Chose minutes or seconds.')

        if self.time_format == 'seconds':
            plt.xlabel('time [s]')
        elif self.time_format == 'minutes':
            plt.xlabel('time [min]')
            x = x/60.0

        if self.phi_plot:
            y =  self.Ddry / y
            plt.ylabel('$\phi$')
        else:
            plt.ylabel('thickness [nm]')

        plt.plot(x,y,'.')
        plt.show()

    def getPath(self, imnum):
        leading_zeros = self.getLeadingZeros(imnum)
        full_path = self.tifpath.replace(len(str(imnum))*'#', leading_zeros*'0' + str(imnum))
        return full_path

    def plotImg(self, imnum):
        qy1 = self.qylim[0]
        qy2 = self.qylim[1]
        qz1 = self.qzlim[0]
        qz2 = self.qzlim[1]

        py1, pz1 = self.getPixels(qy1, qz1)
        py2, pz2 = self.getPixels(qy2, qz2)

        a = self.gisaxs_num.index(imnum)
        thickness = self.thickness[a]
        samx = self.samx[a]
        time_s = self.time[a]
        time_m = time_s/60
        time_string = 'ERROR: time_format set incorrectly. Chose minutes or seconds.'

        if self.time_format == 'minutes':
            time_string = ' | t=' + "%0.1f" % time_m + ' min'
        elif self.time_format == 'seconds':
            time_string = ' | t=' + "%0.1f" % time_s + ' s'

        im = mpimg.imread(self.getPath(imnum))
        # we need to flip the image to get the right pixels
        # and flip again to display correctly
        im = np.flipud(np.flipud(im)[pz1:pz2, py1:py2])

        plt.imshow(im, norm=mpcl.LogNorm(), clim=self.clim, extent=[qy1, qy2, qz1, qz2])
        plt.title('d=' + str(thickness) + 'nm | samx=' + str(samx) + time_string)
        plt.xlabel('q_y [nm^(-1)]')
        plt.ylabel('q_z [nm^(-1]]')

    def showImg(self, imnum):
        self.plotImg(imnum)
        plt.show()

    # given a coordinate qy, qz, find and return the corresponding pixels
    def getPixels(self, qy, qz):
        alpha_f = np.arcsin(qz*self.wavelength/2/np.pi - np.sin(self.alpha_i))
        d_z = self.SD*np.tan(alpha_f + self.alpha_i)
        py = d_z/self.pixelsize + self.db_y

        phi = np.arcsin(qy*self.wavelength/2/np.pi/np.cos(alpha_f))
        d_y = self.SD*np.tan(phi)
        px = d_y/self.pixelsize + self.db_x

        return int(px), int(py)

    # syncs up the data and puts the result in 4 seperate lists of equal length
    def syncData(self):
        for i, t in enumerate(self.gisaxsData[:,1]):
            a = np.argmin(np.abs(t - self.wliData[:,0]))
            self.time.append(t)
            self.thickness.append(self.wliData[a,1])
            self.gisaxs_num.append(int(self.gisaxsData[i,0]))
            self.samx.append(self.gisaxsData[i,2])

    def saveImages(self, range=False):
        if not range:
            range = self.gisaxs_num

        for imnum in range:
            self.plotImg(imnum)
            leading_zeros = self.getLeadingZeros(imnum)
            plt.savefig(leading_zeros*'0' + str(imnum) + '.png', bbox_inches='tight', dpi=self.dpi)
            plt.clf()

    def DpdakPlot(self, dpdak_file, row=0, yrange=False, yname=False, show_wli=True):
        peakfit_data = np.genfromtxt(dpdak_file, skip_header=4)

        fig, ax1 = plt.subplots()
        y = peakfit_data[row,:]

        x = np.array(self.time)
        ax1.set_xlabel('ERROR: time_format set incorrectly. Chose minutes or seconds.')
        if self.time_format == 'seconds':
            ax1.set_xlabel('time [s]')
        elif self.time_format == 'minutes':
            ax1.set_xlabel('time [min]')
            x = x / 60.0

        ax1.plot(x, y, '.')

        if yrange:
            ax1.set_ylim(yrange)

        if yname:
            ax1.set_ylabel(yname, color='b')

        for tl in ax1.get_yticklabels():
            tl.set_color('b')


        if show_wli:
            ax2 = ax1.twinx()
            y2 = self.thickness

            if self.phi_plot:
                y2 = self.Ddry / y2
                ax2.set_ylabel('$\phi$', color='g')
            else:
                ax2.set_ylabel('thickness [nm]', color='g')

            for tl in ax2.get_yticklabels():
                tl.set_color('g')

            ax2.plot(x, y2, '-g')

        plt.show()
