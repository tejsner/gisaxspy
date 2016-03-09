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
    wli_offset = 0

    phi_plot = True

    def __init__(self, wli_file, gisaxs_file, tifpath, offset):
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

    def showWLI(self):
        x = self.wliData[:,0]
        y = self.wliData[:,1]

        if self.phi_plot:
            plt.plot(x, self.Ddry/y)
            plt.ylabel('$\phi$')
        else:
            plt.plot(x/60, y)
            plt.ylabel('thickness [nm]')

        plt.xlabel('time [s]')
        plt.show()

    def getPath(self, imnum):
        num_len = self.tifpath.count('#')
        num_str = str(imnum)
        leading_zeros = num_len - len(num_str)
        full_path = self.tifpath.replace(num_len*'#', leading_zeros*'0' + num_str)
        return full_path

    def showImg(self, imnum):
        qy1 = self.qylim[0]
        qy2 = self.qylim[1]
        qz1 = self.qzlim[0]
        qz2 = self.qzlim[1]

        py1, pz1 = self.getPixels(qy1, qz1)
        py2, pz2 = self.getPixels(qy2, qz2)

        a = self.gisaxs_num.index(imnum)
        thickness = self.thickness[a]
        samx = self.samx[a]

        im = mpimg.imread(self.getPath(imnum))
        # we need to flip the image to get the right pixels
        # and flip again to display correctly
        im = np.flipud(np.flipud(im)[pz1:pz2, py1:py2])

        plt.imshow(im, norm=mpcl.LogNorm(), clim=self.clim, extent=[qy1, qy2, qz1, qz2])
        plt.title('d=' + str(thickness) + 'nm | samx=' + str(samx))
        plt.xlabel('q_y [nm^(-1)]')
        plt.ylabel('q_z [nm^(-1]]')
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

    def syncData(self):
        for i, t in enumerate(self.gisaxsData[:,1]):
            a = np.argmin(np.abs(t - self.wliData[:,0]))
            self.time.append(t)
            self.thickness.append(self.wliData[a,1])
            self.gisaxs_num.append(int(self.gisaxsData[i,0]))
            self.samx.append(self.gisaxsData[i,2])


test = gisaxsExperiment('test_wli.csv', 'test_gisaxs.csv', 'C:\Scattering\CHESS_2015_DATA\corr\pspb05_###.tif', 0)
test.showImg(180)
