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
    db_x = 534
    db_y = 53
    clim = (10.0,3e3)
    pixel_mask = [100, 100, 900, 600]

    phi_plot = True

    def __init__(self, wli_file, gisaxs_file, tifpath):
        self.wliData = np.genfromtxt(wli_file, delimiter=',')
        self.gisaxsData = np.loadtxt(gisaxs_file, delimiter=',')
        self.tifpath = tifpath
        self.Ddry = min(self.wliData[:,1])

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

    def __getpath(self, imnum):
        num_len = self.tifpath.count('#')
        num_str = str(imnum)
        leading_zeros = num_len - len(num_str)
        full_path = self.tifpath.replace(num_len*'#', leading_zeros*'0' + num_str)
        return full_path

    def showImg(self, imnum):
        im = mpimg.imread(self.__getpath(imnum))
        plt.imshow(im, norm=mpcl.LogNorm(), clim=self.clim)
        plt.show()

test = gisaxsExperiment('test_wli.csv', 'test_gisaxs.csv', 'C:\Scattering\CHESS_2015_DATA\corr\pspb05_###.tif')
test.showImg(18)
test.showWLI()
