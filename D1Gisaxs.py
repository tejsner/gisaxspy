from numpy import *
from pylab import *

# conversions (base is nanometer)
micrometer = 1e3
nanometer = 1
milimeter = 1e6
angstrom = 0.1
degree = 2*pi/360

# sync GISAXS and thickness data
# This function takes 4 variables
# t_data is a string with the relative path to the wli thickness data
# g_data is a string with the relative path to the gisaxs data
# t_index and g_index are the indices that syncronise the times of the two datasets
# It returns a dataset with gisaxs measurement number, thickness, time, samx (coloumn 0,1,2,3 respectively)
def gisaxs_thickness(t_data,g_data,t_index,g_index,d1_2013=False,d1_2015=False):
    gisaxs_data = genfromtxt(g_data)

    # data is formatted differently depending on when the experiment was performed
    if (d1_2013 | d1_2015):
        thickness_data = genfromtxt(t_data, delimiter=',',skip_header=9)
        if d1_2015:
            thickness_data = transpose(array([thickness_data[:,0],thickness_data[:,1],
                                              thickness_data[:,2],thickness_data[:,6]]))
        if d1_2013:
            thickness_data = transpose(array([thickness_data[:,0],thickness_data[:,4],
                                              thickness_data[:,5],thickness_data[:,9]]))
    else:
        thickness_data = genfromtxt(t_data, delimiter=',',skip_header=11)

    # measurement_time is the time for the GISAXS measurements starting at t=0 seconds.
    # the index of the vector signifies the measurement number
    measurement_time = gisaxs_data[:,1] - gisaxs_data[0,1]

    # at what index does thickness and gisaxs start?
    thickness_start = t_index
    gisaxs_start = g_index

    # these are the wli thickness measurements after being synced
    wli_thickness_time = thickness_data[thickness_start:-1,3] - thickness_data[thickness_start,3]
    wli_thickness = thickness_data[thickness_start:-1,1]

    # initialize thickness as a function of measurement time
    thickness = zeros((size(measurement_time),4))

    # argmin finds the index of the minimum value in a vector
    for i in arange(size(measurement_time)):
        thickness_arg = argmin(abs(wli_thickness_time - measurement_time[i]))
        thickness[i,1] = wli_thickness[int(thickness_arg)]
        thickness[i,0] = i + gisaxs_start
        thickness[i,2] = measurement_time[i]

    # make a 4th coloumn containing the samx value
    samx_data = open(g_data,'r')
    samx = [];
    for line in samx_data:
        if line[:4] == "#P1 ":
            if d1_2015:
                samx.append(float(line.split(' ')[4]))
            else:
                samx.append(float(line.split(' ')[3]))

    samx = asarray(samx)
    samx_data.close()
    thickness[:,3] = samx

    return thickness

# takes a thickness data file as created in gisaxs_thickness and plots it
# shows both GISAXS image number and thickness on two seperate x-axes
# thanks to http://stackoverflow.com/questions/10514315/how-to-add-a-second-x-axis-in-matplotlib
# requires the data to be continuous (e.g. one thickness data file).
# Otherwise the GISAXS map x-axis is broken
def plot_gisaxs_thickness(thickness,index_phase=False,plot_phi=False, save_figure=False):

    # tight axis
    xlim(0,thickness[-1,2]/60)

    # initial and max thickness found from data
    d_initial = min(thickness[:,1])
    d_max = max(thickness[:,1])

    if plot_phi:
        y_axis = d_initial/thickness[:,1]
        ylabel("$\phi$")
    else:
        y_axis = thickness[:,1]
        ylabel("Thickness [nm]")

    plot(thickness[:,2]/60,y_axis,'.')
    xlabel("Time [minutes]")

    if index_phase:
        y_axis_limit = ylim()
        for i in index_phase:
            plot((thickness[i,2]/60, thickness[i,2]/60), y_axis_limit, 'r-')
        ylim(y_axis_limit)

    disp('d_dry = ' + str(d_initial) + ' nm')
    disp('d_max = ' + str(d_max) + ' nm')
    disp('phi_min = ' + str(d_initial/d_max))

    if save_figure:
        savefig(save_figure + '.png', bbox_inches='tight', dpi=160)
    else:
        show()

# show a gisaxs image given a string dictating the directory, filename, experiment details, and data
# experiment is a vector with details of the experiment: [SD, pixel_size, db_x, db_y, alpha_i, wavelength]
# p_extent is the extent of the pixel area [x1, y1, x2, y2]
# imgtitle is the title of the plot
# colorlim defines the intensity limits of the plot (changes the colors of the plot quite a lot)
def show_gisaxs_image(dir, filename, experiment, p_extent=[0,0,-1,-1], imgtitle='',
                      colorlim=(10.0,3e3), save_figure=False):
    [x1, y1, x2, y2] = p_extent
    im = imread(dir + filename)
    im = flipud(im)
    im = im[y1:y2,x1:x2]
    im = flipud(im)

    [qy1, qz1] = gisaxs_get_q(x1,y1,experiment)
    [qy2, qz2] = gisaxs_get_q(x2,y2,experiment)

    imgplot = imshow(im,norm=matplotlib.colors.LogNorm(),clim=colorlim,
                     extent=[qy1,qy2,qz1,qz2], aspect='auto')
    colorbar()

    xlabel('q_y [nm^(-1)]')
    ylabel('q_z [nm^(-1)]')
    title(imgtitle)

    if save_figure:
        savefig(save_figure + '.png', bbox_inches='tight', dpi=160)

    show()

# uses thickness data to return a string detailing the thickness and time for a specified GISAXS map
# well suited for making titles for show_gisaxs_image
# num: image number
# data: thickness data created by gisaxs_thickness()
# data_offset: the offset (index) that sets t=0
# state: a string describing the phase of the experiment (typically 'swelling' or 'drying')
def get_gisaxs_image_title(num, data, data_offset=0, state='experiment'):
    index = where(data[:,0] == num)
    time = (data[index,2] - data[data_offset,2])/60
    samx = data[index,3]

    title = str(num) + '.tif' + ' | ' + "%.1f" % data[index,1] + ' nm' \
            + ' | ' + "%.1f" % time + ' min (' + str(state) + ') | samx=' + "%.1f" % samx

    return title

# shows a 2d color plot from dpdak data
# data is 2d color plot data from a dpdak export (x-axis and data)
# axis is either "qz" or "qy" (this only changes the label on the plot)
def show_2d_colorplot(cp_data, th_data, axis='qz', colorlim=(10.0,5e2),samx_values='all'):
    colorplot_data = genfromtxt(cp_data, skip_header=4)
    num_inputs = size(colorplot_data,1)/2
    x_axis = colorplot_data[:,0]
    im = colorplot_data[:,num_inputs:-1]

    #disp(size(im[0,:]))
    #disp(size(th_data[:,3]))

    if samx_values == 'all':
        cond = ones(size(th_data[:,3]), dtype=bool)
    else:
        cond = zeros(size(th_data[:,3]), dtype=bool)
        for i in samx_values:
            cond = cond | (th_data[:,3] == i)

    indices = transpose(where(cond))
    indices = indices[:-1,0]
    im = im[:,indices]

    time_frames = size(im,1)

    # Flip the data upside down. This is needed due to the way matrices are ordered.
    imgplot = imshow(flipud(im),norm=matplotlib.colors.LogNorm(),clim=colorlim,
                     extent=[1,time_frames,x_axis[0],x_axis[-1]], aspect='auto',)
    if axis == 'qz':
        ylabel('q_z [nm^(-1)]')
    if axis == 'qy':
        ylabel('q_y [nm^{-1}]')

    title('integration along ' + axis)
    xlabel('time [frames]')
    colorbar()
    show()

# show phi and peakfit data on the same graph.
# when defining the model in DPDAK, the model has to be Lorenzian or Gaussian
# and be defined in the first line of the model box.
# That way, the first three coloums of peakfit_data are pos, height, fwhm
# To obtain the data from DPDAK use DB Export and ONLY select Peak Fit->fit_param
def show_phi_peakfit(data, peakfit_data, fwhm_cond=[0.01,0.05], samx_values='all', axis='q_z',
                     peakfit_axis='fwhm', pos_cond=[0.32,0.40], save_figure=False):
    peakfit = genfromtxt(peakfit_data, skip_header=4)

    # conditions for fit
    cond1 = peakfit[2,:] > fwhm_cond[0] #fwhm
    cond2 = peakfit[2,:] < fwhm_cond[1] #fwhm

    cond3 = peakfit[0,:] > pos_cond[0] # pos
    cond4 = peakfit[0,:] < pos_cond[1] # pos


    if samx_values=='all':
        indices = where(cond1 & cond2 & cond3 & cond4)
    else:
        cond5 = zeros(size(data[:,3]), dtype=bool)
        for i in samx_values:
           cond5 = cond5 | (data[:,3] == i)
        indices = where(cond1 & cond2 & cond3 & cond4 & cond5)

    # draw everything
    fig, ax1 = subplots()
    t = transpose(data[indices,2])/60
    fwhm = transpose(peakfit[2,indices])
    pos = transpose(peakfit[0,indices])
    height = transpose(peakfit[1,indices])

    # output when we first see a peak that meet our conditions
    disp('first point with condition fwhm in range ' + str(fwhm_cond) + ' at t=' + str(data[indices[0][0],2]/60)
         + ', phi=' + str(min(data[:,1])/data[indices[0][0],1]))

    if peakfit_axis == 'fwhm':
        ax1.plot(t, fwhm, 'bo')
        ax1.set_xlabel('time (minutes)')
        if axis == 'q_z':
            ax1.set_ylabel('w_z  [nm^(-1)]', color='b')

        if axis == 'q_y':
            ax1.set_ylabel('w_y  [nm^(-1)]', color='b')

    if peakfit_axis == 'pos':
        ax1.plot(t, pos, 'bo')
        ax1.set_xlabel('time (minutes)')
        ax1.set_ylabel(axis + " [nm^(-1)]", color='b')

    if peakfit_axis == 'height':
        ax1.plot(t, height, 'bo')
        ax1.set_xlabel('time (minutes)')
        ax1.set_ylabel('peak intensity (arb. units)', color='b')

    if peakfit_axis == 'domain':
        if axis == 'q_z':
            disp('WARNING: Domain spacing should be found along q_y')

        disp('hi')
        ax1.plot(t,2*pi/pos,'bo')
        ax1.set_xlabel('time (minutes)')
        ax1.set_ylabel('lattice constant [nm]', color='b')

    # tight x-axis
    ax1.set_xlim(0,data[-1,2]/60)

    # Make the y-axis label and tick labels match the line color.
    for tl in ax1.get_yticklabels():
        tl.set_color('b')

    ax2 = ax1.twinx()
    t2 = data[:,2]/60
    phi = min(data[:,1])/data[:,1]
    ax2.plot(t2, phi, 'g.')
    ax2.set_ylabel('$\phi$', color='g')
    for tl in ax2.get_yticklabels():
        tl.set_color('g')

    if save_figure:
        savefig(save_figure + '.png', bbox_inches='tight', dpi=160)
    else:
        show()

# get (qy,qz) at a pixel_x, pixel_y
# experiment vector contains details of experiment
# experiment = [SD, pixel_size, db_x, db_y, alpha_i, wavelength]
def gisaxs_get_q(pixel_x, pixel_y, experiment):
    [SD, pixel_size, db_x, db_y, alpha_i, wavelength] = experiment

    distance_z = (pixel_y-db_y)*pixel_size
    alpha_f = arctan(distance_z/SD) - alpha_i
    qz = 2*pi/wavelength*(sin(alpha_i) + sin(alpha_f))

    distance_y = (pixel_x-db_x)*pixel_size
    psi = arctan(distance_y/SD)
    qy = 2*pi/wavelength*sin(psi)*cos(alpha_f)
    return [qy, qz]

def load_gisaxs_data(filename):
    datafile = load(filename)
    data = datafile['arr_0']
    datafile.close()
    return data

""" Functions to use D1 data with GenX
    reformat_data_for_genx() reformats the raw specuser data to a readable format for GenX
    furthermore the x-axis is converted to q_z

    display_genx_data simply takes the exports from GenX and displays them while outputting
    fitting parameters to the console. In the Case of SLD it is converted to units of Angstroms """
def display_genx_data(data_file, value_file, sld_file):
    # some values we need
    f_C = 6.01131
    f_H = 0.999979
    f_O = 8.03236
    f_Si = 14.1795
    r_e = 2.81794033e-5

    data = genfromtxt(data_file)
    values = genfromtxt(value_file, delimiter='\t')
    sld = genfromtxt(sld_file)
    x = data[:,0]*10
    I = data[:,2]
    Isim = data[:,1]

    # plot reflectivity and fit
    figure(figsize=(8,6), dpi=160)
    semilogy(x,I,'.', label='data', markersize=10)
    semilogy(x,Isim,'r-',label='model', linewidth=2)
    xlabel('q_z [nm^{-1}]',fontsize='x-large')
    ylabel('R',fontsize='x-large')
    legend(loc='upper right', fontsize='x-large')
    autoscale(enable=True,axis=u'both',tight=True)
    savefig(data_file + '_ref.png', bbox_inches='tight', dpi=160)
    show()

    # plot sld
    figure(figsize=(8,6), dpi=160)
    plot(sld[:,0]*0.1,sld[:,1], linewidth=2)
    xlabel('z [nm]',fontsize='x-large')
    ylabel('SLD [$r_e/ \AA^3}$]',fontsize='x-large')
    autoscale(enable=True,axis=u'both',tight=True)
    savefig(data_file + '_sld.png', bbox_inches='tight', dpi=160)
    show()

    f_polymer = 6*f_C + 7*f_H
    f_SiOx = 1*f_Si + 2*f_O
    polymer_thickness = values[1,1]
    polymer_sld = values[2,1]*f_polymer*r_e
    polymer_sigma = values[3,1]
    siox_thickness = values[4,1]
    siox_sld = values[5,1]*f_SiOx*r_e
    siox_sigma = values[6,1]
    si_sld = values[7,1]*f_Si*r_e

    disp('polymer thickness: ' + str(polymer_thickness) + ' AA')
    disp('polymer SLD: ' + str(polymer_sld) + ' AA^(-2)')
    disp('polymer roughness: ' + str(polymer_sigma))
    disp('SiOx thickness: ' + str(siox_thickness) + ' AA')
    disp('SiOx SLD: ' + str(siox_sld) + ' AA^(-2)')
    disp('SiOx roughness: ' + str(siox_sigma))
    disp('Si SLD: ' + str(si_sld) + ' AA^(-2)')

    # if we want to compare sld profiles of different measurements
    return sld

def reformat_data_for_genx(datafile,Iback,Idirect,seconds,wavelength):
    c4r1 = genfromtxt(datafile)
    angle = c4r1[:,0]
    Imon = c4r1[:,4]/seconds
    Idet = c4r1[:,5]/seconds
    q_z = 4*pi*sin(angle*degree)/wavelength

    R = (Idet/Imon - Iback)/(Idirect)
    #semilogy(q_z,R)
    #xlabel('angle (degrees)')
    #ylabel('intensity (arb. u.)')

    savetxt(datafile + '.out', transpose([q_z, R]))
