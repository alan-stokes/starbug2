from astropy import units as u
from astropy.io import fits
from astropy.table import Table
import numpy as np

# Set values for file names and filters in the data
galaxy_target = 'NGC6822'
observed_filters = ['F115W', 'F200W', 'F356W', 'F444W', 'F770W', 'F1000W', 'F1500W', 'F2100W']

in_dir = './'
out_dir = './'

input_starbug_mag_file = "/Users/olivia.jones/Science/NGC6822_JWST/NGC6822_YSOs/NGC6822_candidateYSOs.fits"
output_file_name = out_dir + galaxy_target + "_YSO_Candidate_Catalog.csv"

# No more user modifications should be needed beyond this point.

# AB - Vega magnitude offsets
offsets = {
    "F070W": ["CLEAR", 0.2497600018978119],
    "F090W": ["CLEAR", 0.5041199922561646],
    "F115W": ["CLEAR", 0.7828400135040283],
    "F140M": ["CLEAR", 1.1191799640655518],
    "F150W": ["CLEAR", 1.2431399822235107],
    "F150W2": ["CLEAR", 1.2230199575424194],
    "F182M": ["CLEAR", 1.5876799821853638],
    "F187N": ["CLEAR", 1.6553699970245361],
    "F200W": ["CLEAR", 1.705970048904419],
    "F210M": ["CLEAR", 1.8102200031280518],
    "F212N": ["CLEAR", 1.8312100172042847],
    "F250M": ["CLEAR", 2.1488699913024902],
    "F277W": ["CLEAR", 2.315419912338257],
    "F300M": ["CLEAR", 2.490180015563965],
    "F322W2": ["CLEAR", 2.570039987564087],
    "F335M": ["CLEAR", 2.71697998046875],
    "F356W": ["CLEAR", 2.823899984359741],
    "F360M": ["CLEAR", 2.870759963989258],
    "F410M": ["CLEAR", 3.107640027999878],
    "F430M": ["CLEAR", 3.212130069732666],
    "F444W": ["CLEAR", 3.2418100833892822],
    "F460M": ["CLEAR", 3.3788299560546875],
    "F480M": ["CLEAR", 3.4422800540924072],
    "F560W": ["CLEAR", 3.76066],
    "F770W": ["CLEAR", 4.38398],
    "F1000W": ["CLEAR", 4.95551],
    "F1130W": ["CLEAR", 5.49349],
    "F1280W": ["CLEAR", 5.24056],
    "F1500W": ["CLEAR", 5.83929],
    "F1800W": ["CLEAR", 6.22752],
    "F2100W": ["CLEAR", 6.53267],
    "F2550W": ["CLEAR", 6.96805]
}


def convert_vegamag_to_mJy(yso_candidates, o_filter):
    # Convert magnitudes from Vega to AB system
    flt_offsets = offsets[o_filter][1]
    ab_mag = + yso_candidates[o_filter] + flt_offsets

    # Convert fluxes [AB magnitude] to flux density [Jansky]
    f_Jy = 10 ** ((23.9 - ab_mag) / 2.5) * 10 ** -6 * u.Jansky

    # Convert flux uncertainties [AB magnitudes] to flux density uncertainty [Jansky]
    error_Filter_data_name = 'e' + o_filter
    ef_Jy = np.log(10) / 2.5 * yso_candidates[error_Filter_data_name] * f_Jy

    # Convert both flux and error values to mJy and no units
    # f_mJy = f_Jy.value * 1000
    # ef_mJy = ef_Jy.value * 1000

    return f_Jy, ef_Jy


def make_yso_fitter_data_file(cat_num, ra, dec, fluxes, output_file_name):
    # Open file for writing sources fluxes in mJy - fitter data format
    yso_cat = open(output_file_name, 'w')
    cntr = 1

    for ii in range(0, len(cat_num)):
        yso_cat.write('Y' + str(cntr) + '_' + str(cat_num[ii]) + ' ' + str(ra[ii]) + ' ' + str(dec[ii]) + ' ')

        # Add in quality flags to file: 0 = nan 1 = data
        for flux_value in fluxes:
            if np.isnan(flux_value[0][ii]):
                yso_cat.write('0 ')
            else:
                yso_cat.write('1 ')

        for flux_value in fluxes:
            yso_cat.write(str(flux_value[0][ii] * 1000) + ' ' + str(flux_value[1][ii] * 1000) + ' ')
        yso_cat.write('\n')

        cntr = cntr + 1
    yso_cat.close()


def main():
    # Read in YSO candidates photometry
    with fits.open(input_starbug_mag_file) as hdu:
        yso_candidates = Table(hdu[1].data)
        fluxes = []

        for o_filter in observed_filters:
            f_jy, ef_jy = convert_vegamag_to_mJy(yso_candidates, o_filter)
            fluxes.append([f_jy.value, ef_jy.value])

            cat_num = yso_candidates['Catalogue_Number']
            ra = yso_candidates['RA']
            dec = yso_candidates['DEC']

        make_yso_fitter_data_file(cat_num, ra, dec, fluxes, output_file_name)


if __name__ == '__main__':
    main()
