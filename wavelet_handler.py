# standard modules
import os
import sys
# Project modules
import scripts.coimbra2024_scripts as coimbra2024
import scripts.wavelet_functions as wavelet_functions

def main(sitename, inputpath, outputpath, datetimerange):
    # Create setup
    configure = coimbra2024.structuredData()

    # Select output file path
    configure.output_path = str(os.path.join(outputpath, 'data', str(sitename)+'_CDWT{}_{}.{}mn.csv'))

    # Select raw input file path
    # e.g.: "<PROJECT FOLDER>/eddypro_output/eddypro_raw_datasets/level_6/"
    configure.raw_kwargs = {'path': str(os.path.join(inputpath, 'eddypro_output', 'eddypro_raw_datasets/level_6'))}

    # Select covariances
    # x*y → Cov(x, y)
    # x*y*z*... → Cov(x, y)|Cov(x, z),Cov(x, ...)
    configure.varstorun = ['co2*co2', 'h2o*h2o', 'ts*ts', 'w*co2*h2o', 'w*h2o*co2', 'w*ts*co2']

    # Select period of interest
    # [START_DATE, END_DATE, FILE_FREQUENCY]
    configure.ymd = [datetimerange.split('-')[0], datetimerange.split('-')[1], '30min']

    # Select wavelet method
    configure.method = 'dwt'

    # RUN WAVELET FLUX PROCESSING
    wavelet_functions.run_wt(**vars(configure), verbosity=5)

    # Merge into a single file
    coimbra2024.concat_into_single_file(
        outputpath, str(sitename)+'_CDWT_full_cospectra.+.30mn.csv', 
        output_path=os.path.join(outputpath, str(sitename)+'_CDWT_full_cospectra_unique.30mn.csv'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('sitename',   type=str)
    parser.add_argument('inputpath',  type=str)
    parser.add_argument('outputpath', type=str)
    parser.add_argument('datetimerange', type=str)
    parser.add_argument('-s', '--sitename',   type=str, required=False)
    parser.add_argument('-i', '--inputpath',  type=str, required=False)
    parser.add_argument('-o', '--outputpath', type=str, required=False)
    parser.add_argument('-d', '--datetimerange', type=str, required=False)
    args = parser.parse_args()

    main(**vars(args))
    #print(vars(args))
