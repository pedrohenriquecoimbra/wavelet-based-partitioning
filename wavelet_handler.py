# standard modules
import os
import sys
# Project modules
import scripts.coimbra2024_scripts as coimbra2024
import scripts.wavelet_functions as wavelet_functions

def main(inputpath, outputpath, datetimerange):
    # CONFIGURE
    coimbra2024.SITES_TO_STUDY = ['Tamu']

    # Create setup
    configure = coimbra2024.structuredData()

    # mother folder (optional)
    configure.mother = 'C:/Users/phherigcoimb/OneDrive/INRAe/thesis-project-1/data/flux/AmeriFlux/Tamu/run_Sun/'

    # Select output file path
    configure.output_path = configure.mother + 'output/data/Tamu_CDWT{}_{}.{}mn.csv'

    # Select raw input file path
    # e.g.: "<PROJECT FOLDER>/eddypro_output/eddypro_raw_datasets/level_6/"
    configure.raw_kwargs = {'path': configure.mother + 'eddypro_output_folder/eddypro_raw_datasets/level_6'}

    # Select covariances
    # x*y → Cov(x, y)
    # x*y*z*... → Cov(x, y)|Cov(x, z),Cov(x, ...)
    configure.varstorun = ['w*co2*h2o', 'w*ts*n2o', 'w*n2o*h2o', 'w*h2o*n2o', 'w*h2o*ts', 'w*ts*h2o']

    # Select period of interest
    # [START_DATE, END_DATE, FILE_FREQUENCY]
    configure.ymd = ['202205250030', '202206020000', '30min']

    # Select wavelet method
    configure.method = 'dwt'

    configure.raw_kwargs.update({'fmt': {'n2o': '4th'}, 
                                'fkwargs': {'date_format': '%Y%m%d-%H%M', 'dt': 0.05, 'file_pattern': '([0-9]{8}-[0-9]{4})_raw_dataset_.*.txt'},
                                'sep': '\s+', 'skiprows': 8, 'na_values': ['NaN', 'nan', -9999]}, 
                                )

    # RUN WAVELET FLUX PROCESSING
    rdw = wavelet_functions.run_wt(**vars(configure), verbosity=5)

    # Merge into a single file
    coimbra2024.concat_into_single_file(
        os.path.dirname(configure.output_path), 'Tamu_CDWT_full_cospectra.+.30mn.csv', 
        output_path=os.path.join(os.path.dirname(configure.output_path), 'Tamu_CDWT_full_cospectra_unique.30mn.csv'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('inputpath', type=str)
    parser.add_argument('outputpath', type=str)
    parser.add_argument('datetimerange', type=str)
    parser.add_argument('-i', '--inputpath', type=str, required=False)
    parser.add_argument('-o', '--outputpath', type=str, required=False)
    parser.add_argument('-d', '--datetimerange', type=str, required=False)
    args = parser.parse_args()
    """
    print("sys.argv", sys.argv)
    args = [a for a in sys.argv if '=' not in a]
    kwargs = dict([a.split('=') for a in sys.argv if '=' in a])
    """
    #main(*args, **kwargs)
    print(vars(args))
