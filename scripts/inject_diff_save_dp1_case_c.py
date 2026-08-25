import pandas as pd
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle
from astropy.io import fits
import astropy.units as u
from astropy.table import Table, vstack
import astropy
from itertools import combinations
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord
import getpass
import sys
import yaml
# import lsst.afw.display as afwDisplay
# import lsst.afw.geom as afwGeom
# import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.geom as geom
from scipy.ndimage import zoom
from lsst.daf.butler import Butler
from lsst.source.injection import ingest_injection_catalog, generate_injection_catalog
from lsst.source.injection import VisitInjectConfig, VisitInjectTask
from lsst.ip.diffim.subtractImages import AlardLuptonSubtractTask, AlardLuptonSubtractConfig
from lsst.meas.algorithms.detection import SourceDetectionTask, SourceDetectionConfig
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
from lsst.meas.deblender import SourceDeblendTask, SourceDeblendConfig
from lsst.meas.base import SingleFrameMeasurementTask, SingleFrameMeasurementConfig
from lsst.meas.base import ForcedMeasurementTask
from lsst.source.injection import CoaddInjectConfig, CoaddInjectTask

import astropy.visualization as aviz
import matplotlib
# matplotlib.use("AGG")
# Force matplotlib defaults
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
import matplotlib.pyplot as plt
from matplotlib import cm
import io
import argparse
def parse_args():
    p = argparse.ArgumentParser("c_inj_double")
    p.add_argument("-path_to_manifest", required=True)          # multi-run YAML
    p.add_argument("-outdir", required=False)                      # output root
    return p.parse_args()
args = parse_args()
# ----- Load collection the coadds

# load butler
butler_1st = Butler("/repo/DP1", collections="u/taceroc/dp1/custom_coadd_window1_60623-60633-t4848-p90-r")
deep_coadd_ref_1st = list(butler_1st.registry.queryDatasets('deep_coadd_predetection'))
butler_2nd = Butler("/repo/DP1", collections="u/taceroc/dp1/custom_coadd_window2_60651-60655-4848-90-r")
deep_coadd_ref_2nd = list(butler_2nd.registry.queryDatasets('deep_coadd_predetection'))


def generate_injection_2nd(butler, deep_coadd_ref, which_fits, return_catalog=True, suffix='', seed='3'):
    
    # ---- Define injection catalog for 1st image
     # I still dont know how to 'remove' the coadds with empty
    coadd = butler.get(deep_coadd_ref[0])
    boxcen = coadd.getBBox().getCenter()
    wcs_b = coadd.getWcs()
    cen = wcs_b.pixelToSky(boxcen)
    radec = SkyCoord(ra=cen[0].asDegrees()*u.deg, dec=cen[1].asDegrees()*u.deg)
    
    print(radec)

    if return_catalog:
        # load metadata of LE
        imsize = coadd.getBBox().getDimensions()[0]*wcs_b.getPixelScale().asDegrees()
        print('Size of calexp in degrees: ', imsize)
        
        inject_size = imsize/2
        print(which_fits)
        hdul = fits.open(which_fits)
        fits_img = hdul[0].data
        print(fits_img.shape)
        fits_img[np.isnan(fits_img)] = 0
        
        # zero_cols = np.where((fits_img == 0).all(axis=0))[0]
        # cols = np.arange(fits_img.shape[-1])
        # mask_cols = [False if x in zero_cols else True for x in cols]
        # fits_img_zero_cols = fits_img[:, mask_cols]
        # zero_row = np.where((fits_img_zero_cols == 0).all(axis=1))[0]
        # rows = np.arange(fits_img_zero_cols.shape[0])
        # mask_rows = [False if x in zero_row else True for x in rows]
        # fits_img = fits_img_zero_cols[mask_rows, :]
        # del fits_img_zero_cols 
        
        mags_plot = -2.5*np.log10(fits_img)-48.6
        mag_mean = np.mean(mags_plot[fits_img>0])
        flux_nJy = 10**(-0.4 * (mags_plot + 48.6))
     
        new_fits = flux_nJy
        new_fits[np.isnan(new_fits)] = 0
        hdul[0].data = new_fits

        hdul.writeto(f'output_{suffix}.fits', overwrite=True)
        hdul.close()
        mags_plot = -2.5*np.log10(new_fits)-48.6
        mags_plot = np.nan_to_num(mags_plot, nan=0.0, posinf=0.0, neginf=0.0)
        print(new_fits.shape)
        ns = imsize/((mags_plot.shape[0])*wcs_b.getPixelScale().asDegrees())
        Ns = int((ns*ns) - 10)
        print(Ns)
        size_img = mags_plot.shape[0]
        scale_fits = False
        if Ns >= 5:
            Ns = 4
            # size_scale = False
            size_scale = imsize/(ns*wcs_b.getPixelScale().asDegrees())
        elif Ns <= 0:
            Ns = 4
            ns = np.sqrt((Ns+1+10))
            size_scale = imsize/(ns*wcs_b.getPixelScale().asDegrees())
            scale_factor =  size_scale/mags_plot.shape[0]
            # scale_factor = 0.5  # 0.5 = half size, 2.0 = double size
            scaled_data = zoom(new_fits, scale_factor, order=1) # order=1 for linear interpolation
            
            # Save the scaled stamp to a new file
            hdul[0].data = scaled_data
            hdul.writeto(f'output_{suffix}.fits', overwrite=True)
            hdul.close()

            print(scaled_data.shape, size_scale, scale_factor)
            scale_fits = scale_factor
            size_img = scaled_data.shape[0]
        else:
            size_scale = imsize/(ns*wcs_b.getPixelScale().asDegrees())
        print("mean surface", np.mean(mags_plot[new_fits>0]))

        
        # if scale_fits != False:
        if size_img >= 700:
            mag = [np.mean(mags_plot[new_fits>0]) if np.mean(mags_plot[new_fits>0])<=18 else 18]
        else:
            mag = [np.mean(mags_plot[new_fits>0]) if np.mean(mags_plot[new_fits>0])<=20 else 20]


        if suffix == '2nd':
            mag =  [mag[0] + 0.5]
        # ---- generate injection of LEs for 1st LE
        my_injection_catalog_LEs = generate_injection_catalog(
            ra_lim=[radec.ra.value-inject_size, radec.ra.value+inject_size],
            dec_lim=[radec.dec.value-inject_size, radec.dec.value+inject_size],
            number=Ns,
            seed=seed,
            source_type= "Stamp",
            mag=mag,
            stamp= [f'output_{suffix}.fits'],
        )
        
        # first_ct = my_injection_catalog_LEs[0]['stamp'].split('ct')[1].split('_')[0]

        mag_source = np.mean(mags_plot[new_fits>0])-1
        if mag_source <= 18:
            mag_source = np.mean(mags_plot[new_fits>0])
    
        return my_injection_catalog_LEs, coadd, wcs_b, radec, inject_size, Ns, mag_source, size_scale, scale_fits
    else:
        return coadd

def generate_injection_1st(butler, deep_coadd_ref, which_fits, return_catalog=True, suffix='', Ns=5, scale_fits=False, seed='3'):
    
    # ---- Define injection catalog for 1st image
     # I still dont know how to 'remove' the coadds with empty
    coadd = butler.get(deep_coadd_ref[0])
    boxcen = coadd.getBBox().getCenter()
    wcs_b = coadd.getWcs()
    cen = wcs_b.pixelToSky(boxcen)
    radec = SkyCoord(ra=cen[0].asDegrees()*u.deg, dec=cen[1].asDegrees()*u.deg)
    
    print(radec)

    if return_catalog:
        # load metadata of LE
        imsize = coadd.getBBox().getDimensions()[0]*wcs_b.getPixelScale().asDegrees()
        print('Size of calexp in degrees: ', imsize)
        
        inject_size = imsize/2
        print(which_fits)
        hdul = fits.open(which_fits)
        fits_img = hdul[0].data
        print(fits_img.shape)
        fits_img[np.isnan(fits_img)] = 0
        
        # zero_cols = np.where((fits_img == 0).all(axis=0))[0]
        # cols = np.arange(fits_img.shape[-1])
        # mask_cols = [False if x in zero_cols else True for x in cols]
        # fits_img_zero_cols = fits_img[:, mask_cols]
        # zero_row = np.where((fits_img_zero_cols == 0).all(axis=1))[0]
        # rows = np.arange(fits_img_zero_cols.shape[0])
        # mask_rows = [False if x in zero_row else True for x in rows]
        # fits_img = fits_img_zero_cols[mask_rows, :]
        # del fits_img_zero_cols 
        
        mags_plot = -2.5*np.log10(fits_img)-48.6
        mag_mean = np.mean(mags_plot[fits_img>0])
        flux_nJy = 10**(-0.4 * (mags_plot + 48.6))
     
        new_fits = flux_nJy
        new_fits[np.isnan(new_fits)] = 0
        hdul[0].data = new_fits

        hdul.writeto(f'output_{suffix}.fits', overwrite=True)
        hdul.close()
        mags_plot = -2.5*np.log10(new_fits)-48.6
        mags_plot = np.nan_to_num(mags_plot, nan=0.0, posinf=0.0, neginf=0.0)
        print(new_fits.shape)
        size_img = new_fits.shape[0]
        if scale_fits != False:
            scaled_data = zoom(new_fits, scale_fits, order=1) # order=1 for linear interpolation
            
            # Save the scaled stamp to a new file
            hdul[0].data = scaled_data
            hdul.writeto(f'output_{suffix}.fits', overwrite=True)
            hdul.close()
            size_img = scaled_data.shape[0]
            print(scaled_data.shape, scale_fits)

        if size_img >= 700:
            mag = [np.mean(mags_plot[new_fits>0]) if np.mean(mags_plot[new_fits>0])<=18 else 18]
        else:
            mag = [np.mean(mags_plot[new_fits>0]) if np.mean(mags_plot[new_fits>0])<=20 else 20]

        if suffix == '2nd':
            mag =  [mag[0] + 0.5]
        # ---- generate injection of LEs for 1st LE
        my_injection_catalog_LEs = generate_injection_catalog(
            ra_lim=[radec.ra.value-inject_size, radec.ra.value+inject_size],
            dec_lim=[radec.dec.value-inject_size, radec.dec.value+inject_size],
            number=Ns,
            seed=seed,
            source_type= "Stamp",
            mag=mag,
            stamp= [f'output_{suffix}.fits'],
        )
        
        # first_ct = my_injection_catalog_LEs[0]['stamp'].split('ct')[1].split('_')[0]

        mag_source = np.mean(mags_plot[new_fits>0])-1
        if mag_source <= 18:
            mag_source = np.mean(mags_plot[new_fits>0])
    
        return my_injection_catalog_LEs, coadd, wcs_b, radec, inject_size, Ns, mag_source
    else:
        return coadd


def do_injections(coadd, injection_catalog):#my_injection_catalog_LEs, my_injection_catalog_source):
    inject_config = CoaddInjectConfig()
    inject_task = CoaddInjectTask(config=inject_config)

    psf = coadd.getPsf()
    photo_calib = coadd.getPhotoCalib()
    wcs_b = coadd.getWcs()
    
    injected_output = inject_task.run(
        injection_catalogs=[*injection_catalog],#[my_injection_catalog_LEs, my_injection_catalog_source],
        input_exposure=coadd.clone(),
        psf=psf,
        photo_calib=photo_calib,
        wcs=wcs_b,
    )
    injected_coadd_1st = injected_output.output_exposure
    injected_catalog = injected_output.output_catalog


    return injected_coadd_1st

        

def do_source_detection_injections(injected_coadd_2nd):

    schema = afwTable.SourceTable.makeMinimalSchema()
    print(schema)
    raerr = schema.addField("coord_raErr", type="F")
    decerr = schema.addField("coord_decErr", type="F")

    schema.addField("sky_source", type="Flag", doc="Sky background value at source position")
    
    config = CharacterizeImageConfig()
    config.psfIterations = 3
    charImageTask = CharacterizeImageTask(config=config)
    del config
    
    config = SourceDetectionConfig()
    config.thresholdValue = 5
    sourceDetectionTask = SourceDetectionTask(schema=schema, config=config)

    config = SourceDeblendConfig()
    sourceDeblendTask = SourceDeblendTask(schema=schema, config=config)
    
    config = SingleFrameMeasurementConfig()
    sourceMeasurementTask = SingleFrameMeasurementTask(schema=schema,
                                                       config=config)
    
    result = charImageTask.run(injected_coadd_2nd)

    try:
        # Try to get sky from the exposure
        sky_value = result.exposure.getSkyBackground().getBackground()
    except:

        image = injected_coadd_2nd.getMaskedImage().getImage()
        sky_value = np.median(image.getArray())


    tab = afwTable.SourceTable.make(schema)
    result = sourceDetectionTask.run(tab, injected_coadd_2nd)
    sources = result.sources

    # Mark all sources as NOT sky sources (sky_source = False)
    for record in sources:
        record.set('sky_source', False)
        
    sourceDeblendTask.run(injected_coadd_2nd, sources)
    sourceMeasurementTask.run(measCat=sources, exposure=injected_coadd_2nd)
    sources = sources.copy(True)

    return sources


def image_subtraction(injected_coadd_1st, injected_coadd_2nd, sources):
    config = AlardLuptonSubtractConfig()
    config.sourceSelector.doRequirePrimary = False
    config.sourceSelector.doSkySources = False
    alTask = AlardLuptonSubtractTask(config=config)
    
    # result = alTask.run(injected_coadd_1st, injected_coadd_2nd, src_catalog[0])
    # template, science
    result = alTask.run(injected_coadd_1st, injected_coadd_2nd, sources)
    return result



def plot_one_image(ax, data, size, scale, name=None):
    """Plot a normalized image on an axis."""
    if name != None:
        norm = aviz.ImageNormalize(
            # focus on a rect of dim 15 at the center of the image.
            data[data.shape[0] // 2 - 7-2:data.shape[0] // 2 + 8-2,
                 data.shape[1] // 2 - 7-2:data.shape[1] // 2 + 8-2],
            interval=aviz.MinMaxInterval(),
            stretch=aviz.AsinhStretch(a=0.1),
        )
    else:
        norm = aviz.ImageNormalize(
            data[data.shape[0] // 2 - 7-2:data.shape[0] // 2 + 8-2,
                 data.shape[1] // 2 - 7-2:data.shape[1] // 2 + 8-2],
            interval=aviz.MinMaxInterval(),
            stretch=aviz.AsinhStretch(a=0.1),
            # stretch=aviz.ZScaleInterval(),
        )
    im = ax.imshow(data, cmap=cm.bone, interpolation="none", norm=norm,
              extent=(0, size, 0, size), origin="lower", aspect="equal")
    x_line = 1
    y_line = 1
    # ax.plot((x_line, x_line + 1.0/scale), (y_line, y_line), color="blue", lw=6)
    # ax.plot((x_line, x_line + 1.0/scale), (y_line, y_line), color="yellow", lw=2)
    ax.axis("off")
    if name is not None:
        ax.set_title(name)

    return im

path_to_manifest = args.path_to_manifest #'runs/runs_082126/manifest.yml'
general_path = '/'.join(path_to_manifest.split('/')[0:-1])    
params =[
    'so_d_ly',
    'dz0_ly',
    'ct_years',
    'a',
    'ay',
    'az',
    'z0ly',
    'angles_deg',
    'wavel',
    'dust_env',
    'composition_s_c']
## ----- load data
def create(manifest_path_ix):
    with open(manifest_path_ix, "r") as f:
        data = yaml.safe_load(f)
    df = pd.DataFrame(data)
    df['outputs'] = df['outputs'].apply(lambda x: x[0])
    df['ct_years'] = df['meta'].apply(lambda x: x['ct_years'])
    for pix in params:
        df[pix] = df['meta'].apply(lambda x: x['params'][pix])
    df.drop(columns=['meta'])

    return df
df_all = pd.DataFrame()
df = create(path_to_manifest)
df_all = pd.concat([df_all, df])

df_all = df_all.sort_values(by=['ct_years'])


df_to_use_to_plot_single_double = df_all.groupby(['so_d_ly',  'a', 'ay',
       'az', 'z0ly', 'wavel', 'dz0_ly']).agg({'ct_years': list, 'angles_deg': list, 'outputs':list})

# Get c) plot t1-t0 for each system, you need two LE, one at t0 and other at t1 for the same system
def create_pairs_dtimeonly(ixy):
    
    sorted_lst = list(enumerate(df_to_use_to_plot_single_double.iloc[ixy]['ct_years']))
    sorted_lst.sort(key=lambda x: x[1])
    
    # Create pairs of adjacent elements
    pairs = []
    index = []
    for i in range(len(sorted_lst) - 1):
        pairs.append((sorted_lst[i][1], sorted_lst[i+1][1]))
        index.append((sorted_lst[i][0], sorted_lst[i+1][0]))
        
    result = [
        {
            "pair": p1, 
            "indices": idx1
        } 
        for p1, idx1 in zip(pairs, index)
    ]

    df_one = pd.DataFrame(index=range(len(result)), columns=['pair_time_ct_years','outputs', 'name_npy'])
    for ix, pp in enumerate(result):
        index = pp['indices']
        
        o_array = np.array(df_to_use_to_plot_single_double.iloc[ixy]['outputs'])
        df_one.loc[ix, 'pair_time_ct_years'] = pp['pair']
        # df_one.loc[ix, 'indices'] = pp['indices']
        df_one.loc[ix, 'outputs'] = o_array[np.array(index)]
        df_one.loc[ix, 'name_npy'] = ''
        
        # print(o_array[np.array(index)])
    return df_one

df_to_use_to_plot_single_double_save = pd.DataFrame()
if args.outdir:
    general_path = args.outdir


numpy_path = os.path.join(general_path, 'injections/cdouble/numpy')
images_path = os.path.join(general_path, 'injections/cdouble/images')

os.makedirs(numpy_path, exist_ok=True)
os.makedirs(images_path, exist_ok=True)


for ix in range(len(df_to_use_to_plot_single_double)):
    # print(ix)
    df_one = create_pairs_dtimeonly(ix)
    df_to_use_to_plot_single_double_save = pd.concat([df_to_use_to_plot_single_double_save, df_one])

df_to_use_to_plot_single_double_save = df_to_use_to_plot_single_double_save.reset_index(drop=True)
df_to_use_to_plot_single_double_save['name_npy'] = ''
df_to_use_to_plot_single_double_save['scaled'] = ''
import random
for ix, row in df_to_use_to_plot_single_double_save.iloc[36:37].iterrows():
    random_seed = random.getrandbits(16)
    print(f"Generated Seed: {random_seed}")
    print(row['outputs'], row['outputs'])
    my_injection_catalog_LEs = []
    which_fits_2nd = row['outputs'][1].replace('arrays', 'fits').replace('surface.npy', 'surface_image.fits')
    my_injection_catalog_LEsix_2nd, coadd_2nd, wcs_b, radec, inject_size, Ns, mag_source, size_scale_2nd, scale_fits = generate_injection_2nd(butler_2nd, 
                                                                                                                              deep_coadd_ref_2nd, 
                                                                                                                              which_fits_2nd,
                                                                                                                             suffix='2nd', seed=f'{random_seed}')
    which_fits_1st = row['outputs'][0].replace('arrays', 'fits').replace('surface.npy', 'surface_image.fits')

    my_injection_catalog_LEsix_1st, coadd_1st, wcs_b, radec, inject_size, Ns, mag_source = generate_injection_1st(butler_1st, 
                                                                                                                          deep_coadd_ref_1st, 
                                                                                                                          which_fits_1st, 
                                                                                                                         suffix='1st', 
                                                                                                                  Ns=Ns, scale_fits=scale_fits, 
                                                                                                                  seed=f'{random_seed}')
    
    df_to_use_to_plot_single_double_save.at[ix, 'scaled'] = [int(size_scale_2nd), int(size_scale_2nd)]
    my_injection_catalog_LEs.append(my_injection_catalog_LEsix_2nd)
    my_injection_catalog_LEs.append(my_injection_catalog_LEsix_1st)
    
    injected_coadd_2nd = do_injections(coadd_2nd, [my_injection_catalog_LEsix_2nd])
    injected_coadd_1st = do_injections(coadd_1st, [my_injection_catalog_LEsix_1st])

    sources = do_source_detection_injections(injected_coadd_2nd)
    subtraction_outputs = image_subtraction(injected_coadd_1st, injected_coadd_2nd, sources)


    path_name_ids = f"{deep_coadd_ref_1st[0].dataId['tract']}_{deep_coadd_ref_1st[0].dataId['patch']}_{deep_coadd_ref_1st[0].dataId['band']}_dp1"
    df_to_use_to_plot_single_double_save.loc[ix, 'coadd_name'] = path_name_ids
    path_name = '_'.join(which_fits_1st.replace('/fits/surface_image.fits', '').split('/'))
    path_name = path_name+'_'+'_'.join(which_fits_2nd.replace('/fits/surface_image.fits', '').split('/'))
        
    numpy_cutouts = {}
    name_npy_paths = []
    for iyy, row in enumerate(my_injection_catalog_LEs[-1]):
        try:
            center = wcs_b.skyToPixel(geom.SpherePoint(row['ra']*geom.degrees, row['dec']*geom.degrees))
            if int(size_scale_2nd)>=600:
                float_value = int(size_scale_2nd)/2
                value = random.random()
                if value>=0.5:
                    new_x = center.x + float_value
                    new_y = center.y #+ float_value
                else:
                    new_x = center.x #- float_value
                    new_y = center.y - float_value
                center = geom.Point2D(new_x, new_y)
            s = 600
            extent = geom.Extent2I(s, s)
            science_cutout = subtraction_outputs.matchedScience.getCutout(center, extent)
            template_cutout = injected_coadd_1st.getCutout(center, extent)
            difference_cutout = subtraction_outputs.difference.getCutout(center, extent)
            dia_source_id = iyy
            # self.numpy_path.mkdir(dia_source_id)
            numpy_cutouts[f"sci_{s}"] = science_cutout.image.array
            numpy_cutouts[f"temp_{s}"] = template_cutout.image.array
            numpy_cutouts[f"diff_{s}"] = difference_cutout.image.array
        
            scale = science_cutout.wcs.getPixelScale(science_cutout.getBBox().getCenter()).asArcseconds()
            fig, axs = plt.subplots(1, 3, constrained_layout=True)
            im1 = plot_one_image(axs[0], template_cutout.image.array, s, scale, "1st inj")
            im2 = plot_one_image(axs[1], science_cutout.image.array, s, scale, "2nd inj")
            imd = plot_one_image(axs[2], difference_cutout.image.array, s, scale, "Difference")
            
            # output = io.BytesIO()
            # plt.show()
            outfile_img = os.path.join(images_path, f"{path_name}_{dia_source_id}_{s}.png")
            plt.savefig(outfile_img, bbox_inches="tight", format="png")
            # output.seek(0)
            plt.close(fig)
            
            for cutout_type, cutout in numpy_cutouts.items():
                outfile = os.path.join(numpy_path, f'{path_name}_{dia_source_id}_{cutout_type}_{s}.npy')
                name_npy_paths.append(outfile)
                np.save(outfile, np.expand_dims(cutout, axis=0))
    
    
            for ax in axs:
                ax.clear()  # Clear the axis
                ax.remove()
        except Exception as ex:
            print('something', ex)
            continue
    
    df_to_use_to_plot_single_double_save.at[ix, 'name_npy'] = name_npy_paths
# df_to_use_to_plot_single_double_save.to_csv(f'a_df_to_use_to_plot_double_onesys_save_runs_{str(days_run[-1])}.csv', index=False)


path_save_csv = os.path.join(general_path, 'injections/cdouble')
if os.path.exists(os.path.join(path_save_csv, 'c_injections_double_onesys_save_runs.csv')):
    import random
    random_seed = random.getrandbits(16)
    df_to_use_to_plot_single_double_save.to_csv(os.path.join(path_save_csv, f'c_injections_double_onesys_save_runs_{random_seed}.csv'), index=False)
else:
    df_to_use_to_plot_single_double_save.to_csv(os.path.join(path_save_csv, 'c_injections_double_onesys_save_runs.csv'), index=False)

    