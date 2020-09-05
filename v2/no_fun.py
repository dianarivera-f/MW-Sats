import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')

from astroquery.gaia import Gaia
tables = Gaia.load_tables(only_names=True)

gal=['Draco','Carina','Ursa Minor','Leo I','Leo II','Sextans','Fornax','Sculptor'] #classical satellite galaxies
num=[393,717,309,464,192,385,5603,2360] #number of stars in the galaxy data
mass=[0.022,0.00151,'NAN',0.02,0.027,0.3,0.054,'NAN'] #billion solar masses
a=pd.read_csv('coor.csv')
b=pd.read_csv('dist.csv')
c=pd.read_csv('err_ra.csv')
e=pd.read_csv('err_dec.csv')
f=pd.read_csv('data.csv')
g=pd.read_csv('data2.csv')

import sat_array
from importlib import reload
reload(sat_array)

d = np.zeros(1, dtype=sat_array.mydtype).view(np.recarray)
d['gal']=gal

draco = pd.read_csv('draco_filtered.csv')
d['draco']['ra']=draco.ra
d['draco']['ra_error']=draco.ra_error
d['draco']['dec']=draco.dec
d['draco']['dec_error']=draco.dec_error
d['draco']['pmra']=draco.pmra
d['draco']['pmra_error']=draco.pmra_error
d['draco']['pmdec']=draco.pmdec
d['draco']['pmdec_error']=draco.pmdec_error
d['draco']['par']=draco.parallax
d['draco']['par_error']=draco.parallax_error
d['draco']['par_over_error']=draco.parallax_over_error
d['draco']['ra_dec_corr']=draco.ra_dec_corr
d['draco']['gflux']=draco.phot_g_mean_flux
d['draco']['gmag']=draco.phot_g_mean_mag
d['draco']['bpflux']=draco.phot_bp_mean_flux
d['draco']['bpmag']=draco.phot_bp_mean_mag
d['draco']['rpflux']=draco.phot_rp_mean_flux
d['draco']['rpmag']=draco.phot_rp_mean_mag
d['draco']['bp_rp']=draco.bp_rp
d['draco']['bp_g']=draco.bp_g
d['draco']['g_rp']=draco.g_rp
d['draco']['center_ra']=a.ra[0]
d['draco']['center_dec']=a.dec[0]
d['draco']['distance']=b.dis[0]
d['draco']['diameter']=b.dia[0]
d['draco']['avm']=b.avm[0]
d['draco']['mass']=mass[0]
d['draco']['num']=num[0]

carina = pd.read_csv('carina_filtered.csv')
d['carina']['ra']=carina.ra
d['carina']['ra_error']=carina.ra_error
d['carina']['dec']=carina.dec
d['carina']['dec_error']=carina.dec_error
d['carina']['pmra']=carina.pmra
d['carina']['pmra_error']=carina.pmra_error
d['carina']['pmdec']=carina.pmdec
d['carina']['pmdec_error']=carina.pmdec_error
d['carina']['par']=carina.parallax
d['carina']['par_error']=carina.parallax_error
d['carina']['par_over_error']=carina.parallax_over_error
d['carina']['ra_dec_corr']=carina.ra_dec_corr
d['carina']['gflux']=carina.phot_g_mean_flux
d['carina']['gmag']=carina.phot_g_mean_mag
d['carina']['bpflux']=carina.phot_bp_mean_flux
d['carina']['bpmag']=carina.phot_bp_mean_mag
d['carina']['rpflux']=carina.phot_rp_mean_flux
d['carina']['rpmag']=carina.phot_rp_mean_mag
d['carina']['bp_rp']=carina.bp_rp
d['carina']['bp_g']=carina.bp_g
d['carina']['g_rp']=carina.g_rp
d['carina']['center_ra']=a.ra[1]
d['carina']['center_dec']=a.dec[1]
d['carina']['distance']=b.dis[1]
d['carina']['diameter']=b.dia[1]
d['carina']['avm']=b.avm[1]
d['carina']['mass']=mass[1]
d['carina']['num']=num[1]

umi = pd.read_csv('ursmin_filtered.csv')
d['umi']['ra']=umi.ra
d['umi']['ra_error']=umi.ra_error
d['umi']['dec']=umi.dec
d['umi']['dec_error']=umi.dec_error
d['umi']['pmra']=umi.pmra
d['umi']['pmra_error']=umi.pmra_error
d['umi']['pmdec']=umi.pmdec
d['umi']['pmdec_error']=umi.pmdec_error
d['umi']['par']=umi.parallax
d['umi']['par_error']=umi.parallax_error
d['umi']['par_over_error']=umi.parallax_over_error
d['umi']['ra_dec_corr']=umi.ra_dec_corr
d['umi']['gflux']=umi.phot_g_mean_flux
d['umi']['gmag']=umi.phot_g_mean_mag
d['umi']['bpflux']=umi.phot_bp_mean_flux
d['umi']['bpmag']=umi.phot_bp_mean_mag
d['umi']['rpflux']=umi.phot_rp_mean_flux
d['umi']['rpmag']=umi.phot_rp_mean_mag
d['umi']['bp_rp']=umi.bp_rp
d['umi']['bp_g']=umi.bp_g
d['umi']['g_rp']=umi.g_rp
d['umi']['center_ra']=a.ra[2]
d['umi']['center_dec']=a.dec[2]
d['umi']['distance']=b.dis[2]
d['umi']['diameter']=b.dia[2]
d['umi']['avm']=b.avm[2]
d['umi']['mass']=mass[2]
d['umi']['num']=num[2]

leo1 = pd.read_csv('leo1_filtered.csv')
d['leo1']['ra']=leo1.ra
d['leo1']['ra_error']=leo1.ra_error
d['leo1']['dec']=leo1.dec
d['leo1']['dec_error']=leo1.dec_error
d['leo1']['pmra']=leo1.pmra
d['leo1']['pmra_error']=leo1.pmra_error
d['leo1']['pmdec']=leo1.pmdec
d['leo1']['pmdec_error']=leo1.pmdec_error
d['leo1']['par']=leo1.parallax
d['leo1']['par_error']=leo1.parallax_error
d['leo1']['par_over_error']=leo1.parallax_over_error
d['leo1']['ra_dec_corr']=leo1.ra_dec_corr
d['leo1']['gflux']=leo1.phot_g_mean_flux
d['leo1']['gmag']=leo1.phot_g_mean_mag
d['leo1']['bpflux']=leo1.phot_bp_mean_flux
d['leo1']['bpmag']=leo1.phot_bp_mean_mag
d['leo1']['rpflux']=leo1.phot_rp_mean_flux
d['leo1']['rpmag']=leo1.phot_rp_mean_mag
d['leo1']['bp_rp']=leo1.bp_rp
d['leo1']['bp_g']=leo1.bp_g
d['leo1']['g_rp']=leo1.g_rp
d['leo1']['center_ra']=a.ra[3]
d['leo1']['center_dec']=a.dec[3]
d['leo1']['distance']=b.dis[3]
d['leo1']['diameter']=b.dia[3]
d['leo1']['avm']=b.avm[3]
d['leo1']['mass']=mass[3]
d['leo1']['num']=num[3]

leo2 = pd.read_csv('leo2_filtered.csv')
d['leo2']['ra']=leo2.ra
d['leo2']['ra_error']=leo2.ra_error
d['leo2']['dec']=leo2.dec
d['leo2']['dec_error']=leo2.dec_error
d['leo2']['pmra']=leo2.pmra
d['leo2']['pmra_error']=leo2.pmra_error
d['leo2']['pmdec']=leo2.pmdec
d['leo2']['pmdec_error']=leo2.pmdec_error
d['leo2']['par']=leo2.parallax
d['leo2']['par_error']=leo2.parallax_error
d['leo2']['par_over_error']=leo2.parallax_over_error
d['leo2']['ra_dec_corr']=leo2.ra_dec_corr
d['leo2']['gflux']=leo2.phot_g_mean_flux
d['leo2']['gmag']=leo2.phot_g_mean_mag
d['leo2']['bpflux']=leo2.phot_bp_mean_flux
d['leo2']['bpmag']=leo2.phot_bp_mean_mag
d['leo2']['rpflux']=leo2.phot_rp_mean_flux
d['leo2']['rpmag']=leo2.phot_rp_mean_mag
d['leo2']['bp_rp']=leo2.bp_rp
d['leo2']['bp_g']=leo2.bp_g
d['leo2']['g_rp']=leo2.g_rp
d['leo2']['center_ra']=a.ra[4]
d['leo2']['center_dec']=a.dec[4]
d['leo2']['distance']=b.dis[4]
d['leo2']['diameter']=b.dia[4]
d['leo2']['avm']=b.avm[4]
d['leo2']['mass']=mass[4]
d['leo2']['num']=num[4]

sextans = pd.read_csv('sextans_filtered.csv')
d['sextans']['ra']=sextans.ra
d['sextans']['ra_error']=sextans.ra_error
d['sextans']['dec']=sextans.dec
d['sextans']['dec_error']=sextans.dec_error
d['sextans']['pmra']=sextans.pmra
d['sextans']['pmra_error']=sextans.pmra_error
d['sextans']['pmdec']=sextans.pmdec
d['sextans']['pmdec_error']=sextans.pmdec_error
d['sextans']['par']=sextans.parallax
d['sextans']['par_error']=sextans.parallax_error
d['sextans']['par_over_error']=sextans.parallax_over_error
d['sextans']['ra_dec_corr']=sextans.ra_dec_corr
d['sextans']['gflux']=sextans.phot_g_mean_flux
d['sextans']['gmag']=sextans.phot_g_mean_mag
d['sextans']['bpflux']=sextans.phot_bp_mean_flux
d['sextans']['bpmag']=sextans.phot_bp_mean_mag
d['sextans']['rpflux']=sextans.phot_rp_mean_flux
d['sextans']['rpmag']=sextans.phot_rp_mean_mag
d['sextans']['bp_rp']=sextans.bp_rp
d['sextans']['bp_g']=sextans.bp_g
d['sextans']['g_rp']=sextans.g_rp
d['sextans']['center_ra']=a.ra[5]
d['sextans']['center_dec']=a.dec[5]
d['sextans']['distance']=b.dis[5]
d['sextans']['diameter']=b.dia[5]
d['sextans']['avm']=b.avm[5]
d['sextans']['mass']=mass[5]
d['sextans']['num']=num[5]

fornax = pd.read_csv('fornax_filtered.csv')
d['fornax']['ra']=fornax.ra
d['fornax']['ra_error']=fornax.ra_error
d['fornax']['dec']=fornax.dec
d['fornax']['dec_error']=fornax.dec_error
d['fornax']['pmra']=fornax.pmra
d['fornax']['pmra_error']=fornax.pmra_error
d['fornax']['pmdec']=fornax.pmdec
d['fornax']['pmdec_error']=fornax.pmdec_error
d['fornax']['par']=fornax.parallax
d['fornax']['par_error']=fornax.parallax_error
d['fornax']['par_over_error']=fornax.parallax_over_error
d['fornax']['ra_dec_corr']=fornax.ra_dec_corr
d['fornax']['gflux']=fornax.phot_g_mean_flux
d['fornax']['gmag']=fornax.phot_g_mean_mag
d['fornax']['bpflux']=fornax.phot_bp_mean_flux
d['fornax']['bpmag']=fornax.phot_bp_mean_mag
d['fornax']['rpflux']=fornax.phot_rp_mean_flux
d['fornax']['rpmag']=fornax.phot_rp_mean_mag
d['fornax']['bp_rp']=fornax.bp_rp
d['fornax']['bp_g']=fornax.bp_g
d['fornax']['g_rp']=fornax.g_rp
d['fornax']['center_ra']=a.ra[6]
d['fornax']['center_dec']=a.dec[6]
d['fornax']['distance']=b.dis[6]
d['fornax']['diameter']=b.dia[6]
d['fornax']['avm']=b.avm[6]
d['fornax']['mass']=mass[6]
d['fornax']['num']=num[6]

sculptor = pd.read_csv('sculptor_filtered.csv')
d['sculptor']['ra']=sculptor.ra
d['sculptor']['ra_error']=sculptor.ra_error
d['sculptor']['dec']=sculptor.dec
d['sculptor']['dec_error']=sculptor.dec_error
d['sculptor']['pmra']=sculptor.pmra
d['sculptor']['pmra_error']=sculptor.pmra_error
d['sculptor']['pmdec']=sculptor.pmdec
d['sculptor']['pmdec_error']=sculptor.pmdec_error
d['sculptor']['par']=sculptor.parallax
d['sculptor']['par_error']=sculptor.parallax_error
d['sculptor']['par_over_error']=sculptor.parallax_over_error
d['sculptor']['ra_dec_corr']=sculptor.ra_dec_corr
d['sculptor']['gflux']=sculptor.phot_g_mean_flux
d['sculptor']['gmag']=sculptor.phot_g_mean_mag
d['sculptor']['bpflux']=sculptor.phot_bp_mean_flux
d['sculptor']['bpmag']=sculptor.phot_bp_mean_mag
d['sculptor']['rpflux']=sculptor.phot_rp_mean_flux
d['sculptor']['rpmag']=sculptor.phot_rp_mean_mag
d['sculptor']['bp_rp']=sculptor.bp_rp
d['sculptor']['bp_g']=sculptor.bp_g
d['sculptor']['g_rp']=sculptor.g_rp
d['sculptor']['center_ra']=a.ra[7]
d['sculptor']['center_dec']=a.dec[7]
d['sculptor']['distance']=b.dis[7]
d['sculptor']['diameter']=b.dia[7]
d['sculptor']['avm']=b.avm[7]
d['sculptor']['mass']=mass[7]
d['sculptor']['num']=num[7]