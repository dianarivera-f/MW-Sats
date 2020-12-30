import numpy as np
import pandas as pd

a=pd.read_csv('coor.csv')
a['center']=list(zip(a['center_ra'], a['center_dec'].astype(float)))
b=pd.read_csv('dist.csv')
c=pd.read_csv('err_ra.csv')
c['center_ra_error']=list(zip(c['plus'], c['minus'].astype(float)))
e=pd.read_csv('err_dec.csv')
e['center_dec_error']=list(zip(e['plus'], e['minus'].astype(float)))

gal=['draco','carina','umi','leo1','leo2','sextans','fornax','sculptor'] #classical satellite galaxies
gal_names=pd.DataFrame({'general names':gal,'gaia names':[sub + '_gaia' for sub in gal],
                        'ap names':[sub + '_ap' for sub in gal],'js names':[sub + '_js' for sub in gal]})

num=[393,717,309,464,192,385,5603,2360] #number of stars in the galaxy data
mass=[0.022,0.00151,'NAN',0.02,0.027,0.3,0.054,'NAN'] #billion solar masses


#[mass(billion solar masses),distance(kpc),]
lmc_params=[10,49.97,]
smc_params=[7,61.7,]

data=pd.DataFrame({'galaxies':gal,'center':a.center,'center_ra_error':c.center_ra_error,'center_dec_error':e.center_dec_error,
                   'distance':b.dis,'diameter':b.dia,'avm':b.avm,'mass':mass,'num':num,'type':b.type})

js=pd.read_csv("js_table.csv")
js['Mv_error']=list(zip(js['dM_V+'], js['dM_V-'].astype(float)))
js['r1_2_error']=list(zip(js['dr_1/2+'], js['dr_1/2-'].astype(float)))
js['dist_error']=list(zip(js['ddist+'], js['ddist-'].astype(float)))
js['Vhel_error']=list(zip(js['dVhel+'], js['dVhel-'].astype(float)))
js['sigma_error']=list(zip(js['dsig+'], js['dsig-'].astype(float)))
js['Fe_H_error']=list(zip(js['d[Fe/H]+'], js['d[Fe/H]-'].astype(float)))
js['sigFe_H_error']=list(zip(js['dsigFe/H+'], js['dsigFe/H-'].astype(float)))
lc=js['dwarf'].str.lower().tolist()
dw=[sub + '_js' for sub in lc]
js['dwarf']=dw
imp=js.drop(columns=['dM_V+', 'dM_V-','dr_1/2+', 'dr_1/2-','ddist+','ddist-','dVhel+','dVhel-','dsig+',
                     'dsig-','d[Fe/H]+','d[Fe/H]-','dsigFe/H+','dsigFe/H-'])
columns=imp.iloc[:,np.r_[1:15]].keys().tolist()

import sat_array
from importlib import reload
reload(sat_array)

d = np.zeros(1, dtype=sat_array.mydtype).view(np.recarray)

d['gals']=gal_names

#smc and lmc data
d['lmc']['mass']=lmc_params[0]
d['lmc']['distance']=lmc_params[1]


d['smc']['mass']=smc_params[0]
d['smc']['distance']=smc_params[1]



#Josh Simon data
for i in range(len(dw)):
    d[dw[i]][columns[0]]=imp.Mv[i]
    d[dw[i]][columns[1]]=imp.Mv_error[i]
    d[dw[i]][columns[2]]=imp.r1_2[i]
    d[dw[i]][columns[3]]=imp.r1_2_error[i]
    d[dw[i]][columns[4]]=imp.dist[i]
    d[dw[i]][columns[5]]=imp.dist_error[i]
    d[dw[i]][columns[6]]=imp.Vhel[i]
    d[dw[i]][columns[7]]=imp.Vhel_error[i]
    d[dw[i]][columns[8]]=imp.sigma[i]
    d[dw[i]][columns[9]]=imp.sigma_error[i]
    d[dw[i]][columns[10]]=imp.Fe_H[i]
    d[dw[i]][columns[11]]=imp.Fe_H_error[i]
    d[dw[i]][columns[12]]=imp.sigFe_H[i]
    d[dw[i]][columns[13]]=imp.sigFe_H_error[i]
    
#general data
norm_names=gal_names['general names'].tolist()
norm_cols=data.keys().tolist()

for i in range(len(gal)):
#(center ra, center dec)
    d[norm_names[i]][norm_cols[1]]=data.center[i]
#center ra error, center dec error
    d[norm_names[i]][norm_cols[2]]=data.center_ra_error[i]
    d[norm_names[i]][norm_cols[3]]=data.center_dec_error[i]
#distance, diameter, AVM
    d[norm_names[i]][norm_cols[4]]=data.distance[i]
    d[norm_names[i]][norm_cols[5]]=data.diameter[i]
    d[norm_names[i]][norm_cols[6]]=data.avm[i]
#mass
    d[norm_names[i]][norm_cols[7]]=data.mass[i]
#num
    d[norm_names[i]][norm_cols[8]]=data.num[i]
#type
    d[norm_names[i]][norm_cols[9]]=data.type[i]




draco = pd.read_csv('draco_filtered.csv')
draco_ap = pd.read_csv('draco_ap.csv')
d['draco_gaia']['ra']=draco.ra
d['draco_gaia']['ra_error']=draco.ra_error
d['draco_gaia']['dec']=draco.dec
d['draco_gaia']['dec_error']=draco.dec_error
d['draco_gaia']['pmra']=draco.pmra
d['draco_gaia']['pmra_error']=draco.pmra_error
d['draco_gaia']['pmdec']=draco.pmdec
d['draco_gaia']['pmdec_error']=draco.pmdec_error
d['draco_gaia']['par']=draco.parallax
d['draco_gaia']['par_error']=draco.parallax_error
d['draco_gaia']['par_over_error']=draco.parallax_over_error
d['draco_gaia']['ra_dec_corr']=draco.ra_dec_corr
d['draco_gaia']['gflux']=draco.phot_g_mean_flux
d['draco_gaia']['gmag']=draco.phot_g_mean_mag
d['draco_gaia']['bpflux']=draco.phot_bp_mean_flux
d['draco_gaia']['bpmag']=draco.phot_bp_mean_mag
d['draco_gaia']['rpflux']=draco.phot_rp_mean_flux
d['draco_gaia']['rpmag']=draco.phot_rp_mean_mag
d['draco_gaia']['bp_rp']=draco.bp_rp
d['draco_gaia']['bp_g']=draco.bp_g
d['draco_gaia']['g_rp']=draco.g_rp
d['draco_ap']['ra']=draco_ap.ra
d['draco_ap']['dec']=draco_ap.dec
d['draco_ap']['Vhel']=draco_ap.Vhel
d['draco_ap']['Vhel_error']=draco_ap.Vhel_error
d['draco_ap']['mem']=draco_ap.mem


carina=pd.read_csv('carina_filtered.csv')
carina_ap=pd.read_csv('carina_ap.csv')
d['carina_gaia']['ra']=carina.ra
d['carina_gaia']['ra_error']=carina.ra_error
d['carina_gaia']['dec']=carina.dec
d['carina_gaia']['dec_error']=carina.dec_error
d['carina_gaia']['pmra']=carina.pmra
d['carina_gaia']['pmra_error']=carina.pmra_error
d['carina_gaia']['pmdec']=carina.pmdec
d['carina_gaia']['pmdec_error']=carina.pmdec_error
d['carina_gaia']['par']=carina.parallax
d['carina_gaia']['par_error']=carina.parallax_error
d['carina_gaia']['par_over_error']=carina.parallax_over_error
d['carina_gaia']['ra_dec_corr']=carina.ra_dec_corr
d['carina_gaia']['gflux']=carina.phot_g_mean_flux
d['carina_gaia']['gmag']=carina.phot_g_mean_mag
d['carina_gaia']['bpflux']=carina.phot_bp_mean_flux
d['carina_gaia']['bpmag']=carina.phot_bp_mean_mag
d['carina_gaia']['rpflux']=carina.phot_rp_mean_flux
d['carina_gaia']['rpmag']=carina.phot_rp_mean_mag
d['carina_gaia']['bp_rp']=carina.bp_rp
d['carina_gaia']['bp_g']=carina.bp_g
d['carina_gaia']['g_rp']=carina.g_rp
d['carina_ap']['ra']=carina_ap.ra
d['carina_ap']['dec']=carina_ap.dec
d['carina_ap']['Vhel']=carina_ap.Vhel
d['carina_ap']['Vhel_error']=carina_ap.Vhel_error
d['carina_ap']['mem']=carina_ap.mem


umi=pd.read_csv('ursmin_filtered.csv')
umi_ap=pd.read_csv('ursaminor_ap.csv')
d['umi_gaia']['ra']=umi.ra
d['umi_gaia']['ra_error']=umi.ra_error
d['umi_gaia']['dec']=umi.dec
d['umi_gaia']['dec_error']=umi.dec_error
d['umi_gaia']['pmra']=umi.pmra
d['umi_gaia']['pmra_error']=umi.pmra_error
d['umi_gaia']['pmdec']=umi.pmdec
d['umi_gaia']['pmdec_error']=umi.pmdec_error
d['umi_gaia']['par']=umi.parallax
d['umi_gaia']['par_error']=umi.parallax_error
d['umi_gaia']['par_over_error']=umi.parallax_over_error
d['umi_gaia']['ra_dec_corr']=umi.ra_dec_corr
d['umi_gaia']['gflux']=umi.phot_g_mean_flux
d['umi_gaia']['gmag']=umi.phot_g_mean_mag
d['umi_gaia']['bpflux']=umi.phot_bp_mean_flux
d['umi_gaia']['bpmag']=umi.phot_bp_mean_mag
d['umi_gaia']['rpflux']=umi.phot_rp_mean_flux
d['umi_gaia']['rpmag']=umi.phot_rp_mean_mag
d['umi_gaia']['bp_rp']=umi.bp_rp
d['umi_gaia']['bp_g']=umi.bp_g
d['umi_gaia']['g_rp']=umi.g_rp
d['umi_ap']['ra']=umi_ap.ra
d['umi_ap']['dec']=umi_ap.dec
d['umi_ap']['Vhel']=umi_ap.Vhel
d['umi_ap']['Vhel_error']=umi_ap.Vhel_error
d['umi_ap']['mem']=umi_ap.mem


leo1=pd.read_csv('leo1_filtered.csv')
leo1_ap=pd.read_csv('leo1_ap.csv')
d['leo1_gaia']['ra']=leo1.ra
d['leo1_gaia']['ra_error']=leo1.ra_error
d['leo1_gaia']['dec']=leo1.dec
d['leo1_gaia']['dec_error']=leo1.dec_error
d['leo1_gaia']['pmra']=leo1.pmra
d['leo1_gaia']['pmra_error']=leo1.pmra_error
d['leo1_gaia']['pmdec']=leo1.pmdec
d['leo1_gaia']['pmdec_error']=leo1.pmdec_error
d['leo1_gaia']['par']=leo1.parallax
d['leo1_gaia']['par_error']=leo1.parallax_error
d['leo1_gaia']['par_over_error']=leo1.parallax_over_error
d['leo1_gaia']['ra_dec_corr']=leo1.ra_dec_corr
d['leo1_gaia']['gflux']=leo1.phot_g_mean_flux
d['leo1_gaia']['gmag']=leo1.phot_g_mean_mag
d['leo1_gaia']['bpflux']=leo1.phot_bp_mean_flux
d['leo1_gaia']['bpmag']=leo1.phot_bp_mean_mag
d['leo1_gaia']['rpflux']=leo1.phot_rp_mean_flux
d['leo1_gaia']['rpmag']=leo1.phot_rp_mean_mag
d['leo1_gaia']['bp_rp']=leo1.bp_rp
d['leo1_gaia']['bp_g']=leo1.bp_g
d['leo1_gaia']['g_rp']=leo1.g_rp
d['leo1_ap']['ra']=leo1_ap.ra
d['leo1_ap']['dec']=leo1_ap.dec
d['leo1_ap']['Vhel']=leo1_ap.Vhel
d['leo1_ap']['Vhel_error']=leo1_ap.Vhel_error
d['leo1_ap']['mem']=leo1_ap.mem


leo2=pd.read_csv('leo2_filtered.csv')
leo2_ap=pd.read_csv('leo2_ap.csv')
d['leo2_gaia']['ra']=leo2.ra
d['leo2_gaia']['ra_error']=leo2.ra_error
d['leo2_gaia']['dec']=leo2.dec
d['leo2_gaia']['dec_error']=leo2.dec_error
d['leo2_gaia']['pmra']=leo2.pmra
d['leo2_gaia']['pmra_error']=leo2.pmra_error
d['leo2_gaia']['pmdec']=leo2.pmdec
d['leo2_gaia']['pmdec_error']=leo2.pmdec_error
d['leo2_gaia']['par']=leo2.parallax
d['leo2_gaia']['par_error']=leo2.parallax_error
d['leo2_gaia']['par_over_error']=leo2.parallax_over_error
d['leo2_gaia']['ra_dec_corr']=leo2.ra_dec_corr
d['leo2_gaia']['gflux']=leo2.phot_g_mean_flux
d['leo2_gaia']['gmag']=leo2.phot_g_mean_mag
d['leo2_gaia']['bpflux']=leo2.phot_bp_mean_flux
d['leo2_gaia']['bpmag']=leo2.phot_bp_mean_mag
d['leo2_gaia']['rpflux']=leo2.phot_rp_mean_flux
d['leo2_gaia']['rpmag']=leo2.phot_rp_mean_mag
d['leo2_gaia']['bp_rp']=leo2.bp_rp
d['leo2_gaia']['bp_g']=leo2.bp_g
d['leo2_gaia']['g_rp']=leo2.g_rp
d['leo2_ap']['ra']=leo2_ap.ra
d['leo2_ap']['dec']=leo2_ap.dec
d['leo2_ap']['Vhel']=leo2_ap.Vhel
d['leo2_ap']['Vhel_error']=leo2_ap.Vhel_error
d['leo2_ap']['Fe_H']=leo2_ap.Fe_H
d['leo2_ap']['Fe_H_error']=leo2_ap.Fe_H_error
d['leo2_ap']['mem']=leo2_ap.mem


sextans=pd.read_csv('sextans_filtered.csv')
sextans_ap=pd.read_csv('sextans_ap.csv')
d['sextans_gaia']['ra']=sextans.ra
d['sextans_gaia']['ra_error']=sextans.ra_error
d['sextans_gaia']['dec']=sextans.dec
d['sextans_gaia']['dec_error']=sextans.dec_error
d['sextans_gaia']['pmra']=sextans.pmra
d['sextans_gaia']['pmra_error']=sextans.pmra_error
d['sextans_gaia']['pmdec']=sextans.pmdec
d['sextans_gaia']['pmdec_error']=sextans.pmdec_error
d['sextans_gaia']['par']=sextans.parallax
d['sextans_gaia']['par_error']=sextans.parallax_error
d['sextans_gaia']['par_over_error']=sextans.parallax_over_error
d['sextans_gaia']['ra_dec_corr']=sextans.ra_dec_corr
d['sextans_gaia']['gflux']=sextans.phot_g_mean_flux
d['sextans_gaia']['gmag']=sextans.phot_g_mean_mag
d['sextans_gaia']['bpflux']=sextans.phot_bp_mean_flux
d['sextans_gaia']['bpmag']=sextans.phot_bp_mean_mag
d['sextans_gaia']['rpflux']=sextans.phot_rp_mean_flux
d['sextans_gaia']['rpmag']=sextans.phot_rp_mean_mag
d['sextans_gaia']['bp_rp']=sextans.bp_rp
d['sextans_gaia']['bp_g']=sextans.bp_g
d['sextans_gaia']['g_rp']=sextans.g_rp
d['sextans_ap']['ra']=sextans_ap.ra
d['sextans_ap']['dec']=sextans_ap.dec
d['sextans_ap']['Vhel']=sextans_ap.Vhel
d['sextans_ap']['Vhel_error']=sextans_ap.Vhel_error
d['sextans_ap']['mem']=sextans_ap.mem


fornax=pd.read_csv('fornax_filtered.csv')
fornax_ap=pd.read_csv('fornax_ap.csv')
d['fornax_gaia']['ra']=fornax.ra
d['fornax_gaia']['ra_error']=fornax.ra_error
d['fornax_gaia']['dec']=fornax.dec
d['fornax_gaia']['dec_error']=fornax.dec_error
d['fornax_gaia']['pmra']=fornax.pmra
d['fornax_gaia']['pmra_error']=fornax.pmra_error
d['fornax_gaia']['pmdec']=fornax.pmdec
d['fornax_gaia']['pmdec_error']=fornax.pmdec_error
d['fornax_gaia']['par']=fornax.parallax
d['fornax_gaia']['par_error']=fornax.parallax_error
d['fornax_gaia']['par_over_error']=fornax.parallax_over_error
d['fornax_gaia']['ra_dec_corr']=fornax.ra_dec_corr
d['fornax_gaia']['gflux']=fornax.phot_g_mean_flux
d['fornax_gaia']['gmag']=fornax.phot_g_mean_mag
d['fornax_gaia']['bpflux']=fornax.phot_bp_mean_flux
d['fornax_gaia']['bpmag']=fornax.phot_bp_mean_mag
d['fornax_gaia']['rpflux']=fornax.phot_rp_mean_flux
d['fornax_gaia']['rpmag']=fornax.phot_rp_mean_mag
d['fornax_gaia']['bp_rp']=fornax.bp_rp
d['fornax_gaia']['bp_g']=fornax.bp_g
d['fornax_gaia']['g_rp']=fornax.g_rp
d['fornax_ap']['ra']=fornax_ap.ra
d['fornax_ap']['dec']=fornax_ap.dec
d['fornax_ap']['Vhel']=fornax_ap.Vhel
d['fornax_ap']['Vhel_error']=fornax_ap.Vhel_error
d['fornax_ap']['mem']=fornax_ap.mem


sculptor=pd.read_csv('sculptor_filtered.csv')
sculptor_ap=pd.read_csv('sculptor_ap.csv')
d['sculptor_gaia']['ra']=sculptor.ra
d['sculptor_gaia']['ra_error']=sculptor.ra_error
d['sculptor_gaia']['dec']=sculptor.dec
d['sculptor_gaia']['dec_error']=sculptor.dec_error
d['sculptor_gaia']['pmra']=sculptor.pmra
d['sculptor_gaia']['pmra_error']=sculptor.pmra_error
d['sculptor_gaia']['pmdec']=sculptor.pmdec
d['sculptor_gaia']['pmdec_error']=sculptor.pmdec_error
d['sculptor_gaia']['par']=sculptor.parallax
d['sculptor_gaia']['par_error']=sculptor.parallax_error
d['sculptor_gaia']['par_over_error']=sculptor.parallax_over_error
d['sculptor_gaia']['ra_dec_corr']=sculptor.ra_dec_corr
d['sculptor_gaia']['gflux']=sculptor.phot_g_mean_flux
d['sculptor_gaia']['gmag']=sculptor.phot_g_mean_mag
d['sculptor_gaia']['bpflux']=sculptor.phot_bp_mean_flux
d['sculptor_gaia']['bpmag']=sculptor.phot_bp_mean_mag
d['sculptor_gaia']['rpflux']=sculptor.phot_rp_mean_flux
d['sculptor_gaia']['rpmag']=sculptor.phot_rp_mean_mag
d['sculptor_gaia']['bp_rp']=sculptor.bp_rp
d['sculptor_gaia']['bp_g']=sculptor.bp_g
d['sculptor_gaia']['g_rp']=sculptor.g_rp
d['sculptor_ap']['ra']=sculptor_ap.ra
d['sculptor_ap']['dec']=sculptor_ap.dec
d['sculptor_ap']['Vhel']=sculptor_ap.Vhel
d['sculptor_ap']['Vhel_error']=sculptor_ap.Vhel_error
d['sculptor_ap']['mem']=sculptor_ap.mem