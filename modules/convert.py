import numpy as np
import pars
import bary_time

'''fit_pars in the module of pars is in the sequence of acceleration, period correction, phase offset, correction of RA, correction of DEC. The units are [bins/hr/hr, bins/hr, bins, radians, radians]'''

'''converts units of bins and hr to second'''
bins = 1. / pars.NPHASEBIN
hr = 3600

'''converts radian to degree and arcsec'''
rad_to_deg = 360/(2*np.pi)
rad_to_arcmin = rad_to_deg * 60
rad_to_arcsec = rad_to_arcmin *60

'''Converts corrections of RA and DEC to arcmin'''
# The cor and cor_err of RA and DEC are in the unit of rad.
RA = pars.RA
RA_cor = pars.fit_pars[3] * rad_to_deg
RA_HMS = str(bary_time.deg_to_HMS(RA + RA_cor))
RA_cor_err = pars.fit_pars_err[3] * rad_to_deg
RA_cor_HMS = str(bary_time.deg_to_HMS(RA_cor_err))
print "New RA: ", RA_HMS, '\+-', RA_cor_HMS

DEC = pars.DEC
DEC_cor = pars.fit_pars[4] * rad_to_deg
DEC_DMS = str(bary_time.deg_to_DMS(DEC + DEC_cor))
DEC_cor_err = pars.fit_pars_err[4] * rad_to_deg
DEC_cor_DMS = str(bary_time.deg_to_DMS(DEC_cor_err))
print "New DEC: ", DEC_DMS, '(DMS)', '\+-', DEC_cor_DMS, '(DMS)'
#print 'Delta RA: ', RA_cor, '\+-', RA_cor_err, '(arcmin)'
#print 'Delta DEC: ', DEC_cor, '\+-', DEC_cor_err, '(arcmin)'

'''Converts accelerations and period to sec. Show epoch and phase offset in Barycentric MJD.'''
accel = pars.fit_pars[0] * (bins / hr / hr) * 2 * (pars.T)**2
accel_err = pars.fit_pars_err[0] * (bins / hr / hr) * (pars.T)**2

period_cor = pars.fit_pars[1] * (bins / hr) * pars.T**2
period_cor_err = pars.fit_pars_err[1] * (bins / hr) * pars.T**2

period_new = pars.T + period_cor
period_new_err =  period_cor_err

epoch = pars.TIME0

phase_offset = pars.fit_pars[2]/ pars.NPHASEBIN 
BAT_phase_offset = epoch + (pars.T * phase_offset)/86400

print 'Period derivative: ', accel, '\+-', accel_err, '(sec 1/sec)'
print 'Period correction: ', period_cor, '\+-', period_cor_err, '(sec)' 
print 'New Period: ', period_new, '\+-', period_new_err, '(sec)'
print 'Epoch: ', epoch, '(MJD)'
print 'Phase offset: ', BAT_phase_offset, '(MJD)'

'''Calculate proper motion'''
'''array in the sequence of [delta_ra, delta_dec, delta_ra_err, delta_dec_err]'''
'''convert unit from rad to arcsec'''
delta_pos_11 = pars.delta_pos_11
delta_pos_15 = pars.delta_pos_15
delta_pos_ra = (delta_pos_15[0] - delta_pos_11[0]) * rad_to_arcsec
print 'delta_pos_ra', delta_pos_ra
delta_pos_ra_err = np.sqrt(delta_pos_15[2]**2 + delta_pos_11[2]**2) * rad_to_arcsec
print 'delta_pos_ra_err', delta_pos_ra_err
delta_pos_dec = (delta_pos_15[1] - delta_pos_11[1]) * rad_to_arcsec
delta_pos_dec_err = np.sqrt(delta_pos_15[3]**2 + delta_pos_11[3]**2) * rad_to_arcsec
years = (pars.mean_mjd_15wp - pars.mean_mjd_11w)/365.25

'''convert unit from arcsec to mas'''
prop_ra = delta_pos_ra / years * 10**3
prop_ra_err = delta_pos_ra_err / years * 10**3
prop_dec = delta_pos_dec / years * 10**3
prop_dec_err = delta_pos_dec_err / years * 10**3
print 'proper motion of RA: ', prop_ra, '\+-', prop_ra_err, '(mas/yr)'
print 'proper motion of DEC: ', prop_dec, '\+-', prop_dec_err, '(mas/yr)'
