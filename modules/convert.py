import numpy as np
import pars
import bary_time

'''fit_pars in the module of pars is in the sequence of acceleration, period correction, phase offset, correction of RA, correction of DEC. The units are [bins/hr/hr, bins/hr, bins, radians, radians]'''

'''converts units of bins and hr to second'''
bins = pars.T / pars.NPHASEBIN
hr = 3600 

'''converts radian to degree and arcsec'''
rad_to_deg = 360/(2*np.pi)
rad_to_arcmin = rad_to_deg * 60
rad_to_arcsec = rad_to_arcmin *60

'''Converts accelerations, period, phase offset to sec'''
accel = pars.fit_pars[0] * bins / hr / hr
accel_err = pars.fit_pars_err[0] * bins / hr / hr

period_cor = pars.fit_pars[1] * bins / hr
period_cor_err = pars.fit_pars_err[1] * bins / hr

period_new = pars.T + period_cor
period_new_err =  period_cor_err

phase_offset = pars.fit_pars[2] * bins

print 'Acceleration: ', accel, '\+-', accel_err, '(1/sec)'
print 'Period correction', period_cor, '\+-', period_cor_err, '(sec)' 
print 'New Period', period_new, '\+-', period_new_err, '(sec)'
#print 'phase offset', phase_offset, '(sec)'


'''Converts corrections of RA and DEC to arcsec'''
RA = pars.RA
DEC = pars.DEC
RA_HMS = bary_time.deg_to_HMS(RA)
DEC_DMS = bary_time.deg_to_DMS(DEC)

RA_cor = pars.fit_pars_err[3] * rad_to_arcsec
RA_cor_err = pars.fit_pars_err[3] * rad_to_arcsec

DEC_cor = pars.fit_pars_err[4] * rad_to_arcsec
DEC_cor_err = pars.fit_pars_err[4] * rad_to_arcsec

print "New RA: ", RA_HMS
print "New DEC: ", DEC_DMS
print 'Delta RA: ', RA_cor, '\+-', RA_cor_err, '(arcsec)'
print 'Delta DEC: ', DEC_cor, '\+-', DEC_cor_err, '(arcsec)' 




