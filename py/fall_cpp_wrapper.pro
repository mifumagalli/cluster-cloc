pro fall_cpp_wrapper, t0=t0,t1=t1, mmin=mmin, mmax=mmax,$
	cmfslope=cmfslope, gamma_min=gamma_min, gamma_max=gamma_max,$
	eta=eta, age_slope=age_slope, start=start, stop_=stop_,$
	step=step, sfr_err=sfr_err, sfr_start=sfr_star, sfr_stop=sfr_stop,$
	length=length, f_c=f_c, grid_out=grid_out, obs_err=obs_err

dir = '~/projects/sfr_l1/calc/'

if ~keyword_set(t0) then t0 = 0.0
if ~keyword_set(t1) then t1 = 0.1
if ~keyword_set(mmin) then mmin = 100
if ~keyword_set(mmax) then mmax = 1e9
if ~keyword_set(cmfslope) then cmfslope = -2.0
if ~keyword_set(gamma_min) then gamma_min = 5.05e18
if ~keyword_set(gamma_max) then gamma_max = 1.198e20
if ~keyword_set(eta) then eta = -0.688
if ~keyword_set(age_slope) then age_slope = -0.9
if ~keyword_set(start) then start = 15.
if ~keyword_set(stop_) then stop_ = 35.
if ~keyword_set(step) then step = 0.0125
if ~keyword_set(length) then length = 1e9
if ~keyword_set(f_c) then f_c = 0.1
if ~keyword_set(obs_err) then obs_err = alog(sqrt(2.512))
if ~keyword_set(sfr_start) then sfr_start = -4
if ~keyword_set(sfr_stop) then sfr_stop = 4
if ~keyword_set(sfr_step) then sfr_step = 0.0125
if ~keyword_set(sfr_err) then sfr_err = 0.5
if ~keyword_set(grid_out) then grid_out='grid.dat'
s=[t0, t1, mmin, mmax, cmfslope, gamma_min, gamma_max, eta, age_slope,$
	start, stop_, step, length, f_c, obs_err, sfr_start, sfr_stop,$
	sfr_step,sfr_err]
s=rstring(s)
s=strjoin(s, ' ')

;spawn, dir+'fall '+s,out
command = dir+'fall '+s+' '+grid_out
print, command
spawn, command


end
