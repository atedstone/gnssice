from scipy import io

plt.figure()
sites = np.arange(1,7)
ax = None
for site in sites:

	data = io.loadmat('lev%svel_2009_spr2013_UTC-2_20130528.mat' %site)


	ax = plt.subplot(6,1,site, sharex=ax)
	plt.plot(data['v_6h'][:,0],data['v_6h'][:,1], label='v_6h', color='#4292C6')
	plt.plot(data['v_24h'][:,0],data['v_24h'][:,1], drawstyle='steps-post', label='v_24h', color='#08519C')
	plt.annotate(xy=(0.05, 0.05), xycoords='axes points', s='lev%s' %site)
	plt.ylim(0,500)
	plt.yticks([0,250,500],[0,250,500])

plt.xlim(121,240)
plt.xlabel('Day of year 2009')


