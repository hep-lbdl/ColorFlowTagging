import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize

def phi_fix(phi):
	if phi > np.pi:
		return phi - 2*np.pi
	return phi
f = np.vectorize(phi_fix)


names = ['qx_qg'] #, 'h_qq', 'h_gg']


def plot_dR_between_jets(n):
	with h5.File('../' + n + '_rot_standard.h5') as s:
		with h5.File('../' + n + '_rot_charged.h5') as c:
			sjp = f(s['meta_variables/jet_phi'][:])
			cjp = f(c['meta_variables/jet_phi'][:])
			# cjp = cjp[:]
			sje = s['meta_variables/jet_eta']
			cje = c['meta_variables/jet_eta']	
			s_im = s['images']			
			c_im = c['images']
			print(len(cjp))
			print(len(sjp))
			idx_s = 0
			for idx, i_cjp in enumerate(cjp):
				if idx % 100 == 0:
					print idx
				ok = False
					
				while not ok:
					d = np.subtract(s_im[idx_s], c_im[idx])
					if idx_s < 1000:
						print "\n\n"
						print idx, idx_s
						print(np.sum(s_im[idx_s]))
						print(np.sum(c_im[idx]))
						print(np.min(d)) # np.hypot(sjp[idx_s] - cjp[idx_s], cje[idx_s] - sje[idx_s]))
					if np.all(d >= 0.0):
						ok = True
					else:
						idx_s += 1

def plot_meta_and_imgs(n):
	n_short = n
	with h5.File('../' + n + '_rot_standard.h5') as s:
		with h5.File('../' + n + '_rot_charged.h5') as c:
			print(s['images'].shape)
			print(c['images'].shape)

			for meta in c['meta_variables'].keys():
				print(meta)
				if not meta == 'jet_phi':
					plt.hist(c['meta_variables/' + meta], bins=50, density=True, histtype='step', color='r', label='Charged')
					plt.hist(s['meta_variables/' + meta], bins=50, density=True, histtype='step', color='b', label='Standard')

				else:
					plt.hist(f(c['meta_variables/jet_phi'][:]), bins=50, density=True, histtype='step', color='r', label='Charged')
					plt.hist(f(s['meta_variables/jet_phi'][:]), bins=50, density=True, histtype='step', color='b', label='Standard')				
				plt.title(n + '  ' + meta)
				plt.legend()
				plt.savefig('imgs/' + n_short + '_' + meta)
				plt.clf()
				plt.cla()

			ms = normalize(np.mean(s['images'][:], axis=0))
			mc = np.mean(c['images'][:], axis=0)
			mc = normalize(mc)
			difference = np.subtract(ms, mc)
			difference = np.log(difference)
			plt.imshow(difference, cmap='plasma')
			plt.title(n + ' - Standard minus Charged')
			plt.colorbar()
			plt.savefig('imgs/' + n_short + '_diff_img')
			plt.clf()
			plt.cla()

			plt.imshow(np.log(ms), cmap='plasma')
			plt.title(n + ' Standard')
			plt.colorbar()
			plt.savefig('imgs/' + n_short + '_standard_img')
			plt.clf()
			plt.cla()

			plt.imshow(np.log(mc), cmap='plasma')
			plt.title(n + ' Charged')
			plt.colorbar()
			plt.savefig('imgs/' + n_short + '_charged_img')
			plt.clf()
			plt.cla()

			meta = 'dR'			
			plt.hist(np.hypot(np.subtract(c['meta_variables/jet_eta'], s['meta_variables/jet_eta']), np.subtract(c['meta_variables/jet_phi'], s['meta_variables/jet_phi'])), bins=200, density=True, histtype='step', color='r', label='Difference', range=[0., 0.1])
			plt.title(n + '  ' + meta)
			plt.legend()
			plt.savefig('imgs/' + n_short + '_' + meta)
			plt.clf()
			plt.cla()



for n in names:
	plot_meta_and_imgs(n)
#	plot_dR_between_jets(n)


