import numpy as np

filename = "my_peaks.txt"
wavelength = 1.5418
peaks = np.loadtxt(filename)
zero_error = 0.0
ratio_tolerance = 0.005

def check_cubic(d_spacings, cub_reflecs, cub_rats, tol=ratio_tolerance,
                suppress_repeats=True, peaks=None):
    ratiosqs = []
    for i0, d0 in enumerate(d_spacings):
        for i1, d1 in enumerate(d_spacings):
            if i0 >= i1:
                continue
            else:
                ratiosq = (d0 / d1)**2
                ratiosqs.append((i0, i1, ratiosq))
    rep_rats = [] #tracking repeated ratios
    for r0 in ratiosqs:
        for r1 in cub_rats:
            if np.abs(1 - r0[2] / r1[2]) < tol:
                if r1[2] in rep_rats:
                    continue
                if suppress_repeats:
                    rep_rats.append(r1[2])
                print "(%f/%f)^2 is %f; could be (%d%d%d)/(%d%d%d) (%f)" \
                        % (d_spacings[r0[0]], d_spacings[r0[1]], r0[2],
                           cub_reflecs[r1[0]][0], cub_reflecs[r1[0]][1],
                           cub_reflecs[r1[0]][2], cub_reflecs[r1[1]][0],
                           cub_reflecs[r1[1]][1], cub_reflecs[r1[1]][2],
                           r1[2])
                if type(peaks) != type(None):
                    print "The above is 2 theta values %f and %f" \
                            % (peaks[r0[0]], peaks[r0[1]])
    if type(peaks) == type(None):
        print "Use peaks=peaks to print out 2 theta values"
    if suppress_repeats:
        print "Use suppress_repeats=False to get other reflection ratios"

def cubic_combs(hkl):
    return [tuple(hkl), (hkl[0], hkl[2], hkl[1]), (hkl[1], hkl[0], hkl[2]),
            (hkl[1], hkl[2], hkl[0]), (hkl[2], hkl[0], hkl[1]),
            (hkl[2], hkl[1], hkl[0])]

def get_cubic_ratios(reflecs):
    ratios = []
    for i0, r0 in enumerate(reflecs):
        for i1, r1 in enumerate(reflecs):
            if i0 >= i1:
                continue
            else:
                ratio = np.sum(np.array(r0, dtype=float)**2) / \
                        np.sum(np.array(r1, dtype=float)**2)
                if ratio < 1:
                    ratio = 1. / ratio
                    ratios.append((i1, i0, ratio))
                else:
                    ratios.append((i0, i1, ratio))
    return ratios

all_reflecs = []
highest_index = 2
for h in range(highest_index + 1):
    for k in range(highest_index + 1):
        for l in range(highest_index + 1):
            if h == 0 and k == 0 and l == 0:
                continue
            else:
                all_reflecs.append((h, k, l))

fcc_reflecs = []
for r in all_reflecs:
    app = False
    if r[0] % 2 == 0 and r[1] % 2 == 0 and r[2] % 2 == 0:
        app = True
    elif r[0] % 2 == 1 and r[1] % 2 == 1 and r[2] % 2 == 1:
        app = True
    for cc in cubic_combs(r):
        if cc in fcc_reflecs:
            app = False
    if app:
        fcc_reflecs.append(r)

bcc_reflecs = []
for r in all_reflecs:
    app = False
    if (r[0] + r[1] + r[2]) % 2 == 0:
        app = True
    for cc in cubic_combs(r):
        if cc in bcc_reflecs:
            app = False
    if app:
        bcc_reflecs.append(r)

pc_reflecs = []
for r in all_reflecs:
    app = True
    for cc in cubic_combs(r):
        if cc in pc_reflecs:
            app = False
    if app:
        pc_reflecs.append(r)

fcc_ratios = get_cubic_ratios(fcc_reflecs)
bcc_ratios = get_cubic_ratios(bcc_reflecs)
pc_ratios = get_cubic_ratios(pc_reflecs)

def tet_combs(hkl):
    return [tuple(hkl), (hkl[1], hkl[0], hkl[2])]

bct_reflecs = []
for r in all_reflecs:
    app = False
    if (r[0] + r[1] + r[2]) % 2 == 0:
        app = True
    for tc in tet_combs(r):
        if tc in bct_reflecs:
            app = False
    if app:
        bct_reflecs.append(r)

pt_reflecs = []
for r in all_reflecs:
    app = True
    for tc in tet_combs(r):
        if tc in pt_reflecs:
            app = False
    if app:
        pt_reflecs.append(r)

d_spacings = wavelength / (2 * np.sin((peaks - zero_error) * np.pi / 360.))


#two_th_sharp = np.array([24.30878, 42.33712, 44.64149, 51.44537, 64.93037,
#                         87.22724])
#two_th_broad = np.array([37.75586, 38.6364, 43.32547, 43.71763, 58.59523,
#                         60.43703, 63.17221, 70.41556, 76.49307, 77.58088,
#                         80.55757])
#
#d_sharp = 1.5418 / (2 * np.sin(two_th_sharp * np.pi / 360))
#d_broad = 1.5418 / (2 * np.sin(two_th_broad * np.pi / 360))
#
##Trying tetragonal a=4.0658, c=8.5528
#bct_ds1 = np.zeros(len(bct_reflecs))
#for i, r in enumerate(bct_reflecs):
#    bct_ds1[i] = 1 / np.sqrt((r[0]**2 + r[1]**2) / 4.0658**2 + \
#                             r[2]**2 / 8.5528**2)
#
#seens = []
#for d in d_sharp:
#    not_indexed = True
#    for i, d1 in enumerate(bct_ds1):
#        if np.abs(d - d1) < 0.02:
#            print '%f observed = %f calc, %s' % (d, d1, str(bct_reflecs[i]))
#            seens.append(i)
#            not_indexed = False
#    if not_indexed:
#        print '%f unindexed' % d
#
#unseens = [d for i, d in enumerate(bct_ds1) if i not in seens]
#unseens.sort(reverse=True)
#for s in unseens:
#    print s

