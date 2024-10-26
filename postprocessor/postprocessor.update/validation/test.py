import numpy as np
from scipy import interpolate

Source = np.asarray([[255,255,255], [246,246,248], [229,230,245], [201,216,236], [101,156,186], [33,157,74], [122,222,0], [245,225,1], [243,142,7], [244,80,37], [248,66,114], [212,61,212], [139,39,139], [71,20,71], [19,9,20], [43,43,43], [86,86,86], [126,126,126], [170,170,170], [213,213,213]], dtype=np.int)

x = np.arange(0, Source.shape[0])

fit = interpolate.interp1d(x, Source, axis=0)

Target = fit(np.linspace(0, Source.shape[0]-1, 256))

print('Target {}'.format(Target))
