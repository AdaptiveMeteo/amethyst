import time
from os.path import join

from dobjects.SounderFOV        import SounderFOV, sounderfov
from dobjects.FirstGuess        import FirstGuess, firstguess
from dobjects.AprioriCovariance import AprioriCovariance , covariance_matrix
from dobjects.Emissivity        import Emissivity, emissivity_model
from mirto.mirto_code_transform import transform


def reader_f(to_compute, obserr, obsnum, L, eigenland, eigensea, co2std, config_vars, var_selection_flags):

    # Read data from disc
    start_read_time = time.time()

    fov = SounderFOV( config_vars ) 
    fov_time = time.time()
    L.log('\n   SounderFov opened in ' +str(fov_time - start_read_time)+' seconds', 4)

    fg = FirstGuess(    config_vars["fg_file"], 
                        eigenland, eigensea, co2std,
                        var_selection=var_selection_flags)
    
    fg_time = time.time()
    L.log('   FirstGuess opened in ' +str(fg_time - fov_time)+' seconds', 4)

    aemiss = Emissivity( config_vars["fg_file"], eigenland, eigensea)
    aemiss_time = time.time()
    L.log('   AemissFile opened in ' +str(aemiss_time - fg_time)+' seconds', 4)

    apriori = AprioriCovariance( config_vars["apriori_file"], aemiss,
                                   eigenland, eigensea,
                                   var_selection=var_selection_flags)
    apriori_time = time.time()
    L.log('   AprioriCovariance opened in ' +str(apriori_time - aemiss_time)+' seconds', 4)

    # Create the initial values
    nlev = fg.xdim(0)[0]
    nselstate = len(fg.LEVEL_SELECTION)
    selchannels = fov.selchannels
    in_val_time = time.time()
    L.log('   Initial values read in ' +str(in_val_time - aemiss_time) +' seconds', 4)

    transformer = transform(obserr, nlev, config_vars)
    transformer_time = time.time()
    L.log('   Transformer object created in ' +str(transformer_time - in_val_time)+' seconds', 4)
    
    # Send initial values to the master process
    to_compute.put((nlev, nselstate, selchannels, transformer))
    while not to_compute.empty():
        time.sleep(0.1)
    # Start to load the data
    for i in obsnum:
        fov_fast = sounderfov(fov, i)
        fg_fast = firstguess(fg, i)
        apriori_fast = covariance_matrix(apriori, i)
        aemiss_fast = emissivity_model(aemiss, i)
        to_compute.put((i, (fov_fast, fg_fast, apriori_fast, aemiss_fast)))
    time.sleep(1.0)
    while not to_compute.empty():
        time.sleep(1.0)
