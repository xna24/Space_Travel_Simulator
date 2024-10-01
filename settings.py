import os, sys
import numpy as np

def is_num(s: str) -> bool: 
    try:
        float(s)
        return True
    except ValueError:
        return False

# if any element of l1 is in l2 return True, otherwise False
def overlaps(l1, l2) -> bool:
    for x in l1:
        for y in l2:
            if y == x: return True
    return False
   
# if l1 is a sub-list of l2 return True, otherwise False
def sublist(l1, l2) -> bool:
    for x in l1:
        if x not in l2: return False
    return True

def settings_read(settings_file):
    settings = dict()

    if os.path.isfile(settings_file):
        fh = open(settings_file, 'r')
        lines = fh.readlines()
        fh.close()
        for line in lines:
            if line == '\n': continue
            if line[0] == '[': continue
            line = line.strip()
            line = line.split(';')[0] # after character ; everything is ignored
            line = line.split('=')
            line = [x.strip() for x in line]
            settings[line[0]] = line[1]
    
    setting_keys = list(settings.keys())
    
    #--------------------------------------------------------------
    # use_latex
    #--------------------------------------------------------------
    if 'LaTeX' in setting_keys:
        use_latex = (settings['LaTeX']=='yes')
    else:
        use_latex = False # by default not to use LaTeX for plot label
    
    #--------------------------------------------------------------
    # star_catalog
    #--------------------------------------------------------------
    if 'StarCatalog' in setting_keys:
        star_catalog = settings['StarCatalog']
    else:
        star_catalog = 'HIP'
    
    #--------------------------------------------------------------
    # pic_width, pic_height
    #--------------------------------------------------------------
    if ('Width' in setting_keys) and ('Height' in setting_keys):
        if is_num(settings['Width']) and is_num(settings['Height']):
            pic_width, pic_height = float(settings['Width']), float(settings['Height'])
        else:
            sys.exit('Width and/or height setting invalid in \'settings.ini\'')
    elif ('Width' not in setting_keys) and ('Height' not in setting_keys):
        pic_width, pic_height = 2405, 1377 # adjusted by hand so output pic has dimension 1920x1080
    else:
        sys.exit('Only one of width, height in \'settings.ini\'')
        
    # aspect_ratio = golden_ratio
    # aspect_ratio = 1920/1080
    aspect_ratio = pic_width/pic_height
    
    #--------------------------------------------------------------
    # HIP_des, des_name
    #--------------------------------------------------------------
    if ('HIP' in setting_keys) and ('Name' in setting_keys):
        if is_num(settings['HIP']):    
            HIP_des = int(settings['HIP'])
            des_name = settings['Name']
            if des_name == '': des_name = 'HIP ' + settings['HIP']
        else:
            sys.exit('Destination HIP number invalid in \'settings.ini\'')
    elif not sublist(['RAdeg', 'DEdeg', 'parallax'], setting_keys):
        # no HIP number, no RAdeg, DEdeg, parallax also, use demo value
        print('Destination info not found. Using: HIP 38549')
        HIP_des = 38594
        des_name = 'HIP 38594'
    else:
        if not ('Name' in setting_keys):
            sys.exit('Name option missing in \'settings.ini\'')
        else:
            HIP_des = 0
            des_name = settings['Name']

    #--------------------------------------------------------------
    # col_names
    #--------------------------------------------------------------
    # column names in the pandas DataFrame containing all the stars
    if star_catalog == 'HIP_TYC':
        col_names = ['RAdeg', 'DEdeg', 'Plx', 'Plx_valid', 'Vmag', 'B-V']
    elif star_catalog == 'HIP':
        col_names = ['lii', 'bii', 'parallax', 'vmag', 'bv_color']
    else:
        if sublist(['StarCatalogFileName','StarCatalogRAColumn', \
        'StarCatalogDEColumn', 'StarCatalogParallaxColumn',\
        'StarCatalogVmagColumn','StarCatalogBVColumn'], setting_keys):
            col_names = [settings['StarCatalogRAColumn'],\
                         settings['StarCatalogDEColumn'],\
                         settings['StarCatalogParallaxColumn'],\
                         settings['StarCatalogVmagColumn'],\
                         settings['StarCatalogBVColumn']]
        else:
            sys.exit('Custom star catalog column names missing or incomplete in \'settings.ini\'')

    #--------------------------------------------------------------
    # star_sz_mult, star_max_size
    #--------------------------------------------------------------
    if ('SizeMultiplier' in setting_keys) and ('MaxSize' in setting_keys) and ('SizePowerLaw' in setting_keys):
        if is_num(settings['SizeMultiplier']) and is_num(settings['MaxSize']) and is_num(settings['SizePowerLaw']):
            star_sz_mult = float(settings['SizeMultiplier'])
            star_max_size = float(settings['MaxSize'])
            star_sz_pw = float(settings['SizePowerLaw'])
        else:
            sys.exit('SizeMultiplier or MaxSize or SizePowerlaw invalid in \'settings.ini\'')
    else:
        star_sz_mult = 5
        star_max_size = 25
        star_sz_pw = 8
        
    #--------------------------------------------------------------
    # azimuth_span, altitude_span
    #--------------------------------------------------------------
    if ('AzimuthSpan' in setting_keys) and ('AltitudeSpan' in setting_keys):
        if is_num(settings['AzimuthSpan']) and is_num(settings['AltitudeSpan']):
            azimuth_span = float(settings['AzimuthSpan'])
            altitude_span = float(settings['AltitudeSpan'])
        elif is_num(settings['AzimuthSpan']) and not is_num(settings['AltitudeSpan']):
            azimuth_span = float(settings['AzimuthSpan'])
            altitude_span = azimuth_span*(1/aspect_ratio)
        elif not is_num(settings['AzimuthSpan']) and is_num(settings['AltitudeSpan']):
            altitude_span = float(settings['AltitudeSpan'])
            azimuth_span = altitude_span*aspect_ratio
        else:
            sys.exit('AzimuthSpan and AltitudeSpan invalid in \'settings\'')
    else:
        azimuth_span = 30
        altitude_span = azimuth_span*(1/aspect_ratio)
        
    #--------------------------------------------------------------
    # ldes, bdes, plx_des if star_catalog == 'custom'
    #--------------------------------------------------------------
    if star_catalog == 'custom':
        if not (is_num(settings['RAdeg']) and is_num(settings['DEdeg']) \
                and is_num(settings['parallax'])):
            sys.exit('Custom destination invalid in \'settings.ini\'')
        else:
            ldes, bdes, plx_des = float(settings['RAdeg']),\
                                  float(settings['DEdeg']),\
                                  float(settings['parallax'])
    else:
        ldes, bdes, plx_des = 0,0,0 # will be replaced by catalog value in main.py

    #--------------------------------------------------------------
    # vs_cutoff
    #--------------------------------------------------------------
    if 'VmagCutoff' in setting_keys:
        if is_num(settings['VmagCutoff']):
            vs_cutoff = float(settings['VmagCutoff'])
        else:
            sys.exit('VmagCutoff invalid in \'settings.ini\'')
    else:
        vs_cutoff = 16

    #--------------------------------------------------------------
    # animated
    #--------------------------------------------------------------
    if os.path.isfile(settings_file) and (not 'Animated' in setting_keys):
        sys.exit('Animated option not present in \'settings.ini\'')
    elif not os.path.isfile(settings_file):
        animated = False
    else:
        animated = (settings['Animated']=='yes')
        
    if animated:
        if ('TripStart' in setting_keys) and ('TripEnd' in setting_keys) and ('TripStep' in setting_keys):
            if is_num(settings['TripStart']) and is_num(settings['TripEnd']) and is_num(settings['TripStep']):
                trip_start = float(settings['TripStart'])
                trip_end = float(settings['TripEnd'])
                trip_step = float(settings['TripStep'])
                if trip_start <= 0.001 or trip_end > 0.999:
                    sys.exit('Trip ratio start/end too close to 0 resp. 1 in \'settings.ini\'')
                elif trip_step <= 0:
                    sys.exit('Trip ratio step invalid in \'settings.ini\'')
                else:
                    trip_ratio_lst = np.arange(trip_start, trip_end, trip_step)
            else:
                sys.exit('Trip ratio info invalid in \'settings.ini\'')
        else:
            trip_ratio_lst = np.arange(0.1,0.9,0.1) # 8 frames
            
        if ('FrameFileHead' in setting_keys) and ('FrameFileTail' in setting_keys):
            file_head = settings['FrameFileHead']
            file_tail = settings['FrameFileTail']
        else:
            file_head = 'frame'
            file_tail = ''
            
        if 'FrameFormat' in setting_keys:
            frame_extn = settings['FrameFormat']
        else:
            frame_extn = 'png'
        
    else:
        trip_ratio_lst = [] # won't be used
        file_head = '' # won't be used
        file_tail = '' # won't be used
        frame_extn = '' # won't be used

    return use_latex, star_catalog, pic_width, pic_height, HIP_des, des_name,\
    col_names, star_sz_mult, star_max_size, star_sz_pw, azimuth_span, altitude_span,\
    vs_cutoff, animated, trip_ratio_lst, file_head, file_tail, frame_extn,\
    ldes, bdes, plx_des