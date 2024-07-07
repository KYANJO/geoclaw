# ----------------------------------------------
# @author:  Brian Kyanjo
# @contact: briankyanjo@u.boisestate.edu
# @date:    2022-10-16
# @version: 1.0
# @desc:    This script is to produce geflood.ini file for the geoflood model
# ------------------------------------------------

# Importing required libraries
from configparser import ConfigParser
 
class GeoFlooddata(object):
    # Several geoflood attributes (ignore for now)

    def __init__(self):
        self.user = {}

        self.minlevel = 0
        self.maxlevel = 0

        self.regrid_interval = 1
        self.refine_threshold = 0.5
        self.coarsen_threshold = 0.5
        self.subcycle = False
        self.output = True
        self.output_gauges = True
        self.cuda = False
        self.verbosity = 'essential'
        self.smooth_refine = False
        self.advance_one_step = False
        self.outstyle_uses_maxlevel = True

        self.mi = 1
        self.mj = 1

        #things we didnt set in fclaw options and geoflood.py
        


    def write(self,rundata):
        geoflood = ConfigParser(allow_no_value=True)
        geoflood.optionxform = str # To maintain original case of comments. 

        clawdata = rundata.clawdata
        geo_data = rundata.geo_data
        refinement_data = rundata.refinement_data
        amrdata = rundata.amrdata

        user = {
        '   # User defined parameters' : None,
        '   cuda' : self.cuda,
        }
        for k in self.user.keys():
            user[f'   {k:}'] = self.user[k]

        geoflood['user'] = user

        geoflood['clawpatch'] = {"\n"
        '   # Grid dimensions' : None,
        '   mx' : clawdata.num_cells[0],

        '   my' : clawdata.num_cells[1], "\n"

        '   # Number of ghost cells': None,
        '   mbc': clawdata.num_ghost, "\n"

        '   # Number of auxiliary variables' : None,
        '   maux': clawdata.num_aux,"\n"

        '   # Number of equations' : None,
        '   meqn': clawdata.num_eqn        
            }

        if clawdata.output_step_interval is None:
            clawdata.output_step_interval = 1


        assert self.verbosity in ['silent','essential','production','info','debug'], \
            "Error (geoflood) : Invalid value for verbosity"


        geoflood['Options'] = {"\n"

        '# Regridding options' : None,
        '   minlevel' : self.minlevel,
        '   maxlevel' : self.maxlevel,"\n"

        "   # Regrid every 'regrid_interval' time steps using threshold values":None, 
        '   regrid_interval' : self.regrid_interval,
        '   refine_threshold' : self.refine_threshold,
        '   coarsen_threshold' : self.coarsen_threshold, "\n"

        "   # Smooth refinement (around finest level)":None,
        '   smooth-refine' : self.smooth_refine,
        '   smooth-level' : self.maxlevel,"\n"


        '# Time stepping' : None, 
        '   # Final time':None,
        '   tfinal': clawdata.tfinal, "\n"     

        '   # Take a fixed time step':None,
        '   use_fixed_dt': not clawdata.dt_variable,"\n"

        "   # Initial time step for 'minlevel'":None,
        '   initial_dt': clawdata.dt_initial,"\n"

        '   # CFL constraints : Timestep will never exceed max_cfl and will ' "\n"
        '   # try not to exceed desired_cfl' : None,
        '   max_cfl': clawdata.cfl_max,
        '   desired_cfl': clawdata.cfl_desired,"\n"

        '   # 1 : Output steps  = tfinal/nout;':None,                  
        '   # 2 : not implemented;':None,                              
        '   # 3 : Take nout steps;  save files every nstep steps.':None,
        '   outstyle': clawdata.output_style,"\n" 

        '   # Used for all three out styles;  has different meaning, though':None,                                     
        '   nout': clawdata.num_output_times,"\n"
        '   # Only used if outstyle is 3':None,
        '   nstep': clawdata.output_step_interval,"\n"

        '   # Advanced time stepping' : None, 
        '   subcycle' : self.subcycle, "\n"

        '# File and console IO' : None,
        '   output' : self.output,
        '   output-gauges' : self.output_gauges,
        '   smooth-refine' : "False",
        '   advance-one-step' : self.advance_one_step,
        '   outstyle-uses-maxlevel' : self.outstyle_uses_maxlevel,
        '   verbosity' : self.verbosity,"\n"


        '# Domain geometry' : None,"\n"
        '   # Domain dimensions [ax,bx]x[ay,by]' : None,
        '   ax': clawdata.lower[0],
        '   bx': clawdata.upper[0],
        '   ay': clawdata.lower[1],
        '   by': clawdata.upper[1],"\n"

        '   # Block dimensions' : None,
        '   mi': self.mi,
        '   mj': self.mj
        }

        #mthbc
        mthbc_in  = [clawdata.bc_lower[0], clawdata.bc_upper[0], 
                     clawdata.bc_lower[1], clawdata.bc_upper[1]]
        mthbc = [0]*4
        mthbc_str = ""
        for i in range(4):
            if mthbc_in[i] in [0,'user']:        
                mthbc[i] = 0
            elif mthbc_in[i] in [1,'extrap']:    
                mthbc[i] = 1
            elif mthbc_in[i] in [2,'periodic']:  
                mthbc[i] = 2
            elif mthbc_in[i] in [3,'wall']:      
                mthbc[i] = 3
            else:
                mthbc_in[i] = mthbc[i]
            mthbc_str += " " + str(mthbc[i])

        #limiters
        lim = [0]*clawdata.num_waves
        lim_str = "" 
        for i in range(clawdata.num_waves):
            if clawdata.limiter[i] in [0,'none']:        
                lim[i] = 0
            elif clawdata.limiter[i] in [1,'minmod']:    
                lim[i] = 1
            elif clawdata.limiter[i] in [2,'superbee']:  
                lim[i] = 2
            elif clawdata.limiter[i] in [3,'vanleer']:   
                lim[i] = 3
            elif clawdata.limiter[i] in [4,'mc']:        
                lim[i] = 4
            else:
                clawdata.limiter[i] = lim[i]
            lim_str += " " + str(lim[i])

        #order 
        ord = [0]*2
        ord_str = ""
        if clawdata.order in [1,'godunov','Godunov']:
            ord[0] = 1
        if clawdata.order in [2,'Lax-Wendroff']:
            ord[0] = 2
        ord_str = str(ord[0])

        #ord_in = clawdata.transverse_waves
        if clawdata.transverse_waves in ['none' , 0]:
            ord[1] = 0
        if clawdata.transverse_waves in ['increment' , 1]:
            ord[1] = 1
        if clawdata.transverse_waves in ['all' , 2]:
            ord[1] = 2
        ord_str += " " + str(ord[1])

        # Source terms splitting:
        #   src_split == 0 or 'none'    ==> no source term (src routine never called)
        #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
        #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.

        if clawdata.source_split in [0, 'none']:
            src_split = 0
        if clawdata.source_split in [1, 'godunov', 'Godunov']:
            src_split = 1
        if clawdata.source_split in [2, 'strang', 'Strang']:
            src_split = 2

        if clawdata.output_format in [1,'ascii']:      # 'ascii' or 'binary' check line 678, "data.py", for explanation
            ascii_out = 'T'
        else:
            ascii_out = 'F'

        #print(clawdata.output_format)
        #print(ascii_out)


        geoflood['geoclaw'] = {
            '   # normal and transverse order': None,
            '   # Order of accuracy:': None,
            '   #   1 => Godunov,': None,  
            '   #   2 => Lax-Wendroff plus limiters': None,

            '   order': ord_str,"\n"

            '   # Location of capacity function in auxiliary array' : None,
            '   mcapa': clawdata.capa_index,"\n"

            '   # Source term splitting' : None,
            '   src_term': src_split,"\n"

            '   # Use an f-waves update (default : True)'
            '   use_fwaves' : clawdata.use_fwaves,"\n"

            '   # Number of waves': None,
            '   mwaves': clawdata.num_waves,"\n"


            "   # mthlim (is a vector in general, with 'mwaves' entries": None,
            '   # List of limiters to use for each wave family:': None,
            '   # Required:  len(limiter) == num_waves': None,
            '   # Some options:': None,
            "   #   0 or 'none'     ==> no limiter (Lax-Wendroff)": None,
            "   #   1 or 'minmod'   ==> minmod": None,
            "   #   2 or 'superbee' ==> superbee": None,
            "   #   3 or 'mc'       ==> MC limiter": None,
            "   #   4 or 'vanleer'  ==> van Leer": None,
            '   mthlim' : lim_str,"\n"

            '   # mthbc (=left,right,bottom,top)' : None,
            '   # Choice of BCs at xlower and xupper:':None,
            '   # 0 => user specified (must modify bcN.f to use this option)':None,
            '   # 1 => extrapolation (non-reflecting outflow)':None,
            '   # 2 => periodic (must specify this at both boundaries)':None,
            '   # 3 => solid wall for systems where q(2) is normal velocity':None,
            '   mthbc' : mthbc_str,"\n"

            '   dry_tolerance_c': geo_data.dry_tolerance,
            '   wave_tolerance_c': refinement_data.wave_tolerance, 
            '   speed_tolerance_c': refinement_data.speed_tolerance, "\n"


            '   # Output' : None,
            '   ascii-out': ascii_out
            }
            
        with open('geoflood.ini','w') as geofloodfile:
            geoflood.write(geofloodfile)