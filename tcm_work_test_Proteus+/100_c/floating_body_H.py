import numpy as np
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import WaveTools as wt
import math

# dependencies for FSI
from proteus.mbd import CouplingFSI as fsi
import pychrono

# general options

opts= Context.Options([
    ("T",30.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.5,"Desired CFL restriction"),
    ("he",0.01,"he relative to Length of domain in x"),
    ("fr",1.0,"Forcing frequency ratio"),
    ("f_n",0.6104,"Natural frequency of body roll in hz"),
    ])

# sim time
T = opts.T
# initial step
dt_init = 0.001
# CFL value
cfl = opts.cfl
# mesh size
he = opts.he
# rate at which values are recorded
sampleRate = opts.dt_output
# for ALE formulation
movingDomain = True
# for added mass stabilization
addedMass = True

# physical options
# water density
rho_0 = 998.2
# water kinematic viscosity
nu_0 = 1.004e-6
# air density
rho_1 = 1.205
# air kinematic viscosity
nu_1 = 1.5e-5
# gravitational acceleration
g = np.array([0., -9.81, 0.])

# body options
fixed = False

# forcing frequency ratio
fr = opts.fr
f_n = opts.f_n

# TLD parameters
tld_w = 0.62
tld_h = 0.049
ic_angle = (0./180.)*math.pi

# channel and body state
water_level = 2.
water_length = 10.

# wave options
wave_period = 1/(f_n*fr)
wave_height = 0.05
wave_direction = np.array([1., 0., 0.])
wave_type = 'Fenton'  #'Linear'
# number of Fourier coefficients
Nf = 8
wave = wt.MonochromaticWaves(period=wave_period,
                             waveHeight=wave_height,
                             mwl=water_level,
                             depth=water_level,
                             g=g,
                             waveDir=wave_direction,
                             waveType=wave_type,
                             Nf=8)
wavelength = wave.wavelength

#  ____                        _
# |  _ \  ___  _ __ ___   __ _(_)_ __
# | | | |/ _ \| '_ ` _ \ / _` | | '_ \
# | |_| | (_) | | | | | | (_| | | | | |
# |____/ \___/|_| |_| |_|\__,_|_|_| |_|
# Domain
# All geometrical options go here (but not mesh options)

domain = Domain.PlanarStraightLineGraphDomain()

# ----- SHAPES ----- #

# TANK (=wave channel)
tank = st.Tank2D(domain, dim=(water_length, 1.5*water_level))

# SPONGE LAYERS
# generation zone: 1 wavelength
# absorption zone: 2 wavelengths
tank.setSponge(x_n=2.)
#tank.setSponge(x_n=1., x_p=1.)
#tank.setSponge(x_n=wavelength, x_p=2*wavelength)

# CAISSON

dims = (1., 0.5)
wall_thickness = 0.5*(dims[0]-tld_w)
bottom_height = 0.25
# define vertices
vertices = np.array([
    [0., 0.],
    [dims[0], 0.],
    [dims[0], dims[1]],
    [dims[0]-wall_thickness, dims[1]],
    [dims[0]-wall_thickness, bottom_height],
    [wall_thickness, bottom_height],
    [wall_thickness, dims[1]],
    [0., dims[1]],
])
# give flags to vertices (1 flag per vertex, here all have the same flag)
vertexFlags = np.array([1 for ii in range(len(vertices))])
# define segments
segments = np.array([[ii-1, ii] for ii in range(1, len(vertices))])
# add last segment
segments = np.append(segments, [[len(vertices)-1, 0]], axis=0)
# give flags to segments (1 flag per segment, here all have the same flag)
segmentFlags = np.array([1 for ii in range(len(segments))])
# define regions inside the body
regions = np.array([[dims[0]/2., bottom_height/2.]]) 
regionFlags = np.array([1])
# define holes inside the body
holes = np.array([[dims[0]/2., bottom_height/2.]]) #
regionFlags = np.array([1])
boundaryTags = {'wall': 1}
# barycenter
barycenter = np.array([dims[0]/2., bottom_height, 0.]) #+++ tsao simplify+++#
caisson = st.CustomShape(
    domain=domain,
    vertices=vertices,
    vertexFlags=vertexFlags,
    segments=segments,
    segmentFlags=segmentFlags,
    regions=regions,
    regionFlags=regionFlags,
    holes=holes,
    boundaryTags=boundaryTags,
    barycenter=barycenter,
)
# translate caisson to middle of the tank
caisson.translate(np.array([0.5*tank.dim[0]-0.5*dims[0], water_level-bottom_height]))
#caisson.rotate(rot = ic_angle)

# PMTLD

#pmtld = st.Rectangle(domain, dim=(tld_w, tld_h), coords=(0., 0.))
# set barycenter in middle of caisson
#pmtld.setBarycenter([0., 0.])
# caisson is considered a hole in the mesh
#pmtld.setHoles([[0., 0.]])
# 2 following lines only for py2gmsh
#pmtld.holes_ind = np.array([0])
#tank.setChildShape(pmtld, 0)
# translate pmtld to where you want
#pmtld.translate(np.array([0.5*tank.dim[0], water_level+0.5*tld_h]))

#   ____ _
#  / ___| |__  _ __ ___  _ __   ___
# | |   | '_ \| '__/ _ \| '_ \ / _ \
# | |___| | | | | | (_) | | | | (_) |
#  \____|_| |_|_|  \___/|_| |_|\___/
# Chrono

# SYSTEM

# create system
system = fsi.ProtChSystem()
# access chrono object
chsystem = system.getChronoObject()
# communicate gravity to system
# can also be set with:
# system.ChSystem.Set_G_acc(pychrono.ChVectorD(g[0], g[1], g[2]))
system.setGravitationalAcceleration(g)
# set maximum time step for system
#system.setTimeStep(1e-4) # default in tutorial
system.setTimeStep(5e-5) # my choice

solver = pychrono.ChSolverMINRES()
chsystem.SetSolver(solver)

# BODY

# set body density
thob = 450.

# create floating body
body = fsi.ProtChBody(system=system)
# give it a name
body.setName(b'my_body')
# attach shape: this automatically adds a body at the barycenter of the caisson shape
body.attachShape(caisson)
# set 2D width (for force calculation)
body.setWidth2D(1.)
# access chrono object
chbody = body.getChronoObject()
# impose constraints
chbody.SetBodyFixed(fixed)
free_x = np.array([1., 1., 0.]) # translational
free_r = np.array([0., 0., 1.]) # rotational
body.setConstraints(free_x=free_x, free_r=free_r)
# access pychrono ChBody
# set mass
# can also be set with:
# body.ChBody.SetMass(14.5)

#---tsao calculate it and simplify it as a rectangule---#
#mb = thob*(dims[0]*bottom_height+2.*wall_thickness*(dims[1]-bottom_height)) 
mb = thob*dims[0]*dims[1]
body.setMass(mb)

mbi = mb*(dims[0]*dims[0]+dims[1]*dims[1])/12.
body.setInertiaXX(np.array([1., 1., mbi]))

# set inertia
# can also be set with:
# body.ChBody.setInertiaXX(pychrono.ChVectorD(1., 1., 0.35))
# body.setInertiaXX(np.array([1., 1., 0.35*body.getMass()/14.5]))

# record values
body.setRecordValues(all_values=True)

#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions

# CAISSON
# set no-slip conditions on caisson
for tag, bc in caisson.BC.items():
    bc.setFreeSlip()

# PMTLD
# atmosphere on top
#pmtld.BC['y+'].setAtmosphere()
# free slip on bottom
#pmtld.BC['y-'].setFreeSlip()
# free slip on the right
#pmtld.BC['x+'].setFreeSlip()
# free slip on left
#pmtld.BC['x-'].setFreeSlip()

# fix in space nodes on the boundaries of the pmtld
#for tag, bc in pmtld.BC.items():
#    bc.setFixedNodes()

# TANK
# atmosphere on top
tank.BC['y+'].setAtmosphere()
# free slip on bottom
tank.BC['y-'].setFreeSlip()
# free slip on the right
tank.BC['x+'].setFreeSlip()
# non material boundaries for sponge interface
#tank.BC['x-'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()

# fix in space nodes on the boundaries of the tank
for tag, bc in tank.BC.items():
    bc.setFixedNodes()

# WAVE AND RELAXATION ZONES

smoothing = he*1.5
dragAlpha = 5*2*np.pi/wave_period/(1.004e-6)
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave,
                                               smoothing=smoothing,
                                               vert_axis=1)
tank.setGenerationZones(x_n=True,
                        waves=wave,
                        smoothing=smoothing,
                        dragAlpha=dragAlpha)
#---no absorption zone---#
#tank.setAbsorptionZones(x_p=True,
#                        dragAlpha=dragAlpha)
#tank.setAbsorptionZones(x_n=True,
#                        dragAlpha=dragAlpha)

#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
smoothing = 1.5*he
nd = domain.nd

x_limit = (0.5*tank.dim[0]-0.5*tld_w, 0.5*tank.dim[0]+0.5*tld_w)
y_limit = (water_level, water_level+tld_h)

class P_IC:
    def uOfXT(self, x, t):
        return 0.0
        if x[0]<=x_limit[0]:
            p_L = 0.0
            phi_L = tank.dim[nd-1] - water_level
            phi = x[nd-1] - water_level
            p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
            return p
        elif x[0]>=x_limit[1]:
            p_L = 0.0
            phi_L = tank.dim[nd-1] - water_level
            phi = x[nd-1] - water_level
            p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
            return p
        else:
            p_L = 0.0
            phi_L = tank.dim[nd-1] - (water_level+tld_h)	#---tsao mod---#
            phi = x[nd-1] - (water_level+tld_h)			#---tsao mod---#
            p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
            return p

class U_IC:
    def uOfXT(self, x, t):
        return 0.0
class V_IC:
    def uOfXT(self, x, t):
        return 0.0
class W_IC:
    def uOfXT(self, x, t):
        return 0.0

class VF_IC:
    def uOfXT(self, x, t):

        if x[0]<=x_limit[0]:
            return smoothedHeaviside(smoothing, x[nd-1]-water_level)
        elif x[0]>=x_limit[1]:
            return smoothedHeaviside(smoothing, x[nd-1]-water_level)
        else:
            return smoothedHeaviside(smoothing, x[nd-1]-(water_level+tld_h)) #---tsao mod---#
#        return smoothedHeaviside(smoothing, x[nd-1]-water_level)

class PHI_IC:
    def uOfXT(self, x, t):
#---tsao mod---#
        if x[0]<=x_limit[0]:
            return x[nd-1] - water_level
        elif x[0]>=x_limit[1]:
            return x[nd-1] - water_level
        else:
            return x[nd-1] - (water_level+tld_h)	#---tsao mod---#
#        return x[nd-1] - water_level

# instanciating the classes for *_p.py files
initialConditions = {'pressure': P_IC(),
                     'vel_u': U_IC(),
                     'vel_v': V_IC(),
                     'vel_w': W_IC(),
                     'vof': VF_IC(),
		     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

#---wave gauge---#
#column_gauge_locations = (((1.0,   0.0, 0.0), (1.0, tank.dim[1], 0.0)),
#			  ((2.0,   0.0, 0.0), (2.0, tank.dim[1], 0.0)),
#			  ((4.0,   0.0, 0.0), (4.0, tank.dim[1], 0.0)),
#			  ((6.0,   0.0, 0.0), (6.0, tank.dim[1], 0.0)),
#                          ((8.0,   0.0, 0.0), (8.0, tank.dim[1], 0.0)))

#tank.attachLineIntegralGauges('vof',
#                              gauges=((('vof',),column_gauge_locations),),
#                              fileName='column_gauges.csv')

#  __  __           _        ___        _   _
# |  \/  | ___  ___| |__    / _ \ _ __ | |_(_) ___  _ __  ___
# | |\/| |/ _ \/ __| '_ \  | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |  | |  __/\__ \ | | | | |_| | |_) | |_| | (_) | | | \__ \
# |_|  |_|\___||___/_| |_|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                |_|


domain.MeshOptions.genMesh = True
domain.MeshOptions.he = he
mesh_fileprefix = 'mesh'
domain.MeshOptions.setOutputFiles(mesh_fileprefix)

st.assembleDomain(domain)


#  _   _                           _
# | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___
# |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|
# | |\  | |_| | | | | | |  __/ |  | | (__\__ \
# |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/
# Numerics

outputStepping = TpFlow.OutputStepping(
    final_time=T,
    dt_init=dt_init,
    # cfl=opts.cfl,
    dt_output=sampleRate,
    nDTout=None,
    dt_fixed=None,
)

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(
    ns_model=None,
    ls_model=None,
    nd=domain.nd,
    cfl=cfl,
    outputStepping=outputStepping,
    structured=False,
    he=he,
    domain=domain,
    initialConditions=initialConditions,
)

# Necessary for moving domains
myTpFlowProblem.movingDomain = movingDomain

params = myTpFlowProblem.Parameters

# PHYSICAL PARAMETERS
params.physical.densityA = rho_0  # water
params.physical.densityB = rho_1  # air
params.physical.kinematicViscosityA = nu_0  # water
params.physical.kinematicViscosityB = nu_1  # air
params.physical.gravity = g
params.physical.surf_tension_coeff = 0.

m = myTpFlowProblem.Parameters.Models

# MODEL PARAMETERS
ind = -1
# first model is mesh motion (if any)
if movingDomain:
    m.moveMeshElastic.index = ind+1
    ind += 1
# navier-stokes
m.rans2p.index = ind+1
ind += 1
# volume of fluid
m.vof.index = ind+1
ind += 1
# level set
m.ncls.index = ind+1
ind += 1
# redistancing
m.rdls.index = ind+1
ind += 1
# mass correction
m.mcorr.index = ind+1
ind += 1
# added mass estimation
if addedMass is True:
    m.addedMass.index = ind+1
    ind += 1

# ADD RELAXATION ZONES TO AUXILIARY VARIABLES
m.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
# ADD SYSTEM TO AUXILIARY VARIABLES
m.rans2p.auxiliaryVariables += [system]
m.rans2p.p.coefficients.NONCONSERVATIVE_FORM=0
if addedMass is True:
    # passed in added_mass_p.py coefficients
    m.addedMass.auxiliaryVariables += [system.ProtChAddedMass]
    max_flag = 0
    max_flag = max(domain.vertexFlags)
    max_flag = max(domain.segmentFlags+[max_flag])
    max_flag = max(domain.facetFlags+[max_flag])
    flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
    for s in system.subcomponents:
        if type(s) is fsi.ProtChBody:
            for i in s.boundaryFlags:
                flags_rigidbody[i] = 1
    m.addedMass.p.coefficients.flags_rigidbody = flags_rigidbody
