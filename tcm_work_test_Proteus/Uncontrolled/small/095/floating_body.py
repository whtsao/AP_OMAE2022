import numpy as np
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import WaveTools as wt

# dependencies for FSI
from proteus.mbd import CouplingFSI as fsi
import pychrono

opts= Context.Options([
    ("final_time",30.0,"Final time for simulation"),
    ("dt_output",0.01,"Time interval to output solution"),
    ("cfl",0.5,"Desired CFL restriction"),
    ("he",0.05,"he relative to Length of domain in x"),
    ("tank_dim", (10., 4.), "Dimensions of the tank"),
    ("fr",1.0 ,"Forcing frequency ratio"),
    ])

# general options
# sim time
T = opts.final_time
# initial step
dt_init = 0.001
# CFL value
cfl = opts.cfl #0.5
# mesh size
he = opts.he #0.05 # should be opts.he so I can change it
# rate at which values are recorded
sampleRate = 0.01
# for ALE formulation
movingDomain = True
# for added mass stabilization
addedMass = True
# forcing frequency ratio
fr = opts.fr

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

# wave options
water_level = 2.0

#---tsao del wave---#
wave_period = 1./(fr*0.6104)
wave_height = 0.1
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

# TANK
tank_dim = opts.tank_dim
tank = st.Tank2D(domain, dim=(tank_dim[0], tank_dim[1])) #(2*wavelength, 2*water_level)) #---tsao mod---#

# SPONGE LAYERS
# generation zone: 1 wavelength
# absorption zone: 2 wavelengths
#tank.setSponge(x_n=1., x_p=1.) #---tsao add both side sponge---#
tank.setSponge(x_n=2.) #---tsao add one side sponge---#

caisson = st.Rectangle(domain, dim=(1.0, 0.5), coords=(0., 0.))
# set barycenter in middle of caisson
caisson.setBarycenter([0., 0.])
# caisson is considered a hole in the mesh
caisson.setHoles([[0., 0.]])
# 2 following lines only for py2gmsh
caisson.holes_ind = np.array([0])
tank.setChildShape(caisson, 0)
# translate caisson to middle of the tank
#caisson.translate(np.array([1*wavelength, water_level]))
#caisson.rotate(rot=0.26178) #---tsao mod 15 degree---#
#caisson.rotate(rot=0.1745329252) #---tsao mod 10 degree---#

caisson.translate(np.array([0.5*tank_dim[0], water_level])) #---tsao mod---#

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
system.setTimeStep(1e-4)

solver = pychrono.ChSolverMINRES()
chsystem.SetSolver(solver)

# BODY

# create floating body
body = fsi.ProtChBody(system=system)
# give it a name
body.setName(b'my_body')
# attach shape: this automatically adds a body at the barycenter of the caisson shape
body.attachShape(caisson)
# set 2D width (for force calculation)
body.setWidth2D(1.0)
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
body.setMass(225)
# set inertia
# can also be set with:
# body.ChBody.setInertiaXX(pychrono.ChVectorD(1., 1., 0.35))
body.setInertiaXX(np.array([1., 1., 23.4375]))
# record values
body.setRecordValues(all_values=True,ang_disp=True)


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
    bc.setNoSlip()

# TANK

# atmosphere on top
tank.BC['y+'].setAtmosphere()
# free slip on bottom
tank.BC['y-'].setFreeSlip()
# free slip on the right
tank.BC['x+'].setFreeSlip()
# free fix on the right
#tank.BC['x+'].setFreeSlip() #---tsao mod---#

# free slip on the left #---tsao mod---#
#tank.BC['x-'].setFreeSlip()
# non material boundaries for sponge interface
tank.BC['sponge'].setNonMaterial()
#tank.BC['x-'].setFreeSlip() #---tsao mod---#

# fix in space nodes on the boundaries of the tank
for tag, bc in tank.BC.items():
    bc.setFixedNodes()

# WAVE AND RELAXATION ZONES

#---tsao mod both side absorption zone---#
smoothing = he*1.5
dragAlpha = 1.e6 #5*2*np.pi/wave_period/(1.004e-6)
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave,
                                               smoothing=smoothing,
                                               vert_axis=1)
tank.setGenerationZones(x_n=True,
                        waves=wave,
                        smoothing=smoothing,
                        dragAlpha=dragAlpha)
#tank.setAbsorptionZones(x_p=True,
#                        dragAlpha=dragAlpha)
#tank.setAbsorptionZones(x_n=True,
#                        dragAlpha=dragAlpha)

#--- tsao add ---#
column_gauge_locations = ((( 1.0,  0.0, 0.0), ( 1.0, tank_dim[1], 0.0)),
                          (( 2.0,  0.0, 0.0), ( 2.0, tank_dim[1], 0.0)),
                          (( 4.0,  0.0, 0.0), ( 4.0, tank_dim[1], 0.0)),
                          (( 6.0,  0.0, 0.0), ( 6.0, tank_dim[1], 0.0)),
                          (( 8.0,  0.0, 0.0), ( 8.0, tank_dim[1], 0.0)))

tank.attachLineIntegralGauges('vof',
                              gauges=((('vof',),column_gauge_locations),),
                              fileName='column_gauges.csv')

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

class P_IC:
    def uOfXT(self, x, t):
        p_L = 0.0
        phi_L = tank.dim[nd-1] - water_level
        phi = x[nd-1] - water_level
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
        return smoothedHeaviside(smoothing, x[nd-1]-water_level)
class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - water_level

# instanciating the classes for *_p.py files
initialConditions = {'pressure': P_IC(),
                     'vel_u': U_IC(),
                     'vel_v': V_IC(),
                     'vel_w': W_IC(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

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
    useSuperlu=True
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
