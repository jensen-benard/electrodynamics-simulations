from pint import UnitRegistry
import numpy as np
import pyvista as pv

ureg = UnitRegistry()


## Setup the simulation
q = 0.1 * ureg.coulomb
qPos1 = np.array([0, 0, 0,]) * ureg.meter

k = 1 / (4 * np.pi *ureg.electric_constant)

shellRadius = 0.1 * ureg.meter
shellRadius = shellRadius.to_base_units().magnitude

# Helper functions
def getFieldVector(position, charge):
    position = position

    norm = np.linalg.norm(position)
    if norm < shellRadius:
        return (np.array([0, 0, 0]) * ureg.newton / ureg.coulomb).to_base_units().magnitude

    magnitude = (k * charge  / norm**3).to_base_units().magnitude
    return (magnitude * position)

## Setup the graph
# DASHBOARD
GRID_SIZE_X = 1 * ureg.meter
GRID_SIZE_Y = 1 * ureg.meter
GRID_SIZE_Z = 1 * ureg.meter

CELL_COUNT_X = 30
CELL_COUNT_Y = 30
CELL_COUNT_Z = 30

## Create the grid
cellSizeX = GRID_SIZE_X / CELL_COUNT_X
cellSizeY = GRID_SIZE_Y / CELL_COUNT_Y
cellSizeZ = GRID_SIZE_Z / CELL_COUNT_Z

axisX = (np.arange(1, CELL_COUNT_X + 1) * cellSizeX).to_base_units().magnitude
axisY = (np.arange(1, CELL_COUNT_Y + 1) * cellSizeY).to_base_units().magnitude
axisZ = (np.arange(1, CELL_COUNT_Z + 1) * cellSizeZ).to_base_units().magnitude

axisX = axisX - np.mean(axisX)
axisY = axisY - np.mean(axisY)
axisZ = axisZ - np.mean(axisZ)

x, y, z = np.meshgrid(axisX, axisY, axisZ, indexing='ij')

grid = pv.StructuredGrid(x, y, z)

# Calculate the electric 
qPos1 = qPos1.to_base_units().magnitude
grid["E"] = np.array([getFieldVector(p - qPos1, q) for p in grid.points])

# Plotting
plotter = pv.Plotter()
plotter.add_arrows(grid.points, grid["E"], mag=0.000000000001)
sphereMesh1 = pv.Sphere(radius=shellRadius, center=qPos1)
plotter.add_mesh(sphereMesh1, color='red', opacity=0.5)
plotter.show_bounds(
    grid='back',          
    ticks='outside',      
    xtitle='X-axis', ytitle='Y-axis', ztitle='Z-axis',
    font_size=15,
    n_xlabels=5,          
    n_ylabels=5,          
    n_zlabels=6,          
    location='outer',
    minor_ticks=True,
    use_3d_text=True      
)

plotter.show()
