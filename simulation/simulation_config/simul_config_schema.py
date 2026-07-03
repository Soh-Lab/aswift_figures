# simulation_config/simul_config_schema.py

from pydantic import BaseModel
from typing import List, Tuple


class ModelsConfig(BaseModel):
    generate_data: bool = True
    solve_data: bool = True
    algos: List[str]
    baseline: str
    peak: str
    x: Tuple[float, float, int]
    replicates: int
    constant_noise: bool = False
    inject_laplacian: bool = False


class ParametersConfig(BaseModel):
    amplitudes: List[float]
    peak_widths: List[float]
    peak_loc: Tuple[float, float, int] = (-0.2, -0.2, 1)
    psnrs_db: List[float]

    bl_slope: Tuple[float, float] = (0., 0.)
    bl_intercept: Tuple[float, float] = (0., 0.)


class LocationsConfig(BaseModel):
    simulation_path: str = ''
    save_path: str = ''


class SimulationConfig(BaseModel):
    model: ModelsConfig
    parameters: ParametersConfig
    locations: LocationsConfig
