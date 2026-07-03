# config/simul_config_schema.py

from pydantic import BaseModel
from typing import List


class ParametersConfig(BaseModel):
    fitting_methods: List['str']
    background_method: 'str'
    peak_method: 'str'
    hz_values: List[int]
    polynomial_penalty: float = 0.0
    polynomial_max_degree: int = 4
    calibration_lower_idx: int
    calibration_upper_idx: int
    bg_buffer: int = 0
    baseline_boundary: float = 0.05
    pre_filter: bool = False
    huber_reweight: bool = True
    use_file0: bool = True
    full_cutoff: float = 1.
    peak_lambda_scale: float = 1.0
    peak_prominence: float = 0.5
    huber_cutoff: float = 2.0
    mad_window: int = 21
    filter_outliers: bool = False
    outliers_cutoff: float = 10.0


class LocationsConfig(BaseModel):
    input_dir: 'str'


class Config(BaseModel):
    parameters: ParametersConfig
    locations: LocationsConfig
