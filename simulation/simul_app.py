import os
import logging
import json
import argparse
import pandas as pd

from simulation.simulation_config import load_config, simul_config
from simulation.peak_simulation import generate_simulation_data, solve_simulation_data


# Set up logger
logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s: %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def parse_arguments():
    """ Args parser for config and visualization flag """

    parser = argparse.ArgumentParser(description='echem peak simulator')
    parser.add_argument('-c', '--config', type=str, help='Path to configuration TOML file', required=True)
    parser.add_argument('--display', default=False, action='store_true', help='Visualize results?')
    parser.add_argument('--save', default=False, action='store_true', help='Save simulated data?')
    return parser.parse_args()


def main():
    """ Data simulation and solvers """
    # Parse arguments
    args = parse_arguments()
    load_config(args.config)
    model = simul_config.model

    # Simulate data
    if model.generate_data:
        logger.info('Generating simulated data')
        df = generate_simulation_data()
    else:
        logger.info('Loading simulated data')
        df = pd.read_feather(simul_config.locations.simulation_path)

    # Solve for simulated peak and baseline
    if model.solve_data:
        logger.info('Fitting simulated data')
        df = solve_simulation_data(df)

    # Save the results
    if args.save:
        feather_path = simul_config.locations.save_path
        df.to_feather(feather_path)
        logger.info(f"Output saved in: {feather_path}")

    # Visualize the result
    if args.display:
        logger.info('Visualizing results')


if __name__ == '__main__':
    main()