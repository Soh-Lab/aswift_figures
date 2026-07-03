from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
import tomli

from simulation.peak_simulation import generate_simulation_data, solve_simulation_data
from simulation.simulation_config import load_config, simul_config
from simulation.statistics import summarize_file


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _resolve_path(path: str | Path, base_dir: Path) -> Path:
    path = Path(path)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def _run_simulation(config_path: Path, save: bool, n_workers: int | None, backend: str) -> pd.DataFrame:
    load_config(str(config_path))
    model = simul_config.model

    if model.generate_data:
        logger.info("Generating simulated data from %s", config_path)
        df = generate_simulation_data()
    else:
        input_path = _resolve_path(simul_config.locations.simulation_path, _repo_root())
        logger.info("Loading simulated data from %s", input_path)
        df = pd.read_feather(input_path)

    if model.solve_data:
        logger.info("Fitting simulated data from %s", config_path)
        df = solve_simulation_data(df, n_workers=n_workers, backend=backend)

    if save:
        output_path = _resolve_path(simul_config.locations.save_path, _repo_root())
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_feather(output_path)
        logger.info("Saved simulation output to %s", output_path)

    return df


def run_workflow(config_path: str | Path, run_simulations: bool = True, run_statistics: bool = True) -> None:
    config_path = Path(config_path)
    with config_path.open("rb") as f:
        workflow_config = tomli.load(f)

    base_dir = config_path.parent
    workflow = workflow_config.get("workflow", {})
    default_backend = workflow.get("backend", "pool")
    default_n_workers = workflow.get("n_workers")

    if run_simulations:
        for item in workflow_config.get("simulations", []):
            if not item.get("enabled", True):
                continue
            simulation_config = _resolve_path(item["config"], base_dir)
            logger.info("Running simulation step: %s", item.get("name", simulation_config.stem))
            _run_simulation(
                simulation_config,
                save=item.get("save", True),
                n_workers=item.get("n_workers", default_n_workers),
                backend=item.get("backend", default_backend),
            )

    if run_statistics:
        for item in workflow_config.get("statistics", []):
            if not item.get("enabled", True):
                continue
            logger.info("Running statistics step: %s", item.get("name", item["source"]))
            summarize_file(
                source_path=_resolve_path(item["source"], base_dir),
                output_path=_resolve_path(item["output"], base_dir),
                methods=item["methods"],
                kind=item.get("kind", "error"),
                bad_threshold=item.get("bad_threshold", 0.27),
                drop_baselines=item.get("drop_baselines", []),
            )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run ASWIFT figure simulation and statistics workflow.")
    parser.add_argument(
        "-c",
        "--config",
        default="simulation/workflow_config/paper_simulations.toml",
        help="Path to workflow TOML file.",
    )
    parser.add_argument("--skip-simulations", action="store_true", help="Only run statistics steps.")
    parser.add_argument("--skip-statistics", action="store_true", help="Only run simulation steps.")
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    run_workflow(
        args.config,
        run_simulations=not args.skip_simulations,
        run_statistics=not args.skip_statistics,
    )


if __name__ == "__main__":
    main()
