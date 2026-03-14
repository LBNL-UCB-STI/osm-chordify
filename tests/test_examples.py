import ast
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

from osm_chordify.utils.geo import name_osm_network


ROOT = Path(__file__).resolve().parents[1]


def _load_module(name: str, relative_path: str):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _parse_imports(relative_path: str):
    source = (ROOT / relative_path).read_text()
    return ast.parse(source)


def test_build_examples_define_string_custom_filters():
    for name, relative_path in (
        ("build_sfbay", "examples/build_sfbay.py"),
        ("build_seattle", "examples/build_seattle.py"),
    ):
        module = _load_module(name, relative_path)
        for layer_name, layer_cfg in module.osm_config["graph_layers"].items():
            if "custom_filter" in layer_cfg:
                assert isinstance(layer_cfg["custom_filter"], str), (
                    f"{relative_path}:{layer_name} custom_filter must be a string"
                )


def test_diagnose_sfbay_tracks_build_sfbay_network_name():
    build_sfbay = _load_module("build_sfbay", "examples/build_sfbay.py")
    diagnose_sfbay = _load_module("diagnose_sfbay", "examples/diagnose_sfbay.py")

    expected_name = name_osm_network(
        build_sfbay.area_config["name"],
        build_sfbay.osm_config["graph_layers"],
        build_sfbay.osm_config["strongly_connected_components"],
    )

    assert diagnose_sfbay.osm_name == expected_name
    assert diagnose_sfbay.utm_epsg == build_sfbay.geo_config["utm_epsg"]
    assert diagnose_sfbay.pbf_path.endswith(f"{expected_name}.osm.pbf")


def test_examples_use_public_package_imports():
    for relative_path in (
        "examples/build_sfbay.py",
        "examples/build_seattle.py",
        "examples/common.py",
    ):
        tree = _parse_imports(relative_path)
        public_imports = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module == "osm_chordify"
            for alias in node.names
        }
        if relative_path == "examples/common.py":
            assert "build_osm_by_pop_density" in public_imports
            assert "create_osm_highway_filter" in public_imports


def test_diagnose_sfbay_runs_as_a_script_without_executing_diagnostics():
    env = os.environ.copy()
    env["OSM_CHORDIFY_SKIP_DIAGNOSE"] = "1"

    result = subprocess.run(
        [sys.executable, str(ROOT / "examples/diagnose_sfbay.py")],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "diagnose_sfbay.py failed to import/run as a script.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )


def test_build_examples_expose_cli_help():
    for script_name in ("build_sfbay.py", "build_seattle.py"):
        result = subprocess.run(
            [sys.executable, str(ROOT / "examples" / script_name), "--help"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        assert result.returncode == 0, (
            f"{script_name} --help failed.\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        assert "--census-api-key" in result.stdout
        assert "--skip-validation" in result.stdout
