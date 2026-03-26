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


def test_build_examples_use_fused_regional_and_cbg_residential_layer_roles():
    for name, relative_path in (
        ("build_sfbay", "examples/build_sfbay.py"),
        ("build_seattle", "examples/build_seattle.py"),
    ):
        module = _load_module(name, relative_path)
        layers = module.osm_config["graph_layers"]

        assert layers["backbone"]["layer_role"] == "backbone"
        assert layers["backbone"]["geo_level"] == "county"
        assert layers["connector"]["layer_role"] == "connector"
        assert layers["connector"]["geo_level"] == "county"

        if "ferry" in layers:
            assert layers["ferry"]["layer_role"] == "ferry"
            assert layers["ferry"]["geo_level"] == "county"

        assert layers["residential"]["layer_role"] == "residential"
        assert layers["residential"]["geo_level"] == "cbg"


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


def test_diagnose_osm_pbf_exposes_help():
    result = subprocess.run(
        [sys.executable, str(ROOT / "examples/diagnose_osm_pbf.py"), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "diagnose_osm_pbf.py --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "--epsg-utm" in result.stdout
    assert "--graph-path" in result.stdout
    assert "--osm-xml" in result.stdout


def test_compare_osm_pbf_exposes_help():
    result = subprocess.run(
        [sys.executable, str(ROOT / "examples/compare_osm_pbf.py"), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "compare_osm_pbf.py --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "--epsg-utm-a" in result.stdout
    assert "--epsg-utm-b" in result.stdout
    assert "--graph-a" in result.stdout
    assert "--graph-b" in result.stdout


def test_map_network_csv_to_osm_pbf_exposes_help():
    result = subprocess.run(
        [sys.executable, str(ROOT / "examples/map_network_csv_to_osm_pbf.py"), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "map_network_csv_to_osm_pbf.py --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "--network-osm-id-col" in result.stdout


def test_download_sfbay_counties_exposes_help():
    result = subprocess.run(
        [sys.executable, str(ROOT / "examples/download_sfbay_counties.py"), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "download_sfbay_counties.py --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "--work-dir" in result.stdout
    assert "--output-path" in result.stdout
    assert "--epsg" in result.stdout


def test_download_sfbay_masks_exposes_help():
    result = subprocess.run(
        [sys.executable, str(ROOT / "examples/download_sfbay_masks.py"), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "download_sfbay_masks.py --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "--land-output-path" in result.stdout
    assert "--whole-area-output-path" in result.stdout
    assert "--buffer-m" in result.stdout


def test_intersect_zones_exposes_help():
    result = subprocess.run(
        [sys.executable, str(ROOT / "examples/intersect_zones.py"), "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "intersect_zones.py --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "--output-epsg" in result.stdout
    assert "--output-path" in result.stdout


def test_main_module_exposes_cli_help():
    env = os.environ.copy()
    src_path = str(ROOT / "src")
    env["PYTHONPATH"] = src_path if not env.get("PYTHONPATH") else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"

    result = subprocess.run(
        [sys.executable, "-m", "osm_chordify.main", "--help"],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "python -m osm_chordify.main --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "build" in result.stdout
    assert "intersect" in result.stdout
    assert "map" in result.stdout
    assert "diagnose" in result.stdout
    assert "diagnose-built" in result.stdout
    assert "compare-pbf" in result.stdout
    assert "map-pbf" in result.stdout


def test_main_module_exposes_compare_pbf_help():
    env = os.environ.copy()
    src_path = str(ROOT / "src")
    env["PYTHONPATH"] = src_path if not env.get("PYTHONPATH") else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"

    result = subprocess.run(
        [sys.executable, "-m", "osm_chordify.main", "compare-pbf", "--help"],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "python -m osm_chordify.main compare-pbf --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "--pbf-a" in result.stdout
    assert "--pbf-b" in result.stdout
    assert "--epsg-utm-a" in result.stdout


def test_main_module_exposes_map_pbf_help():
    env = os.environ.copy()
    src_path = str(ROOT / "src")
    env["PYTHONPATH"] = src_path if not env.get("PYTHONPATH") else f"{src_path}{os.pathsep}{env['PYTHONPATH']}"

    result = subprocess.run(
        [sys.executable, "-m", "osm_chordify.main", "map-pbf", "--help"],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "python -m osm_chordify.main map-pbf --help failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "--network-csv-path" in result.stdout
    assert "--osm-pbf-path" in result.stdout
    assert "--output-path" in result.stdout


def test_diagnose_osm_pbf_runs_as_a_script_without_executing_diagnostics(tmp_path):
    env = os.environ.copy()
    env["OSM_CHORDIFY_SKIP_DIAGNOSE"] = "1"
    pbf_path = tmp_path / "sample.osm.pbf"
    pbf_path.write_bytes(b"placeholder")

    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "examples/diagnose_osm_pbf.py"),
            str(pbf_path),
            "--epsg-utm",
            "26910",
        ],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, (
        "diagnose_osm_pbf.py failed to import/run as a script.\n"
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
