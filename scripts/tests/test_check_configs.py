#!/usr/bin/env python3
"""
test_check_configs.py

Unit tests for configuration file syntax validator.
"""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

SCRIPT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

import check_configs


class TestCheckConfigs(unittest.TestCase):
    """Test TOML, JSON, and YAML configuration validation."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_valid_toml(self) -> None:
        p = self.root / "valid.toml"
        p.write_text("[section]\nkey = 'value'\n", encoding="utf-8")
        errors = check_configs.check_toml_file(p)
        self.assertEqual(errors, [])

    def test_invalid_toml(self) -> None:
        p = self.root / "invalid.toml"
        p.write_text("[section\nkey = 'missing bracket'\n", encoding="utf-8")
        errors = check_configs.check_toml_file(p)
        self.assertTrue(len(errors) > 0)
        self.assertIn("[TOML]", errors[0])

    def test_valid_json(self) -> None:
        p = self.root / "valid.json"
        p.write_text('{"key": "value", "list": [1, 2, 3]}', encoding="utf-8")
        errors = check_configs.check_json_file(p)
        self.assertEqual(errors, [])

    def test_invalid_json(self) -> None:
        p = self.root / "invalid.json"
        p.write_text('{"key": "unclosed', encoding="utf-8")
        errors = check_configs.check_json_file(p)
        self.assertTrue(len(errors) > 0)
        self.assertIn("[JSON]", errors[0])

    @unittest.skipIf(check_configs.yaml is None, "PyYAML not installed")
    def test_valid_yaml(self) -> None:
        p = self.root / "valid.yml"
        p.write_text("name: ci\non: [push]\njobs:\n  test:\n    runs-on: ubuntu-latest\n", encoding="utf-8")
        errors = check_configs.check_yaml_file(p)
        self.assertEqual(errors, [])

    @unittest.skipIf(check_configs.yaml is None, "PyYAML not installed")
    def test_invalid_yaml(self) -> None:
        p = self.root / "invalid.yml"
        p.write_text("key: [unclosed list\n  bad indentation", encoding="utf-8")
        errors = check_configs.check_yaml_file(p)
        self.assertTrue(len(errors) > 0)
        self.assertIn("[YAML]", errors[0])

    def test_should_skip(self) -> None:
        self.assertTrue(check_configs.should_skip(Path("build/temp.json")))
        self.assertTrue(check_configs.should_skip(Path(".venv/lib/python.toml")))
        self.assertTrue(check_configs.should_skip(Path(".git/config")))
        self.assertFalse(check_configs.should_skip(Path("mise.toml")))
        self.assertFalse(check_configs.should_skip(Path(".github/workflows/ci.yml")))


if __name__ == "__main__":
    unittest.main()
