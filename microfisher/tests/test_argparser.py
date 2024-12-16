
import os
from unittest.mock import patch, MagicMock
from src.microfisher.microfisher import parse_output_dir


def test_parse_output_dir_abs_output():
    args = MagicMock()
    args.out_dir = "/tmp/test_dir"
    args.workspace = "work_dir"

    with patch("os.makedirs") as mock_makedirs:
        args = parse_output_dir(args)
        # mock_makedirs.assert_called_once_with("test_dir")
        assert args.out_dir == os.path.join("/tmp/test_dir")


def test_parse_output_dir_in_workspace():
    args = MagicMock()
    args.out_dir = "test_dir"
    args.workspace = "work_dir"

    with patch("os.makedirs") as mock_makedirs:
        args = parse_output_dir(args)
        assert args.out_dir == os.path.join("work_dir", "test_dir")


def test_parse_output_dir_abs_workspace():
    args = MagicMock()
    args.out_dir = "test_dir"
    args.workspace = "/tmp/work_dir"

    with patch("os.makedirs") as mock_makedirs:
        args = parse_output_dir(args)
        assert args.out_dir == os.path.join("/tmp/work_dir", "test_dir")


def test_parse_output_dir_abs_both():
    args = MagicMock()
    args.out_dir = "/tmp/test_dir"
    args.workspace = "/tmp/work_dir"

    with patch("os.makedirs") as mock_makedirs:
        args = parse_output_dir(args)
        assert args.out_dir == os.path.join("/tmp/test_dir")


def test_parse_output_dir_file_exists():
    args = MagicMock()
    args.out_dir = "test_dir"
    args.workspace = "work_dir"

    with patch("os.path.realpath", return_value="test_dir"), \
         patch("os.makedirs", side_effect=FileExistsError):
        parse_output_dir(args)
        # No exception should be raised


def test_parse_output_dir_permission_error():
    args = MagicMock()
    args.out_dir = "test_dir"
    args.workspace = "work_dir"

    with patch("os.path.realpath", return_value="test_dir"), \
         patch("os.makedirs", side_effect=PermissionError), \
         patch("sys.exit") as mock_exit:
        parse_output_dir(args)
        mock_exit.assert_called_once_with(-1)
