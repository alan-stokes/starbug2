from typing import Final

from starbug2.bin.main import starbug_internal_main
from starbug2.constants import ExitStates
from starbug2.star_bug_config import StarBugMainConfig
from tests.generic import TEST_PATH
import os


TEST_NGC_FITS: Final[str] = str(
    os.path.join(str(TEST_PATH), "ngc6822_F770W_i2d.fits"))

class TestSystemResults:

    def test_detection_on_proper_fits_file(self, capsys):
        config: StarBugMainConfig = StarBugMainConfig()
        config.custom_filter = 'F770W'
        config.fits_images = [TEST_NGC_FITS]
        config.do_star_detection = True
        config.verbose_logs = True
        exit_state: int = starbug_internal_main(config)
        assert exit_state == ExitStates.EXIT_SUCCESS

        captured = capsys.readouterr()
        lines = captured.out.splitlines()

        plain_line: str | None = None
        background_line: str | None = None
        con_vol_line: str | None = None
        cleaning: str | None = None
        total_line: str | None = None

        for line in lines:
            if "[PLAIN] pass:" in line:
                plain_line = line
            if "[BGD2D] pass:" in line:
                background_line = line
            if "[CONVL] pass" in line:
                con_vol_line = line
            if "cleaning" in line:
                cleaning = line
            if "Total" in line:
                total_line = line

        assert (
            plain_line is not None and background_line is not None
            and con_vol_line is not None and cleaning is not None
            and total_line is not None)

        plain_count: int | None = None
        background_count: int | None = None
        con_vol_count: int | None = None
        cleaning_count: int | None = None
        total_line_count: int | None = None


        count_part: str = plain_line.split("pass: ")[1]
        plain_count: int = int(count_part.split()[0])

        count_part = background_line.split("pass: ")[1]
        background_count: int = int(count_part.split()[0])

        count_part = con_vol_line.split("pass: ")[1]
        con_vol_count: int = int(count_part.split()[0])

        count_part = cleaning.split("cleaning ")[1]
        cleaning_count: int = int(count_part.split()[0])

        count_part = total_line.split("Total: ")[1]
        total_line_count: int = int(count_part.split()[0])

        expected_plain: int = 4461
        expected_background: int = 4508
        expected_con_vol: int = 21854
        expected_cleaning: int = 3433
        expected_total: int = 18421

        assert (expected_plain * 0.8) <= plain_count <= (expected_plain * 1.2)
        assert ((expected_background * 0.8) <= background_count
                <= (expected_background * 1.2))
        assert ((expected_con_vol * 0.8) <= con_vol_count
                <= (expected_con_vol * 1.2))
        assert ((expected_cleaning * 0.8) <= cleaning_count
                <= (expected_cleaning * 1.2))
        assert ((expected_total * 0.8) <= total_line_count
                <= (expected_total * 1.2))


