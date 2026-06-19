from starbug2.bin.main import starbug_internal_main
from starbug2.constants import ExitStates
from starbug2.star_bug_config import StarBugMainConfig
from tests.generic import TEST_NGC_FITS


class TestSystemResults:

    def test_detection_on_proper_fits_file(self, capsys):
        config: StarBugMainConfig = StarBugMainConfig()
        config.custom_filter = 'F770W'
        config.fits_images = [TEST_NGC_FITS]
        config.do_star_detection = True
        config.verbose_logs = True
        exit_state: ExitStates = starbug_internal_main(config)
        assert exit_state == ExitStates.EXIT_SUCCESS

        captured = capsys.readouterr()
        assert "[PLAIN] pass: 4461 sources" in captured.out
        assert "[BGD2D] pass: 4508 sources" in captured.out
        assert "[CONVL] pass: 21854 sources" in captured.out

        assert "cleaning 3433 unlikely point sources" in captured.out
        assert "Total: 18421 sources" in captured.out
