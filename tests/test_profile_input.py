"""
Tests whether VARGRAM output data is the same when equivalent inputs are given.
Nextclade CLI must be installed.
"""

from vargram import vargram

class TestFileInput():

    def setup_method(self):
        """Obtain output data for equivalent inputs."""

        # Using local data
        vg = vargram(seq='tests/test_data/batches', 
                     ref='tests/test_data/sc2_wuhan_2019.fasta', 
                     gene='tests/test_data/sc2.gff')
        vg.profile(threshold=1, ytype='counts')
        vg.save('tests/test_data/vargram_analysis.csv')
        self.local_ref = vg.stat()

        # Using VARGRAM output
        vg = vargram(data='tests/test_data/vargram_analysis.csv',
                     format='vargram')
        vg.profile()
        self.stat = vg.stat()

        # Using Nextclade reference
        vg = vargram(seq='tests/test_data/batches', 
                     ref='sars-cov-2', ytype='counts')
        vg.profile(threshold=1, ytype='counts')
        self.online_ref = vg.stat()

        # Using Nextclade CLI analysis output
        vg = vargram(data='tests/test_data/batches_analysis_cli.tsv')
        vg.profile(threshold=1, ytype='counts')
        self.analysis_file = vg.stat()

    def test_stat(self):
        """Output using local data must equal to output when VARGRAM analysis output is provided as input."""
        assert self.local_ref.equals(self.stat)

    def test_ref(self):
        """Output using local reference must equal output using downloaded reference."""
        assert self.local_ref.equals(self.online_ref)
    
    def test_analysis(self):
        """Output using local data must equal output using Nextclade CLI analysis output as input."""
        shared_columns = ['gene', 'mutation', 'position', 'type', 'sum']
        sh_local_ref = self.local_ref[shared_columns]
        sh_analysis_file = self.analysis_file[shared_columns]
        assert sh_local_ref.equals(sh_analysis_file)

    